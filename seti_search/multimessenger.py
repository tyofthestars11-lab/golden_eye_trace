"""The strict radio-time gate: a coincidence needs a clock, not a date.

The neutrino lane can only answer its question if the radio side supplies two
things: a position and a time with a resolution. It currently supplies neither,
and the failure this module prevents is the quiet substitution of one for the
other - treating ``observation_date`` as a timestamp, and then reporting a
coincidence whose window is a day long because the file is named by day.

So the status is graded and it is *strict*: an observation with a date and no
time of day is ``MISSING``, not ``PARTIAL``, because a date is not a clock. The
distinction matters arithmetically: a 24-hour window in a 4.68 degree cone at
100 alerts a year carries a chance probability of 4.6e-4, which is 300 times the
background of the survey's own 256 s integration. Calling that "a timestamp"
would buy a coincidence with a box nobody measured.

The vocabulary is closed, like every other one in this repository:

``radio_time_status``
    ``MISSING``    no usable event time: absent, or a date with no time of day
    ``UTC_START``  a timezone-aware start and a positive duration, but no position
    ``FULL``       start, duration and a position with an uncertainty

``multimessenger_status``
    ``NOT TESTABLE``  the question cannot be asked, and the reason is named
    ``TESTABLE``      both ingredients exist; the query has not been run
    ``TESTED``        the query was run against a catalog

Two rules are enforced rather than trusted:

1. **A coincidence may not be reported for a ``NOT TESTABLE`` record.**
   ``enforce`` raises, so a crossmatch that produced ``COINCIDENT`` while the
   time was missing is a bug that stops the run instead of a row in a table.
2. **``TESTED`` needs a position as well as a time.** A time match with no cone
   is not a coincidence; it is a coincidence with the whole sky.

Nothing here can raise a disposition, and nothing here changes a radio verdict.
A ``FULL`` record is testable, not confirmed, and the neutrino lane still may
never restore a record the radio pipeline excluded.
"""

from __future__ import annotations

import datetime as dt
import math

from .sensitivity import background_events

__all__ = [
    "Assessment",
    "MULTIMESSENGER_STATUSES",
    "MultimessengerNotTestable",
    "RADIO_TIME_STATUSES",
    "assess",
    "enforce",
    "require_testable",
]

#: Closed vocabularies. A status outside them is a caller's bug.
RADIO_TIME_STATUSES = ("MISSING", "UTC_START", "FULL")
MULTIMESSENGER_STATUSES = ("NOT TESTABLE", "TESTABLE", "TESTED")

#: The reference rate and cone used to price what a missing timestamp costs.
#: They are the same numbers the power analysis uses, so the "a day is not a
#: timestamp" argument is quoted in the units the rest of the lane works in.
REFERENCE_RATE_PER_YEAR = 100.0
REFERENCE_CONE_DEG = 4.68
DAY_S = 86400.0


class MultimessengerNotTestable(Exception):
    """Raised when a coincidence is claimed for a record without a clock."""


def _parse_utc(raw, what: str) -> dt.datetime:
    if not isinstance(raw, str) or not raw.strip():
        raise ValueError(f"{what} must be an ISO 8601 string")
    try:
        stamp = dt.datetime.fromisoformat(raw.strip().replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError(f"{what} {raw!r} is not ISO 8601: {exc}")
    if stamp.tzinfo is None or stamp.utcoffset() is None:
        raise ValueError(f"{what} {raw!r} is naive: a wall-clock time with no "
                         "zone is not a timestamp")
    return stamp.astimezone(dt.timezone.utc)


def _position(record) -> tuple[float | None, float | None, float | None]:
    metadata = getattr(record, "metadata", {}) or {}
    ra = metadata.get("ra_deg")
    dec = metadata.get("dec_deg")
    err = metadata.get("radio_localisation_deg")
    if ra is None or dec is None:
        return None, None, None
    return float(ra), float(dec), (None if err is None else float(err))


def _start(record) -> dt.datetime | None:
    metadata = getattr(record, "metadata", {}) or {}
    raw = metadata.get("start_time_utc")
    if not raw:
        return None
    try:
        return _parse_utc(raw, "start_time_utc")
    except ValueError:
        return None


class Assessment:
    """What the radio record can and cannot support, with the reason named."""

    __slots__ = ("radio_time_status", "multimessenger_status", "window_s",
                 "pointing_status", "reasons", "start_utc")

    def __init__(self, *, radio_time_status: str, multimessenger_status: str,
                 window_s: float | None, pointing_status: str,
                 reasons: tuple[str, ...],
                 start_utc: dt.datetime | None) -> None:
        if radio_time_status not in RADIO_TIME_STATUSES:
            raise ValueError(f"radio_time_status {radio_time_status!r} is not "
                             f"one of {RADIO_TIME_STATUSES}")
        if multimessenger_status not in MULTIMESSENGER_STATUSES:
            raise ValueError(f"multimessenger_status {multimessenger_status!r} "
                             f"is not one of {MULTIMESSENGER_STATUSES}")
        if multimessenger_status != "NOT TESTABLE" and not reasons:
            raise ValueError("a testable status still needs its ingredients "
                             "recorded")
        self.radio_time_status = radio_time_status
        self.multimessenger_status = multimessenger_status
        self.window_s = window_s
        self.pointing_status = pointing_status
        self.reasons = reasons
        self.start_utc = start_utc

    @property
    def testable(self) -> bool:
        return self.multimessenger_status in ("TESTABLE", "TESTED")

    def to_dict(self) -> dict:
        return {
            "radio_time_status": self.radio_time_status,
            "multimessenger_status": self.multimessenger_status,
            "window_s": self.window_s,
            "pointing_status": self.pointing_status,
            "reasons": list(self.reasons),
            "start_utc": (None if self.start_utc is None
                          else self.start_utc.isoformat()),
        }


def assess(record, sigma: float = 3.0) -> Assessment:
    """Grade one radio record against the two ingredients a coincidence needs.

    ``sigma`` multiplies the combined timing uncertainty to set the window, and
    the window is reported so a reader can price it rather than trust it.
    """
    start = _start(record)
    duration = float(getattr(record, "duration_s", 0.0) or 0.0)
    ra, dec, localisation = _position(record)
    raw_start = (getattr(record, "metadata", {}) or {}).get("start_time_utc")

    reasons: list[str] = []
    if start is None:
        cost = background_events(REFERENCE_RATE_PER_YEAR, REFERENCE_CONE_DEG,
                                 DAY_S)
        if raw_start:
            reasons.append(
                f"start_time_utc {raw_start!r} is present but unusable: it is "
                "not ISO 8601 with a timezone, and a wall-clock time with no "
                "zone is not a timestamp")
        else:
            reasons.append(
                "no timezone-aware start_time_utc in the record metadata: the "
                "files are named by date, and a date is not a clock")
        reasons.append(
            f"substituting a 24 h window would carry a chance probability of "
            f"{cost:.2e} in a {REFERENCE_CONE_DEG:.2f} deg cone at "
            f"{REFERENCE_RATE_PER_YEAR:.0f} alerts/yr, about "
            f"{cost / background_events(REFERENCE_RATE_PER_YEAR, REFERENCE_CONE_DEG, 256.0):.0f}x "
            "the survey's own 256 s integration")
        return Assessment(radio_time_status="MISSING",
                          multimessenger_status="NOT TESTABLE",
                          window_s=None,
                          pointing_status=("PRESENT" if ra is not None
                                           else "MISSING"),
                          reasons=tuple(reasons), start_utc=None)

    window_s = sigma * 2.0 * duration if duration > 0 else None
    if window_s is None or not math.isfinite(window_s) or window_s <= 0:
        reasons.append("the record has a start but no positive duration, so the "
                       "window has no end")
        return Assessment(radio_time_status="UTC_START",
                          multimessenger_status="NOT TESTABLE",
                          window_s=None,
                          pointing_status=("PRESENT" if ra is not None
                                           else "MISSING"),
                          reasons=tuple(reasons), start_utc=start)

    if ra is None or dec is None:
        reasons.append("a start time and a duration exist, but the record has no "
                       "position: a time match with no cone is a coincidence "
                       "with the whole sky")
        return Assessment(radio_time_status="UTC_START",
                          multimessenger_status="NOT TESTABLE",
                          window_s=window_s, pointing_status="MISSING",
                          reasons=tuple(reasons), start_utc=start)

    if localisation is None:
        reasons.append("the position has no uncertainty, so the cone is set by "
                       "the neutrino's angular error alone")
    reasons.append(f"window {window_s:.0f} s from a {duration:.0f} s "
                   f"integration at {sigma:g} sigma")
    return Assessment(radio_time_status="FULL",
                      multimessenger_status="TESTABLE",
                      window_s=window_s,
                      pointing_status="PRESENT",
                      reasons=tuple(reasons), start_utc=start)


def require_testable(assessment: Assessment, what: str = "record") -> None:
    """Refuse to proceed unless the record can actually be tested."""
    if not assessment.testable:
        detail = "; ".join(assessment.reasons) or "no reason recorded"
        raise MultimessengerNotTestable(
            f"{what} is {assessment.multimessenger_status} "
            f"(radio_time_status={assessment.radio_time_status}): {detail}")


def enforce(assessment: Assessment, coincidence, what: str = "record") -> None:
    """Refuse a reported coincidence for a record that cannot be tested.

    This is the strict half of the gate. It is called on every result, so the
    failure mode it catches - a table row saying ``COINCIDENT`` next to a
    ``MISSING`` clock - stops the run rather than being noticed later.
    """
    status = getattr(coincidence, "status", None)
    if status == "COINCIDENT" and assessment.multimessenger_status == "NOT TESTABLE":
        detail = "; ".join(assessment.reasons) or "no reason recorded"
        raise MultimessengerNotTestable(
            f"a coincidence was reported for a record that is NOT TESTABLE "
            f"(radio_time_status={assessment.radio_time_status}): {detail}")
