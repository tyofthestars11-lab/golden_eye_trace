"""The neutrino lane: a second messenger that can add, never rescue.

A radio transient has one obvious independent check that is not another radio
telescope: a high-energy neutrino from the same place at the same time. This
module states what a neutrino event is, what it takes for it to be coincident
with a radio record, and what the coincidence is allowed to change.

It is allowed to change very little, and the asymmetry is the whole point:

* A **coincidence** adds a counterpart to an already-valid record. It is a
  reason to look again with a pointing telescope, not a disposition. It cannot
  upgrade a candidate, because the radio gates that failed are checks on the
  radio, and a neutrino does not perform them.
* **Absence proves nothing.** A neutrino detector has a sky coverage, an energy
  threshold, an angular resolution of a degree or worse, a duty cycle of order
  unity only in the best case, and a flux that may simply be below threshold for
  a given source. ``Coincidence.absence_is_informative`` is therefore always
  ``False``, and it is a field rather than a comment so that a reader cannot
  build a claim on an empty query without meeting it. Absence of a neutrino is
  not evidence that a radio feature is artificial.
* **A coincidence cannot restore an invalid record.** If the radio pipeline
  determined that a signal was manufactured by its own mask, a neutrino in the
  same part of the sky is a coincidence between an artifact and a real particle,
  which is not evidence about either. ``guard_support`` raises rather than
  returning a flag, so the mistake cannot be made quietly.

Every statistic here is reported beside its null. A coincidence has a chance
probability computed from the exposure - the event rate per unit solid angle per
unit time - so "coincident" without an exposure is an ``UNINFORMATIVE`` result
and not a detection. When no exposure is supplied the probability is ``None``,
and ``None`` is printed as ``no null available``, never as zero.
"""

from __future__ import annotations

import datetime as dt
import math
from dataclasses import asdict, dataclass, field

__all__ = [
    "ALERT_STATUSES",
    "Coincidence",
    "DetectorExposure",
    "ENERGY_CLASSES",
    "EVENT_TYPES",
    "NeutrinoAsEvidence",
    "NeutrinoEvent",
    "angular_separation_deg",
    "chance_probability",
    "guard_support",
    "match",
    "match_catalog",
    "solid_angle_sr",
]

#: Closed vocabularies, for the same reason the disposition list is closed: a
#: value outside them is a typo in a catalog, not a new kind of event.
ENERGY_CLASSES = ("shower", "track", "cascade", "unknown")
EVENT_TYPES = ("through-going", "starting", "contained", "double-bang",
               "cascade", "unknown")
ALERT_STATUSES = ("GOLD", "BRONZE", "archival", "test")

#: A coincidence is never a discovery on its own, so the statuses say what was
#: found and handed back the judgement to the radio record it belongs to.
MATCH_STATUSES = (
    "COINCIDENT",          # inside the time and angle box, with a null
    "NO EVENT IN WINDOW",  # queried, nothing there - uninformative, see below
    "UNINFORMATIVE",       # the query could not be made at all
)


class NeutrinoAsEvidence(Exception):
    """Raised when a neutrino is asked to do more than it can.

    Two things are refused: using a coincidence to raise a disposition, and
    using any neutrino result to restore a record the radio pipeline invalidated.
    """


def _utc(value: dt.datetime, what: str) -> dt.datetime:
    if not isinstance(value, dt.datetime):
        raise ValueError(f"{what} must be a datetime, got {type(value).__name__}")
    if value.tzinfo is None or value.utcoffset() is None:
        # A naive timestamp in a multi-messenger table is an hour of ambiguity
        # disguised as precision. Refusing it is cheaper than guessing the zone.
        raise ValueError(f"{what} must be timezone-aware UTC, not naive")
    return value.astimezone(dt.timezone.utc)


def solid_angle_sr(half_angle_deg: float) -> float:
    """Solid angle of a cone of the given half-angle, in steradians."""
    if half_angle_deg <= 0:
        raise ValueError("a cone needs a positive half-angle")
    theta = math.radians(half_angle_deg)
    return 2.0 * math.pi * (1.0 - math.cos(theta))


def angular_separation_deg(ra1: float, dec1: float,
                           ra2: float, dec2: float) -> float:
    """Great-circle separation in degrees, with RA wrap handled.

    The haversine form is used because the flat-sky difference of two RA values
    is wrong near the poles and near RA = 0, which is exactly where a convenient
    approximation would be used without checking.
    """
    for ra, dec in ((ra1, dec1), (ra2, dec2)):
        if not (-90.0 <= dec <= 90.0):
            raise ValueError(f"declination {dec} is outside [-90, 90]")
    phi1, phi2 = math.radians(ra1), math.radians(ra2)
    dphi = phi2 - phi1
    d1, d2 = math.radians(dec1), math.radians(dec2)
    h = (math.sin((d2 - d1) / 2.0) ** 2
         + math.cos(d1) * math.cos(d2) * math.sin(dphi / 2.0) ** 2)
    return math.degrees(2.0 * math.asin(min(1.0, math.sqrt(h))))


@dataclass(frozen=True)
class NeutrinoEvent:
    """One public neutrino event, as published by its detector.

    ``data_release`` and ``public_source_url`` are required because an event
    without them cannot be checked by anyone else, and an uncheckable event is
    not evidence. A synthetic fixture must say so in ``data_release``; the
    crossmatch refuses to treat it as a real coincidence claim.
    """

    event_id: str
    detector: str
    event_time_utc: dt.datetime
    time_uncertainty_s: float
    ra_deg: float
    dec_deg: float
    angular_error_deg: float
    energy_estimate_gev: float
    energy_class: str
    event_type: str
    alert_status: str
    public_source_url: str
    data_release: str

    def __post_init__(self) -> None:
        for name in ("event_id", "detector", "public_source_url", "data_release"):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"a neutrino event needs {name}")
        _utc(self.event_time_utc, "event_time_utc")
        if self.time_uncertainty_s < 0:
            raise ValueError("time_uncertainty_s cannot be negative")
        if not (0.0 <= self.ra_deg < 360.0):
            raise ValueError(f"ra_deg {self.ra_deg} is outside [0, 360)")
        if not (-90.0 <= self.dec_deg <= 90.0):
            raise ValueError(f"dec_deg {self.dec_deg} is outside [-90, 90]")
        if self.angular_error_deg <= 0:
            # A zero-degree error would make every event trivially coincident
            # with everything, which is the one failure mode a cone cannot have.
            raise ValueError("angular_error_deg must be positive")
        if self.energy_estimate_gev <= 0:
            raise ValueError("energy_estimate_gev must be positive")
        if self.energy_class not in ENERGY_CLASSES:
            raise ValueError(f"energy_class {self.energy_class!r} not one of "
                             f"{ENERGY_CLASSES}")
        if self.event_type not in EVENT_TYPES:
            raise ValueError(f"event_type {self.event_type!r} not one of "
                             f"{EVENT_TYPES}")
        if self.alert_status not in ALERT_STATUSES:
            raise ValueError(f"alert_status {self.alert_status!r} not one of "
                             f"{ALERT_STATUSES}")

    @property
    def synthetic(self) -> bool:
        return (self.alert_status == "test"
                or "synthetic" in self.data_release.lower())

    def to_dict(self) -> dict:
        out = asdict(self)
        out["event_time_utc"] = self.event_time_utc.isoformat()
        return out

    @classmethod
    def from_dict(cls, row: dict) -> "NeutrinoEvent":
        payload = dict(row)
        raw = payload.get("event_time_utc")
        if isinstance(raw, str):
            payload["event_time_utc"] = dt.datetime.fromisoformat(
                raw.replace("Z", "+00:00"))
        return cls(**payload)


@dataclass(frozen=True)
class DetectorExposure:
    """What turns a coincidence into a probability.

    ``density_per_sr_s`` is the event rate per steradian per second in the part
    of the sky the detector could see at the relevant energy. Everything else
    here exists to make that number's assumptions explicit, because a chance
    probability computed from an unstated exposure is a number with no null.
    """

    livetime_s: float
    n_events: int
    sky_fraction: float = 1.0
    energy_threshold_gev: float | None = None
    note: str = ""

    def __post_init__(self) -> None:
        if self.livetime_s <= 0:
            raise ValueError("livetime_s must be positive")
        if self.n_events < 0:
            raise ValueError("n_events cannot be negative")
        if not (0.0 < self.sky_fraction <= 1.0):
            raise ValueError("sky_fraction is a fraction of 4 pi, in (0, 1]")

    @property
    def density_per_sr_s(self) -> float:
        return self.n_events / (self.livetime_s * 4.0 * math.pi * self.sky_fraction)


def chance_probability(half_angle_deg: float, time_window_s: float,
                       exposure: DetectorExposure | None) -> float | None:
    """Probability that one random event lands inside the box, or ``None``.

    Poisson: the expected count is the density times the space-time volume,
    ``mu = lambda * Omega * (2 * dt)``, and ``P(>= 1) = 1 - exp(-mu)``. With no
    exposure the answer is ``None`` and must be reported as such: the absence of
    a null is not a small null.
    """
    if exposure is None:
        return None
    if time_window_s < 0:
        raise ValueError("time_window_s cannot be negative")
    mu = exposure.density_per_sr_s * solid_angle_sr(half_angle_deg) * (2.0 * time_window_s)
    return 1.0 - math.exp(-mu)


@dataclass(frozen=True)
class Coincidence:
    """The outcome of asking the neutrino catalog about one radio record."""

    status: str
    event: NeutrinoEvent | None = None
    separation_deg: float | None = None
    dt_s: float | None = None
    allowed_separation_deg: float | None = None
    allowed_dt_s: float | None = None
    chance_probability: float | None = None
    #: Always ``False``. A field rather than a comment, so that a caller reading
    #: the dataclass meets the rule before forming a claim from an empty query.
    absence_is_informative: bool = False
    #: What limited the query. An empty tuple with a non-COINCIDENT status is a bug.
    limits: tuple[str, ...] = ()
    note: str = ""

    def __post_init__(self) -> None:
        if self.status not in MATCH_STATUSES:
            raise ValueError(f"status {self.status!r} not one of {MATCH_STATUSES}")
        if self.absence_is_informative:
            raise ValueError(
                "absence of a neutrino is never informative about a radio "
                "feature: a detector has a threshold, a coverage and a duty "
                "cycle, and this module will not represent otherwise")
        if self.status == "COINCIDENT":
            if self.event is None:
                raise ValueError("a COINCIDENT result needs the event")
            if self.chance_probability is None:
                raise ValueError(
                    "a COINCIDENT result needs its null: supply a "
                    "DetectorExposure so the chance probability is computed")
        elif self.status == "NO EVENT IN WINDOW" and not self.limits:
            raise ValueError(
                "an empty query must name the limits that make it uninformative")

    @property
    def is_counterpart_candidate(self) -> bool:
        return (self.status == "COINCIDENT"
                and self.chance_probability is not None
                and self.chance_probability < 0.01)

    def to_dict(self) -> dict:
        out = asdict(self)
        if self.event is not None:
            out["event"] = self.event.to_dict()
        out["limits"] = list(self.limits)
        return out


def match(event: NeutrinoEvent, *,
          ra_deg: float | None,
          dec_deg: float | None,
          time_utc: dt.datetime | None,
          radio_time_uncertainty_s: float = 0.0,
          radio_localisation_deg: float | None = None,
          exposure: DetectorExposure | None = None,
          sigma: float = 3.0) -> Coincidence:
    """Test one event against one radio position and time.

    The box is ``sigma`` times the combined uncertainty in each dimension, added
    in quadrature rather than taking the larger of the two: two independent
    one-degree errors make a 1.41 degree error, and using the larger would
    silently buy sensitivity that the data does not have.
    """
    limits: list[str] = []
    if ra_deg is None or dec_deg is None:
        return Coincidence(status="UNINFORMATIVE",
                           limits=("the radio record has no pointing: a "
                                   "single-beam product with no WCS cannot be "
                                   "crossmatched in angle",))
    if time_utc is None:
        return Coincidence(status="UNINFORMATIVE",
                           limits=("the radio record has no UTC timestamp, so "
                                   "no time window can be drawn",))
    stamp = _utc(time_utc, "time_utc")
    if radio_time_uncertainty_s < 0:
        raise ValueError("radio_time_uncertainty_s cannot be negative")

    sep = angular_separation_deg(ra_deg, dec_deg, event.ra_deg, event.dec_deg)
    dt_s = abs((event.event_time_utc - stamp).total_seconds())

    if radio_localisation_deg is None:
        angle_sigma = event.angular_error_deg
        limits.append("the radio position has no uncertainty of its own; only "
                      "the neutrino angular error is used")
    else:
        angle_sigma = math.hypot(event.angular_error_deg, radio_localisation_deg)
    time_sigma = math.hypot(event.time_uncertainty_s, radio_time_uncertainty_s)
    allowed_angle = sigma * angle_sigma
    allowed_dt = sigma * time_sigma

    allowed_angle = min(allowed_angle, 180.0)
    inside = sep <= allowed_angle and dt_s <= allowed_dt
    if not inside:
        return Coincidence(
            status="NO EVENT IN WINDOW",
            separation_deg=sep, dt_s=dt_s,
            allowed_separation_deg=allowed_angle, allowed_dt_s=allowed_dt,
            limits=tuple(limits + [
                "an event outside the box is not evidence about the record; the "
                "detector's threshold and coverage still apply"]),
            note=(f"{event.event_id} is {sep:.2f} deg and {dt_s:.0f} s away, "
                  f"outside {allowed_angle:.2f} deg and {allowed_dt:.0f} s"))

    p = chance_probability(allowed_angle, allowed_dt, exposure)
    if p is None:
        # Coincident, but with no exposure there is no null, and without a null
        # "coincident" is a statement about a catalog row, not about the sky.
        return Coincidence(
            status="UNINFORMATIVE",
            event=event, separation_deg=sep, dt_s=dt_s,
            allowed_separation_deg=allowed_angle, allowed_dt_s=allowed_dt,
            limits=tuple(limits + [
                "no detector exposure was supplied, so no chance probability "
                "could be computed"]),
            note=f"{event.event_id} is inside the box but has no null")
    return Coincidence(status="COINCIDENT", event=event,
                       separation_deg=sep, dt_s=dt_s,
                       allowed_separation_deg=allowed_angle,
                       allowed_dt_s=allowed_dt,
                       chance_probability=p,
                       limits=tuple(limits))


def match_catalog(events, **kwargs) -> Coincidence:
    """Best (highest chance-probability-of-being-real) result over a catalog.

    The best event is the one with the smallest chance probability, not the
    closest one: the closest event in a crowded catalog is often a coincidence,
    and picking it would be a selection on the null.
    """
    best: Coincidence | None = None
    for event in events:
        result = match(event, **kwargs)
        if best is None:
            best = result
            continue
        if result.status == "COINCIDENT" and (
                best.status != "COINCIDENT"
                or (result.chance_probability or 1.0)
                < (best.chance_probability or 1.0)):
            best = result
        elif (best.status != "COINCIDENT"
              and result.status == "UNINFORMATIVE" and best.status == "NO EVENT IN WINDOW"):
            best = result
    if best is None:
        return Coincidence(status="NO EVENT IN WINDOW",
                           limits=("the catalog is empty",))
    return best


def guard_support(record, coincidence: Coincidence, *,
                  claimed: str | None = None) -> dict:
    """Refuse the two ways a neutrino could be used dishonestly.

    ``claimed`` is the disposition the caller wants the coincidence to support.
    Passing it is how the caller states a promotion, and a promotion is refused:
    a coincidence may be attached to a record's metadata, never to its rank.

    A record whose disposition is in ``INVALID_CANDIDATE`` refuses every
    coincidence, because the radio signal was manufactured by the processing and
    a real particle nearby does not make it real.
    """
    from .signal_record import INVALID_CANDIDATE  # local: keeps the cycle out

    if record.disposition in INVALID_CANDIDATE:
        raise NeutrinoAsEvidence(
            f"{record.observation_date} {record.frequency_hz / 1e6:.4f} MHz is "
            f"{record.disposition}: it is excluded from discovery claims and a "
            "neutrino coincidence cannot restore it, because the radio event it "
            "would be coincident with does not exist in the data")
    if claimed is not None and claimed != record.disposition:
        raise NeutrinoAsEvidence(
            f"a neutrino coincidence cannot move a record from "
            f"{record.disposition} to {claimed}: the neighbourhood of a real "
            "particle is not a measurement of the radio source")
    if coincidence.status != "COINCIDENT":
        raise NeutrinoAsEvidence(
            f"no coincidence to attach: status is {coincidence.status}, and "
            "absence of a neutrino is not evidence about the radio record")
    return {"neutrino_event_id": coincidence.event.event_id,
            "neutrino_chance_probability": coincidence.chance_probability,
            "neutrino_separation_deg": coincidence.separation_deg,
            "neutrino_dt_s": coincidence.dt_s,
            "neutrino_data_release": coincidence.event.data_release,
            "neutrino_synthetic": coincidence.event.synthetic}
