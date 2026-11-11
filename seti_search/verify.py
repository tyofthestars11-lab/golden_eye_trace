"""ON/OFF cadence verification.

A single file can only produce candidates: the strongest thing in it may be a
satellite, a microwave oven or a loose cable. The test that separates a sky
signal from local interference is a cadence - point at the target (ON), point
away (OFF), and require the signal to appear in ON and vanish in OFF while its
drift stays consistent with the sky position. Interference is not switched off
by where the telescope points; a source is.

What this module does, for each ON candidate:

1. carry the observation forward to the OFF epoch with
   :func:`seti_search.cadence.predict_observation`, which accounts for the
   Doppler difference between the two epochs and for the source's own drift
   accumulating across the gap;
2. look for an OFF candidate at that predicted frequency;
3. if one is there, reject the candidate - it is present when the telescope is
   pointing somewhere else, which no sky source can be. Report whether the OFF
   hit's drift also matched the prediction, which distinguishes the same
   transmitter seen twice from two unrelated things in the same channel;
4. if none is there, confirm the candidate as *cadence-consistent*: it behaved
   like a source. That is still not a detection, it is the absence of the
   objection that kills nearly every candidate.

The check is only as strong as the cadence. If the ON and OFF scans are close
together in time, the predicted drift barely moves and step 3 cannot tell a
celestial drift from a stationary interferer; each verdict carries
``drift_discrimination_hz_s`` so that weakness is visible rather than hidden.
Likewise a much shallower OFF scan will fail to show real signals, so compare
sensitivities before believing a confirmation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .cadence import Observer, SkyPosition, predict_observation
from .search import Candidate

__all__ = ["Verdict", "verify_cadence"]

CONFIRMED = "confirmed"
REJECTED = "rejected"


@dataclass
class Verdict:
    """The outcome for one ON candidate."""

    on: Candidate
    status: str
    reason: str
    off_hit: Candidate | None = None
    predicted_freq_mhz: float | None = None
    predicted_drift_hz_s: float | None = None
    freq_window_mhz: float = 0.0
    drift_discrimination_hz_s: float | None = None

    @property
    def confirmed(self) -> bool:
        return self.status == CONFIRMED


@dataclass
class CadenceResult:
    verdicts: list[Verdict] = field(default_factory=list)
    predicted: bool = False
    notes: list[str] = field(default_factory=list)

    @property
    def confirmed(self) -> list[Verdict]:
        return [v for v in self.verdicts if v.confirmed]

    @property
    def rejected(self) -> list[Verdict]:
        return [v for v in self.verdicts if not v.confirmed]


def verify_cadence(
    on_hits: list[Candidate],
    off_hits: list[Candidate],
    *,
    on_meta: dict,
    off_meta: dict,
    observer: Observer | None = None,
    sky: SkyPosition | None = None,
    on_mjd: float | None = None,
    off_mjd: float | None = None,
    drift_tolerance_hz_s: float = 0.05,
    freq_tolerance_hz: float | None = None,
    drift_uncertainty_hz_s: float = 0.05,
) -> CadenceResult:
    """Decide which ON candidates survive the OFF scan.

    Args:
        on_hits, off_hits: candidates from each file, from ``search_file``.
        on_meta, off_meta: the matching metadata dicts (used for channel width).
        observer, sky, on_mjd, off_mjd: everything needed to predict the drift.
            If any is missing the drift is not predicted and only frequency
            coincidence is tested; ``CadenceResult.notes`` says so.
        drift_tolerance_hz_s: how far an OFF hit's drift may sit from the
            prediction and still count as the same transmitter.
        freq_tolerance_hz: half-width of the frequency window. Default grows
            with the gap between scans, since the source's own drift (known only
            to within ``drift_uncertainty_hz_s``) accumulates across it.
    """
    result = CadenceResult()
    df_off_hz = abs(float(off_meta["df_mhz"])) * 1e6
    can_predict = (
        observer is not None
        and sky is not None
        and on_mjd is not None
        and off_mjd is not None
    )
    result.predicted = bool(can_predict)

    if not can_predict:
        missing = [
            name
            for name, value in (
                ("observer (site latitude/longitude)", observer),
                ("sky position (RA/Dec)", sky),
                ("ON start time", on_mjd),
                ("OFF start time", off_mjd),
            )
            if value is None
        ]
        result.notes.append(
            "drift prediction unavailable (missing " + ", ".join(missing) + "); "
            "comparing frequencies only"
        )
        gap_s = 0.0
    else:
        gap_s = abs(float(off_mjd) - float(on_mjd)) * 86400.0

    default_window_hz = max(
        4.0 * df_off_hz,
        drift_uncertainty_hz_s * gap_s + 3.0 * df_off_hz,
    )
    window_hz = float(freq_tolerance_hz) if freq_tolerance_hz else default_window_hz
    window_mhz = window_hz / 1e6
    # a stationary interferer does not move between scans, so its window is a
    # couple of channels wide however long the gap is
    static_window_hz = (
        float(freq_tolerance_hz) if freq_tolerance_hz else max(3.0 * df_off_hz, 3.0)
    )
    static_window_mhz = static_window_hz / 1e6

    for hit in sorted(on_hits, key=lambda h: -h.snr):
        if can_predict:
            predicted_hz, predicted_drift = predict_observation(
                hit.freq_mhz * 1e6,
                hit.drift_hz_s,
                observer,
                on_mjd,
                off_mjd,
                sky,
            )
            discrimination = abs(predicted_drift - hit.drift_hz_s)
        else:
            predicted_hz, predicted_drift, discrimination = hit.freq_mhz * 1e6, None, None

        predicted_mhz = predicted_hz / 1e6
        # Two things could sit in the OFF scan, and they get different windows.
        # The same sky source, where the prediction puts it: that window has to
        # absorb the prediction's own uncertainty, which grows with the gap. Or
        # terrestrial interference, which does not move with the sky and stays
        # where the ON candidate was, within a couple of channels. Either one
        # being present means the candidate was not switched off by pointing
        # away - which is what catches a carrier whose drift, taken at face
        # value, would imply a barycentric drift of its own.
        anchors = [(hit.freq_mhz, static_window_mhz)]
        if can_predict:
            anchors.append((predicted_mhz, window_mhz))
        matches = [
            other
            for other in off_hits
            if any(abs(other.freq_mhz - anchor) <= width for anchor, width in anchors)
        ]

        if not matches:
            reason = (
                "absent from the OFF scan at its predicted frequency"
                if can_predict
                else "absent from the OFF scan at the same frequency"
            )
            result.verdicts.append(
                Verdict(
                    on=hit,
                    status=CONFIRMED,
                    reason=reason,
                    predicted_freq_mhz=predicted_mhz,
                    predicted_drift_hz_s=predicted_drift,
                    freq_window_mhz=window_mhz,
                    drift_discrimination_hz_s=discrimination,
                )
            )
            continue

        best = min(
            matches,
            key=lambda other: min(abs(other.freq_mhz - a) for a, _ in anchors),
        )
        if predicted_drift is not None:
            drift_ok = abs(best.drift_hz_s - predicted_drift) <= drift_tolerance_hz_s
            if drift_ok:
                reason = (
                    "present in the OFF scan at the same frequency and drift: "
                    "not switched off by pointing away, so not a sky signal"
                )
            else:
                reason = (
                    "present in the OFF scan at the same frequency but a different "
                    "drift: either the same transmitter seen twice with an unstable "
                    "drift, or unrelated interference in that channel"
                )
        else:
            reason = (
                "present in the OFF scan at the same frequency (no drift check: "
                "time or coordinates unavailable)"
            )
        result.verdicts.append(
            Verdict(
                on=hit,
                status=REJECTED,
                reason=reason,
                off_hit=best,
                predicted_freq_mhz=predicted_mhz,
                predicted_drift_hz_s=predicted_drift,
                freq_window_mhz=window_mhz,
                drift_discrimination_hz_s=discrimination,
            )
        )

    if can_predict and result.verdicts:
        weakest = min(
            (v.drift_discrimination_hz_s or 0.0) for v in result.verdicts
        )
        if weakest < drift_tolerance_hz_s:
            result.notes.append(
                f"the scans predict drifts differing by only {weakest:.3f} Hz/s, "
                f"less than the {drift_tolerance_hz_s:g} Hz/s tolerance: a longer "
                "cadence would make the drift test decisive"
            )
    return result
