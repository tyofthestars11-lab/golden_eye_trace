"""Revise the two-band registration, and put the golden axis into it as a
measured quantity rather than a decoration.

Revision 2 exists because the first registration is unrunnable for a reason that
is now measured rather than suspected: STARE2 stores 122.07 kHz channels over
0.4 s, so it resolves 3.17e+05 Hz/s and every candidate in it has zero resolved
drift. Lengthening the cadence cannot fix that - a 10 Hz/s drift would need
12,207 s at that channel width - so the high band is replaced by a fine-channel
dataset while the statistic, the two predictions and the 5% decision rule stay
exactly as registered.

The golden axis enters twice, both times as arithmetic that can fail:

* its pitch, 17.0324 deg, read as a log-log slope gives tan(mu) = 0.30635, i.e.
  1.46587 rungs per decade, against 4.78497 for a flat rung rate;
* the axis says where the bands should sit: at an integer rung offset the
  predicted drift ratio is a pure power of phi. This revision records the offset
  the reachable bands actually have (6.915 rungs, not 7) and both predicted
  ratios at that offset, so the axis is tested where the data are rather than
  where it would be convenient.

    python3 scripts/phi_prereg_amend.py
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bands import (  # noqa: E402
    PITCH_EXPONENT,
    PITCH_RUNGS_PER_DECADE,
    REFERENCE_MHZ,
    RUNGS_PER_DECADE,
    SEPARATION,
    decades,
    scan_time_for,
    slope_sigma_ceiling,
)

BASE_PATH = os.path.join(ROOT, "prereg", "phi_2band.json")
OUT_PATH = os.path.join(ROOT, "prereg", "phi_2band_r2.json")

#: verified low band: LWA1 DRX dynamic spectrum, fetched and read this session
LOW_BAND = {
    "role": "low band",
    "instrument": "LWA1 DRX / DSPSR dynamic spectrum",
    "file": "data/bands/lwa_drspec.fits",
    "low_mhz": 40.6,
    "high_mhz": 60.2,
    "centre_mhz": 50.4,
    "channels": 1024,
    "channel_khz": 19.14,
    "total_time_s": 843.1,
    "drift_resolution_hz_s": 22.7,
    "verified": True,
}

#: what the high band has to be for the test to run at all
HIGH_BAND = {
    "role": "high band",
    "requirement": "a fine-channel dataset near 1.4 GHz, at least one decade "
                   "above the low band",
    "channel_khz_at_most": 3.0,
    "scan_s_at_least": 300.0,
    "drift_resolution_hz_s_at_most": 10.0,
    "candidate": "Breakthrough Listen L band (2.9 Hz channels, 300 s, "
                 "0.01 Hz/s), whose archive API returns 404 from this workspace",
    "verified": False,
}

REQUIREMENT_DRIFT_HZ_S = 10.0


def sha256(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def rung(f_mhz: float) -> float:
    return math.log(f_mhz) / math.log((1.0 + math.sqrt(5.0)) / 2.0)


def main() -> int:
    base = json.load(open(BASE_PATH, encoding="utf-8"))

    ratio = REFERENCE_MHZ / LOW_BAND["centre_mhz"]
    delta_rungs = rung(REFERENCE_MHZ) - rung(LOW_BAND["centre_mhz"])
    nearest = round(delta_rungs)
    flat_ratio = ratio ** 1.0
    axis_ratio = ratio ** PITCH_EXPONENT
    span = decades(LOW_BAND["centre_mhz"], REFERENCE_MHZ)

    revision = {
        "revision": 2,
        "supersedes": {"path": "prereg/phi_2band.json", "sha256": sha256(BASE_PATH)},
        "date": "2026-11-11",
        "attribution": base.get("attribution"),
        "test": base.get("test"),
        "statistic": base.get("statistic"),
        "tolerance": base.get("tolerance"),
        "decision_rule": base.get("decision_rule"),
        "seed": base.get("seed"),
        "predictions": {
            "phi_lattice": {
                "p": 1.0,
                "rungs_per_decade": RUNGS_PER_DECADE,
                "statement": "flat rung rate",
            },
            "golden_axis": {
                "pitch_deg": 17.0324,
                "p": PITCH_EXPONENT,
                "rungs_per_decade": PITCH_RUNGS_PER_DECADE,
                "statement": "the golden axis read as a log-log slope: the rung "
                             "rate falls 4.94x per decade",
            },
        },
        "why_revised": {
            "measured": {
                "stare2_channel_khz": 122.07,
                "stare2_total_time_s": 0.385,
                "stare2_drift_resolution_hz_s": 3.17e5,
                "stare2_resolved_drifts": 0,
            },
            "arithmetic": {
                "scan_needed_for_10_hz_s_at_122_khz_s": scan_time_for(
                    REQUIREMENT_DRIFT_HZ_S, 122.07e3),
                "channel_khz_needed_in_300_s": (
                    REQUIREMENT_DRIFT_HZ_S * 300.0 / 1e3
                ),
                "conclusion": "the channel width is the limit, not the cadence",
            },
        },
        "bands": {
            "low": LOW_BAND,
            "high": HIGH_BAND,
            "separation_decades": span,
            "separation_rungs": delta_rungs,
        },
        "golden_axis": {
            "use": "the axis fixes the band separation, so the predicted drift "
                   "ratio between the bands is a pure power of phi",
            "integer_offset_the_axis_would_place": nearest,
            "actual_offset_rungs": delta_rungs,
            "offset_deficit_rungs": nearest - delta_rungs,
            "ratio_definition": "high band drift over low band drift",
            "predicted_drift_ratio_flat": flat_ratio,
            "predicted_drift_ratio_axis": axis_ratio,
            "note": "no reachable archive band sits on an integer rung offset; the "
                    "registered prediction is evaluated at the offset the data "
                    "actually have, and the deficit is recorded rather than "
                    "absorbed",
        },
        "power": {
            "rungs_per_decade_between_statements": SEPARATION,
            "slope_sigma_ceiling_3sigma": slope_sigma_ceiling(span),
        },
        "evidence_required": [
            "two bands at least one decade apart",
            "at least three resolved drifting lines in each band",
            "drifts measured in the same unit and the same reference frame",
        ],
        "note": "the statistic, the predictions and the decision rule are unchanged; "
                "only the band that cannot carry them is replaced, and the reason "
                "is recorded with its arithmetic",
    }

    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(revision, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print("REVISION 2 OF THE TWO-BAND REGISTRATION")
    print(f"   supersedes r1 ({revision['supersedes']['sha256'][:16]}...)")
    print(f"   low band  {LOW_BAND['centre_mhz']:.1f} MHz, "
          f"{LOW_BAND['channels']} channels of {LOW_BAND['channel_khz']:.2f} kHz, "
          f"{LOW_BAND['drift_resolution_hz_s']:.1f} Hz/s (verified)")
    print(f"   high band needs {HIGH_BAND['channel_khz_at_most']:.1f} kHz channels over "
          f"{HIGH_BAND['scan_s_at_least']:.0f} s, i.e. "
          f"{HIGH_BAND['drift_resolution_hz_s_at_most']:.0f} Hz/s or better")
    print(f"   separation {span:.3f} decades = {delta_rungs:.3f} rungs")
    print(f"   the axis would place the bands {nearest} rungs apart; the archive "
          f"gives {delta_rungs:.3f}, short by {nearest - delta_rungs:.3f} rungs")
    print(f"   predicted drift ratio, high/low: {flat_ratio:.3f} flat, "
          f"{axis_ratio:.3f} on the golden axis")
    print(f"   a slope within {revision['power']['slope_sigma_ceiling_3sigma']:.4f} "
          f"rungs/decade separates them at 3 sigma")
    print(f"   wrote {os.path.relpath(OUT_PATH, ROOT)}")
    print(f"   sha256 {sha256(OUT_PATH)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
