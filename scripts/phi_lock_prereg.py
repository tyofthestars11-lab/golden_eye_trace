"""Pre-registration for the decisive test of the phi phase lock on drift lines.

Written before any confirmatory data are examined. `GEOMETRY.md` lists the five
conditions a lock test needs to be decisive; this file fixes them in advance and
hashes itself, so the decision rule cannot be chosen after seeing which band
happens to contain something.

The conditions, and how each is fixed here:

1. `tsamp` comes from the file header, i.e. the instrument. No dedoppler-summed
   blocks, no averaging chosen because it would put a drift on the lock.
2. The drift window runs to Nyquist. The lock sits at 76.4% of it, so a window
   capped at a few Hz/s cannot reach the lock at all.
3. The integer advance is bounded at `|k| = 1`. `k = 2` is outside Nyquist, so
   one physical lock exists; allowing aliased `k` would let `k <= 16` cover every
   rate, which is the failure mode recorded in `GEOMETRY.md` section 4.
4. Sample count: at least 987 samples for eleven Fibonacci lags is required
   before the echo contrast may be quoted. The rate test needs only the drift and
   can run on fewer samples; the two are reported separately and never combined.
5. The chance rate is stated and applied as a binomial, not per candidate alone.

The exploratory search already run on LWA `drspec` (2026-11-11, +/-7454 Hz/s, 11
candidates, 0 resolved drifts) is recorded as **exploratory** and excluded from
the confirmatory count: it was run before this pre-registration existed, so it
cannot count as a prediction that came true.

    python3 scripts/phi_lock_prereg.py [--fil data/lwa_drspec.fil]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.filterbank import Filterbank  # noqa: E402
from seti_search.rungs import PHI  # noqa: E402

THETA = 1.0 / (PHI * PHI)
TOLERANCE = 0.05          # of theta, on the lock distance
K_MAX = 1                 # k = 2 lies outside Nyquist
MIN_SAMPLES_ECHO = 987    # eleven Fibonacci lags
EXPLORATORY = {
    "date": "2026-11-11",
    "file": "data/lwa_drspec.fil",
    "window_hz_s": 7453.7,
    "candidates": 11,
    "resolved_drifts": 0,
    "status": "exploratory, excluded from the confirmatory count",
}


def chance_per_candidate() -> float:
    """Probability a drift drawn uniformly in Nyquist lands within tolerance*theta."""
    return 2.0 * TOLERANCE * THETA / 0.5


def build(channel_hz: float, tsamp_s: float, nchans: int, nints: int,
          source: str) -> dict:
    nyquist = 0.5 * channel_hz / tsamp_s
    theta_drift = THETA * channel_hz / tsamp_s
    resolution = channel_hz / (nints * tsamp_s)
    p0 = chance_per_candidate()
    return {
        "attribution": "framework: Tyree Jones (tyofthestarz), golden_eye_trace",
        "revision": 1,
        "theta": THETA,
        "tolerance_of_theta": TOLERANCE,
        "band": {
            "source": source, "channel_hz": channel_hz, "tsamp_s": tsamp_s,
            "nchans": nchans, "nints": nints,
            "nyquist_hz_s": nyquist, "drift_resolution_hz_s": resolution,
            "theta_drift_hz_s": theta_drift,
            "theta_fraction_of_nyquist": theta_drift / nyquist,
        },
        "conditions": [
            "tsamp is taken from the instrument header and never re-binned",
            "the drift window is searched to Nyquist, not capped below the lock",
            "the integer advance is bounded at |k| = 1, k = 2 being outside Nyquist",
            "the echo contrast is quoted only when nints >= 987 (eleven lags)",
            "the chance rate is applied as a binomial over candidates",
        ],
        "statistic": {
            "rate": "nu = drift * tsamp / channel, cycles per sample",
            "lock_distance": "|nu/theta - round(nu/theta)|, evaluated only for resolved drifts",
            "resolved_drift_requires": "|drift| > drift_resolution_hz_s",
            "echo": "mean |autocorrelation| at Fibonacci lags minus mean over all "
                    "non-Fibonacci lags, as in scripts/phi_echo.py",
            "k_bound": K_MAX,
        },
        "predictions": {
            "framework_locked_drift_hz_s": theta_drift,
            "framework_window_hz_s": [theta_drift * (1 - TOLERANCE),
                                      theta_drift * (1 + TOLERANCE)],
            "chance_per_candidate": p0,
            "binomial_null": "P(locked >= m | N evaluable) with p = chance_per_candidate",
        },
        "decision_rule": (
            "A candidate is locked when its drift is resolved, |k| = 1, and the "
            "lock distance is within 0.05. Report the count against the binomial "
            "null. If no candidate has a resolved drift the verdict is NOT "
            "MEASURED: the rate is not quoted, a null is not claimed, and a "
            "failure to resolve is never reported as a non-detection."
        ),
        "minimum_samples_for_echo": MIN_SAMPLES_ECHO,
        "echo_available_here": nints >= MIN_SAMPLES_ECHO,
        "multiplicity": {
            "look_elsewhere": "1 - (1 - p)^N over evaluable candidates"},
        "exploratory_excluded": EXPLORATORY,
        "seed": 1729,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fil", default=None,
                    help="band to pre-register; defaults to the LWA dynamic spectrum")
    args = ap.parse_args()

    if args.fil:
        fb = Filterbank.open(args.fil)
        channel_hz = abs(fb.foff) * 1e6
        payload = build(channel_hz, fb.tsamp, fb.nchans, fb.nints,
                        os.path.basename(args.fil))
    else:
        payload = build(19_140.6, 1.284, 1024, 709, "lwa_drspec (nominal)")

    canonical = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    payload["sha256"] = digest
    out = os.path.join(ROOT, "prereg", "phi_lock.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    with open(os.path.join(ROOT, "prereg", "phi_lock.sha256"), "w",
              encoding="utf-8") as fh:
        fh.write(digest + "\n")

    band = payload["band"]
    print("PRE-REGISTRATION: the phi phase lock on drift lines")
    print(f"   band {band['source']}: {band['nchans']} channels, "
          f"{band['nints']} samples, tsamp {band['tsamp_s']} s")
    print(f"   Nyquist drift {band['nyquist_hz_s']:,.0f} Hz/s, "
          f"resolution {band['drift_resolution_hz_s']:.1f} Hz/s")
    print(f"   locked drift {band['theta_drift_hz_s']:,.1f} Hz/s = "
          f"{100 * band['theta_fraction_of_nyquist']:.1f}% of Nyquist")
    print(f"   framework window {payload['predictions']['framework_window_hz_s'][0]:,.0f}"
          f" to {payload['predictions']['framework_window_hz_s'][1]:,.0f} Hz/s")
    print(f"   chance per candidate {payload['predictions']['chance_per_candidate']:.4f}")
    print(f"   echo available here: {payload['echo_available_here']} "
          f"(needs {MIN_SAMPLES_ECHO} samples, has {band['nints']})")
    print(f"   sha256 {digest}")
    print(f"   wrote prereg/phi_lock.json and prereg/phi_lock.sha256")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
