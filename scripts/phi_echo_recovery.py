"""Recover injected phi-locked chirps, and hunt the LWA band at full drift range.

Two questions the echo harness left open. Does the lock actually recover a signal
that is phi-locked, through the search rather than only in a synthetic track? And
does the LWA band contain any drifting carrier at all once the drift window is
opened to Nyquist instead of the 8 Hz/s the earlier runs used?

Part 1 injects chirps into a simulated file at the rate theta and at arbitrary
rates, runs the real search, and checks that the recovered rate lands on the
lattice only for theta. Part 2 checks the echo contrast on a track at the length
real data has, where ten Fibonacci lags make it a test rather than a coin.

The rate deserves a warning, because it is a physical statement and not a
formality. A phi lock means theta = 0.382 channels of drift per integration, which
is 76% of Nyquist. For the LWA dynamic spectrum (19.14 kHz channels, 1.284 s
integrations) that is 5,695 Hz/s: a chirp that crosses 270 channels in 709
integrations. A beacon drifting at a realistic 1 Hz/s has a rate of 6.7e-5
channels per sample, five thousand times below the lock. So the lock tests fast
chirps, not slow beacons, and no integration choice available here changes that:
theta is a fraction of a channel per sample by construction.

    python3 scripts/phi_echo_recovery.py
"""

from __future__ import annotations

import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

from phi_echo import (  # noqa: E402
    CALIBRATION_N,
    LOCK_TOLERANCE,
    THETA,
    autocorr,
    drift_lock,
    echo_contrast,
    lag_sets,
    null_contrast,
    p_value,
)
from seti_search.filterbank import Filterbank, InjectedSignal, simulate  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

SEED = 1729
CHANNEL_MHZ = 0.0191406   # 19.14 kHz, the LWA dynamic-spectrum channel
TSAMP_S = 0.05


def main() -> int:
    rng = np.random.default_rng(SEED)
    channel_hz = CHANNEL_MHZ * 1e6
    print("RECOVERY OF PHI-LOCKED CHIRPS, AND THE LWA DRIFT WINDOW")
    print(f"   theta = {THETA:.6f} cycles per sample; "
          f"one sample = {TSAMP_S} s, channel = {channel_hz / 1e3:.4f} kHz")
    print(f"   theta is {THETA * channel_hz / TSAMP_S:,.0f} Hz/s here; "
          f"Nyquist drift is {0.5 * channel_hz / TSAMP_S:,.0f} Hz/s")
    print()

    # --- part 1: inject, search, recover the rate ------------------------
    print("   A. INJECTION THROUGH THE SEARCH")
    # +theta and -theta both stay inside Nyquist; 1 - theta = 0.618 aliases past
    # it and the dedoppler would clamp the injection, so the conjugate is
    # injected negatively instead
    rates = {"theta": THETA, "minus theta": -THETA, "arbitrary 0.31": 0.31}
    results: dict = {"theta_hz_s": THETA, "part1": {}, "part2": {}}
    print(f"   {'injected':<16}{'drift Hz/s':>12}{'found MHz':>12}"
          f"{'found Hz/s':>12}{'rate':>8}{'dist':>8}  verdict")
    for name, nu in rates.items():
        drift = nu * channel_hz / TSAMP_S
        path = os.path.join(ROOT, "data", f"inject_{name.replace(' ', '_').replace('*', '')}.fil")
        # start where a chirp of this sign stays inside the 50.1-55.0 MHz band
        start_mhz = 51.0 if drift > 0 else 53.0
        simulate(path, nchans=256, nints=200, tsamp=TSAMP_S,
                 fch1=55.0, foff=-CHANNEL_MHZ,
                 signals=(InjectedSignal(freq_mhz=start_mhz, drift_hz_s=drift,
                                         amplitude=3.0, width_channels=1.5),),
                 glitch=False, seed=SEED)
        cfg = SearchConfig(boxcar=3, max_drift_hz_s=0.5 * channel_hz / TSAMP_S,
                           snr_threshold=4.0, max_hits=10)
        hits, _meta = search_file(path, cfg, chan_chunk=2048, max_ints=None)
        if not hits:
            print(f"   {name:<16}{drift:>12.0f}{'-':>12}{'-':>12}{'-':>8}{'-':>8}"
                  f"  NOT RECOVERED")
            results["part1"][name] = {"recovered": False}
            continue
        best = min(hits, key=lambda h: abs(h.drift_hz_s - drift))
        got_nu, distance = drift_lock(best.drift_hz_s, channel_hz, TSAMP_S)
        locked = distance <= LOCK_TOLERANCE
        print(f"   {name:<16}{drift:>12.0f}{best.freq_mhz:>12.4f}"
              f"{best.drift_hz_s:>12.0f}{got_nu:>8.3f}{distance:>8.3f}  "
              f"{'locked' if locked else 'off lattice'}")
        results["part1"][name] = {
            "recovered": True, "injected_hz_s": drift,
            "found_hz_s": best.drift_hz_s, "rate": got_nu,
            "distance": distance, "locked": bool(locked),
        }
    print()

    # --- part 2: echo contrast at the length real tracks have -------------
    print(f"   B. ECHO CONTRAST AT N={CALIBRATION_N} (ten Fibonacci lags)")
    n = CALIBRATION_N
    fibs, others = lag_sets(n)
    null = null_contrast(n, draws=500, rng=rng)
    t = np.arange(n)
    echo_rows = {}
    print(f"   {'track':<20}{'contrast':>10}{'p':>8}  verdict")
    for name, nu in (("theta", THETA), ("phi*theta", 1.0 - THETA),
                     ("arbitrary 0.31", 0.31)):
        x = np.cos(2 * math.pi * nu * t) + 0.5 * rng.normal(size=n)
        contrast = echo_contrast(x, fibs, others)
        p = p_value(contrast, null)
        echo_rows[name] = {"contrast": contrast, "p": p}
        phi_case = name in ("theta", "phi*theta")
        ok = (p <= 0.02) if phi_case else (p >= 0.05)
        print(f"   {name:<20}{contrast:>+10.4f}{p:>8.4f}  "
              f"{'phi detected' if ok and phi_case else ('goes through' if ok else 'MISSED')}")
    results["part2"] = echo_rows
    print(f"   null over 500 random rates: {null.mean():+.3f} +/- {null.std():.3f}")
    print()

    # --- part 3: the LWA band at the full drift range --------------------
    path = os.path.join(ROOT, "data", "lwa_drspec.fil")
    if os.path.exists(path):
        fb = Filterbank.open(path)
        channel_hz = abs(fb.foff) * 1e6
        nyquist = 0.5 * channel_hz / fb.tsamp
        print("   C. THE LWA WATERFALL, DRIFT WINDOW OPENED TO NYQUIST")
        print(f"   {fb.nchans} channels, {fb.nints} rows, tsamp {fb.tsamp:.4f} s, "
              f"resolution {channel_hz / fb.total_time:.1f} Hz/s")
        print(f"   search window +/-{nyquist:,.0f} Hz/s, theta = "
              f"{THETA * channel_hz / fb.tsamp:,.0f} Hz/s")
        cfg = SearchConfig(boxcar=3, max_drift_hz_s=nyquist, snr_threshold=4.0,
                           max_hits=100)
        hits, _meta = search_file(path, cfg, chan_chunk=2048, max_ints=None)
        drifts = [h for h in hits if abs(h.drift_hz_s) > channel_hz / fb.total_time]
        print(f"   {len(hits)} candidates at 4 sigma, {len(drifts)} with resolved drift")
        rows = []
        for hit in sorted(hits, key=lambda h: -abs(h.drift_hz_s))[:8]:
            nu, distance = drift_lock(hit.drift_hz_s, channel_hz, fb.tsamp)
            rows.append({"freq_mhz": hit.freq_mhz, "drift_hz_s": hit.drift_hz_s,
                         "snr": hit.snr, "rate": nu, "distance": distance})
            print(f"      {hit.freq_mhz:9.4f} MHz  drift {hit.drift_hz_s:>12.1f} "
                  f"Hz/s  snr {hit.snr:5.2f}  rate {nu:.4f}  dist {distance:.3f}")
        results["part3"] = {"candidates": len(hits), "resolved": len(drifts),
                            "nyquist_hz_s": nyquist, "rows": rows}
        chance = 2.0 * LOCK_TOLERANCE * THETA / 0.5
        print(f"   {len(drifts)} resolved drifts; chance a random rate is within "
              f"{LOCK_TOLERANCE} theta of an integer: {chance:.4f}")
        results["part3"]["chance_rate"] = chance

    out = os.path.join(ROOT, "prereg", "phi_echo_recovery.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"\n   wrote {os.path.relpath(out, ROOT)}")

    recovered = all(v.get("recovered") for v in results["part1"].values())
    theta_locked = results["part1"].get("theta", {}).get("locked", False)
    arbitrary_locked = results["part1"].get("arbitrary 0.31", {}).get("locked", False)
    if not (recovered and theta_locked and not arbitrary_locked):
        print("   RECOVERY FAILED: theta must be recovered and locked, and the")
        print("   arbitrary rate must not be")
        return 1
    print(f"   recovery holds: theta locked, arbitrary rate "
          f"{'locked (BAD)' if arbitrary_locked else 'not locked'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
