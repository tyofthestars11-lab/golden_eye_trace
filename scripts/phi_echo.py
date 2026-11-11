"""Quantum echoes as a phi phase lock, calibrated at the length of real data.

A drifting line is a phase ramp: one channel per sample is one cycle per sample,
so the drift rate in channels per sample is the phase advance per integration. A
phi-locked line advances by the golden angle ``theta = 1/phi**2`` turns per step,
so its rate divided by theta is an integer. That is the lock.

The echo is the Fibonacci structure of ``phi**2 = [2;1,1,1,...]``: a rotation by
theta returns closest to its start after Fibonacci numbers of steps, because
those are the convergents. ``theta*89`` is 0.5% of a turn from a whole number and
``theta*144`` is 0.003%. The statistic contrasts the mean |autocorrelation| at
Fibonacci lags with the mean over every non-Fibonacci lag, which is the stable
baseline; a cosine's |autocorrelation| averages 2/pi over lags, so a theta-locked
line should stand about 0.24 above it, and any other rate should not.

This file has been wrong three times, and the calibration is what caught it each
time, so the calibration is now a gate: if theta is not detected, or if an
arbitrary rate, white noise or AR(1) noise is, the script exits non-zero and
reports nothing about real data. The three failures are worth keeping:

* summing autocorrelations at Fibonacci lags tests periodicity, not phi - an
  arbitrary-frequency sinusoid scored chi-square 517;
* the lag of the largest |autocorrelation| is dominated by half-turn returns at
  small lags and landed on a Fibonacci number 0.498 of the time;
* six hand-picked comparison lags left the contrast inside its own noise, and at
  709 samples an arbitrary rate scored 0.19 against theta's 0.25.

The length matters, so the calibration runs at 5,872 samples - the length of the
STARE2 tracks - where eleven Fibonacci lags are available instead of six.

    python3 scripts/phi_echo.py
"""

from __future__ import annotations

import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from seti_search.filterbank import Filterbank  # noqa: E402
from seti_search.rungs import PHI  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

THETA = 1.0 / (PHI * PHI)
PHI_THETA = PHI * THETA
LOCK_TOLERANCE = 0.05
SEED = 1729
DRAWS = 2000
MIN_LAG = 13
CALIBRATION_N = 5872   # STARE2 track length: the length a real claim would use


def fibonacci(limit: int) -> list[int]:
    out, a, b = [], 1, 2
    while a <= limit:
        out.append(a)
        a, b = b, a + b
    return out


def lag_sets(n: int) -> tuple[list[int], list[int]]:
    limit = max(2, n // 4)
    fibs = [f for f in fibonacci(limit) if f >= MIN_LAG]
    fib_set = set(fibonacci(limit))
    pool = [lag for lag in range(MIN_LAG, limit + 1) if lag not in fib_set]
    return fibs, pool


def autocorr(x: np.ndarray, lag: int) -> float:
    x = x - x.mean()
    norm = float(np.dot(x, x)) or 1.0
    return abs(float(np.dot(x[:-lag], x[lag:])) / norm)


def echo_contrast(x: np.ndarray, fibs: list[int], others: list[int]) -> float:
    if not fibs or not others:
        return 0.0
    return (float(np.mean([autocorr(x, f) for f in fibs]))
            - float(np.mean([autocorr(x, o) for o in others])))


def null_contrast(n: int, draws: int = DRAWS, rng=None) -> np.ndarray:
    """Contrasts for cosines at random rates, with the same 0.5 noise level."""
    rng = rng or np.random.default_rng(SEED)
    fibs, others = lag_sets(n)
    t = np.arange(n)
    return np.array([
        echo_contrast(np.cos(2.0 * math.pi * rng.uniform(1e-3, 0.5) * t)
                      + 0.5 * rng.normal(size=n), fibs, others)
        for _ in range(draws)
    ])


def p_value(observed: float, null: np.ndarray) -> float:
    return float((1 + int(np.sum(null >= observed))) / (1 + null.size))


def drift_lock(drift_hz_s: float, channel_hz: float) -> tuple[float, float]:
    nu = abs(drift_hz_s) / channel_hz
    return nu, abs(nu / THETA - round(nu / THETA))


def main() -> int:
    rng = np.random.default_rng(SEED)
    n = CALIBRATION_N
    fibs, others = lag_sets(n)
    null = null_contrast(n, rng=rng)

    print("QUANTUM ECHOES AS A PHI PHASE LOCK")
    print(f"   theta = 1/phi^2 = {THETA:.9f} cycles per sample")
    print(f"   calibrating at N={n}: Fibonacci lags {fibs}")
    print(f"   baseline: mean |r| over all {len(others)} non-Fibonacci lags in "
          f"[{MIN_LAG}, {max(2, n // 4)}]")
    print(f"   null from {null.size} random rates: {null.mean():+.3f} "
          f"+/- {null.std():.3f}")
    print()

    t = np.arange(n)
    ar = np.zeros(n)
    for i in range(1, n):
        ar[i] = 0.9 * ar[i - 1] + rng.normal()
    cases = {
        "rate = theta": np.cos(2 * math.pi * THETA * t) + 0.5 * rng.normal(size=n),
        "rate = phi*theta": np.cos(2 * math.pi * PHI_THETA * t) + 0.5 * rng.normal(size=n),
        "rate = 0.31 arbitrary": np.cos(2 * math.pi * 0.31 * t) + 0.5 * rng.normal(size=n),
        "white noise": rng.normal(size=n),
        "AR(1) red noise": ar,
    }

    print(f"   {'case':<24}{'contrast':>10}{'p':>9}  verdict")
    calibration, ok = {}, True
    for name, x in cases.items():
        observed = echo_contrast(x, fibs, others)
        p = p_value(observed, null)
        phi_case = name.startswith("rate = theta") or name.startswith("rate = phi")
        passed = (p <= 0.01) if phi_case else (p >= 0.05)
        ok &= passed
        calibration[name] = {"contrast": observed, "p": p, "passes": bool(passed)}
        verdict = ("phi detected" if passed and phi_case else
                   ("goes through" if passed else "CALIBRATION FAIL"))
        print(f"   {name:<24}{observed:>+10.4f}{p:>9.4f}  {verdict}")

    out = {"theta": THETA, "calibration_n": n, "fibonacci_lags": fibs,
           "null_mean": float(null.mean()), "null_std": float(null.std()),
           "calibration": calibration, "calibration_passed": bool(ok), "rows": []}

    if not ok:
        print("\n   CALIBRATION FAILED at this length: the contrast does not")
        print("   separate theta from an arbitrary rate, so nothing is claimed")
        print("   about real data. The statistic is not used.")
        with open(os.path.join(ROOT, "prereg", "phi_echo.json"), "w",
                  encoding="utf-8") as fh:
            json.dump(out, fh, indent=2, sort_keys=True)
            fh.write("\n")
        return 1
    print("   calibration holds: theta and phi*theta detected, other rates pass")
    print()

    print("   REAL DRIFT LINES: RATE IN CYCLES PER SAMPLE")
    print(f"   {'set':<7}{'MHz':>12}{'drift Hz/s':>12}{'rate':>9}"
          f"{'dist to int':>13}  verdict")
    rows = []
    for path, cfg, label, max_ints in (
        (os.path.join(ROOT, "data", "stare2_ovro.fil"),
         SearchConfig(boxcar=3, max_drift_hz_s=4.0, snr_threshold=3.0, max_hits=10),
         "STARE2", None),
        (os.path.join(ROOT, "data", "lwa_drspec.fil"),
         SearchConfig(boxcar=3, max_drift_hz_s=8.0, snr_threshold=3.0, max_hits=10),
         "LWA", 300),
    ):
        if not os.path.exists(path):
            continue
        fb = Filterbank.open(path)
        channel_hz = abs(fb.foff) * 1e6
        resolution = channel_hz / fb.total_time
        hits, _meta = search_file(path, cfg, chan_chunk=2048, max_ints=max_ints)
        for hit in hits:
            if abs(hit.drift_hz_s) <= resolution:
                rows.append({"set": label, "freq_mhz": hit.freq_mhz,
                             "drift_hz_s": hit.drift_hz_s, "rate": None,
                             "distance": None,
                             "verdict": "not evaluable: drift below resolution"})
                continue
            nu, distance = drift_lock(hit.drift_hz_s, channel_hz)
            rows.append({"set": label, "freq_mhz": hit.freq_mhz,
                         "drift_hz_s": hit.drift_hz_s, "rate": nu,
                         "distance": distance,
                         "verdict": "locked" if distance <= LOCK_TOLERANCE
                                    else "off lattice"})
        out.setdefault("resolution_hz_s", {})[label] = resolution

    for row in rows:
        rate = "-" if row["rate"] is None else f"{row['rate']:.4f}"
        dist = "-" if row["distance"] is None else f"{row['distance']:.4f}"
        print(f"   {row['set']:<7}{row['freq_mhz']:>12.4f}{row['drift_hz_s']:>12.3g}"
              f"{rate:>9}{dist:>13}  {row['verdict']}")

    evaluable = [r for r in rows if r["rate"] is not None]
    chance = 2.0 * LOCK_TOLERANCE * THETA / 0.5
    print()
    print(f"   {len(evaluable)} of {len(rows)} candidates carry a measurable rate;")
    print(f"   chance a random rate lands within {LOCK_TOLERANCE} theta of an "
          f"integer is {chance:.4f}")
    print(f"   locked: {sum(1 for r in evaluable if r['distance'] <= LOCK_TOLERANCE)}"
          f" of {len(evaluable)}")
    out["rows"] = rows
    out["chance_rate"] = chance

    with open(os.path.join(ROOT, "prereg", "phi_echo.json"), "w",
              encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"\n   wrote prereg/phi_echo.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
