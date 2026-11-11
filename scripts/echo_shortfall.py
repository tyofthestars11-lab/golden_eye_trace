"""How much echo sensitivity does 709 samples lose, and what dwell would fix it?

`phi_echo.lag_sets` caps Fibonacci lags at `n // 4` so every autocorrelation still
averages over many pairs. That cap ties the sample count to the number of lags, and
the lag count sets the null spread: the contrast is a mean over `m` lags, so its
uncertainty falls as `1/sqrt(m)` while the signal stays put.

Detecting theta is not the requirement, though — the calibration requires theta to
be *separated* from an arbitrary rate, because a sta- tistic that scores every
chirp highly is not measuring the lock. This script therefore reports both p
values per sample count: p(theta), which says the contrast is unusual, and
p(0.31), which says whether a non-locked rate would have passed too. A count is
usable only when the first is at most 0.01 and the second is at least 0.05, which
is exactly the gate `phi_echo.py` runs at its calibration length.

Dwell is quoted at the LWA cadence, 1.284 s per sample, so the requirement reads
as observing time rather than sample count.

    python3 scripts/echo_shortfall.py
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
    SEED,
    THETA,
    echo_contrast,
    lag_sets,
    p_value,
)

TSAMP_S = 1.284      # LWA drspec row duration
DRAWS = 2000     # same null size as the phi_echo calibration, so the two agree
CASE_DRAWS = 12
ARBITRARY = 0.31


def draws(n: int, fibs: list[int], others: list[int], rng,
          nu: float | None, count: int) -> np.ndarray:
    t = np.arange(n)
    out = []
    for _ in range(count):
        rate = rng.uniform(1e-3, 0.5) if nu is None else nu
        x = np.cos(2.0 * math.pi * rate * t) + 0.5 * rng.normal(size=n)
        out.append(echo_contrast(x, fibs, others))
    return np.array(out)


def main() -> int:
    rng = np.random.default_rng(SEED)
    counts = [177, 246, 354, 709, 987, 1974, 3948, 5872]
    rows = []
    print("ECHO SENSITIVITY VERSUS SAMPLE COUNT")
    print(f"   theta = {THETA:.6f}; noise 0.5 against a unit-amplitude carrier")
    print(f"   lags are Fibonacci numbers in [13, n/4]; null {DRAWS} random "
          f"rates, cases {CASE_DRAWS} draws each")
    print()
    print(f"   {'N':>6}{'lags':>6}{'max':>6}{'null sd':>9}{'contrast':>10}"
          f"{'p theta':>9}{'p 0.31':>8}{'dwell h':>9}  usable")
    for n in counts:
        fibs, others = lag_sets(n)
        if not fibs or not others:
            continue
        null = draws(n, fibs, others, rng, None, DRAWS)
        c_theta = float(np.median(draws(n, fibs, others, rng, THETA, CASE_DRAWS)))
        c_arb = float(np.median(draws(n, fibs, others, rng, ARBITRARY, CASE_DRAWS)))
        p_theta = p_value(c_theta, null)
        p_arb = p_value(c_arb, null)
        usable = p_theta <= 0.01 and p_arb >= 0.05
        borderline = abs(p_arb - 0.05) <= 0.01
        dwell_h = n * TSAMP_S / 3600.0
        rows.append({"n": n, "n_lags": len(fibs), "max_lag": max(fibs),
                     "null_mean": float(null.mean()), "null_sd": float(null.std()),
                     "contrast_theta": c_theta, "contrast_arbitrary": c_arb,
                     "p_theta": p_theta, "p_arbitrary": p_arb,
                     "dwell_h": dwell_h, "usable": usable,
                     "borderline": borderline})
        mark = "no" if not usable else ("edge" if borderline else "yes")
        print(f"   {n:>6}{len(fibs):>6}{max(fibs):>6}{null.std():>9.3f}"
              f"{c_theta:>10.3f}{p_theta:>9.4f}{p_arb:>8.4f}{dwell_h:>9.2f}  {mark}")

    usable = [r for r in rows if r["usable"] and not r["borderline"]]
    first = usable[0] if usable else None
    current = next((r for r in rows if r["n"] == 709), None)
    print()
    if current:
        verdict = "passes" if current["usable"] else "fails"
        edge = " on the boundary" if current["borderline"] else ""
        print(f"   the 709-sample LWA file has {current['n_lags']} lags: "
              f"p(theta) = {current['p_theta']:.4f}, p(0.31) = "
              f"{current['p_arbitrary']:.4f}, so it {verdict} the gate"
              f"{edge} (margin {current['p_arbitrary'] - 0.05:+.4f})")
    if first:
        extra = max(0, first["n"] - 709) * TSAMP_S / 60.0
        print(f"   first count clear of the boundary: N = {first['n']} "
              f"({first['n_lags']} lags, p(0.31) = {first['p_arbitrary']:.4f}), "
              f"dwell {first['dwell_h']:.2f} h")
        print(f"   from 709 samples that is {extra:.0f} more minutes of "
              f"continuous integration at 1.284 s")
        print(f"   the calibrated configuration is N = 5872, 10 lags, "
              f"{5872 * TSAMP_S / 60:.0f} min")
    else:
        print("   no count here clears the boundary: the contrast never separates")

    out = {"tsamp_s": TSAMP_S, "theta": THETA, "arbitrary_rate": ARBITRARY,
           "rows": rows, "first_usable": first,
           "gate": "p(theta) <= 0.01 and p(arbitrary) >= 0.05"}
    path = os.path.join(ROOT, "prereg", "echo_shortfall.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"\n   wrote prereg/echo_shortfall.json")
    return 0 if first else 1


if __name__ == "__main__":
    raise SystemExit(main())
