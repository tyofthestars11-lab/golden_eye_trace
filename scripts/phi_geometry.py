"""Which instrument geometry could put a drifting beacon on the phi phase lock.

The lock says a drifting line advances by ``theta = 1/phi**2`` turns per sample,
where a sample is one integration. The rate in cycles per sample is

    nu = drift * tsamp / channel,

so a drift sits on the lock when ``nu = theta``, and the geometry that puts it
there is

    tsamp = theta * channel / drift.

This script answers the two questions that follow. For a realistic drift, what
integration time does the lock demand, and is that time something an instrument
can give? And if the integration time is left free, as the analysis binning is,
what is left of the claim?

It also measures the lattice's coverage of a logarithmic axis, which matters more
than any single hit: at 5% tolerance successive rungs are 1.2018 apart, so the
tolerance bands cover a fraction of log-drift space that is worth knowing before a
hit is called a hit.

    python3 scripts/phi_geometry.py
"""

from __future__ import annotations

import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.rungs import PHI  # noqa: E402

THETA = 1.0 / (PHI * PHI)
TOLERANCE = 0.05
# phi_echo.lag_sets caps lags at n//4, so ten Fibonacci lags (the number the
# echo is calibrated at) need n//4 >= 987, i.e. n >= 3948. At 987 samples only
# seven lags are available, which was an error in the earlier note.
ECHO_N = 3948
ECHO_N_LONG = 5872    # the length of real tracks used in the echo harness


def tsamp_for(channel_hz: float, drift_hz_s: float) -> float:
    return THETA * channel_hz / drift_hz_s


def alias_points(kmax: int) -> list[float]:
    """Lock rates a discrete-time measurement aliases onto, inside Nyquist.

    A rate is defined modulo one cycle per sample, so an advance of k theta is
    indistinguishable from k theta mod 1. theta = 2 - phi is irrational, so those
    folded points are equidistributed: the more integers are allowed, the denser
    the lock becomes in the window and the less it constrains.
    """
    folded = [min(p, 1.0 - p) for p in ((k * THETA) % 1.0 for k in range(1, kmax + 1))]
    return sorted(folded)


def main() -> int:
    out: dict = {"theta": THETA, "tolerance": TOLERANCE,
                 "frequency_ratio_per_rung": PHI}

    print("GEOMETRY OF THE PHI PHASE LOCK")
    print(f"   theta = {THETA:.6f} cycles per sample")
    print(f"   lock condition: drift * tsamp / channel = theta")
    print(f"   so tsamp = theta * channel / drift")
    print()

    print("   A. INTEGRATION TIME THE LOCK DEMANDS (seconds)")
    channels = [("LWA DRX", 19_140.0), ("STARE2", 122_070.0),
                ("3 kHz channel", 3_000.0), ("700 Hz channel", 700.0),
                ("100 Hz channel", 100.0)]
    drifts = [0.01, 0.1, 1.0, 10.0, 100.0, 1_000.0]
    header = "   " + f"{'channel':<16}" + "".join(f"{d:>12g}" for d in drifts)
    print(header)
    print("   " + f"{'':<16}" + "".join(f"{'Hz/s':>12}" for _ in drifts))
    table = {}
    for name, ch in channels:
        row = {d: tsamp_for(ch, d) for d in drifts}
        table[name] = row
        print(f"   {name:<16}" + "".join(f"{row[d]:>12.4g}" for d in drifts))
    out["tsamp_for_lock_s"] = table
    print()

    print("   B. THE TWO BANDS ON HAND")
    bands = [
        ("LWA DRX", 19_140.0, 1.284, 7_453.7, 19.1, 11),
        ("STARE2", 122_070.0, 0.4, 152_588.0, 0.0, 3),
    ]
    print(f"   {'band':<10}{'channel Hz':>12}{'tsamp s':>10}"
          f"{'theta Hz/s':>12}{'Nyquist':>11}{'% Nyq':>8}"
          f"{'max |D|':>10}{'ratio':>9}")
    out["bands"] = {}
    for name, ch, tsamp, nyq, observed, n in bands:
        d_theta = THETA * ch / tsamp
        ratio = (d_theta / observed) if observed else None
        shown = f"{ratio:>8.0f}x" if ratio else f"{'n/a':>9}"
        print(f"   {name:<10}{ch:>12,.0f}{tsamp:>10.3f}{d_theta:>12,.0f}"
              f"{nyq:>11,.0f}{100 * d_theta / nyq:>7.1f}%{observed:>10.1f}{shown}")
        out["bands"][name] = {"channel_hz": ch, "tsamp_s": tsamp,
                              "theta_hz_s": d_theta, "nyquist_hz_s": nyq,
                              "fraction_of_nyquist": d_theta / nyq,
                              "max_observed_hz_s": observed,
                              "observed_ratio": ratio}
    print()

    print("   C. WHAT A SLOW BEACON WOULD NEED")
    dwell = {}
    for label, n in (("echo (10 lags)", ECHO_N), ("echo (real length)", ECHO_N_LONG)):
        for drift in (1.0, 0.1):
            t = tsamp_for(19_140.0, drift)
            total = n * t
            bw = THETA * n * 19_140.0
            dwell[f"{label} D={drift:g}"] = {
                "tsamp_s": t, "dwell_s": total, "dwell_days": total / 86_400.0,
                "bandwidth_hz": bw,
            }
            print(f"   {label:<20} D={drift:>4g} Hz/s  tsamp {t:>9,.0f} s  "
                  f"dwell {total / 86_400:>7.1f} d  bandwidth {bw / 1e6:>6.2f} MHz")
    out["slow_beacon"] = dwell
    print()

    per_value = 2.0 * TOLERANCE * THETA / 0.5
    need = 2.0 * TOLERANCE * THETA
    print("   D. HOW MANY LOCK POINTS FIT INSIDE NYQUIST")
    print(f"   the drift lattice is additive: nu = k*theta, spacing {THETA:.4f}")
    print(f"   Nyquist is 0.5, so k = 1 at {THETA:.4f} is {100 * THETA / 0.5:.1f}% of it")
    print(f"   and k = 2 at {2 * THETA:.4f} is already outside: one physical rate")
    print(f"   chance a random drift falls within {100 * TOLERANCE:.0f}% of it: "
          f"{per_value:.4f}")
    kmax, worst = None, None
    for k in range(1, 500):
        pts = [0.0] + alias_points(k) + [0.5]
        gap = max(b - a for a, b in zip(pts, pts[1:]))
        if kmax is None and gap <= need:
            kmax = k
        worst = (k, gap)
    print()
    print("   but a rate is defined modulo one cycle per sample, so k*theta folds")
    print("   back into the window, and theta is irrational: those points fill it.")
    if kmax:
        print(f"   allowing k up to {kmax} already leaves no gap wider than "
              f"{need:.4f} in [0, 0.5],")
        print(f"   so every drift matches some integer advance within {100 * TOLERANCE:.0f}%")
        print("   unless k is bounded by something other than the data")
    out["alias"] = {"per_value_chance": per_value, "tolerance_in_rate": need,
                    "k_that_covers_nyquist": kmax,
                    "points_at_kmax": alias_points(kmax) if kmax else None}
    print()

    print("   E. THE FREE KNOB")
    for drift in (1.0, 10.0):
        print(f"   a {drift:g} Hz/s beacon in a 19.14 kHz channel sits on the lock")
        print(f"      when tsamp = {tsamp_for(19_140.0, drift):,.0f} s "
              f"({tsamp_for(19_140.0, drift) / 3600:.2f} h)")
    print("   the analysis may choose that tsamp by dedoppler-summing blocks,")
    print("   so the lock constrains the binning, not the source, unless the")
    print("   sample time is fixed by the instrument before the data are seen.")
    out["free_knob"] = {"note": "tsamp set by the analyst places any drift on theta",
                        "lwa_1hz_tsamp_s": tsamp_for(19_140.0, 1.0)}

    print("   F. CHANNEL OVER CADENCE: WHAT THE LOCK NEEDS AND WHAT EXISTS")
    need = 1.0 / THETA      # channel Hz per second of sample, for a 1 Hz/s drift
    print(f"   a 1 Hz/s chirp locks at channel/tsamp = drift/theta = "
          f"{need:.2f} Hz/s")
    archives = [
        ("LWA1 drspec", 19_140.6, 1.284, "open, read here"),
        ("STARE2", 122_070.0, 0.4, "open, one integration"),
        ("CAMRAS Dwingeloo", 136_719.0, 468.11429e-6, "open, direct GET"),
        ("OVRO-LWA spec_fits", 100_547.0, 0.256, "open S3, tested here"),
        ("BL L band (nominal)", 2.9, 300.0, "metadata only, 404 here"),
    ]
    print(f"   {'archive':<22}{'channel Hz':>12}{'tsamp s':>11}"
          f"{'channel/tsamp':>15}{'theta drift':>14}  access")
    for name, ch, t, access in archives:
        ratio = ch / t
        print(f"   {name:<22}{ch:>12,.1f}{t:>11.4g}{ratio:>15,.1f}"
              f"{THETA * ratio:>14,.4g}  {access}")
    out["channel_over_tsamp"] = {
        "needed_for_1hz_s": need,
        "archives": {name: {"channel_hz": ch, "tsamp_s": t,
                            "ratio_hz_s": ch / t,
                            "theta_drift_hz_s": THETA * ch / t,
                            "access": access}
                     for name, ch, t, access in archives},
    }
    print()

    path = os.path.join(ROOT, "prereg", "phi_geometry.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"\n   wrote {os.path.relpath(path, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
