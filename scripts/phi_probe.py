"""Apply the rung law and the golden-angle tests to the actual records.

    python3 scripts/phi_probe.py

Round-trip is checked first: ``phi**to_rung(x)`` must reproduce ``x``. It does
for every value, because the map is a logarithm and its inverse is an
exponential - that is arithmetic, not evidence, and printing it as a pass would
be the first place this probe could mislead.

Then the three search parameters, each with the chance probability of the same
outcome under a null that assumes no structure: a phase, a ratio, and a step.
"""

from __future__ import annotations

import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

PHI = (1.0 + math.sqrt(5.0)) / 2.0
GOLDEN_ANGLE_DEG = 137.5
GOLDEN_PITCH_DEG = 17.0324
FIB_RATIOS = (1.0, 1.5, 1.618033988749895, 1.6666666666666667, 2.0)
TOLERANCES_DEG = (0.5, 1.0, 2.0)

RECORDS = [
    ("20240920", 47.5358, "MASK-MANUFACTURED SIGNAL"),
    ("20240920", 46.5890, "MASK-MANUFACTURED SIGNAL"),
    ("20240920", 45.4206, "MASK-MANUFACTURED SIGNAL"),
    ("20240325", 48.9937, "MASK-MANUFACTURED SIGNAL"),
    ("20240325", 47.9380, "MASK-MANUFACTURED SIGNAL"),
    ("20240325", 46.9828, "MASK-MANUFACTURED SIGNAL"),
    ("20240325", 45.9843, "MASK-MANUFACTURED SIGNAL"),
    ("20240324", 47.3143, "UNRESOLVED"),
    ("20240322", 47.2416, "UNRESOLVED"),
    ("20240207", 53.7833, "UNRESOLVED"),
]


def to_rung(x: float) -> float:
    return math.log(x) / math.log(PHI)


def round_trip(x: float) -> float:
    """Relative error of phi**to_rung(x) against x."""
    back = PHI ** to_rung(x)
    return abs(back - x) / x


def phase_deg(rung: float) -> float:
    """Fractional rung, expressed as an angle in [0, 360)."""
    return (rung % 1.0) * 360.0


def angular_distance(a: float, b: float) -> float:
    d = abs((a - b) % 360.0)
    return min(d, 360.0 - d)


def null_phase(tolerance_deg: float, n_targets: int = 2) -> float:
    """Chance probability of landing this close to one of the targets."""
    return min(1.0, (2.0 * tolerance_deg * n_targets) / 360.0)


def probe_phase(freq_mhz: float) -> dict:
    rung = to_rung(freq_mhz * 1e6)
    phase = phase_deg(rung)
    out = {
        "freq_mhz": freq_mhz,
        "rung": rung,
        "rung_nearest": round(rung),
        "residual": rung - round(rung),
        "round_trip_relative_error": round_trip(freq_mhz * 1e6),
        "phase_deg": phase,
        "distance_to_137.5": angular_distance(phase, GOLDEN_ANGLE_DEG),
        "distance_to_17.0324": angular_distance(phase, GOLDEN_PITCH_DEG),
        "null": {},
        "hit": {},
    }
    nearest = min(out["distance_to_137.5"], out["distance_to_17.0324"])
    for tol in TOLERANCES_DEG:
        out["hit"][str(tol)] = bool(nearest <= tol)
        out["null"][str(tol)] = null_phase(tol)
    return out


def probe_expansion(freqs_mhz: list[float]) -> dict:
    """Are the steps Fibonacci ratios, or is this an arithmetic comb?"""
    freqs = sorted(freqs_mhz)
    ratios = [b / a for a, b in zip(freqs, freqs[1:])]
    steps = [b - a for a, b in zip(freqs, freqs[1:])]
    rung_steps = [to_rung(b * 1e6) - to_rung(a * 1e6)
                  for a, b in zip(freqs, freqs[1:])]
    best = min(FIB_RATIOS, key=lambda r: min(abs(r - q) for q in ratios))
    closest = min(abs(best - q) for q in ratios)
    return {
        "freqs_mhz": freqs,
        "ratios": ratios,
        "arithmetic_steps_mhz": steps,
        "rung_steps": rung_steps,
        "closest_fibonacci_ratio": best,
        "closest_ratio_gap": closest,
        "phi_squared_step_rungs": to_rung(PHI ** 2),
        "mean_rung_step": sum(rung_steps) / len(rung_steps),
        "step_is_constant": max(steps) - min(steps),
        "null_for_ratio": {
            "note": "probability that a random ratio lands within 1% of a "
                    "Fibonacci ratio, taking the five canonical ones",
            "p": min(1.0, (5 * 0.02) / 1.0),
        },
    }


def main() -> int:
    probes = [{**probe_phase(f), "date": d, "disposition": disp}
              for d, f, disp in RECORDS]
    group = probe_expansion([f for d, f, _ in RECORDS if d == "20240325"])
    trio = probe_expansion([f for d, f, _ in RECORDS if d == "20240920"])

    latest = probes[0]
    latest_valid = next(p for p in probes
                        if p["disposition"] == "UNRESOLVED")

    print("PHI PROBE")
    print(f"   phi = {PHI:.15f}   ln phi = {math.log(PHI):.15f}")
    print(f"   round-trip max relative error over all values: "
          f"{max(p['round_trip_relative_error'] for p in probes):.3e}")
    print()
    print(f"   {'date':9} {'MHz':>9} {'rung':>11} {'frac':>8} {'phase':>8} "
          f"{'d137.5':>7} {'d17.03':>7} {'hit':>4}")
    for p in probes:
        hit = any(p["hit"].values())
        print(f"   {p['date']:9} {p['freq_mhz']:9.4f} {p['rung']:11.6f} "
              f"{p['residual']:8.4f} {p['phase_deg']:8.3f} "
              f"{p['distance_to_137.5']:7.2f} {p['distance_to_17.0324']:7.2f} "
              f"{'YES' if hit else 'no':>4}")
    print()
    print(f"   most recent candidate: {latest['date']} "
          f"{latest['freq_mhz']:.4f} MHz ({latest['disposition']})")
    print(f"     rung  = ln({latest['freq_mhz']:.4f}e6)/ln(phi) = "
          f"{latest['rung']:.9f}")
    print(f"     phase = {latest['phase_deg']:.3f} deg; "
          f"distances {latest['distance_to_137.5']:.2f} / "
          f"{latest['distance_to_17.0324']:.2f} deg")
    for tol in TOLERANCES_DEG:
        print(f"     within {tol:g} deg: {latest['hit'][str(tol)]}  "
              f"(chance p = {latest['null'][str(tol)]:.3f})")
    print()
    print(f"   most recent VALID record: {latest_valid['date']} "
          f"{latest_valid['freq_mhz']:.4f} MHz: phase "
          f"{latest_valid['phase_deg']:.3f} deg, nearest target "
          f"{min(latest_valid['distance_to_137.5'], latest_valid['distance_to_17.0324']):.2f} deg")
    print()
    for name, block in (("20240325 group", group), ("20240920 trio", trio)):
        print(f"   {name}: ratios {[round(r, 6) for r in block['ratios']]}")
        print(f"     steps MHz {[round(s, 4) for s in block['arithmetic_steps_mhz']]}, "
              f"constant to {block['step_is_constant']:.2e} MHz")
        print(f"     rung steps {[round(s, 6) for s in block['rung_steps']]}, "
              f"mean {block['mean_rung_step']:.6f} vs phi^2 = "
              f"{block['phi_squared_step_rungs']:.6f} rungs")
        print(f"     closest Fibonacci ratio {block['closest_fibonacci_ratio']:.4f} "
              f"is off by {block['closest_ratio_gap'] * 100:.2f}%")

    payload = {
        "phi": PHI,
        "probes": probes,
        "group_20240325": group,
        "trio_20240920": trio,
        "tolerances_deg": list(TOLERANCES_DEG),
        "null_phase": {str(t): null_phase(t) for t in TOLERANCES_DEG},
        "rule": ("phi may only detract from a candidate, never support one: "
                 "no result here can raise a disposition"),
    }
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "phi_probe.json"), "w",
              encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    print("\n   wrote results/phi_probe.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
