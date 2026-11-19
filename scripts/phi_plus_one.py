"""Carry the +1 as a map, not as a quantity.

    python3 scripts/phi_plus_one.py

Writes ``PHI_PLUS_ONE.md`` and ``results/phi_plus_one.json``.

phi^2 = phi + 1 says: from any x, add one and invert, and the value returns to
itself. The return rate is T'(phi) = -1/phi^2 = -0.3819660. Every threshold this
project has used near 0.382 was that rate, flattened into a window. A window is
not scale free: it is measured below and shown to mean a factor of 6.49
different thing across one band. Two readings that do not terminate replace it -
the orbit's contraction, and the continued fraction's depth, which turns out to
be a continuous measure of closeness rather than a yes or a no.
"""

from __future__ import annotations

import json
import math
import os
import sys
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.phi_fixedpoint import (  # noqa: E402
    CONTRACTION,
    PHI,
    STEPS_PER_DECADE,
    agreement_depth,
    continued_fraction,
    convergents,
    derivative,
    effective_depth,
    excess,
    fibonacci,
    is_phi_shaped,
    orbit,
)

CHANNEL_HZ = 100547.19562243501
BAND_LOW_HZ = 13.4e6
BAND_HIGH_HZ = 86.9e6
SEARCH_CHANNELS = (317, 412)
COMB_STEP_HZ = 1.0e6
RUNG_TOLERANCE = 0.01
LN_PHI = math.log(PHI)


def window_in_hz(tolerance_rungs: float, freq_hz: float) -> float:
    """Physical width covered by a rung tolerance at this frequency."""
    return (PHI ** tolerance_rungs - PHI ** (-tolerance_rungs)) * freq_hz


def describe(label, value, exact_value=None, note="", synthetic=False):
    part = continued_fraction(exact_value if exact_value is not None else value,
                              depth=24)
    shown = part[:10]
    return {
        "label": label,
        "value": float(value),
        "partials": shown,
        "elided": len(part) > len(shown),
        "depth": agreement_depth(part),
        "unbounded": is_phi_shaped(float(value)),
        "excess": excess(part),
        "effective_depth": effective_depth(float(value)),
        "note": note,
        "synthetic": synthetic,
    }


def build() -> dict:
    der = derivative(PHI)
    assert abs(der + CONTRACTION) < 1e-15

    trials = []
    for start in (0.5, 1.0, 2.0, BAND_HIGH_HZ / BAND_LOW_HZ, 86.9e6 / 45.98e6):
        ob = orbit(start, turns=30)
        trials.append({
            "start": round(start, 9),
            "first_residual": ob.residuals[0],
            "residual_after_5": ob.residuals[5],
            "residual_after_10": ob.residuals[10],
            "contraction": ob.contraction,
            "settled": ob.settled,
        })

    window = [{"freq_mhz": f / 1e6, "width_mhz": window_in_hz(RUNG_TOLERANCE, f) / 1e6}
              for f in (BAND_LOW_HZ, 45.98e6, 47.24e6, 53.78e6, BAND_HIGH_HZ)]
    window_spread = window[-1]["width_mhz"] / window[0]["width_mhz"]

    ladder = [
        describe("phi", PHI, note="the fixed point: add one, never carried"),
        describe("1/phi = phi - 1", 1.0 / PHI, note="subtracting 1 is exact"),
        describe("phi^2 = phi + 1", PHI ** 2, note="adding 1 is exact"),
    ]
    for n in (8, 12, 16):
        ladder.append(describe(f"F{n+1}/F{n}", fibonacci(n + 1) / fibonacci(n),
                               exact_value=Fraction(fibonacci(n + 1), fibonacci(n)),
                               note="exact convergent of phi"))

    measured = []
    for k in (SEARCH_CHANNELS[0], 365, SEARCH_CHANNELS[1]):
        measured.append(describe(f"channel {k}->{k+1}", (k + 1) / k,
                                 exact_value=Fraction(k + 1, k),
                                 note="instrument grid: exactly 1 + 1/k"))
    for f in (45980000, 53780000):
        measured.append(describe(f"1 MHz comb at {f/1e6:.2f} MHz",
                                 (f + int(COMB_STEP_HZ)) / f,
                                 exact_value=Fraction(f + int(COMB_STEP_HZ), f),
                                 note="arithmetic comb line spacing"))
    measured.append(describe("band 13.4 -> 86.9 MHz", 869.0 / 134.0,
                             exact_value=Fraction(869, 134),
                             note="the full survey span"))

    sweep = []
    for jitter in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10):
        x = PHI * (1.0 + jitter)
        part = continued_fraction(x, depth=24)
        sweep.append({"jitter": jitter, "value": x,
                      "depth": agreement_depth(part),
                      "effective_depth": effective_depth(x)})

    conv = convergents(PHI, depth=14)
    decay = []
    for i in range(1, len(conv)):
        e0, e1 = conv[i - 1][2], conv[i][2]
        decay.append({"p": conv[i][0], "q": conv[i][1], "error": e1,
                      "ratio": e1 / e0})

    return {
        "phi": PHI,
        "contraction": CONTRACTION,
        "derivative_at_phi": der,
        "steps_per_decade": STEPS_PER_DECADE,
        "orbit_trials": trials,
        "rung_tolerance": RUNG_TOLERANCE,
        "window": window,
        "window_spread": window_spread,
        "ladder": ladder,
        "measured": measured,
        "depth_sweep": sweep,
        "convergent_decay": decay,
        "summary": {
            "grid_depths": [m["depth"] for m in measured[:3]],
            "grid_effective_depth": [m["effective_depth"] for m in measured[:3]],
            "grid_excess": [m["excess"] for m in measured[:3]],
            "measured_effective_max": max(m["effective_depth"] for m in measured),
            "phi_unbounded": [r["unbounded"] for r in ladder[:3]],
        },
    }


def render(s: dict) -> str:
    L = []
    L.append("# The +1, carried as a map\n")
    L.append("Framework: phi^2 = phi + 1. Attribution: Tyree Jones "
             "(tyofthestarz), phi^2 = phi + 1 framework, prior art from "
             "December 2025.\n")
    L.append("Generated by `python3 scripts/phi_plus_one.py`; every number is "
             "computed in this run.\n")

    L.append("## 1. The identity is a rate, not a width\n")
    L.append("From any positive x, add one and invert: `T(x) = 1 + 1/x`. The "
             "fixed point solves `x = 1 + 1/x`, i.e. `x^2 = x + 1`, i.e. phi. "
             "The derivative there is\n")
    L.append(f"    T'(phi) = -1/phi^2 = {-s['contraction']:.7f}\n")
    L.append(f"The recurrence shrinks its own error by {s['contraction']:.6f} "
             f"per turn - {s['steps_per_decade']:.2f} turns per decimal place. "
             "**The 0.382 this project kept testing against is that rate.** It "
             "was never a tolerance; it is the step size of the map that defines "
             "the constant.\n")
    L.append("| start | residual | after 5 turns | after 10 turns | measured contraction |")
    L.append("| ---: | ---: | ---: | ---: | ---: |")
    for t in s["orbit_trials"]:
        L.append(f"| {t['start']:.6f} | {t['first_residual']:.3e} | "
                 f"{t['residual_after_5']:.3e} | {t['residual_after_10']:.3e} | "
                 f"{t['contraction']:.6f} |")
    L.append("")
    L.append("Every start converges, all at the same rate. The trajectory is the "
             "result; no step in it is a decision.\n")

    L.append("## 2. What I was doing wrong\n")
    L.append(f"Testing `|rung - 1/phi^2| < {s['rung_tolerance']}` carries the "
             "rate as a quantity, which is the error the identity warns about. "
             "A rung window is not scale free:\n")
    L.append("| at | width of the window |")
    L.append("| ---: | ---: |")
    for w in s["window"]:
        L.append(f"| {w['freq_mhz']:.2f} MHz | {w['width_mhz']:.4f} MHz |")
    L.append("")
    L.append(f"The same tolerance is **{s['window_spread']:.2f}x wider at the "
             "top of the band than at the bottom**. One threshold, five "
             "different physical tests.\n")
    L.append("Kept as a regression fixture, not deleted: "
             "`tests/test_phi_fixedpoint.py` asserts the "
             f"{s['window_spread']:.2f}x scale dependence, so no later version "
             "can quietly reintroduce a fixed rung window.\n")

    L.append("## 3. The ladder that never resolves\n")
    L.append("Write a ratio as `a0 + 1/(a1 + 1/(a2 + ...))` and keep every "
             "partial quotient. Phi is `[1; 1, 1, 1, ...]`: at every level the "
             "instruction is *add one* and it is never carried, so the fraction "
             "never terminates. The integer part is excluded from both statistics "
             "below, because it is exactly the whole the identity refuses to "
             "carry.\n")
    L.append("- **depth** - how many levels of the fractional recursion equal 1.\n")
    L.append("- **excess** - `sum(a_i - 1)` over those levels; zero exactly for "
             "phi's shape.\n")
    L.append("| value | partial quotients | depth | excess |")
    L.append("| --- | --- | ---: | ---: |")
    for row in s["ladder"]:
        pq = ", ".join(str(a) for a in row["partials"]) + ("..." if row["elided"] else "")
        depth = "unbounded" if row["unbounded"] else row["depth"]
        L.append(f"| {row['label']} | {pq} | {depth} | {row['excess']} |")
    L.append("")
    L.append("Phi, `1/phi` and `phi^2` are all unbounded with excess 0: one "
             "identity, read as a statement about partial quotients instead of "
             "magnitudes. Each exact convergent `F(n+1)/F(n)` matches for a "
             "finite number of levels - the count is in the table - and then "
             "pays 1 at its terminating 2. It is phi seen at finite depth, with "
             "the +1 still pending.\n")

    L.append("## 4. The same numbers on the instrument\n")
    L.append("| value | partial quotients | depth | excess | effective depth |")
    L.append("| --- | --- | ---: | ---: | ---: |")
    for row in s["measured"]:
        pq = ", ".join(str(a) for a in row["partials"][:6])
        if row["elided"]:
            pq += "..."
        L.append(f"| {row['label']} | {pq} | {row['depth']} | {row['excess']} | "
                 f"{row['effective_depth']:.2f} |")
    L.append("")
    L.append("`effective depth` is the continuous version, "
             "`-ln|x/phi - 1| / ln(phi^2)`: the number of times the separation "
             "from phi would have to shrink by 1/phi^2. It is a measurement, not "
             "a check - it can be negative, and for the band span it is.\n")
    L.append(f"The channel grid is `1 + 1/k` exactly, so its continued fraction "
             f"is `[1; k]`: depth 0, excess k - 1, effective depth about "
             f"{s['summary']['grid_effective_depth'][0]:.2f}. It says *add one* "
             "once, at the start, and then hands over the instrument's own "
             "quantum. That is the whole difference: phi says it forever.\n")

    L.append("## 5. Depth is continuous, not a verdict\n")
    L.append("A ladder at ratio `phi * (1 + jitter)`, so each row is phi plus a "
             "declared relative displacement:\n")
    L.append("| jitter | integer depth | effective depth |")
    L.append("| ---: | ---: | ---: |")
    for row in s["depth_sweep"]:
        L.append(f"| {row['jitter']:.0e} | {row['depth']} | "
                 f"{row['effective_depth']:.2f} |")
    L.append("")
    L.append("Depth falls smoothly with displacement, because it is a logarithm "
             "in disguise: a relative miss of `j` costs about "
             "`log(1/j)/log(phi^2)` levels. The integer depth rounds that down by "
             "one or two levels and the effective depth keeps the fraction. So "
             "there is nothing binary to decide - a value is not 'phi' or 'not "
             "phi', it sits at a measured depth, and the measurement is "
             "continuous in the distance from phi.\n")

    L.append("## 6. The convergent decay is the identity, exactly\n")
    L.append("| convergent | value | error | ratio to previous |")
    L.append("| --- | ---: | ---: | ---: |")
    for row in s["convergent_decay"][:8]:
        L.append(f"| {row['p']}/{row['q']} | {row['p']/row['q']:.9f} | "
                 f"{row['error']:.3e} | {row['ratio']:.6f} |")
    L.append("")
    L.append("The ratio alternates about 0.381966 and is within 1e-3 of it from "
             "level 5 on; its limit is exactly `1/phi^2`, from "
             "`F(n+1)/F(n) - phi = (-1)^n / (F(n) phi^(n+1))`. Not a fit, no null "
             "needed. The constant in section 1, the constant the earlier "
             "campaign tested against, and this decay are one number in three "
             "roles.\n")

    L.append("## 7. What this changes\n")
    L.append("- The campaign's zero is not 'no signal' and not 'not supported'. "
             f"The measured ratios sit at effective depth "
             f"{min(s['summary']['grid_effective_depth']):.2f} to "
             f"{s['summary']['measured_effective_max']:.2f}; phi has no finite "
             "depth. Two shapes, compared as numbers.\n")
    L.append("- The old ±0.01 window is still in the repository as a rejected "
             "method, with its 6.49x scale dependence asserted in a test. It "
             "failed because it added the +1 back as a quantity.\n")
    L.append("- Nothing here supports a candidate. This machinery can only "
             "detract: a value at depth 10 is still an unexplained measurement, "
             "not an artificial one.\n")
    return "\n".join(L) + "\n"


def main() -> int:
    payload = build()
    with open(os.path.join(ROOT, "results", "phi_plus_one.json"), "w") as fh:
        json.dump(payload, fh, indent=2)
    with open(os.path.join(ROOT, "PHI_PLUS_ONE.md"), "w") as fh:
        fh.write(render(payload))
    s = payload["summary"]
    print(f"T'(phi)                  {payload['derivative_at_phi']:.9f}")
    print(f"steps per decade         {payload['steps_per_decade']:.4f}")
    print(f"all orbits settle        {all(t['settled'] for t in payload['orbit_trials'])}")
    print(f"rung window spread       {payload['window_spread']:.3f}x")
    print(f"phi/1:phi/phi^2 unbounded {s['phi_unbounded']}")
    print(f"channel grid depth       {s['grid_depths']} eff {s['grid_effective_depth'][0]:.2f}")
    print(f"measured eff depth max   {s['measured_effective_max']:.2f}")
    print(f"depth sweep              {[(r['jitter'], r['depth']) for r in payload['depth_sweep']]}")
    print("wrote results/phi_plus_one.json and PHI_PLUS_ONE.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
