"""Nitrogen-vacancy centres, read through the rung lattice with its null attached.

An NV centre is a much better-defined object than a radio candidate: its
transition frequencies are tabulated, its linewidths are known, and its
geometry is exact. That makes it the right place to ask what the framework's
address arithmetic would actually have to do to say anything about a solid.

Five sections:

* the constants, with the uncertainties that come with them;
* the rung address of every constant expressed in six units - the same
  physical quantity, six different answers about lattice membership;
* the lattice geometry: every angle between low-index diamond directions is
  an arccosine of a rational or a simple radical, and the golden angle is not
  one of them;
* the pitch, applied to the one place it could bite (D against the hydrogen
  line);
* the tolerance budget: what the 0.05-rung cone is in MHz at 2.87 GHz, against
  the 2 MHz the instrument actually resolves.

    python3 scripts/nv_check.py
"""

from __future__ import annotations

import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.rungs import LN_PHI, PHI, rung  # noqa: E402

MU_DEG = 17.0324
SLOPE = math.tan(math.radians(MU_DEG))
GOLDEN_ANGLE_DEG = 360.0 / (PHI * PHI)
TOL = 0.05
H_PLANCK = 6.62607015e-34
E_CHARGE = 1.602176634e-19
C_LIGHT = 299792458.0
NU_HI_HZ = 1420.405751786e6

# (label, frequency in Hz, stated uncertainty in Hz, note)
CONSTANTS = [
    ("D, NV- ground-state ZFS", 2.870e9, 1e6, "T-dependent, dD/dT = -74.2 kHz/K"),
    ("14N A_parallel", 2.16e6, 0.02e6, "hyperfine, axial"),
    ("14N A_perp", 2.70e6, 0.05e6, "hyperfine, transverse"),
    ("14N P (quadrupole)", 5.01e6, 0.05e6, "nuclear electric quadrupole"),
    ("15N A_parallel", 3.03e6, 0.05e6, "isotope shift of the same term"),
    ("15N A_perp", 3.65e6, 0.05e6, "isotope shift of the same term"),
    ("13C nearest-neighbour", 130.0e6, 1.0e6, "position-dependent by family"),
    ("NV- ZPL, 637 nm", C_LIGHT / 637.0e-9, 0.3e12, "1.945 eV, 15687 cm^-1"),
    ("NV0 ZPL, 575 nm", C_LIGHT / 575.0e-9, 0.4e12, "2.156 eV, 17389 cm^-1"),
    ("NV- IR line, 1042 nm", C_LIGHT / 1042.0e-9, 5e11, "1.190 eV"),
]

UNITS = [
    ("Hz", lambda v: v, "{:.6e}"),
    ("MHz", lambda v: v / 1e6, "{:.6f}"),
    ("GHz", lambda v: v / 1e9, "{:.6f}"),
    ("eV", lambda v: v * H_PLANCK / E_CHARGE, "{:.9f}"),
    ("cm^-1", lambda v: v / (C_LIGHT * 100.0), "{:.6f}"),
    ("nm", lambda v: C_LIGHT / v * 1e9, "{:.4f}"),
]


def phase_of(value: float) -> float:
    return (math.log(value) / LN_PHI) % 1.0


def dist_of(value: float) -> float:
    return abs(phase_of(value) - round(phase_of(value)))


def binomial_tail(trials: int, observed: int, p: float) -> float:
    return sum(math.comb(trials, k) * p ** k * (1.0 - p) ** (trials - k)
               for k in range(observed, trials + 1))


def section_constants() -> None:
    print("A. THE SPECTROSCOPY (values and what they are uncertain by)")
    print(f"   {'quantity':>26}{'Hz':>14}{'+/-':>12}  note")
    for label, hz, unc, note in CONSTANTS:
        print(f"   {label:>26}{hz:>14.6e}{unc:>12.1e}  {note}")
    print("   one line of D moves 74 kHz per kelvin, which is 2.6e-5 of it;")
    print("   the hyperfine terms are known to about 1% of themselves.")
    print()


def section_units() -> None:
    print("B. THE RUNG ADDRESS IN SIX UNITS: THE SAME PHYSICS, SIX VERDICTS")
    print(f"   tolerance {TOL} rung; 'on' means within {TOL} of an integer rung")
    print(f"   {'quantity':>26}{'unit':>7}{'value':>19}{'phase':>9}{'dist':>8}  on?")
    hits = 0
    tests = 0
    per_unit = {name: 0 for name, _, _ in UNITS}
    for label, hz, _unc, _note in CONSTANTS:
        for name, convert, fmt in UNITS:
            value = convert(hz)
            dist = dist_of(value)
            on = dist <= TOL
            hits += on
            tests += 1
            per_unit[name] += on
            print(f"   {label:>26}{name:>7}"
                  f"{fmt.format(value):>19}{phase_of(value):>9.4f}"
                  f"{dist:>8.4f}  {'yes' if on else 'no'}")
    chance = 2.0 * TOL
    expected = tests * chance
    print(f"   {hits} of {tests} land on a rung; a coin expects {expected:.1f}")
    print(f"   binomial P(>= {hits} of {tests} | {chance:.2f}) = "
          f"{binomial_tail(tests, hits, chance):.3f}")
    print("   per unit: " + ", ".join(f"{name} {count}"
                                     for name, count in per_unit.items()))
    print("   the counts differ unit to unit because log(1e6)/log(phi) = "
          f"{math.log(1e6) / LN_PHI:.4f} shifts")
    print("   every phase by a fixed amount. Nothing about the defect changed;")
    print("   only the ruler did. So 'is this frequency on a rung' is not a")
    print("   property of the NV centre until the framework names its unit.")
    print()


def _directions(limit: int = 4):
    seen = {}
    for x in range(-limit, limit + 1):
        for y in range(-limit, limit + 1):
            for z in range(-limit, limit + 1):
                if (x, y, z) == (0, 0, 0):
                    continue
                norm = math.sqrt(x * x + y * y + z * z)
                key = tuple(round(component / norm, 6) for component in (x, y, z))
                if key not in seen:
                    seen[key] = (x, y, z)
    return list(seen.values())


def section_geometry() -> None:
    print("C. THE LATTICE GEOMETRY: WHERE 137.5078 DEG WOULD HAVE TO COME FROM")
    vectors = _directions(4)
    angles = set()
    for i, a in enumerate(vectors):
        na = math.sqrt(sum(c * c for c in a))
        for b in vectors[i + 1:]:
            nb = math.sqrt(sum(c * c for c in b))
            cosine = sum(p * q for p, q in zip(a, b)) / (na * nb)
            cosine = max(-1.0, min(1.0, cosine))
            angles.add(round(math.degrees(math.acos(cosine)), 4))
    ordered = sorted(angles)
    crowded = min(ordered, key=lambda value: abs(value - GOLDEN_ANGLE_DEG))
    canonical = sorted({round(math.degrees(math.acos(c)), 4) for c in
                        (0.0, 0.5, -0.5, 1 / math.sqrt(2), -1 / math.sqrt(2),
                         1 / math.sqrt(3), -1 / math.sqrt(3), 1 / 3, -1 / 3,
                         math.sqrt(2 / 3), -math.sqrt(2 / 3), 2 / 3, -2 / 3,
                         1.0, -1.0)})
    crystallographic = min(canonical, key=lambda v: abs(v - GOLDEN_ANGLE_DEG))
    print(f"   angles between integer-direction axes, |hkl| <= 4: {len(ordered)} "
          "distinct values")
    print(f"   nearest of those, high-index included: {crowded:.4f} deg, off by "
          f"{abs(crowded - GOLDEN_ANGLE_DEG):.4f} deg - with")
    print(f"   {len(ordered)} exact angles between 0 and 180 deg, some pair lands")
    print("   within hundredths of a degree of any target. That is a statement")
    print("   about how many pairs there are, not about the target.")
    print(f"   nearest of the angles a bond or axis pair actually makes: "
          f"{crystallographic:.4f} deg,")
    print(f"   off by {abs(crystallographic - GOLDEN_ANGLE_DEG):.4f} deg, and that gap "
          "is exact: the")
    print("   crystal's angles are arccosines of 0, +/-1/2, +/-1/sqrt(2), "
          "+/-1/sqrt(3), +/-1/3,")
    print("   +/-sqrt(2/3) and +/-2/3. 135 deg has cos = -1/sqrt(2); the golden")
    print("   angle has no such cosine, so it is not a lattice angle at any index.")
    print(f"   the tetrahedral angle, the NV axis to a bond: "
          f"{math.degrees(math.acos(-1.0 / 3.0)):.4f} deg (cos = -1/3 exactly)")
    print(f"   cos(golden angle) = {math.cos(math.radians(GOLDEN_ANGLE_DEG)):.6f}, "
          "which is arccos of a")
    print("   fifth-degree algebraic number, not of any of those rationals.")
    print("   a periodic lattice also admits no five-fold rotation at all")
    print("   (crystallographic restriction), and phi is the five-fold divergence,")
    print("   so diamond's structure cannot carry it at any index.")
    print()


def section_pitch() -> None:
    print("D. THE PITCH, APPLIED TO THE ONE RATIO IT CAN REACH")
    ratio = 2.870e9 / NU_HI_HZ
    print(f"   D / nu(HI) = {ratio:.6f}")
    forms = [
        ("phi", PHI),
        ("phi^2 - 1", PHI ** 2 - 1.0),
        ("2", 2.0),
        ("sqrt(phi^3)", math.sqrt(PHI ** 3)),
        ("1 + 1/phi^2", 1.0 + 1.0 / PHI ** 2),
    ]
    for name, value in sorted(forms, key=lambda item: abs(item[1] - ratio)):
        print(f"   {name:>12} = {value:.6f}   error {100.0 * (value / ratio - 1.0):+.3f}%")
    print("   the nearest of them misses by one percent and the next by nearly")
    print("   two, and the list was chosen after seeing the number, which is")
    print("   exactly how such coincidences get made. D is a spin-spin coupling of the NV")
    print("   electronic ground state; the hydrogen line is a hyperfine")
    print("   transition of an atom. Nothing couples them, and no baseline T")
    print("   enters, so the pitch has nothing to act on here.")
    print()


def section_tolerance() -> None:
    print("E. WHAT NV COULD AND COULD NOT SETTLE")
    d_hz = 2.870e9
    cone_hz = d_hz * (PHI ** TOL - 1.0)
    linewidth = 2.0e6
    print(f"   the {TOL}-rung tolerance at 2.87 GHz is {cone_hz / 1e6:.1f} MHz wide,")
    print(f"   i.e. any line within +/-{cone_hz / 2e6:.0f} MHz counts as on-lattice.")
    print(f"   a room-temperature ensemble ODMR line is about "
          f"{linewidth / 1e6:.0f} MHz wide, so the instrument pins")
    print(f"   the frequency a factor {cone_hz / linewidth:.0f} finer than the cone. The claim is")
    print("   therefore falsifiable - a line 100 MHz off its rung would be")
    print("   caught - but uninformative: one frequency in ten satisfies it,")
    print("   so passing it costs nothing.")
    tight_hz = d_hz * (PHI ** 1e-3 - 1.0)
    print(f"   at 1e-3 rungs the cone is {tight_hz / 1e6:.2f} MHz, comparable to an")
    print("   ensemble line and a thousand times wider than a single centre")
    print("   resolves, so the tight version is within reach of a bench. There")
    print(f"   the temperature term bites: 74.2 kHz/K needs about "
          f"{tight_hz / 74.2e3:.0f} K of control.")
    print("   the local data do not reach 2.87 GHz at all - every file here tops")
    print("   out at 1530 MHz - and an NV ensemble is a magnetometer reading its")
    print("   own spin splitting, not a radio receiver, so no filterbank in this")
    print("   project can test a claim about NV centres. That needs a bench.")
    print()


def main() -> int:
    print("NV CENTRES IN DIAMOND, READ AGAINST THE RUNG LATTICE AND ITS NULL")
    print()
    section_constants()
    section_units()
    section_geometry()
    section_pitch()
    section_tolerance()
    print("SUMMARY")
    print(f"   address: membership flips between units; the coin rate is "
          f"{2 * TOL:.2f} per test")
    print("   geometry: 137.5078 deg is not an angle of a cubic lattice")
    print("   pitch: D/nu(HI) is near 2, and near-2 was picked after the fact")
    print("   instrument: the cone is ~35x coarser than an ODMR line, so a pass")
    print("               costs 10%; at 1e-3 it is within reach and empty")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
