"""The tightened claim: one named unit, 1e-3 rungs. Testable, and true?

Tightening the rung tolerance is the only way the framework becomes a statement
about an instrument instead of an arithmetic exercise. Two numbers decide it:

    content(t) = 2t                 a uniform phase passes with probability 2t
    checkable(t) <=> t > tau        tau = ln(1 + error/f) / ln(phi)

The second is the one the framework never quotes. A claim says a frequency sits
inside a cone of width f(phi^t - 1). A measurement with error delta can refute
that only while the cone is wider than delta; once the cone is narrower than the
error bar, every honest measurement lands inside it and the statement cannot be
false. So the claim is checkable at the tolerances above tau, and tau is a floor
on how tight anyone can test it - which is why a 1%-known frequency cannot be
tested below 0.019 rungs, and an exactly known one can be tested at any
tolerance above zero.

Content and checkability therefore pull in opposite directions: smaller t means
a pass is rarer (more informative) but demands a better measurement. The
interesting region is tau < t << 1, and that is where this script looks, with
the unit fixed to Hz.

    python3 scripts/phi_tight.py
"""

from __future__ import annotations

import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.rungs import LN_PHI, PHI  # noqa: E402

TOL = 1e-3
FRAMEWORK_TOL = 0.05

# (label, frequency in Hz, measurement error in Hz, note)
ENTRIES = [
    ("Cs 133 hyperfine", 9192631770.0, 0.0, "the SI second, exact"),
    ("H 21 cm line", 1420405751.77, 1.0, "about 1e-9 relative"),
    ("Rb 87 hyperfine", 6834682610.904, 1.0, "about 1 Hz"),
    ("NV- D", 2870000000.0, 1.0e3, "single-centre ODMR, kHz"),
    ("14N A_parallel", 2.16e6, 2.0e4, "1% of itself"),
    ("13C nearest neighbour", 130.0e6, 1.0e6, "family spread"),
]


def phase(value: float) -> float:
    return (math.log(value) / LN_PHI) % 1.0


def distance(value: float) -> float:
    return abs(phase(value) - round(phase(value)))


def tau(hz: float, err: float) -> float:
    """Smallest tolerance at which this measurement could refute an on-rung value."""
    if err <= 0.0:
        return 0.0
    return math.log1p(err / hz) / LN_PHI


def section_structure() -> None:
    print("A. THE TWO NUMBERS, AND WHY THE FRAMEWORK QUOTES ONLY ONE")
    print(f"   content   = 2t = {2 * TOL:g} at t = {TOL:g} rungs")
    print("             (uniform phase in a band: 0.2% pass, so a hit here")
    print("              would be a 1-in-500 event instead of a 1-in-10)")
    print("   threshold = tau = ln(1 + error/f) / ln(phi)")
    print("             a claim is checkable only above tau: below it the")
    print("             cone is narrower than the error bar and no measurement")
    print("             can land outside the cone.")
    print("   so a 1%-known frequency cannot even be tested below ~0.019 rungs,")
    print("   while an exactly known one can be tested arbitrarily tight.")
    print("   content falls as t falls; checkability needs t > tau. Both at")
    print("   once means tau < t << 1, which is the only region worth quoting.")
    print()


def section_table(tol: float, title: str) -> None:
    print(title)
    print(f"   {'quantity':>24}{'Hz':>17}{'+/- Hz':>11}{'tau':>11}"
          f"{'phase':>8}{'dist':>9}{'pass':>6}{'checkable':>11}")
    passes = 0
    checkable = 0
    for label, hz, err, _note in ENTRIES:
        d = distance(hz)
        on = d <= tol
        limit = tau(hz, err)
        can = limit < tol
        passes += on
        checkable += can
        print(f"   {label:>24}{hz:>17.6e}{err:>11.1e}{limit:>11.2e}"
              f"{phase(hz):>8.4f}{d:>9.4f}{'yes' if on else 'no':>6}"
              f"{'yes' if can else 'no':>11}")
    chance = 2.0 * tol
    print(f"   {passes} of {len(ENTRIES)} pass, where uniform phase gives "
          f"{chance:.4f} each ({len(ENTRIES) * chance:.4f} expected); "
          f"{checkable} of {len(ENTRIES)} are checkable at t = {tol:g}")
    print()


def section_verdict() -> None:
    print("C. THE WINDOW")
    least = max(tau(hz, err) for _l, hz, err, _n in ENTRIES)
    print(f"   the coarsest entry in this table cannot be tested below "
          f"{least:.4f} rungs;")
    print(f"   the framework asks for {FRAMEWORK_TOL}, which is "
          f"{FRAMEWORK_TOL / least:.1f}x coarser than that floor - so at 0.05")
    print(f"   everything here is checkable, and every pass is a "
          f"{100 * 2 * FRAMEWORK_TOL:g}% event.")
    print("   tightening to 1e-3 keeps 4 of the 6 checkable and drops the pass")
    print("   rate to 0.2%; there, nothing passes. Tightening further keeps")
    print("   only the exact and the 1e-9 entries, and still nothing passes.")
    print("   an informative pass rate is available at every tolerance this")
    print("   table can test. What is missing is a pass, not a measurement.")
    print()


def main() -> int:
    print("THE TIGHTENED CLAIM: UNIT = Hz, TOLERANCE = 1e-3 RUNGS")
    print()
    section_structure()
    section_table(TOL, "B. AT THE TIGHTENED TOLERANCE")
    section_table(FRAMEWORK_TOL, "B'. THE SAME ENTRIES AT THE FRAMEWORK'S 0.05")
    section_verdict()
    print("SUMMARY")
    print(f"   unit fixed to Hz, tolerance {TOL:g}: "
          f"{sum(distance(h) <= TOL for _l, h, _e, _n in ENTRIES)} of "
          f"{len(ENTRIES)} pass, "
          f"{sum(tau(h, e) < TOL for _l, h, e, _n in ENTRIES)} of "
          f"{len(ENTRIES)} checkable")
    print(f"   at the framework's {FRAMEWORK_TOL}, all six are checkable and a")
    print("   pass is worth 10%: checkable, and uninformative")
    print("   verdict: empty at the tolerance where a pass would mean something")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
