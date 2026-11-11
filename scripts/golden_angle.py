"""Apply the golden angle where the framework puts it, and give each use a null.

The angle is exact, not approximate: ``360 / 137.50776405003785... = phi**2``, so
a turn is exactly phi^2 golden angles and the constant is fixed by phi^2 = phi+1.
Everything below uses that exact value, and every number is reported next to the
chance rate for a random direction or a random frequency.

Three uses, because the framework makes three different claims with it:

A.  as a sky address - the angle is a rotation, so a source at longitude
    ``137.5078 deg`` is a special place. Galactic and ecliptic longitude are
    both tested, on directions whose coordinates are known exactly: the LWA beam
    (B0329+54, from the observation's own .inf) and Sgr A*.
B.  as a rung step - the angle is ``1/phi**2 = 0.381966`` of a turn, so one
    golden-angle step in the rung coordinate is ``0.381966`` rungs, i.e. a
    frequency ratio of ``phi ** 0.381966 = 1.20182``. Band separations and
    candidate gaps are tested against multiples of that step.
C.  as geometry - ``cos(137.5078 deg) = -0.737369``, and because 360/GOLDEN is
    irrational the angle closes no circle: no integer count of these steps
    returns to the start. That is the structural reason it tiles a disc the way
    phyllotaxis does.

    python3 scripts/golden_angle.py
"""

from __future__ import annotations

import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import numpy as np  # noqa: E402

from seti_search.filterbank import Filterbank  # noqa: E402
from seti_search.rungs import GOLDEN_ANGLE_DEG, LN_PHI, PHI, rung  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

GOLDEN = GOLDEN_ANGLE_DEG
TURN_RATIO = 360.0 / GOLDEN
STEP_RUNGS = GOLDEN / 360.0
STEP_RATIO = PHI ** STEP_RUNGS

TARGETS = [
    ("B0329+54 (LWA beam 3)", 53.2472916667, 54.5786666667),
    ("Sgr A*", 266.41683, -29.00781),
]


def ecliptic_lon(ra_deg: float, dec_deg: float, obl_deg: float = 23.4392911) -> float:
    a, d, e = math.radians(ra_deg), math.radians(dec_deg), math.radians(obl_deg)
    lam = math.atan2(math.sin(a) * math.cos(e) + math.tan(d) * math.sin(e),
                     math.cos(a))
    return math.degrees(lam) % 360.0


def galactic_lon(ra_deg: float, dec_deg: float) -> float:
    a_ngp, d_ngp, l_ncp = 192.85948, 27.12825, 122.93192
    a, d = math.radians(ra_deg), math.radians(dec_deg)
    an, dn = math.radians(a_ngp), math.radians(d_ngp)
    y = math.cos(d) * math.sin(a - an)
    x = math.sin(d) * math.cos(dn) - math.cos(d) * math.sin(dn) * math.cos(a - an)
    return (l_ncp - math.degrees(math.atan2(y, x))) % 360.0


def signed_offset(lon_deg: float) -> float:
    return (lon_deg - GOLDEN + 180.0) % 360.0 - 180.0


def chance_rate(delta_deg: float) -> float:
    """Chance a random longitude falls this close to the angle."""
    return min(1.0, 2.0 * abs(delta_deg) / 360.0)


def main() -> int:
    print("THE GOLDEN ANGLE, APPLIED")
    print(f"   golden angle            {GOLDEN:.10f} deg")
    print(f"   360 / angle             {TURN_RATIO:.10f}  (= phi^2 = {PHI ** 2:.10f})")
    print(f"   exactly additive        {abs(TURN_RATIO - PHI ** 2) < 1e-12}")
    print(f"   cos(angle)              {math.cos(math.radians(GOLDEN)):.10f}")
    print()

    print("A. AS A SKY ADDRESS")
    for name, ra, dec in TARGETS:
        lam = ecliptic_lon(ra, dec)
        lon = galactic_lon(ra, dec)
        d_lam = signed_offset(lam)
        d_lon = signed_offset(lon)
        ra_off = signed_offset(ra)
        print(f"   {name}")
        print(f"      RA {ra:9.4f}  offset from the angle {ra_off:+9.3f} deg"
              f"  (chance {chance_rate(ra_off):.3f})")
        print(f"      ecliptic longitude {lam:9.4f}  offset {d_lam:+9.3f} deg"
              f"  (chance {chance_rate(d_lam):.3f})")
        print(f"      galactic longitude {lon:9.4f}  offset {d_lon:+9.3f} deg"
              f"  (chance {chance_rate(d_lon):.3f})")
    print()

    print("B. AS A RUNG STEP")
    print(f"   one step = {STEP_RUNGS:.6f} rungs, a ratio of {STEP_RATIO:.6f}")
    band_a = 1405.122  # STARE2 OVRO centre
    band_b = 50.390    # LWA DRX dynamic spectrum centre
    gap_rungs = rung(band_a) - rung(band_b)
    steps = gap_rungs / STEP_RUNGS
    nearest = round(steps)
    print(f"   band separation {gap_rungs:.4f} rungs = {steps:.4f} steps;")
    print(f"   nearest integer {nearest} steps = {nearest * STEP_RUNGS:.4f} rungs, "
          f"off by {abs(steps - nearest) * STEP_RUNGS:.4f} rungs")
    print(f"   chance of landing that close: "
          f"{2 * 0.05 / STEP_RUNGS:.3f} at a 0.05-rung tolerance "
          f"(the lattice spacing is {STEP_RUNGS:.3f} rungs)")

    path = os.path.join(ROOT, "data", "stare2_ovro.fil")
    if os.path.exists(path):
        cfg = SearchConfig(boxcar=3, max_drift_hz_s=4.0, snr_threshold=3.0,
                           max_hits=20)
        hits, _meta = search_file(path, cfg, chan_chunk=2048, max_ints=None)
        freqs = sorted(h.freq_mhz for h in hits)
        print(f"   candidate gaps, {len(freqs)} candidates in "
              f"{os.path.basename(path)}:")
        for lo, hi in zip(freqs[:-1], freqs[1:]):
            step = (rung(hi) - rung(lo)) / STEP_RUNGS
            print(f"      {hi - lo:9.3f} MHz = {rung(hi) - rung(lo):.4f} rungs "
                  f"= {step:.4f} steps, off {abs(step - round(step)) * STEP_RUNGS:.4f} rungs")
    print()

    print("C. AS GEOMETRY")
    print("   360/angle = phi^2 is irrational, so no integer number of golden-angle")
    print("   steps closes the circle: the sequence never repeats, which is what")
    print("   makes the phyllotaxis packing work. The angle is a real constant;")
    print("   whether a given source or band sits on it is the question the")
    print("   numbers above answer, and each answer carries its chance rate.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
