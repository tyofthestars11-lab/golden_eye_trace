"""Do the address arithmetic in Z[phi], not in degrees.

The golden angle is ``phi**-2`` of a turn, and ``phi**-2 = 2 - phi`` is an
algebraic integer, so the lattice of golden-angle multiples is ``n*(2-phi)``: an
element of Z[phi], written here as an integer pair ``(a, b)`` meaning
``a + b*phi``. Degrees appear only when printing.

Distances are measured two ways, and the second is the one that matters:

* additive - ``|u - n*s|`` where ``u`` is the address and ``s = 2 - phi``;
* multiplicative - ``|log_phi(u / (n*s))|``, the phi-natural distance between two
  addresses, since addresses compose by multiplication.

The null has to be computed in the same metric as the distance, or a change of
metric looks like a change of physics. For a random address uniform over a lattice
spacing, the additive chance is ``2d/s``; the multiplicative one is
``2*d*ln(phi)*u/s``, and the extra factor ``u`` is exactly what cancels the
shrinking that a logarithmic distance gives at large addresses. With it, the two
readings agree - which is the useful thing this file establishes, because the
multiplicative number alone looks three times more interesting than it is.

    python3 scripts/phi_metric.py
"""

from __future__ import annotations

import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.rungs import LN_PHI, PHI, rung  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

#: one golden angle, in turns, as an exact element of Z[phi]: 2 - phi
GOLDEN_TURN = (2, -1)
GOLDEN_DEG = 360.0 * (2 - PHI)

TARGETS = [
    ("B0329+54 (LWA beam 3)", 53.2472916667, 54.5786666667),
    ("Sgr A*", 266.41683, -29.00781),
]


def fib(n: int) -> int:
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def phi_pow(n: int) -> tuple[int, int]:
    """phi**n as (a, b) with a + b*phi, for positive and negative n."""
    if n == 0:
        return (1, 0)
    if n > 0:
        return (fib(n), fib(n - 1))
    m = -n
    sign = -1 if m % 2 else 1
    return (sign * fib(m + 1), -sign * fib(m))


def to_float(pair: tuple[int, int]) -> float:
    return pair[0] + pair[1] * PHI


def fmt(pair: tuple[int, int]) -> str:
    a, b = pair
    if b == 0:
        return f"{a}"
    coeff = "" if abs(b) == 1 else str(abs(b))
    return f"{a} {'+' if b > 0 else '-'} {coeff}phi"


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


def measure(address: float, step: float) -> dict:
    """Nearest lattice point, both distances, and each metric's own chance."""
    index = int(round(address / step))
    target = index * step
    additive = abs(address - target)
    p_add = min(1.0, 2.0 * additive / step)
    if index == 0:
        return {"index": 0, "additive": additive, "p_add": p_add,
                "multiplicative": None, "p_mul": 1.0}
    multiplicative = abs(math.log(address / target)) / LN_PHI
    # the null in the multiplicative metric carries the address scale; without
    # it the logarithmic distance is compared against an additive null and the
    # result is three times too small
    p_mul = min(1.0, 2.0 * multiplicative * LN_PHI * address / step)
    return {"index": index, "additive": additive, "p_add": p_add,
            "multiplicative": multiplicative, "p_mul": p_mul}


def main() -> int:
    step = to_float(GOLDEN_TURN)
    print("IN Z[phi]: THE GOLDEN ANGLE AS AN ALGEBRAIC INTEGER")
    print(f"   phi**-2                    = {fmt(GOLDEN_TURN)} = {step:.12f} turns")
    print(f"   one turn in golden angles  = phi**2 = {fmt(phi_pow(2))} "
          f"= {360.0 / GOLDEN_DEG:.12f}")
    print(f"   in degrees                   {GOLDEN_DEG:.10f}")
    print("   the lattice is n*(2-phi) turns, and no integer n closes the circle,")
    print("   because 1/(2-phi) = phi**2 is irrational")
    print()

    print("A. SKY ADDRESSES AS MULTIPLES OF THE GOLDEN ANGLE")
    chances = []
    for name, ra, dec in TARGETS:
        for label, lon in (("RA", ra),
                           ("ecliptic", ecliptic_lon(ra, dec)),
                           ("galactic", galactic_lon(ra, dec))):
            address = lon / 360.0
            m = measure(address, step)
            chances.append(m["p_add"])
            mult = ("-" if m["multiplicative"] is None
                    else f"{m['multiplicative']:.4f} rungs (p={m['p_mul']:.3f})")
            print(f"   {name:<22} {label:<9} {lon:8.4f} deg = {address:.6f} turns "
                  f"= {address / step:.5f} turns of the angle")
            print(f"      nearest {m['index']} steps; additive {m['additive']:.6f} "
                  f"turns (p={m['p_add']:.3f}); multiplicative {mult}")
    best = min(chances)
    look = 1.0 - (1.0 - best) ** len(chances)
    print(f"   {len(chances)} tests; best additive chance {best:.3f}, "
          f"look-elsewhere over them {look:.3f}")
    print()

    print("B. BAND SEPARATION IN THE SAME LATTICE")
    band_a, band_b = 1405.122, 50.390
    separation = rung(band_a) - rung(band_b)
    m = measure(separation, step)
    exact = m["index"] * step
    print(f"   {band_a:.3f} MHz / {band_b:.3f} MHz = {separation:.4f} rungs")
    print(f"   nearest lattice point {m['index']} steps = {exact:.4f} rungs "
          f"= {fmt((m['index'] * GOLDEN_TURN[0], m['index'] * GOLDEN_TURN[1]))} rungs")
    print(f"   additive shortfall {m['additive']:.4f} rungs  p={m['p_add']:.3f}")
    print(f"   multiplicative     {m['multiplicative']:.4f} rungs  "
          f"p={m['p_mul']:.3f}")
    print(f"   the multiplicative distance is {m['multiplicative'] / m['additive']:.2f}x "
          f"the additive one,")
    print(f"   and its chance is {m['p_mul'] / m['p_add']:.2f}x the additive chance:"
          f" the same factor, so the metric is not doing the work")
    print()

    path = os.path.join(ROOT, "data", "stare2_ovro.fil")
    if os.path.exists(path):
        cfg = SearchConfig(boxcar=3, max_drift_hz_s=4.0, snr_threshold=3.0,
                           max_hits=20)
        hits, _meta = search_file(path, cfg, chan_chunk=2048, max_ints=None)
        freqs = sorted(h.freq_mhz for h in hits)
        print("C. CANDIDATE GAPS, SAME LATTICE")
        print(f"   {'MHz':>10}{'rungs':>9}{'steps':>8}{'additive':>10}"
              f"{'p_add':>7}{'p_mul':>7}")
        for lo, hi in zip(freqs[:-1], freqs[1:]):
            gap = rung(hi) - rung(lo)
            m = measure(gap, step)
            print(f"   {hi - lo:>10.3f}{gap:>9.4f}{gap / step:>8.4f}"
                  f"{m['additive']:>10.4f}{m['p_add']:>7.3f}{m['p_mul']:>7.3f}")
        print("   these chances are what random pairs give; the gaps carry no")
        print("   lattice information beyond their own arithmetic")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
