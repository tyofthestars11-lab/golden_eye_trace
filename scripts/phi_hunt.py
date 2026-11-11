"""Hunt at the tightened tolerance: any channel within 1e-3 rungs of a rung?

At 0.05 rungs a pass is a 10% event, so counting passes says nothing. At 1e-3
rungs a pass is a 0.2% event, which is small enough that a handful of them
would mean something - so this script looks at every channel in the band, not
just the candidates, and reports the count against the count uniform phase
predicts.

It also reports the scale of the thing being tested. At 1.4 GHz the 1e-3-rung
cone is 683 kHz against a 122 kHz channel, so the file's own channel grid is
5.6 times finer than the tolerance: the test is not being decided by where the
channels happen to fall.

    python3 scripts/phi_hunt.py data/stare2_ovro.fil data/demo.fil
"""

from __future__ import annotations

import argparse
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from seti_search.filterbank import Filterbank  # noqa: E402
from seti_search.rungs import LN_PHI, PHI  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

TIGHT = 1e-3
LOOSE = 0.05


def distances(hz: np.ndarray) -> np.ndarray:
    values = np.asarray(hz, dtype=float)
    grid = np.log(values) / LN_PHI
    return np.abs(grid - np.round(grid))


def coverage(low_hz: float, high_hz: float, tol: float) -> float:
    """Fraction of this band that lies inside a tolerance cone, by geometry.

    The channels are evenly spaced in frequency, so the fraction of channels
    inside the cones is the fraction of the band inside them. Unlike N*2t this
    is correct for a band narrower than a rung, where most placements never
    come near a rung at all.
    """
    lo = math.log(low_hz) / LN_PHI
    hi = math.log(high_hz) / LN_PHI
    total = 0.0
    for n in range(math.floor(lo - tol), math.ceil(hi + tol) + 1):
        lower = max(low_hz, PHI ** (n - tol))
        upper = min(high_hz, PHI ** (n + tol))
        if upper > lower:
            total += upper - lower
    return total / (high_hz - low_hz)


def poisson_tail(observed: int, expected: float) -> float:
    if expected <= 0.0:
        return 0.0
    total = 0.0
    term = math.exp(-expected)
    for k in range(0, observed):
        total += term
        term *= expected / (k + 1)
    return max(0.0, min(1.0, 1.0 - total))


def hunt(path: str) -> None:
    fb = Filterbank.open(path)
    freqs_mhz = fb.frequencies()
    hz = freqs_mhz.astype(float) * 1e6
    low, high = float(hz.min()), float(hz.max())
    dist = distances(hz)

    cone_tight = low * (PHI ** TIGHT - 1.0)
    cone_loose = low * (PHI ** LOOSE - 1.0)
    channel = abs(fb.foff) * 1e6

    print(f"   {os.path.basename(path)}: {fb.nchans} channels, "
          f"{low / 1e6:.3f}-{high / 1e6:.3f} MHz, channel {channel / 1e3:.3f} kHz")
    print(f"   the {TIGHT:g}-rung cone is {cone_tight / 1e3:.1f} kHz wide, "
          f"{cone_tight / channel:.1f} channels;")
    print(f"   the {LOOSE:g}-rung cone is {cone_loose / 1e3:.1f} kHz, "
          f"{cone_loose / channel:.1f} channels")

    span = math.log(high / low) / LN_PHI
    print(f"   the band spans {span:.4f} rungs, so a whole-rung null (N * 2t)")
    print("   would overstate what can be found here; the null used below is the")
    print("   fraction of the band the cones actually cover.")
    for tol, label in ((TIGHT, "tightened"), (LOOSE, "framework")):
        inside = int(np.sum(dist <= tol))
        cover = coverage(low, high, tol)
        expected = len(dist) * cover
        p_value = poisson_tail(inside, expected) if expected > 0 else 0.0
        print(f"   {label}: {inside} of {len(dist)} channels pass; the cones")
        print(f"      cover {100 * cover:.2f}% of the band, so {expected:.3f} were")
        if cover == 0.0:
            print("      possible: the band never reaches a rung, so there is")
            print("      nothing here to pass or fail")
        else:
            print(f"      possible, P(>= {inside}) = {p_value:.3f}")

    cfg = SearchConfig(boxcar=3, max_drift_hz_s=4.0, snr_threshold=3.0, max_hits=200)
    hits, _meta = search_file(path, cfg, chan_chunk=2048, max_ints=None)
    print(f"   {len(hits)} candidates at {cfg.snr_threshold:g} sigma:")
    print(f"   {'MHz':>14}{'SNR':>7}{'dist':>10}{'tight':>8}{'loose':>8}")
    tight_hits = 0
    for hit in sorted(hits, key=lambda h: -h.snr)[:10]:
        d = float(distances(np.array([hit.freq_mhz * 1e6]))[0])
        on = d <= TIGHT
        tight_hits += on
        print(f"   {hit.freq_mhz:>14.6f}{hit.snr:>7.1f}{d:>10.5f}"
              f"{'yes' if on else 'no':>8}{'yes' if d <= LOOSE else 'no':>8}")
    cover = coverage(low, high, TIGHT)
    expected_hits = len(hits) * cover
    print(f"   {tight_hits} of {len(hits)} candidates pass the tightened test; "
          f"{expected_hits:.4f} were possible"
          + (" (no channel is within reach)" if cover == 0.0 else ""))
    print()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fils", nargs="*",
                        default=[os.path.join(ROOT, "data", "stare2_ovro.fil")])
    args = parser.parse_args()

    print("HUNT AT 1e-3 RUNGS: EVERY CHANNEL, AGAINST THE BAND'S OWN RATE")
    print()
    for path in args.fils:
        hunt(path)
    print("SUMMARY")
    print(f"   at {TIGHT:g} rungs a pass costs about {2 * TIGHT:.1e} per channel")
    print("   where the band covers the cones, and nothing where it does not.")
    print("   the wide band above misses the nearest rung by 0.051 rungs, so")
    print("   even the framework's 0.05 cone finds no channel to test; the")
    print("   narrow bands sit inside a single gap end to end. Neither is empty")
    print("   evidence about the lattice: both are files that never meet it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
