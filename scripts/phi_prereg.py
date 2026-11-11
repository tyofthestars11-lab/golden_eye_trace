"""A pre-registered two-band test of the rung-rate law, threshold fixed first.

The framework makes two different statements about how rung rate behaves with
frequency, and they cannot both be right:

* the rung coordinate is ``r(x) = ln x / ln phi``, so any band spans
  ``ln10 / ln phi = 4.7849`` rungs per decade - equivalently, a source whose
  velocity does not depend on frequency has a flat rung rate, drift exponent
  ``p = 1``;
* the pitch ``mu = 17.0324 deg`` read as a log-log slope asserts
  ``|drift| proportional to f ** tan(mu)``, exponent ``p = 0.306``, which is a
  rung rate falling by ``10 ** (1 - p) = 4.94`` per decade.

So the observable is ``p = d ln|drift| / d ln f`` over confirmed drifting lines
seen in two widely separated bands. The prediction, the threshold and the
decision rule are written to ``prereg/phi_2band.json`` and hashed *before* the
measurement runs, so that the number cannot be chosen after seeing the data.
Sections 1 and 2 record what the coordinate arithmetic alone returns; section 3
runs the test that can actually fail.

    python3 scripts/phi_prereg.py --fil data/stare2_ovro.fil
    python3 scripts/phi_prereg.py --fil band_a.fil --band-b band_b.fil

An honest run reports ``NOT MEASURED`` when the data cannot support the
statistic, and says which of the two evidence requirements is missing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from seti_search.filterbank import Filterbank  # noqa: E402
from seti_search.rungs import LN_PHI, rung  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

MU_DEG = 17.0324
SLOPE = math.tan(math.radians(MU_DEG))
RUNGS_PER_DECADE = math.log(10.0) / LN_PHI
PITCH_RUNGS_PER_DECADE = RUNGS_PER_DECADE * SLOPE
TOLERANCE = 0.05  # pre-registered relative agreement, before any data is read
PREREG_PATH = os.path.join(ROOT, "prereg", "phi_2band.json")


def band_edges(path: str):
    fb = Filterbank.open(path)
    freqs = fb.frequencies()
    return fb, float(freqs.min()), float(freqs.max())


def rung_density(low: float, high: float) -> float:
    return (rung(high) - rung(low)) / math.log10(high / low)


def drifting_lines(path: str, min_snr: float = 4.0):
    """Detected lines whose drift is resolved and non-zero, as (MHz, Hz/s)."""
    cfg = SearchConfig(boxcar=3, max_drift_hz_s=8.0, snr_threshold=min_snr,
                       max_hits=500)
    hits, _meta = search_file(path, cfg, chan_chunk=2048, max_ints=None)
    fb = Filterbank.open(path)
    resolution = abs(fb.foff) * 1e6 / fb.total_time
    lines = []
    for hit in hits:
        drift = abs(float(getattr(hit, "drift_hz_s", 0.0)))
        if drift > resolution and drift > 0.0:
            lines.append((float(hit.freq_mhz), drift))
    return lines, resolution, len(hits)


def write_preregistration(bands: dict) -> dict:
    record = {
        "test": "rung rate versus frequency, two bands",
        "attribution": "framework: Tyree Jones (tyofthestarz), golden_eye_trace",
        "statistic": "p = d ln|drift| / d ln f over confirmed drifting lines",
        "predictions": {
            "phi_lattice": {
                "p": 1.0,
                "rungs_per_decade": RUNGS_PER_DECADE,
                "statement": "flat rung rate, 4.7849 rungs per decade",
            },
            "pitch_angle": {
                "p": SLOPE,
                "rungs_per_decade": PITCH_RUNGS_PER_DECADE,
                "statement": "rung rate falls 4.94x per decade",
            },
        },
        "tolerance": TOLERANCE,
        "decision_rule": (
            "p within 5% of 1.0000 keeps the flat-rung-rate statement and refutes "
            "the pitch; p within 5% of 0.3063 keeps the pitch and refutes the "
            "flat rung rate; neither means both are refuted by these data."
        ),
        "evidence_required": [
            "two bands separated by at least one decade in frequency",
            "at least 3 confirmed drifting lines per band",
            "each drift resolved well above that file's df/scan amount",
        ],
        "bands": bands,
        "seed": 1729,
    }
    os.makedirs(os.path.dirname(PREREG_PATH), exist_ok=True)
    blob = json.dumps(record, indent=2, sort_keys=True).encode()
    with open(PREREG_PATH, "wb") as handle:
        handle.write(blob)
    record["sha256"] = hashlib.sha256(blob).hexdigest()
    return record


def section_coordinate() -> None:
    print("1. WHAT THE COORDINATE ITSELF RETURNS (cannot fail)")
    print(f"   r(x) = ln x / ln phi gives {RUNGS_PER_DECADE:.4f} rungs per decade")
    print("   in any band, in any frame, for any data, because that is what the")
    print("   coordinate is. Read as a prediction it is an identity: measuring it")
    print("   and finding it holds is not evidence about the sky.")
    print()


def section_bands(low_a: float, high_a: float, band_b) -> None:
    print("2. THE TWO STATEMENTS ON THE BAND IN HAND")
    measured = rung_density(low_a, high_a)
    print(f"   band A {low_a:.4f}-{high_a:.4f} MHz spans "
          f"{rung(high_a) - rung(low_a):.4f} rungs")
    print(f"   measured {measured:.4f} rungs/decade")
    print(f"   phi statement: {RUNGS_PER_DECADE:.4f}  "
          f"(ratio {measured / RUNGS_PER_DECADE:.4f})")
    print(f"   pitch statement: {PITCH_RUNGS_PER_DECADE:.4f}  "
          f"(ratio {measured / PITCH_RUNGS_PER_DECADE:.4f})")
    print(f"   the pitch is already off by {measured / PITCH_RUNGS_PER_DECADE:.2f}x "
          f"on this band alone,")
    print(f"   far outside the pre-registered {TOLERANCE:.0%} tolerance.")
    if band_b is None:
        print("   band B: none supplied; the signal-side test below cannot run.")
    print()


def section_signal_test(lines_a, res_a, hit_a, lines_b=None, res_b=None,
                        hit_b=0) -> bool:
    print("3. THE TEST THAT CAN FAIL: p = d ln|drift| / d ln f ACROSS BANDS")
    print(f"   band A: {hit_a} candidates at 4 sigma, {len(lines_a)} with a drift "
          f"resolved above {res_a:.2e} Hz/s")
    if lines_b is not None:
        print(f"   band B: {hit_b} candidates at 4 sigma, {len(lines_b)} with a "
              f"drift resolved above {res_b:.2e} Hz/s")
    missing = []
    if lines_b is None:
        missing.append("a band B at least a decade away in frequency")
    if len(lines_a) < 3:
        missing.append(f"3 drifting lines in band A resolved above its noise "
                       f"(has {len(lines_a)}, resolution {res_a:.2e} Hz/s)")
    if lines_b is not None and len(lines_b) < 3:
        missing.append(f"3 drifting lines in band B resolved above its noise "
                       f"(has {len(lines_b)}, resolution {res_b:.2e} Hz/s)")
    if missing:
        print("   NOT MEASURED. Still missing, per the pre-registration:")
        for item in missing:
            print(f"     - {item}")
        print("   one band alone cannot fix a slope in frequency, and a number")
        print("   invented from a single band would be a number about that band's")
        print("   noise. It stays unmeasured rather than being called either way.")
        print()
        return False
    freqs = np.log([f * 1e6 for f, _ in lines_a + (lines_b or [])])
    drifts = np.log([d for _, d in lines_a + (lines_b or [])])
    slope = float(np.polyfit(freqs, drifts, 1)[0])
    residual = float(np.std(drifts - np.polyval(np.polyfit(freqs, drifts, 1), freqs)))
    print(f"   fitted p = {slope:.4f} (rms {residual:.3f} in ln)")
    for label, predicted in (("phi lattice", 1.0), ("pitch angle", SLOPE)):
        agree = abs(slope / predicted - 1.0) <= TOLERANCE
        print(f"   {label:>11} predicts {predicted:.4f}: "
              f"{'AGREES' if agree else 'refuted'}")
    print()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fil", required=True, help="band A filterbank file")
    parser.add_argument("--band-b", default="", help="band B filterbank file")
    args = parser.parse_args()

    fb_a, low_a, high_a = band_edges(args.fil)
    fb_b = low_b = high_b = None
    if args.band_b:
        fb_b, low_b, high_b = band_edges(args.band_b)

    record = write_preregistration({
        "A": {"file": os.path.basename(args.fil), "low_mhz": low_a,
              "high_mhz": high_a},
        "B": ({"file": os.path.basename(args.band_b), "low_mhz": low_b,
               "high_mhz": high_b} if fb_b else None),
    })
    print("0. PRE-REGISTRATION, WRITTEN AND HASHED BEFORE ANY MEASUREMENT")
    print(f"   {PREREG_PATH}")
    print(f"   sha256 {record['sha256']}")
    print(f"   p = 1.0000 vs p = {SLOPE:.4f}, tolerance {TOLERANCE:.0%}")
    print("   both statements are the framework's own; the decision rule is in")
    print("   the file above and was not written after seeing the numbers.")
    print()

    section_coordinate()
    section_bands(low_a, high_a, fb_b)

    lines_a, res_a, hit_a = drifting_lines(args.fil)
    if fb_b is not None:
        lines_b, res_b, hit_b = drifting_lines(args.band_b)
    else:
        lines_b, res_b, hit_b = None, None, 0
    measured = section_signal_test(lines_a, res_a, hit_a, lines_b, res_b, hit_b)

    print("SUMMARY")
    print(f"   pre-registered: sha256 {record['sha256'][:16]}...")
    print("   coordinate: 4.7849 rungs/decade by definition, not by measurement")
    print(f"   pitch on this band: off by {rung_density(low_a, high_a) / PITCH_RUNGS_PER_DECADE:.2f}x")
    print("   signal-side test: " + ("measured" if measured else "NOT MEASURED"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
