"""Run the pre-registered two-band drift-exponent test, or say exactly which of
its evidence requirements is missing.

The statistic is ``p = d ln|drift| / d ln f`` over confirmed drifting lines in
two bands. Two predictions compete: a flat rung rate, ``p = 1`` (4.7849 rungs
per decade), and the pitch read as a slope, ``p = 0.3063`` (1.4659 rungs per
decade). The decision rule and the threshold live in
``prereg/phi_2band.json``, written and hashed before any of this runs.

This script exists so that the requirement cannot be quietly relaxed: no second
band, no resolved drifts, or fewer than three lines in either band all produce
``NOT MEASURED`` with the missing quantity printed, never a slope.

    python3 scripts/phi_twoband.py --band-a data/stare2_ovro.fil
    python3 scripts/phi_twoband.py --band-a a.fil --band-b b.fil --min-snr 5
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

from bands import (  # noqa: E402
    MIN_DECADES,
    PITCH_RUNGS_PER_DECADE,
    REFERENCE_MHZ,
    RUNGS_PER_DECADE,
    SEPARATION,
    decades,
    scan_time_for,
    slope_sigma_ceiling,
)
from seti_search.filterbank import Filterbank  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

MIN_LINES = 3
PREREG_PATH = os.path.join(ROOT, "prereg", "phi_2band.json")


def band_summary(path: str, min_snr: float, max_ints: int | None) -> dict:
    """Band edges, drift resolution, and the drifting lines found in it."""
    fb = Filterbank.open(path)
    freqs = fb.frequencies()
    channel_hz = abs(fb.foff) * 1e6
    total_time = fb.total_time if math.isfinite(fb.tsamp) and fb.tsamp > 0 else 0.0
    resolution = channel_hz / total_time if total_time > 0 else float("inf")

    base = {
        "path": path,
        "low_mhz": float(freqs.min()),
        "high_mhz": float(freqs.max()),
        "centre_mhz": float(0.5 * (freqs.min() + freqs.max())),
        "channel_khz": channel_hz / 1e3,
        "total_time_s": total_time,
        "resolution_hz_s": resolution,
        "nchans": int(fb.nchans),
    }
    if fb.nchans < 2:
        # a dedispersed time series has no frequency axis at all: whatever it
        # holds, it cannot carry a measurement of drift against frequency
        return {**base, "hits": 0, "lines": [], "truncated": False,
                "n_ints_used": 0, "single_channel": True}

    cfg = SearchConfig(boxcar=3, max_drift_hz_s=8.0, snr_threshold=min_snr,
                       max_hits=500)
    hits, meta = search_file(path, cfg, chan_chunk=2048, max_ints=max_ints)

    lines = [
        (float(hit.freq_mhz), abs(float(hit.drift_hz_s)), float(hit.snr))
        for hit in hits
        if abs(float(hit.drift_hz_s)) > 2.0 * resolution
    ]
    return {
        **base,
        "hits": len(hits),
        "lines": lines,
        "truncated": bool(meta.get("truncated")),
        "n_ints_used": int(meta.get("n_ints", 0)),
        "single_channel": False,
    }


def bootstrap_slope(lines: list[tuple[float, float]], draws: int = 2000) -> tuple[float, float]:
    """Slope of ln|drift| on ln f, with a 1-sigma bootstrap interval."""
    freq = np.log(np.array([f for f, _ in lines], dtype=float))
    drift = np.log(np.array([d for _, d in lines], dtype=float))

    def slope(idx):
        x, y = freq[idx], drift[idx]
        if np.ptp(x) == 0:
            return float("nan")
        return float(np.polyfit(x, y, 1)[0])

    point = slope(np.arange(freq.size))
    rng = np.random.default_rng(1729)
    values = []
    for _ in range(draws):
        idx = rng.integers(0, freq.size, freq.size)
        candidate = slope(idx)
        if math.isfinite(candidate):
            values.append(candidate)
    if not values:
        return point, float("nan")
    return point, float(np.std(values))


def print_band(label: str, info: dict) -> None:
    print(f"   {label}  {os.path.basename(info['path'])}")
    print(f"      {info['low_mhz']:.3f}-{info['high_mhz']:.3f} MHz, "
          f"channel {info['channel_khz']:.4g} kHz")
    if info["total_time_s"]:
        print(f"      {info['total_time_s']:.1f} s, drift resolution "
              f"{info['resolution_hz_s']:.3g} Hz/s")
    else:
        print("      header has no tsamp: drift resolution unknown")
    if info.get("single_channel"):
        print("      single channel: no frequency axis, so no drift can be")
        print("      measured against frequency in this file")
        return
    print(f"      {info['hits']} candidates at threshold, "
          f"{len(info['lines'])} with resolved non-zero drift"
          + (f" (first {info['n_ints_used']} integrations)" if info["truncated"] else ""))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--band-a", required=True)
    parser.add_argument("--band-b")
    parser.add_argument("--min-snr", type=float, default=5.0)
    parser.add_argument("--max-ints", type=int, default=4096,
                        help="integrations per band; 0 means all")
    parser.add_argument("--min-decades", type=float, default=MIN_DECADES)
    parser.add_argument("--min-lines", type=int, default=MIN_LINES)
    args = parser.parse_args()

    max_ints = None if args.max_ints == 0 else args.max_ints

    print("PRE-REGISTERED TWO-BAND DRIFT-EXPONENT TEST")
    if os.path.exists(PREREG_PATH):
        with open(PREREG_PATH, "rb") as fh:
            digest = hashlib.sha256(fh.read()).hexdigest()
        print(f"   registration sha256 {digest[:32]}...")
    else:
        print("   WARNING: prereg/phi_2band.json is missing; run phi_prereg.py")
    print()

    a = band_summary(args.band_a, args.min_snr, max_ints)
    print_band("A", a)
    print()

    missing = []
    if not args.band_b:
        missing.append("a second band (--band-b)")
        span = 0.0
    else:
        b = band_summary(args.band_b, args.min_snr, max_ints)
        print_band("B", b)
        print()
        span = decades(a["centre_mhz"], b["centre_mhz"])
        print(f"   separation {span:.2f} decades "
              f"(requirement {args.min_decades:g})")
        if span < args.min_decades:
            need_low = a["centre_mhz"] / 10 ** args.min_decades
            need_high = a["centre_mhz"] * 10 ** args.min_decades
            missing.append(
                f"{args.min_decades - span:.2f} more decades "
                f"(needs <= {need_low:.0f} MHz or >= {need_high:.0f} MHz)"
            )
        for label, info in (("A", a), ("B", b)):
            if len(info["lines"]) < args.min_lines:
                missing.append(
                    f"band {label} has {len(info['lines'])} drifting lines, "
                    f"needs {args.min_lines}"
                )

    print()
    print(f"   the two predictions are {SEPARATION:.4f} rungs/decade apart, so at "
          f"{span:.2f} decades")
    print(f"   a slope known to {slope_sigma_ceiling(max(span, 1e-9)):.3f} "
          "rungs/decade separates them at 3 sigma")
    print()

    if missing:
        print("NOT MEASURED")
        for item in missing:
            print(f"   missing: {item}")
        if a["total_time_s"] and a["resolution_hz_s"] > 1.0:
            need = scan_time_for(10.0, a["channel_khz"] * 1e3)
            print(f"   band A resolves {a['resolution_hz_s']:.3g} Hz/s: at its "
                  f"{a['channel_khz']:.4g} kHz channels a 10 Hz/s drift needs")
            print(f"   a {need:,.0f} s scan ({need / 3600:.1f} h) and it has "
                  f"{a['total_time_s']:.3g} s, so lengthening the cadence alone")
            print("   cannot make band A usable; the channel width is the limit")
        print("   no slope is quoted: at this separation it would describe the data")
        print("   that produced it, which is the failure mode the registration exists")
        print("   to prevent.")
        return 0

    lines = a["lines"] + b["lines"]
    slope_ln, sigma_ln = bootstrap_slope(lines)
    slope = slope_ln * RUNGS_PER_DECADE
    sigma = sigma_ln * RUNGS_PER_DECADE

    print(f"   {len(lines)} lines, slope {slope:.4f} +/- {sigma:.4f} rungs/decade")
    print(f"   flat rung rate      {RUNGS_PER_DECADE:.4f}  "
          f"({abs(slope - RUNGS_PER_DECADE) / sigma:.2f} sigma away)")
    print(f"   pitch as a slope    {PITCH_RUNGS_PER_DECADE:.4f}  "
          f"({abs(slope - PITCH_RUNGS_PER_DECADE) / sigma:.2f} sigma away)")
    print()
    flat_ok = abs(slope - RUNGS_PER_DECADE) <= 0.05 * RUNGS_PER_DECADE
    pitch_ok = abs(slope - PITCH_RUNGS_PER_DECADE) <= 0.05 * PITCH_RUNGS_PER_DECADE
    if flat_ok and not pitch_ok:
        print("MEASURED: the flat rung rate survives the pre-registered 5% rule;")
        print("   the pitch does not.")
    elif pitch_ok and not flat_ok:
        print("MEASURED: the pitch survives the pre-registered 5% rule; the flat")
        print("   rate does not.")
    elif flat_ok and pitch_ok:
        print("MEASURED: both survive at this separation, so the test did not")
        print("   discriminate; a wider lever arm is needed.")
    else:
        print("MEASURED: neither statement survives; report the slope as it stands.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
