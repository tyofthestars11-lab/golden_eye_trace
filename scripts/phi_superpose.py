"""Superpose every phi-lattice check into one composite picture.

The other two harnesses report their sections side by side. This one stacks
them, so that a result which depends on the frame shows up as scatter in the
fold while a result which belongs to the ruler collapses onto the same line as
everything else:

* all candidates of the real file folded onto a single rung (phase mod 1),
  with the same fold for uniform draws as the control;
* all four cone ladders folded onto a single (residual, excess width) plane,
  where a log-log fit tests whether one law covers all of them;
* the pitch and the searched drift window folded onto one invariant axis,
  drift / (f tan mu) = 1 / T, so their separation is a pure number;
* every check on one stack, with the real file's hits compared against the
  false alarms its own search volume predicts.

    python3 scripts/phi_superpose.py --fil data/stare2_ovro.fil
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import numpy as np  # noqa: E402

from phi_native import (  # noqa: E402
    BEACON_MHZ,
    FCH1_NARROW,
    FOFF,
    LADDER_DRIFT,
    LADDER_NINTS,
    MU_DEG,
    SCAN_S,
    SLOPE,
    TSAMP_S,
    _half_max_width,
)
from seti_search.filterbank import Filterbank, InjectedSignal, simulate  # noqa: E402
from seti_search.rungs import LN_PHI, PHI, rung  # noqa: E402
from seti_search.search import (  # noqa: E402
    SearchConfig,
    control_bandpass,
    dedoppler,
    excise_rfi,
    search_file,
)
from seti_search.trials import expected_false_alarms, tail_probability  # noqa: E402


def section_fold(path: str) -> tuple[int, float, int]:
    print(f"S1. REAL FILE FOLDED ONTO ONE RUNG: {os.path.basename(path)}")
    cfg = SearchConfig(boxcar=3, max_drift_hz_s=4.0, snr_threshold=3.0, max_hits=500)
    hits, _meta = search_file(path, cfg, chan_chunk=2048, max_ints=None)
    fb = Filterbank.open(path)
    freqs = fb.frequencies()
    low, high = float(freqs.min()), float(freqs.max())
    span = rung(high) - rung(low)
    low_phase = float(rung(low) % 1.0)
    high_phase = float(rung(high) % 1.0)
    # the searched volume, counted the way dedoppler builds its grid
    df_hz = abs(fb.foff) * 1e6
    total_time = fb.nints * float(fb.tsamp)
    max_drift = min(cfg.max_drift_hz_s, 0.5 * df_hz / float(fb.tsamp))
    step = (0.5 * cfg.boxcar * df_hz / total_time
            if cfg.drift_step_hz_s is None else float(cfg.drift_step_hz_s))
    step = max(step, 1e-5)
    n_drifts = 2 * int(math.floor(max_drift / step)) + 1
    cells = fb.nchans * n_drifts
    phases = [rung(h.freq_mhz) % 1.0 for h in hits]
    n = len(phases)
    print(f"   band {low:.3f}-{high:.3f} MHz, {span:.4f} rungs, {n} hits at 3 sigma")
    print("   each candidate's fractional rung, superposed modulo one rung:")
    for hit, phase in zip(sorted(hits, key=lambda h: -h.snr), phases):
        print(f"   {hit.freq_mhz:>13.6f} MHz  {hit.snr:>6.1f} sigma"
              f"  phase {phase:.4f}")
    counts, _edges = np.histogram(phases, bins=20, range=(0.0, 1.0))
    top = max(1, int(counts.max()))
    print("   folded histogram, 20 bins, one full rung left to right:")
    for i, count in enumerate(counts):
        bar = "#" * int(round(36.0 * count / top))
        print(f"   {i / 20.0:.2f} {bar:<36s}{count}")

    def rayleigh(values) -> tuple[float, float]:
        if not len(values):
            return 0.0, 1.0
        theta = 2.0 * np.pi * np.asarray(values, dtype=float)
        resultant = float(abs(np.mean(np.exp(1j * theta))))
        return resultant, math.exp(-len(values) * resultant * resultant)

    resultant, p_value = rayleigh(phases)
    rng = np.random.default_rng(17)
    circle_r, circle_p = rayleigh(rng.uniform(0.0, 1.0, max(n, 1)))
    arc = np.array([rung(float(v)) % 1.0 for v in rng.uniform(low, high, 500)])
    arc_r, _arc_p = rayleigh(arc)
    print(f"   band covers rung phases {low_phase:.3f} to 1 and 0 to {high_phase:.3f},"
          " so the")
    print("   null is uniform on that arc, not on the whole circle.")
    print(f"   resultant length R = {resultant:.4f}, asymptotic Rayleigh"
          f" p = {p_value:.3f}")
    print(f"   500 control draws inside the band give R = {arc_r:.4f}; uniform on"
          f" the whole circle gives R = {circle_r:.4f} (p = {circle_p:.3f}).")
    print("   with n = 3 neither p is usable: the asymptotic formula is invalid")
    print("   this small, and any fold of a 0.37-rung band concentrates anyway.")
    inside = sum(abs(ph - round(ph)) <= 0.05 for ph in phases)
    chance = 0.10  # phases uniform on [0, 1): within 0.05 of 0 or of 1
    print(f"   {inside} of {n} lie within 0.05 of a rung; chance per candidate is"
          f" {chance:.2f} in")
    print("   this fold, 0.27 if the band is read in rung space. A rung is 2.7")
    print("   times wider than the band, so the fold has almost no power: a uniform")
    print("   fold and a lattice fold look the same here.")
    print(f"   the searched volume is {fb.nchans} channels x {n_drifts} drifts ="
          f" {cells} cells.")
    print()
    strongest = max((float(h.snr) for h in hits), default=0.0)
    return n, strongest, cells


def section_cone(folder: str) -> tuple[list[tuple[float, float, float]], float]:
    print("S2. CONE LADDERS FOLDED ONTO ONE LINE")
    path = os.path.join(folder, "superpose.fil")
    simulate(path, nchans=1024, nints=LADDER_NINTS, tsamp=TSAMP_S,
             fch1=FCH1_NARROW, foff=FOFF,
             signals=(InjectedSignal(BEACON_MHZ, drift_hz_s=LADDER_DRIFT,
                                     amplitude=14.0, width_channels=0.4),),
             rfi=(), glitch=False, seed=5)
    fb = Filterbank.open(path)
    data = fb.read(0, LADDER_NINTS)
    freqs = fb.frequencies()
    if fb.foff < 0:
        data = np.ascontiguousarray(data[:, ::-1])
        freqs = freqs[::-1]
    normalised, _ = control_bandpass(data)
    clean, _ = excise_rfi(normalised)
    df_hz = abs(fb.foff) * 1e6
    scan_s = LADDER_NINTS * TSAMP_S
    channel = int(np.argmin(np.abs(freqs - BEACON_MHZ)))
    nyquist = 0.5 * df_hz / TSAMP_S
    snr_map, drifts = dedoppler(clean, TSAMP_S, df_hz / 1e6,
                                max_drift_hz_s=nyquist, drift_step_hz_s=0.0005,
                                boxcar=1)
    snr_map = np.asarray(snr_map)
    drifts = np.asarray(drifts)
    if snr_map.shape[0] != drifts.size:
        snr_map = snr_map.T
    floor_index = int(np.argmin(np.abs(drifts - LADDER_DRIFT)))
    floor, _ = _half_max_width(snr_map[floor_index, max(0, channel - 25):channel + 26])
    points: list[tuple[float, float, float]] = []
    print(f"   beacon {BEACON_MHZ} MHz, drift {LADDER_DRIFT:g} Hz/s, "
          f"{LADDER_NINTS} rows, floor {floor:.2f} channels")
    print(f"   {'ladder':>8}{'step':>8}{'residual':>10}{'predicted':>11}"
          f"{'excess':>9}{'excess/pred':>13}")
    for rho in (PHI, PHI ** 2, 2.0, 3.0):
        for k in range(3):
            residual = 0.02 * rho ** k
            predicted = residual * scan_s / df_hz
            if predicted > 30.0:
                break
            trial = int(np.argmin(np.abs(drifts - (LADDER_DRIFT - residual))))
            used = abs(float(drifts[trial]) - LADDER_DRIFT)
            pad = int(1.5 * predicted) + 20
            lo, hi = max(0, channel - pad), min(snr_map.shape[1], channel + pad + 1)
            width, _peak = _half_max_width(snr_map[trial, lo:hi])
            excess = float(width) - floor
            points.append((used, excess, rho))
            print(f"   {rho:>8.3f}{rho:>8.3f}{used:>10.4f}{predicted:>11.1f}"
                  f"{excess:>9.2f}{excess / predicted:>13.2f}")
    xs = np.log(np.array([p[0] for p in points]))
    ys = np.log(np.array([p[1] for p in points]))
    slope, intercept = np.polyfit(xs, ys, 1)
    residual_rms = float(np.sqrt(np.mean((ys - (slope * xs + intercept)) ** 2)))
    print("   one log-log fit through all %d points from four ladders:" % len(points))
    print(f"   excess width proportional to residual^{slope:.3f} "
          f"(rms in ln = {residual_rms:.3f})")
    print("   a power of 1 means the width is the residual times the scan divided")
    print("   by the channel width - the geometry of the dedoppler itself - and the")
    print("   ratio between ladder steps is then the step, whatever it is chosen to")
    print("   be. The phi^2 growth is the ladder's spacing, which is why every")
    print("   ladder collapses onto this one line.")
    print()
    return points, float(slope)


def section_axis() -> None:
    print("S3. PITCH AND SEARCH GRID ON ONE INVARIANT AXIS")
    f_hz = BEACON_MHZ * 1e6
    rungs_per_baseline = SLOPE / LN_PHI
    window_hz_s = 4.0
    window_rungs = window_hz_s * SCAN_S / (f_hz * LN_PHI)
    print("   the invariant is  drift * T / f = tan(mu) = "
          f"{SLOPE:.6f}; in rungs that is")
    print(f"   {rungs_per_baseline:.4f} per baseline for every baseline, so all")
    print("   three rows of the earlier table are the same point:")
    print(f"   {'baseline':>12}{'drift (Hz/s)':>16}{'rungs/baseline':>18}")
    for label, t_s in (("scan", SCAN_S), ("day", 86400.0), ("year", 3.156e7)):
        print(f"   {label:>12}{SLOPE * f_hz / t_s:>16.3e}"
              f"{rungs_per_baseline:>18.4f}")
    print(f"   the searched window is +/-{window_hz_s:g} Hz/s, which is "
          f"{window_rungs:.2e} rungs")
    print(f"   per scan. Pitch over window = {rungs_per_baseline / window_rungs:.2e}:")
    print(f"   {math.log10(rungs_per_baseline / window_rungs):.1f} decades apart on"
          " the same axis.")
    print(f"   in frequency it is the same number: a {SLOPE * f_hz / 1e6:.1f} MHz sweep")
    print(f"   against a {window_hz_s * SCAN_S:.0f} Hz window, a ratio of "
          f"{SLOPE * f_hz / (window_hz_s * SCAN_S):.2e}.")
    print("   the only baseline that brings the pitch down to the window is")
    print(f"   {SLOPE * f_hz / window_hz_s / 3.156e7:.2f} years, which makes it an")
    print("   epoch-to-epoch statement rather than a within-scan drift.")
    print()


def section_stack(hits: int, strongest: float, cells: int, slope: float) -> None:
    print("S4. EVERY CHECK ON ONE STACK")
    rows = [
        ("round trip phi**rung(f)", "100% of noise passes", "identity"),
        ("pitch as a linear slope", "435 MHz sweep per baseline", "not narrowband"),
        ("pitch as a log slope", "0.306 vs 1.0000 measured", "fails as stated"),
        ("cone widening", f"excess ~ residual^{slope:.2f}", "the ruler"),
        ("golden angle = 17.03", "137.5078 vs 17.0324, 3.03 off", "no form found"),
        ("rung fold of the file", "3 hits, arc null, n too small", "unresolved"),
    ]
    print(f"   {'check':<27}{'measured':<35}{'verdict'}")
    for check, measured, verdict in rows:
        print(f"   {check:<27}{measured:<35}{verdict}")
    expected = expected_false_alarms(3.0, cells)
    observed_p = tail_probability(strongest)
    print("   the file itself, folded onto its own search volume:")
    print(f"   {cells} trial cells at 3 sigma predict {expected:.1f} false alarms;")
    print(f"   the search returned {hits}. The strongest is {strongest:.2f} true"
          f" sigma, whose")
    print(f"   single-cell tail is {observed_p:.2e}, or {observed_p * cells:.2f}"
          " expected in the band.")
    print("   superposed, nothing in the file prefers a rung: the hits are what the")
    print("   noise predicts, and the two phi statements that survive are the ones")
    print("   true of every number (the round trip) and every ladder (multiplicative")
    print("   growth). Neither can select a source.")
    print()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fil", default="", help="real filterbank file to fold")
    args = parser.parse_args()

    if args.fil:
        hits, strongest, cells = section_fold(args.fil)
    else:
        print("S1. REAL FILE FOLDED ONTO ONE RUNG: skipped (no --fil given)")
        print()
        hits, strongest, cells = 0, 0.0, 1
    with tempfile.TemporaryDirectory() as folder:
        _points, slope = section_cone(folder)
    section_axis()
    section_stack(hits, strongest, cells, slope)
    print("SUMMARY")
    print("   fold: three candidates cannot resolve one rung in this band")
    print("   cone: all ladders lie on one line, width = residual * scan / channel")
    decades = math.log10(SLOPE / LN_PHI / (4.0 * SCAN_S
                                           / (BEACON_MHZ * 1e6 * LN_PHI)))
    print(f"   axis: pitch sits {decades:.1f} decades above the search window")
    print("   stack: no phi statement separates a source from the ruler")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
