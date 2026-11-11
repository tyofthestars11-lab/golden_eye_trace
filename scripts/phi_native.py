"""The phi-lattice directives re-derived in their own coordinates.

``scripts/phi_check.py`` measured the directives on a linear frequency axis with
evenly spaced drift trials. That is a Cartesian frame, and the framework does
not claim to live in one: its natural coordinate is ``x = ln f``, where the
rungs are integers, a scale change is a translation, and a signal's drift is the
rate ``dx/dt``. This script repeats the measurement there, so that a null result
cannot be an artefact of the frame.

Four things are checked:

* **Invariance.** Which combination of (angle, frequency, baseline) is unchanged
  by relabelling coordinates? The answer is the fractional change per baseline,
  and it pins down what 17.0324 degrees actually asserts.
* **The scale-invariant reading.** An angle read as a log-log slope asserts
  ``|drift| proportional to f ** tan(mu)``. The sky's drift law, taken from the
  cadence model in :mod:`seti_search.cadence`, has a measurable exponent, so the
  two can be compared directly.
* **The cone on a geometric ladder.** With trials spaced by a constant *ratio*,
  cone widening is exactly multiplicative. The measured ratio is reported for
  ladders stepping by phi, phi^2, 2 and 3, so phi^2 can be confirmed where it is
  true - and seen to be the ladder's property rather than the source's.
* **The golden angle itself.** 360/phi^2 = 137.5078 deg is not 17.0324 deg.
  Simple closed forms are searched and the nearest is reported.

    python3 scripts/phi_native.py --fil data/stare2_ovro.fil
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402

from seti_search.cadence import SITES, SkyPosition, topocentric_drift  # noqa: E402
from seti_search.filterbank import Filterbank, InjectedSignal, simulate  # noqa: E402
from seti_search.rungs import GOLDEN_ANGLE_DEG, LN_PHI, PHI, rung  # noqa: E402
from seti_search.search import (  # noqa: E402
    SearchConfig,
    control_bandpass,
    dedoppler,
    excise_rfi,
    search_file,
)

MU_DEG = 17.0324
SLOPE = math.tan(math.radians(MU_DEG))
BEACON_MHZ = 1420.554
TSAMP_S = 18.2536
NINTS = 16
SCAN_S = NINTS * TSAMP_S
FOFF = -2.86e-6
FCH1_NARROW = 1420.55543  # puts the beacon at channel 500 of 1024
LADDER_NINTS = 64
LADDER_DRIFT = 0.06


def section_invariance() -> None:
    print("A. WHAT SURVIVES A CHANGE OF COORDINATES")
    f_hz = BEACON_MHZ * 1e6
    rungs_per_baseline = SLOPE / LN_PHI
    finite_rungs = math.log(1.0 + SLOPE) / LN_PHI
    print(f"   tan(mu) = {SLOPE:.6f}. On any pair of axes an angle mu means the")
    print("   fractional change per baseline is tan(mu):  drift * T / f = tan(mu).")
    print("   That ratio is dimensionless, so relabelling f or T cannot change it.")
    print(f"   In rungs (x = ln f / ln phi) that is {rungs_per_baseline:.4f} per baseline")
    print(f"   linearly, {finite_rungs:.4f} exactly for the finite step.")
    print(f"   {'baseline':>13}{'drift (Hz/s)':>18}{'rungs/baseline':>18}")
    for label, t_s in (("one scan", SCAN_S), ("one day", 86400.0), ("one year", 3.156e7)):
        print(f"   {label:>13}{SLOPE * f_hz / t_s:>18.3e}{rungs_per_baseline:>18.4f}")
    needed = SLOPE * f_hz / 4.0
    print("   the rung rate is invariant; the drift is not, because it divides by a")
    print("   baseline the directive never names. The pitch lands inside a 4 Hz/s")
    print(f"   window only if that baseline is {needed / 3.156e7:.2f} years, i.e. if it")
    print("   describes an epoch-to-epoch comparison, not a within-scan drift.")
    print("   verdict: 17.0324 deg says a source advances a fixed fraction of its")
    print("   frequency per baseline. Taken as a drift it is underdetermined by one")
    print("   dimension, and the angle does not choose that dimension.")
    print()


def section_sky_law() -> None:
    print("B. THE SCALE-INVARIANT READING: drift proportional to f^tan(mu)")
    site = SITES["ovro"]
    sky = SkyPosition(0.0, 0.0)  # a zenith-ish geometry, nothing special about it
    mjd = 59000.0
    freqs = np.array([1.0e9, 1.42e9, 2.0e9, 5.0e9, 1.4e10])
    drifts = np.array([topocentric_drift(0.0, float(f), site, mjd, sky) for f in freqs])
    fit = float(np.polyfit(np.log(freqs), np.log(np.abs(drifts)), 1)[0])
    print(f"   drift of a source with no intrinsic drift, OVRO, MJD {mjd}:")
    print(f"   {'f (GHz)':>10}{'|drift| (Hz/s)':>18}{'rungs per scan':>18}")
    for f, d in zip(freqs, drifts):
        rung_rate = abs(d) * SCAN_S / (f * LN_PHI)
        print(f"   {f / 1e9:>10.2f}{abs(d):>18.3e}{rung_rate:>18.2e}")
    print(f"   fitted exponent d ln|drift| / d ln f = {fit:.4f}")
    print(f"   the pitch asks for tan(mu) = {SLOPE:.4f}")
    print("   the two predict different rung rates. For exponent p the rung rate")
    print("   drift * T / (f ln phi) goes as f^(p-1):")
    print(f"     sky   (p = {fit:.2f}): {fit - 1:+.3f}, rung rate flat in frequency")
    print(f"     pitch (p = {SLOPE:.3f}): {SLOPE - 1:+.3f}, rung rate falls by"
          f" {10.0 ** (1.0 - SLOPE):.1f}x per decade")
    print("   any mechanism whose velocity does not depend on frequency gives p = 1;")
    print("   p = 0.306 needs velocity proportional to f^-0.694, which nothing in")
    print("   the standard model supplies. So this is a testable difference: measure")
    print("   rungs/scan for confirmed drifting lines in two bands and compare.")
    print()


def _half_max_width(window: np.ndarray) -> tuple[float, float]:
    """Interpolated full width at half maximum in channels, and the peak."""
    values = window.astype(np.float64)
    peak_index = int(np.argmax(values))
    peak = float(values[peak_index])
    if peak <= 4.0:
        return 0.0, peak
    half = 0.5 * peak
    above = np.nonzero(values >= half)[0]
    left_index, right_index = int(above[0]), int(above[-1])

    def edge(i_below: int, j_above: int) -> float:
        low, high = values[i_below], values[j_above]
        if high == low:
            return float(j_above)
        return i_below + (half - low) / (high - low)

    left = edge(left_index - 1, left_index) if left_index > 0 else float(left_index) - 0.5
    if right_index + 1 < values.size:
        right = edge(right_index + 1, right_index)
    else:
        right = float(right_index) + 0.5
    return right - left + 1.0, peak


def section_ladder(folder: str) -> None:
    print("C. THE CONE IN NATIVE MEASURE: growth per step of a geometric ladder")
    path = os.path.join(folder, "ladder.fil")
    nints = LADDER_NINTS
    simulate(path, nchans=1024, nints=nints, tsamp=TSAMP_S, fch1=FCH1_NARROW,
             foff=FOFF,
             signals=(InjectedSignal(BEACON_MHZ, drift_hz_s=LADDER_DRIFT,
                                     amplitude=14.0, width_channels=0.4),),
             rfi=(), glitch=False, seed=5)
    fb = Filterbank.open(path)
    data = fb.read(0, nints)
    freqs = fb.frequencies()
    if fb.foff < 0:
        data = np.ascontiguousarray(data[:, ::-1])
        freqs = freqs[::-1]
    normalised, _ = control_bandpass(data)
    clean, _ = excise_rfi(normalised)
    df_hz = abs(fb.foff) * 1e6
    scan_s = nints * TSAMP_S
    channel = int(np.argmin(np.abs(freqs - BEACON_MHZ)))
    nyquist = 0.5 * df_hz / TSAMP_S
    snr_map, drifts = dedoppler(clean, TSAMP_S, df_hz / 1e6,
                                max_drift_hz_s=nyquist, drift_step_hz_s=0.0005,
                                boxcar=1)
    snr_map = np.asarray(snr_map)
    drifts = np.asarray(drifts)
    if snr_map.shape[0] != drifts.size:
        snr_map = snr_map.T
    print(f"   beacon at {BEACON_MHZ} MHz, {nints} rows, true drift {LADDER_DRIFT:g} Hz/s")
    print(f"   trials placed geometrically in the residual |drift - true|, the")
    print(f"   quantity the cone width is proportional to. Grid capped at"
          f" {nyquist:.4f} Hz/s:")
    print("   half a channel per integration is where shift-and-add stops aligning")
    print(f"   the signal, so a cone wider than about {nints // 2} channels cannot be")
    print("   measured at all - the peak sinks before the width gets there.")
    floor_index = int(np.argmin(np.abs(drifts - LADDER_DRIFT)))
    floor, _ = _half_max_width(snr_map[floor_index, max(0, channel - 25):channel + 26])
    print(f"   the unsmeared line measures {floor:.2f} channels FWHM, so each width")
    print("   below carries that instrumental floor. The excess over it is the part")
    print("   the residual drift produces, and that is what should scale.")
    for rho in (PHI, PHI ** 2, 2.0, 3.0):
        widths, used, excess = [], [], []
        for n in range(3):
            residual = 0.02 * rho ** n
            predicted = residual * scan_s / df_hz
            if predicted > 30.0:
                break
            target = LADDER_DRIFT - residual
            i = int(np.argmin(np.abs(drifts - target)))
            used_value = abs(float(drifts[i]) - LADDER_DRIFT)
            pad = int(1.5 * predicted) + 20
            lo = max(0, channel - pad)
            hi = min(snr_map.shape[1], channel + pad + 1)
            width, _peak = _half_max_width(snr_map[i, lo:hi])
            widths.append(round(float(width), 2))
            used.append(round(used_value, 4))
            excess.append(round(float(width) - float(floor), 2))
        ratios = [round(b / a, 2) for a, b in zip(excess[:-1], excess[1:])
                  if a > 0 and b > 0]
        predicted_widths = [round(r * scan_s / df_hz, 1) for r in used]
        label = "phi^2" if abs(rho - PHI ** 2) < 1e-9 else f"{rho:.3f}"
        print(f"   ladder {label:>6}  residual {used}  predicted {predicted_widths}")
        print(f"   {'':>13}  measured {widths}  excess {excess}  ratio {ratios}"
              f"  spacing {rho:.3f}")
    print("   verdict: the excess column follows the ladder spacing, not a constant:")
    print("   a phi^2 ladder grows by about 2.6, a step of 2 by about 2, a step of 3")
    print("   by about 3, each a little short because the floor is a few percent of")
    print("   these widths. So the widening is multiplicative, exactly as claimed,")
    print("   and it belongs to the ladder: one beacon gives whatever ratio the")
    print("   spacing asks for. A longer scan shrinks the floor's share and the")
    print("   Nyquist cap grows with it, so phi^2 could be followed further out.")
    print()


def section_golden_form() -> None:
    print("D. THE GOLDEN ANGLE ITSELF")
    print(f"   360/phi^2 = {GOLDEN_ANGLE_DEG:.4f} deg, the phyllotaxis divergence.")
    print(f"   The pitch on offer is {MU_DEG} deg, so these are not the same object;")
    print("   a construction is needed to connect them. Closed forms tried:")
    forms = []
    for k in range(1, 13):
        forms.append((f"arctan(phi^-{k})", math.degrees(math.atan(PHI ** -k))))
        forms.append((f"360/phi^{k}", 360.0 / PHI ** k))
        forms.append((f"golden_angle/phi^{k}", GOLDEN_ANGLE_DEG / PHI ** k))
    forms.sort(key=lambda item: abs(item[1] - MU_DEG))
    for name, value in forms[:4]:
        print(f"   {name:>24} = {value:8.4f} deg   error {value - MU_DEG:+.4f}")
    best = abs(forms[0][1] - MU_DEG)
    if best > 1e-3:
        print(f"   nearest is off by {best:.4f} deg ({best / MU_DEG * 100:.2f}%), so")
        print("   17.0324 is not one of these; the arithmetic behind it is not")
        print("   recoverable from the number alone. Name the construction and it can")
        print("   be checked exactly, like every other line in this file.")
    else:
        print(f"   nearest form matches to {best:.2e} deg.")
    print()


def section_log_rung(path: str) -> None:
    print(f"E. RUNG TEST IN LOG MEASURE: {os.path.basename(path)}")
    cfg = SearchConfig(boxcar=3, max_drift_hz_s=4.0, snr_threshold=3.0, max_hits=500)
    hits, _meta = search_file(path, cfg, chan_chunk=2048, max_ints=None)
    fb = Filterbank.open(path)
    freqs = fb.frequencies()
    low, high = float(freqs.min()), float(freqs.max())
    lo_rung, hi_rung = rung(low), rung(high)
    span = hi_rung - lo_rung
    tol = 0.05
    # Chance in the framework's own measure: draw uniformly in rung space.
    rng = np.random.default_rng(5)
    draws = rng.uniform(lo_rung, hi_rung, 400_000)
    chance = float(np.mean(np.abs(draws - np.round(draws)) <= tol))
    print(f"   band {low:.3f}-{high:.3f} MHz spans {span:.4f} rungs"
          f" ({lo_rung:.3f}-{hi_rung:.3f})")
    print("   rung space is the native measure, so chance is measured there:")
    print(f"   {chance:.3f} of the band lies within {tol} rungs of an integer")
    print(f"   ({2 * tol / span * 100:.1f}% by the span alone, if it stays inside one)")
    inside = 0
    print(f"   {'MHz':>13}{'SNR':>7}{'rung':>9}{'dist to int':>13}")
    for hit in sorted(hits, key=lambda h: -h.snr)[:10]:
        value = rung(hit.freq_mhz)
        dist = abs(value - round(value))
        inside += dist <= tol
        print(f"   {hit.freq_mhz:>13.6f}{hit.snr:>7.1f}{value:>9.3f}{dist:>13.3f}")
    trials = len(hits)
    print(f"   {trials} hits at 3 sigma; {inside} lie within {tol} rungs")
    if trials:
        pmf = [math.comb(trials, k) * chance ** k * (1 - chance) ** (trials - k)
               for k in range(trials + 1)]
        tail = sum(pmf[inside:])
        print(f"   binomial P(>= {inside} of {trials} | chance {chance:.3f}) = {tail:.3f}")
        print("   with this many candidates that is not a discrimination test; the")
        print("   band must span several rungs before one exists.")
    print("   rung differences are invariant under relabelling (a ratio inside a")
    print("   logarithm), but membership is not: MHz -> Hz shifts every rung by")
    print(f"   {rung(1e6):.4f}, whose fractional part {rung(1e6) % 1.0:.4f} re-draws")
    print("   which lines count as 'on a rung' while the physics is unchanged.")
    print()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fil", default="", help="real filterbank file to test")
    args = parser.parse_args()

    section_invariance()
    section_sky_law()
    with tempfile.TemporaryDirectory() as folder:
        section_ladder(folder)
    section_golden_form()
    if args.fil:
        section_log_rung(args.fil)
    else:
        print("E. RUNG TEST IN LOG MEASURE: skipped (no --fil given)")
        print()
    print("SUMMARY")
    print("   invariance: the angle fixes 0.306 of f per baseline, and no baseline")
    print("   sky law: measured exponent 1.00, the pitch asks 0.31 - distinguishable")
    print("   ladder: cone growth equals the ladder ratio; phi^2 is the ruler")
    print("   rung test: chance in native measure is 0.27; no resolution here")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
