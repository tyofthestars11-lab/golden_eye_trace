"""Test the three phi-lattice directives as stated, and report what they return.

The proposal is that (1) a rung test ``rung(f) = ln f / ln phi`` identifies an
intelligent address, (2) a "golden pitch angle" of about 17.03 degrees is the
drift a real signal must have, and (3) a genuine signal's smear cone widens by
phi^2 = 2.618 per step instead of at random. Each of those is a quantitative
claim about measurable quantities, so each is checked here rather than argued.

Three of the checks are properties of arithmetic and of the detector, so they
need no observations at all:

* the round-trip ``phi ** (ln f / ln phi)`` is ``f`` for *every* positive number,
  by construction. A test that every input passes has no discriminating power;
  the pass rate on pure noise is the number that matters.
* an angle between two axes fixes a drift once both axes have consistent units,
  and the drift it fixes is a sweep of tan(mu) * f = 435 MHz per baseline. That
  is not a narrowband line; it is spread across the whole band. The only
  scale-invariant reading of the same angle is a log-log slope (drift growing as
  f^tan(mu)), which scripts/phi_native.py tests against the sky.
* a filterbank record is real total power. Squaring a complex field destroys its
  phase, so there is no phase left in the file for a phase-lock filter to use.

The remaining two are measured from the pipeline: the width of a beacon's
dedoppler cone as the search grid changes, and the rung drift of whatever the
real file actually contains against the same statistic for random frequencies.

    python3 scripts/phi_check.py --fil data/stare2_ovro.fil
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402

from seti_search.filterbank import Filterbank, InjectedSignal, simulate  # noqa: E402
from seti_search.rungs import (  # noqa: E402
    LN_PHI,
    PHI,
    rung,
    rung_drift,
    rung_uniformity,
)
from seti_search.search import (  # noqa: E402
    SearchConfig,
    control_bandpass,
    dedoppler,
    excise_rfi,
    search_file,
)

MU_DEG = 17.0324
BEACON_MHZ = 1420.554
CHANNEL_HZ = 2.86
TSAMP_S = 18.2536
NINTS = 16
FCH1 = 1420.56
FOFF = -2.86e-6
SCAN_S = NINTS * TSAMP_S


def section_roundtrip() -> None:
    print("1. ROUND TRIP: rung(f) = ln f / ln phi, then phi ** rung(f)")
    rng = np.random.default_rng(7)
    n = 200_000
    freqs = rng.uniform(1280.0, 1530.0, n)
    noise = np.abs(rng.normal(size=n)) * 1500.0 + 1e-12
    for label, values in (("observing band", freqs), ("pure noise", noise)):
        back = np.array([PHI ** rung(float(v)) for v in values])
        rel = float(np.max(np.abs(back - values) / values))
        hits = int(np.count_nonzero(np.abs(back - values) / values < 1e-12))
        print(f"   {label:14s} n={values.size:>7,}  max relative error {rel:.2e}"
              f"  verified {hits / values.size:.4f}")
    print("   verdict: the identity exp(ln x) = x holds for every positive number,")
    print("   so the round trip passes 100% of signal AND 100% of noise. A test with")
    print("   a 100% pass rate on noise rejects nothing and cannot act as a filter.")
    print("   It is a definition, not a measurement.\n")


def section_angle() -> None:
    print(f"2. GOLDEN PITCH: mu = {MU_DEG} deg as a slope on a linear axis")
    slope = math.tan(math.radians(MU_DEG))
    f_hz = BEACON_MHZ * 1e6
    print(f"   tan(mu) = {slope:.6f};  f = {f_hz:.5g} Hz")
    print("   drift = tan(mu) * f / T, with f in Hz and T in seconds:")
    print(f"   {'T (s)':>10}{'drift (Hz/s)':>20}{'sweep per baseline (Hz)':>26}")
    for label, t_s in (("1", 1.0), (f"{TSAMP_S:g}", TSAMP_S), (f"{SCAN_S:g}", SCAN_S)):
        drift = slope * f_hz / t_s
        print(f"   {label:>10}{drift:>20.3e}{drift * t_s:>26.3e}")
    sweep = slope * f_hz
    print(f"   the sweep is tan(mu) * f = {sweep / 1e6:.1f} MHz per baseline at every")
    print("   baseline, and the drift is 1.5e6 Hz/s or more. A line sweeping 435 MHz")
    print("   per scan is not narrowband: it fills the band by construction.")
    smallest = math.degrees(math.atan(4.0 * SCAN_S / f_hz))
    print(f"   staying inside the +/-4 Hz/s searched here needs an angle of at most"
          f" {smallest:.2e} deg,")
    print(f"   which is {MU_DEG / smallest:.2e} times smaller than 17.03 deg.")
    print("   verdict: with consistent units the pitch fixes a drift of order 1e6 Hz/s,")
    print("   five orders above any narrowband search - it is not a line the dedoppler")
    print("   could find. The scale-invariant reading (drift proportional to")
    print("   f^tan(mu), a log-log slope) is tested against the sky in")
    print("   scripts/phi_native.py.\n")


def section_phase(sample_path: str) -> None:
    print("3. PHASE LOCK: what a filterbank file actually contains")
    fb = Filterbank.open(sample_path)
    block = fb.read(0, 2)
    print(f"   {os.path.basename(sample_path)}: dtype {block.dtype},"
          f" complex={np.iscomplexobj(block)}, {fb.nchans} ch x {fb.nints} rows")
    n = 4096
    t = np.arange(n)
    tone = np.exp(2j * np.pi * 0.1 * t)
    rng = np.random.default_rng(11)
    phase_per_block = np.repeat(rng.uniform(0.0, 2 * np.pi, n // 64), 64)
    scrambled = np.exp(2j * np.pi * 0.1 * t + 1j * phase_per_block)
    locked_power = np.abs(tone) ** 2
    scrambled_power = np.abs(scrambled) ** 2
    print(f"   complex fields differ by up to {np.max(np.abs(tone - scrambled)):.3f}")
    print(f"   after square-law detection the total power differs by"
          f" {np.max(np.abs(locked_power - scrambled_power)):.2e}")
    print("   verdict: a filterbank stores |E|^2. Phase-coherent structure is gone")
    print("   before the file is written, and the dedoppler sum is incoherent by")
    print("   design, so a 'phase-lock filter' has no phase to lock to.\n")


def section_smear(folder: str) -> None:
    print("4. SMEAR CONE: does a beacon widen by phi^2 = %.3f per drift step?" % (PHI ** 2))
    path = os.path.join(folder, "smear.fil")
    true_drift = 0.1
    # The beacon has to move: a stationary line is indistinguishable from the
    # bandpass and gets divided out by the control spectrum before search.
    # 4096 channels at 2.86 Hz spans 11.7 kHz below fch1, so the beacon at
    # 1420.554 MHz (channel 2098) is inside the band; a 2048-channel band would
    # stop short of it and this section would measure pure noise.
    simulate(path, nchans=4096, nints=NINTS, tsamp=TSAMP_S, fch1=FCH1, foff=FOFF,
             signals=(InjectedSignal(BEACON_MHZ, drift_hz_s=true_drift,
                                     amplitude=14.0, width_channels=0.4),),
             rfi=(), glitch=False, seed=5)
    fb = Filterbank.open(path)
    data = fb.read(0, NINTS)
    freqs = fb.frequencies()
    if fb.foff < 0:
        data = np.ascontiguousarray(data[:, ::-1])
        freqs = freqs[::-1]
    normalised, _ = control_bandpass(data)
    clean, _ = excise_rfi(normalised)
    df_hz = abs(fb.foff) * 1e6
    channel = int(np.argmin(np.abs(freqs - BEACON_MHZ)))
    print(f"   beacon at {BEACON_MHZ} MHz, {true_drift:g} Hz/s, {NINTS} x {TSAMP_S:g} s")
    print(f"   {'step':>6}{'|drift|':>9}{'residual':>10}{'width ch':>10}"
          f"{'cells >5s':>11}{'peak':>7}{'ratio':>8}")
    sequences = {}
    for step in (0.05, 0.03):
        snr_map, drifts = dedoppler(clean, TSAMP_S, df_hz / 1e6,
                                    max_drift_hz_s=0.6, drift_step_hz_s=step,
                                    boxcar=3)
        snr_map = np.asarray(snr_map)
        drifts = np.asarray(drifts)
        if snr_map.shape[0] != drifts.size:
            snr_map = snr_map.T
        rows = []
        for i, drift in enumerate(drifts):
            residual = abs(float(drift) - true_drift)
            if residual == 0.0 or residual > 0.3:
                continue
            width = residual * SCAN_S / df_hz
            pad = int(2.0 * width) + 15
            lo = max(0, channel - pad)
            hi = min(snr_map.shape[1], channel + pad + 1)
            window = snr_map[i, lo:hi]
            cells = int(np.count_nonzero(window >= 5.0))
            rows.append((residual, abs(float(drift)), width, cells, float(window.max())))
        rows.sort()
        ratios = []
        previous = None
        for residual, drift_abs, width, cells, peak in rows:
            ratio = "" if not previous else f"{cells / previous:.2f}"
            if previous:
                ratios.append(cells / previous)
            print(f"   {step:>6.2f}{drift_abs:>9.2f}{residual:>10.2f}"
                  f"{width:>10.1f}{cells:>11d}{peak:>7.1f}{ratio:>8}")
            previous = cells or previous
        sequences[step] = ratios
    print(f"   step 0.05 ratios: {[round(r, 2) for r in sequences[0.05][:6]]}")
    print(f"   step 0.03 ratios: {[round(r, 2) for r in sequences[0.03][:6]]}")
    print("   verdict: the cone width is |drift - true drift| x scan / channel, and the")
    print("   cell count follows it. The ratio between successive trials is set by the")
    print("   trial spacing and the boxcar, so it changes when the grid changes while")
    print("   the signal does not. There is no fixed phi^2, and a ratio cannot be")
    print("   evidence anyway: any value is available by choosing a step.\n")


def section_real(path: str) -> None:
    print(f"5. RUNG TEST ON THE REAL FILE: {os.path.basename(path)}")
    cfg = SearchConfig(boxcar=3, max_drift_hz_s=4.0, snr_threshold=3.0, max_hits=500)
    hits, meta = search_file(path, cfg, chan_chunk=2048, max_ints=None)
    fb = Filterbank.open(path)
    freqs = fb.frequencies()
    low, high = float(freqs.min()), float(freqs.max())
    rng = np.random.default_rng(3)
    random_mhz = rng.uniform(low, high, 50_000)
    random_drifts = np.array([rung_drift(float(f)) for f in random_mhz])
    tops = sorted(hits, key=lambda h: -h.snr)[:10]
    span = rung(high) - rung(low)
    df_hz = abs(fb.foff) * 1e6
    channel_rungs = math.log(1.0 + df_hz / (low * 1e6)) / LN_PHI
    per_rung_mhz = (PHI - 1.0) * 0.5 * (low + high)
    print(f"   {len(hits)} hits at 3 sigma; band {low:.3f}-{high:.3f} MHz")
    print(f"   the band spans {span:.3f} rungs; one rung is {per_rung_mhz:.1f} MHz;"
          f" one channel is")
    print(f"   {channel_rungs:.2e} rungs - every channel in the band has the same"
          f" nearest rung")
    print(f"   ({round(rung(low))}), so the address carries no information about"
          " frequency inside the band")
    print(f"   {'MHz':>13}{'SNR':>7}{'rung':>9}{'drift':>8}{'d rung':>9}")
    if tops:
        base = rung(tops[0].freq_mhz)
        for hit in tops:
            d_mhz = rung_drift(hit.freq_mhz)
            print(f"   {hit.freq_mhz:>13.6f}{hit.snr:>7.1f}"
                  f"{rung(hit.freq_mhz):>9.3f}{d_mhz:>8.3f}"
                  f"{rung(hit.freq_mhz) - base:>+9.3f}")
    print(f"   random frequencies in this band with drift <= 0.05:"
          f" {np.mean(random_drifts <= 0.05):.3f}, not 0.10, because 50,000 draws")
    print("   spread over 0.37 rungs cannot be uniform in a quantity they barely span")
    if len(hits) >= 4:
        test = rung_uniformity([h.freq_mhz for h in hits])
        print(f"   KS uniformity of all {test['n']} fractional rungs:"
              f" D={test['statistic']:.3f}, p={test['p_value']:.3f}")
        print("   (uniform is the null for a smooth population; a small p would mean")
        print("    the search grid is quantized, not the sky)")
    shift = rung(1e6)
    print(f"   changing the unit MHz -> Hz shifts every rung by {shift:.4f}, which is"
          " not an integer,")
    print("   so with a band wider than one rung the membership test would change with"
          " the label.")
    print("   verdict: within one 250 MHz band the rung test has no resolution at all:")
    print("   every frequency is at the same address, so it can neither select a")
    print("   candidate nor reject one. The round trip is satisfied by all of them")
    print("   because it is satisfied by every number.\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fil", default="", help="real filterbank file to test")
    args = parser.parse_args()

    section_roundtrip()
    section_angle()
    with tempfile.TemporaryDirectory() as folder:
        sample = args.fil or os.path.join(folder, "phase.fil")
        if not args.fil:
            simulate(sample, nchans=1024, nints=8, tsamp=1.0, fch1=FCH1,
                     foff=FOFF, signals=(), rfi=(), glitch=False, seed=1)
        section_phase(sample)
        section_smear(folder)
    if args.fil:
        section_real(args.fil)
    else:
        print("5. RUNG TEST ON THE REAL FILE: skipped (no --fil given)")
    print("SUMMARY")
    print("   round trip: passes 100% of noise - definition, not a filter")
    print("   pitch angle: implies a 435 MHz sweep per baseline - not narrowband")
    print("   phase lock: no phase exists in total-power data - not applicable")
    print("   smear cone: geometry of the search grid - not a signal property")
    print("   rung test: band is narrower than one rung - the test has no resolution")
    print("   None of these can promote a candidate; at most they can demote one.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
