"""Does the quoted 5-sigma flux limit mean what the limits block says it means?

A limit is only worth quoting if a source at that flux actually behaves the way
the number claims. This script puts a synthetic beacon at a known flux through
the whole pipeline and measures three things:

1. **flux transfer** - inject well above the limit, recover the hit, convert its
   SNR back to Jy, and compare with what was injected. If the conversion carries
   an offset, every quoted limit inherits it.
2. **the limit itself** - a source sitting exactly at the 5-sigma limit should be
   detected about half the time, and one at the 95% completeness limit about
   95% of the time. That is what `limit_jy`'s two numbers mean, and it is a
   statement about the *distribution* of recovered SNR, so it needs many noise
   realisations, not one.
3. **the size of any offset** - because the answer is a rate, not a single
   number, it comes with the binomial error of the trials used.

The file is synthetic (exponential noise, SIGPROC header, Breakthrough
Listen-like fine channelisation: 2.86 Hz channels, 18.2536 s integrations), not
real observing data. It exercises the statistics of the pipeline; it cannot
confirm any behaviour of a real receiver.

    python3 scripts/limit_check.py --trials 60

Two conventions have to be stated, because they are where this can go wrong:

* An amplitude in the simulator is not a flux. It becomes one only once a line
  width is assumed, and the limits assume an *unresolved* line. The beacon is
  therefore injected with ``width_channels=0.4`` - a 0.94-channel FWHM, one
  channel to the pipeline - and not with the simulator's default 1.0, which is a
  2.35-channel FWHM line the pipeline would credit with 1.16x the flux.
* The flux of a given amplitude is measured, not assumed: the transfer from
  amplitude to reported Jy is calibrated once at high SNR, and the trials are
  then injected at the amplitude that this transfer attributes to the nominal
  flux. Otherwise the nominal flux and the injected flux differ by whatever the
  profile discretisation does to it, and the measured rate answers the wrong
  question.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402

from seti_search.filterbank import InjectedSignal, simulate  # noqa: E402
from seti_search.flux import (  # noqa: E402
    COMPLETENESS_MARGIN,
    limit_jy,
    noise_per_sample_jy,
)
from seti_search.run import scan_sensitivity, spectrum_flux  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

FCH1 = 1420.56
FOFF = -2.86e-6
NINTS = 16
TSAMP = 18.2536
BOXCAR = 3
MAX_DRIFT = 4.0
SEFD = 10.0
THRESHOLD = 5.0
SIGNAL_MHZ = 1420.554
DRIFT_HZ_S = 0.05
DF_HZ = abs(FOFF) * 1e6
WIDTH_CHANNELS = 0.4
DRIFT_TOLERANCE = 0.05  # the CLI default; must exceed one drift step


def config() -> SearchConfig:
    return SearchConfig(
        boxcar=BOXCAR,
        max_drift_hz_s=MAX_DRIFT,
        snr_threshold=THRESHOLD,
        max_hits=50,
    )


def search(path: str, nchans: int):
    return search_file(path, config(), chan_chunk=nchans, max_ints=None)


def limits(meta: dict) -> tuple[float, float]:
    common = dict(
        sensitivity=scan_sensitivity(meta),
        sefd_jy=SEFD,
        df_hz=DF_HZ,
        tsamp_s=TSAMP,
        boxcar=BOXCAR,
        line_width_channels=1.0,
    )
    return (
        limit_jy(threshold=THRESHOLD, margin_sigma=0.0, **common),
        limit_jy(threshold=THRESHOLD, margin_sigma=COMPLETENESS_MARGIN, **common),
    )


def simulate_file(folder: str, amplitude: float, seed: int, nchans: int) -> str:
    path = os.path.join(folder, "beacon.fil")
    simulate(
        path,
        nchans=nchans,
        nints=NINTS,
        tsamp=TSAMP,
        fch1=FCH1,
        foff=FOFF,
        signals=(
            InjectedSignal(
                SIGNAL_MHZ,
                drift_hz_s=DRIFT_HZ_S,
                amplitude=amplitude,
                width_channels=WIDTH_CHANNELS,
            ),
        )
        if amplitude > 0
        else (),
        rfi=(),
        glitch=False,
        seed=seed,
    )
    return path


def find_beacon(hits):
    """The hit the injected beacon should produce, or None if it was missed.

    The drift window is the CLI's own default tolerance. It has to exceed the
    search's drift step (0.015 Hz/s here): a real signal's measured drift sits up
    to a step down from the true rate, and a window narrower than that reports
    missed detections that the pipeline actually made.
    """
    df_mhz = DF_HZ / 1e6
    close = [
        hit
        for hit in hits
        if abs(hit.freq_mhz - SIGNAL_MHZ) <= 3.0 * df_mhz
        and abs(hit.drift_hz_s - DRIFT_HZ_S) < DRIFT_TOLERANCE
    ]
    return max(close, key=lambda hit: hit.snr) if close else None


def recovered_flux(snr: float, meta: dict) -> float:
    """Jy for a one-channel line, through the same path the CLI uses."""
    flux = spectrum_flux(snr, meta, SEFD, 1.0)
    if flux is None:
        raise RuntimeError("flux conversion gave no value")
    return float(flux[0])


def profile_flux() -> float:
    """Sampled flux of the injected profile, in units of the peak amplitude.

    A flat one-channel line has 1.0; the Gaussian beacon leaks a little into its
    neighbours, and that is a real difference in injected flux.
    """
    offsets = np.arange(-4, 5, dtype=np.float64)
    return float(np.sum(np.exp(-0.5 * (offsets / WIDTH_CHANNELS) ** 2)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--trials", type=int, default=60,
                        help="noise realisations per injected flux (default 60)")
    parser.add_argument("--nchans", type=int, default=4096,
                        help="simulated band width (default 4096)")
    args = parser.parse_args()
    nchans = int(args.nchans)

    with tempfile.TemporaryDirectory() as folder:
        # 1. the limits this configuration produces, from a no-signal file
        noise_path = simulate_file(folder, 0.0, seed=100, nchans=nchans)
        _noise_hits, noise_meta = search(noise_path, nchans)
        limit50, limit95 = limits(noise_meta)
        calibration = float(noise_meta.get("snr_calibration", 1.0))
        print(f"band {nchans:,} ch @ {DF_HZ:.3g} Hz, {NINTS} x {TSAMP:g} s,"
              f" SEFD {SEFD:g} Jy, boxcar {BOXCAR}, threshold {THRESHOLD:g} sigma")
        print(f"noise calibration on this synthetic file: {calibration:.3f}")
        print(f"channels blanked as RFI: {noise_meta['blanked_channels']:,}"
              f" of {noise_meta['nchans']:,} (the simulator's exponential noise is"
              " heavier-tailed than real radiometer noise)")
        print(f"quoted limit: {limit50:.3g} Jy (50% complete),"
              f" {limit95:.3g} Jy (95% complete),"
              f" unresolved line {WIDTH_CHANNELS:g} ch FWHM"
              f" {WIDTH_CHANNELS * 2.3548:.2f} ch")

        # 2. transfer, measured far enough above the limit that SNR is not an
        #    issue, but below the 20-sigma carrier cut that would rightly blank
        #    the beacon as a terrestrial tone
        strong_flux = 5.0 * limit50
        probe_amplitude = 6.0
        path = simulate_file(folder, probe_amplitude, seed=101, nchans=nchans)
        hits, meta = search(path, nchans)
        hit = find_beacon(hits)
        if hit is None:
            print(f"FAILED: a beacon at amplitude {probe_amplitude:g} was not"
                  " detected at all")
            return 1
        flux_at_probe = recovered_flux(hit.snr, meta)
        transfer = flux_at_probe / probe_amplitude  # Jy per unit amplitude
        gain = hit.snr / probe_amplitude  # SNR per unit amplitude
        flat_flux = probe_amplitude * profile_flux() * noise_per_sample_jy(
            SEFD, DF_HZ, TSAMP
        )
        print(f"\ntransfer at amplitude {probe_amplitude:g}:"
              f" SNR {hit.snr:.1f}, reported {flux_at_probe:.3g} Jy,"
              f" injected {flat_flux:.3g} Jy"
              f" (ratio {flux_at_probe / flat_flux:.3f})")
        print(f"  the conversion is conservative by"
              f" {100.0 * (flux_at_probe / flat_flux - 1.0):+.0f}% for this"
              " 0.94-channel line, which is what assuming an unresolved line"
              " costs")

        # 3. detection rate at and around the limit, injecting the amplitude the
        #    transfer attributes to each nominal flux
        print(f"\n{'flux (Jy)':>10}{'SNR in':>9}{'predicted':>11}{'measured':>10}"
              f"{'recovered':>11}")
        for ratio in (0.6, 1.0, limit95 / limit50, 2.0):
            flux = ratio * limit50
            amplitude = flux / transfer
            snr_in = gain * amplitude
            predicted = 0.5 * math.erfc((THRESHOLD - snr_in) / math.sqrt(2.0))
            found = 0
            fluxes = []
            for trial in range(args.trials):
                path = simulate_file(
                    folder, amplitude, seed=1000 + trial, nchans=nchans
                )
                hits, trial_meta = search(path, nchans)
                beacon = find_beacon(hits)
                if beacon is not None:
                    found += 1
                    fluxes.append(recovered_flux(beacon.snr, trial_meta))
            measured = found / args.trials
            error = math.sqrt(max(measured * (1.0 - measured), 1e-12) / args.trials)
            mean = float(np.mean(fluxes)) if fluxes else float("nan")
            print(f"{flux:>10.3g}{snr_in:>9.2f}"
                  f"{predicted:>11.3f}{measured:>10.3f}"
                  f"{mean:>11.3g}"
                  f"   (+/-{error:.3f})")
        print("\nmeasured is the fraction of trials the beacon cleared the"
              " threshold within 3 channels")
        print(f"and {DRIFT_TOLERANCE:g} Hz/s; predicted is the Gaussian chance"
              " that a source at that SNR in lands")
        print("above the threshold. The quoted 50% and 95% points are the first"
              " two and third rows.")
        print("Recovered flux is averaged over detections only, so it is biased"
              " high near the limit.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
