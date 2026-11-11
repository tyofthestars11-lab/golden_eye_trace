"""Converting amplitudes to flux density.

The unit test here that matters is the pair: two scans of the same source on
receiveds 10x apart in SEFD must agree on the flux, because each scan's
amplitude is measured in its own noise. The absolute value is checked only
loosely, because the boxcar mismatch depends on a line width the pipeline does
not measure.
"""

import tempfile
import unittest

from seti_search.cadence import SITES, SkyPosition, doppler_factor, topocentric_drift
from seti_search.cli import _resolve_quantity, _scan_sefd
from seti_search.filterbank import InjectedSignal, RFILine, simulate
from seti_search.flux import (
    COMPLETENESS_MARGIN,
    K_M2_PER_JY,
    boxcar_scale,
    combine_fluxes,
    flux_density_jy,
    limit_jy,
    noise_per_sample_jy,
    run_limit,
    sefd_from_gain,
    sefd_from_tsys,
)
from seti_search.run import Scan, link_run, spectrum_flux
from seti_search.search import SearchConfig, search_file

FCH1 = 1420.4
FOFF = -2.86e-6
NCHANS = 4096
TSAMP = 1.0
MJD0 = 58967.23
OBSERVER = SITES["ovro"]
SKY = SkyPosition(75.0, 20.0)
SIGNAL_CHANNEL = 1200

DF_HZ = abs(FOFF) * 1e6
TOPO_HZ_AT_MJD0 = (FCH1 + FOFF * SIGNAL_CHANNEL) * 1e6
# the source's barycentric frequency, which is what a run links on
BARY_HZ = TOPO_HZ_AT_MJD0 / doppler_factor(OBSERVER, MJD0, SKY)


class UnitTests(unittest.TestCase):
    def test_sefd_from_area(self):
        # SEFD = 2k Tsys / A_eff; 2k/Jy is 2761.3 K m^2 / Jy
        self.assertAlmostEqual(K_M2_PER_JY, 2761.298, places=3)
        self.assertAlmostEqual(sefd_from_tsys(20.0, 100.0), 552.26, places=2)
        # twice the collecting area, half the SEFD
        self.assertAlmostEqual(
            sefd_from_tsys(20.0, 200.0), sefd_from_tsys(20.0, 100.0) / 2, places=6
        )

    def test_sefd_from_gain(self):
        self.assertAlmostEqual(sefd_from_gain(20.0, 0.1), 200.0)
        with self.assertRaises(ValueError):
            sefd_from_gain(20.0, 0.0)

    def test_noise_per_sample(self):
        self.assertAlmostEqual(noise_per_sample_jy(100.0, 1.0, 1.0), 100.0)
        # four times the bandwidth halves the noise
        self.assertAlmostEqual(noise_per_sample_jy(100.0, 4.0, 1.0), 50.0)
        # four times the integration time halves it too
        self.assertAlmostEqual(noise_per_sample_jy(100.0, 1.0, 4.0), 50.0)

    def test_boxcar_scale(self):
        self.assertAlmostEqual(boxcar_scale(3, 1.0), 3.0)
        self.assertAlmostEqual(boxcar_scale(3, 3.0), 1.0)
        # a line wider than the boxcar gains nothing further
        self.assertAlmostEqual(boxcar_scale(3, 10.0), 1.0)
        self.assertAlmostEqual(boxcar_scale(1, 1.0), 1.0)


class FluxMathTests(unittest.TestCase):
    def test_unresolved_line_pays_the_boxcar_penalty(self):
        # amplitude 2 in noise units, sensitivity 10, SEFD 100 Jy, 1 Hz, 1 s:
        # scale = 3 * 100 Jy, and the amplitude estimate is good to 1/10
        flux, error = flux_density_jy(
            amplitude=2.0, sensitivity=10.0, sefd_jy=100.0,
            df_hz=1.0, tsamp_s=1.0, boxcar=3, line_width_channels=1.0,
        )
        self.assertAlmostEqual(flux, 600.0)
        self.assertAlmostEqual(error, 30.0)

    def test_matched_line_has_no_penalty(self):
        flux, _ = flux_density_jy(
            amplitude=2.0, sensitivity=10.0, sefd_jy=100.0,
            df_hz=1.0, tsamp_s=1.0, boxcar=3, line_width_channels=3.0,
        )
        self.assertAlmostEqual(flux, 200.0)

    def test_sensitivity_sets_the_uncertainty_not_the_flux(self):
        common = dict(
            amplitude=2.0, sefd_jy=100.0, df_hz=1.0, tsamp_s=1.0, boxcar=3
        )
        weak, weak_err = flux_density_jy(sensitivity=10.0, **common)
        strong, strong_err = flux_density_jy(sensitivity=40.0, **common)
        self.assertAlmostEqual(weak, strong)
        self.assertAlmostEqual(strong_err, weak_err / 4)

    def test_combine_is_inverse_variance(self):
        combined = combine_fluxes([(100.0, 10.0, "a"), (100.0, 1.0, "b")])
        self.assertAlmostEqual(combined.jy, 100.0)
        self.assertAlmostEqual(combined.uncertainty_jy, 10.0 / (101.0**0.5))
        self.assertEqual(combined.n_scans, 2)
        self.assertEqual(combined.per_scan, {"a": 100.0, "b": 100.0})

    def test_combine_lets_the_better_scan_win(self):
        combined = combine_fluxes([(50.0, 1.0, "good"), (150.0, 100.0, "poor")])
        self.assertLess(abs(combined.jy - 50.0), 0.1)

    def test_combine_without_usable_estimates(self):
        self.assertIsNone(combine_fluxes([]))
        self.assertIsNone(combine_fluxes([(float("nan"), 0.0, "x")]))


class LimitTests(unittest.TestCase):
    """What a null result rules out, and why the weakest scan sets it."""

    COMMON = dict(sensitivity=10.0, sefd_jy=100.0, df_hz=1.0, tsamp_s=1.0, boxcar=3)

    def test_limit_is_where_the_threshold_sits(self):
        # a source at the threshold must have exactly the flux the limit quotes,
        # or the limit and the detections would not describe the same band
        limit = limit_jy(threshold=10.0, **self.COMMON)
        flux, _ = flux_density_jy(amplitude=10.0 / 10.0, **self.COMMON)
        self.assertAlmostEqual(limit, flux)

    def test_limit_scales_with_the_threshold_and_the_margin(self):
        base = limit_jy(threshold=10.0, **self.COMMON)
        stricter = limit_jy(threshold=10.0, margin_sigma=COMPLETENESS_MARGIN, **self.COMMON)
        self.assertAlmostEqual(stricter / base, (10.0 + COMPLETENESS_MARGIN) / 10.0)
        self.assertAlmostEqual(limit_jy(threshold=20.0, **self.COMMON), 2 * base)
        # assuming a wider line lowers the limit by the boxcar factor
        self.assertAlmostEqual(
            limit_jy(threshold=10.0, line_width_channels=3.0, **self.COMMON),
            base / 3.0,
        )

    def test_run_limit_is_the_weakest_scan_that_counts(self):
        limits = {"deep": 1.0, "middle": 3.0, "shallow": 10.0}
        self.assertEqual(run_limit(limits, 1), (1.0, "deep"))
        self.assertEqual(run_limit(limits, 2), (3.0, "middle"))
        self.assertEqual(run_limit(limits, 3), (10.0, "shallow"))
        # asking for more scans than exist rules nothing out
        self.assertIsNone(run_limit(limits, 4))
        self.assertIsNone(run_limit({}, 1))


class LimitCsvTests(unittest.TestCase):
    """A null result has to survive on disk, not just on screen."""

    def setUp(self):
        from types import SimpleNamespace
        from seti_search.cli import (
            LIMIT_FIELDS,
            _limit_row,
            _limits_path,
            _write_limits,
        )

        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.fields = LIMIT_FIELDS
        self.row_for = _limit_row
        self.path_for = _limits_path
        self.write = _write_limits
        self.args = SimpleNamespace(
            snr=10.0, line_width=1.0, max_drift=4.0, drift_step=None, false_alarm=None
        )
        self.meta = {
            "n_ints": 200,
            "boxcar": 3,
            "df_mhz": -1.0e-6,
            "tsamp": 1.0,
            "nchans": 1024,
            "blanked_channels": 24,
            "fch1": 1420.5,
            "foff": -1.0e-3,
        }
        self.path = f"{self.tmp.name}/run.linked.csv"

    def test_limits_sit_beside_the_results(self):
        self.assertEqual(self.path_for("run.linked.csv"), "run.linked.limits.csv")
        self.assertEqual(self.path_for(f"{self.tmp.name}/x"), f"{self.tmp.name}/x.limits.csv")

    def test_row_carries_every_column_and_its_own_inputs(self):
        row = self.row_for(self.args, "a.fil", self.meta, 10.0)
        self.assertEqual(list(row), self.fields)
        self.assertEqual(row["kind"], "file")
        self.assertEqual(row["sefd_jy"], "10")
        self.assertEqual(row["threshold_sigma"], "10")
        self.assertEqual(row["df_hz"], "1")
        self.assertAlmostEqual(float(row["searched_fraction"]), 0.9765625, places=3)
        self.assertAlmostEqual(float(row["freq_max_mhz"]), 1420.5)
        self.assertAlmostEqual(float(row["freq_min_mhz"]), 1420.5 - 1.023)
        # 200 integrations at 1 s and 1 Hz channels: drift steps of 7.5 mHz/s over
        # the 0.5 Hz/s Nyquist limit give 133 trials per channel, whole band
        self.assertEqual(int(row["trials_per_channel"]), 133)
        self.assertEqual(int(row["total_trials"]), 133 * 1024)
        # 10 sigma is about 1 false alarm in 10^23 cells
        self.assertLess(float(row["expected_false_alarms"]), 1e-15)
        self.assertAlmostEqual(float(row["false_alarm_log10"]), -23.118, places=2)

    def test_written_row_reads_back_as_the_same_limit(self):
        import csv

        from seti_search.run import scan_sensitivity

        row = self.row_for(self.args, "a.fil", self.meta, 10.0)
        path = self.path_for(self.path)
        self.write(path, [row])
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        self.assertEqual(list(rows[0]), self.fields)
        self.assertEqual(float(rows[0]["limit_jy"]), float(row["limit_jy"]))
        # and that number is the one the limit maths gives, not a rounded copy
        self.assertAlmostEqual(
            float(rows[0]["limit_jy"]),
            limit_jy(
                threshold=10.0,
                sensitivity=scan_sensitivity(self.meta),
                sefd_jy=10.0,
                df_hz=1.0,
                tsamp_s=1.0,
                boxcar=3,
            ),
            # the file keeps 6 significant figures, so compare at that width
            places=4,
        )

    def test_empty_record_is_still_a_record(self):
        path = self.path_for(self.path)
        self.write(path, [])
        with open(path) as handle:
            lines = handle.read().splitlines()
        # the header alone: a run with no calibrated scans leaves no file, but a
        # calibrated run that found nothing leaves this
        self.assertEqual(lines, [",".join(self.fields)])
        self.assertIn("limit_jy_95", lines[0])


class CommandLineTests(unittest.TestCase):
    def test_quantity_flags(self):
        self.assertEqual(_resolve_quantity([]), (None, {}))
        self.assertEqual(_resolve_quantity(["3.5"]), (3.5, {}))
        default, by_path = _resolve_quantity(["3.5", "a.fil=9", "b.fil=12"])
        self.assertEqual(default, 3.5)
        self.assertEqual(by_path, {"a.fil": 9.0, "b.fil": 12.0})
        with self.assertRaises(SystemExit):
            _resolve_quantity(["abc"])

    def test_scan_sefd_prefers_what_is_given(self):
        self.assertEqual(_scan_sefd("a.fil", (5.0, {}), (None, {}), (None, {})), 5.0)
        self.assertEqual(
            _scan_sefd("a.fil", (None, {"a.fil": 7.0}), (None, {}), (None, {})), 7.0
        )
        # Tsys and gain stand in when no SEFD is given
        self.assertEqual(
            _scan_sefd("a.fil", (None, {}), (20.0, {}), (0.1, {})), 200.0
        )
        self.assertIsNone(_scan_sefd("a.fil", (None, {}), (None, {}), (None, {})))
        self.assertIsNone(_scan_sefd("a.fil", (None, {}), (20.0, {}), (None, {})))


class ReceiverCombinationTests(unittest.TestCase):
    """Two receivers, one source: the flux must agree, the better one must win."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def build_scan(self, name, mjd, *, amplitude, sefd_jy, nints=200, seed=1):
        topo_hz = (FCH1 + FOFF * SIGNAL_CHANNEL) * 1e6
        bary_hz = topo_hz / doppler_factor(OBSERVER, MJD0, SKY)
        seen_hz = bary_hz * doppler_factor(OBSERVER, mjd, SKY)
        drift = topocentric_drift(0.0, seen_hz, OBSERVER, mjd, SKY)
        path = f"{self.tmp.name}/{name}"
        simulate(
            path,
            nchans=NCHANS,
            nints=nints,
            tsamp=TSAMP,
            fch1=FCH1,
            foff=FOFF,
            signals=(InjectedSignal(seen_hz / 1e6, drift_hz_s=drift, amplitude=amplitude),),
            rfi=(RFILine(channel=3000, amplitude=3.0),),
            glitch=False,
            seed=seed,
            mjd=mjd,
        )
        config = SearchConfig(snr_threshold=10.0, max_drift_hz_s=1.4, max_hits=20)
        hits, meta = search_file(path, config, chan_chunk=2048)
        return Scan(path=path, kind="on", mjd=mjd, hits=hits, meta=meta, sefd_jy=sefd_jy)

    def test_same_source_on_two_receivers(self):
        # A 10 Jy receiver sees this source at per-sample SNR 20, a 100 Jy one at
        # 2 - the same flux, ten times less signal-to-noise.
        good = self.build_scan("good.fil", MJD0, amplitude=20.0, sefd_jy=10.0, seed=3)
        poor = self.build_scan(
            "poor.fil", MJD0 + 0.02, amplitude=2.0, sefd_jy=100.0, seed=4
        )
        result = link_run([good, poor], observer=OBSERVER, sky=SKY, min_scans=2)
        # the stationary carrier links too; pick the source by its frequency
        event = min(result.events, key=lambda e: abs(e.freq_mhz * 1e6 - BARY_HZ))
        self.assertLess(
            abs(event.freq_mhz * 1e6 - BARY_HZ), 2e3, "the source should link"
        )
        self.assertIsNotNone(event.flux, "both scans carry a SEFD")
        self.assertEqual(event.flux.n_scans, 2)

        from_good = event.flux.per_scan[good.path]
        from_poor = event.flux.per_scan[poor.path]
        # each scan's amplitude is in its own noise, so dividing out its own
        # SEFD must land on the same physical flux
        self.assertLess(abs(from_poor / from_good - 1.0), 0.4, (from_good, from_poor))

        # and being compared in the same units, the better receiver dominates
        self.assertLess(abs(event.flux.jy / from_good - 1.0), 0.15)
        self.assertLess(event.flux.uncertainty_jy, abs(from_good) * 0.05)

        # absolute scale, loosely: line width is assumed, not measured
        truth = 20.0 * 10.0 / (DF_HZ * TSAMP) ** 0.5
        self.assertLess(abs(event.flux.jy / truth - 1.0), 1.5, (event.flux.jy, truth))

    def test_single_file_flux_matches_the_run(self):
        # the same hit certified two ways: once by the single-file helper, once
        # inside a linked event, must give the same number
        scan = self.build_scan("solo.fil", MJD0, amplitude=20.0, sefd_jy=10.0, seed=7)
        other = self.build_scan(
            "solo2.fil", MJD0 + 0.02, amplitude=20.0, sefd_jy=10.0, seed=8
        )
        result = link_run([scan, other], observer=OBSERVER, sky=SKY, min_scans=2)
        event = min(result.events, key=lambda e: abs(e.freq_mhz * 1e6 - BARY_HZ))
        member = next(m for m in event.members if m.scan == scan.path)

        single = spectrum_flux(member.candidate.snr, scan.meta, scan.sefd_jy)
        self.assertIsNotNone(single)
        self.assertAlmostEqual(single[0], event.flux.per_scan[scan.path], places=6)
        self.assertGreater(single[1], 0.0)

        # no SEFD, no flux - never a default guess
        self.assertIsNone(spectrum_flux(member.candidate.snr, scan.meta, None))

    def test_uncalibrated_scan_leaves_flux_unset(self):
        good = self.build_scan("cal.fil", MJD0, amplitude=20.0, sefd_jy=10.0, seed=5)
        other = self.build_scan("uncal.fil", MJD0 + 0.02, amplitude=20.0, sefd_jy=None, seed=6)
        self.assertIsNone(other.sefd_jy)
        result = link_run([good, other], observer=OBSERVER, sky=SKY, min_scans=2)
        event = min(result.events, key=lambda e: abs(e.freq_mhz * 1e6 - BARY_HZ))
        self.assertEqual(event.flux.n_scans, 1)
        self.assertEqual(set(event.flux.per_scan), {good.path})


if __name__ == "__main__":
    unittest.main()
