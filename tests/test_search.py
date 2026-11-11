"""Injection-recovery tests: does the search find what we buried, with the right drift?"""

import os
import tempfile
import unittest

import numpy as np

from seti_search.filterbank import InjectedSignal, RFILine, simulate
from seti_search.search import SearchConfig, calibrate_map, map_tail, search_file
from seti_search.trials import drift_trials

FCH1 = 1420.4
FOFF = -2.86e-6
NCHANS = 8192
NINTS = 200
TSAMP = 1.0


def chan_freq(channel: int) -> float:
    return FCH1 + FOFF * channel


class CalibrationTests(unittest.TestCase):
    """The printed sigma has to be the sigma of this file's own map."""

    def test_map_with_inflated_spread_is_rescaled(self):
        rng = np.random.default_rng(0)
        snr_map = rng.standard_normal((3, 2000)) * 1.4
        calibrated, spread = calibrate_map(snr_map)
        self.assertAlmostEqual(spread, 1.4, delta=0.06)
        self.assertAlmostEqual(float(calibrated.std()), 1.0, places=6)
        self.assertAlmostEqual(
            float(calibrated.max()), float(snr_map.max()) / spread, places=5
        )

    def test_blanked_channels_do_not_bias_the_calibration(self):
        rng = np.random.default_rng(1)
        snr_map = rng.standard_normal((2, 1000)) * 1.3
        live = np.ones(1000, dtype=bool)
        live[:500] = False
        snr_map[:, :500] = 0.0
        _, spread = calibrate_map(snr_map, live)
        self.assertAlmostEqual(spread, 1.3, delta=0.09)

    def test_a_real_signal_does_not_set_the_noise_level(self):
        """A 200-sigma peak must not drag the calibration up with it."""
        rng = np.random.default_rng(2)
        snr_map = rng.standard_normal((2, 2048)) * 1.2
        snr_map[0, 700] = 240.0
        _, spread = calibrate_map(snr_map)
        self.assertAlmostEqual(spread, 1.2, delta=0.06)

    def test_too_small_a_map_is_left_alone(self):
        snr_map = np.ones((1, 4), dtype=np.float32) * 3.0
        calibrated, spread = calibrate_map(snr_map)
        self.assertEqual(spread, 1.0)
        self.assertTrue(np.array_equal(calibrated, snr_map))


class InjectionRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def run_search(self, path: str, **overrides):
        config = SearchConfig(
            snr_threshold=overrides.pop("snr_threshold", 10.0),
            max_drift_hz_s=overrides.pop("max_drift_hz_s", 1.4),
            max_hits=overrides.pop("max_hits", 20),
        )
        return search_file(path, config, chan_chunk=4096)

    def test_recovers_drifting_signals(self):
        signals = (
            InjectedSignal(chan_freq(1000), drift_hz_s=0.5, amplitude=1.0),
            InjectedSignal(chan_freq(5000), drift_hz_s=-1.0, amplitude=0.9),
            InjectedSignal(chan_freq(7000), drift_hz_s=0.05, amplitude=1.2),
        )
        path = os.path.join(self.tmp.name, "injected.fil")
        simulate(
            path,
            nchans=NCHANS,
            nints=NINTS,
            tsamp=TSAMP,
            fch1=FCH1,
            foff=FOFF,
            signals=signals,
            rfi=(RFILine(channel=3000, amplitude=40.0),),
            seed=1,
        )

        hits, meta = self.run_search(path)

        self.assertEqual(
            len(hits), 3,
            msg="expected exactly the three injected signals, got "
            + ", ".join(f"{h.freq_mhz:.6f}@{h.drift_hz_s:+.3f}" for h in hits),
        )

        channel_width = abs(FOFF)
        for signal in signals:
            match = min(hits, key=lambda h: abs(h.freq_mhz - signal.freq_mhz))
            self.assertLess(
                abs(match.freq_mhz - signal.freq_mhz), 2 * channel_width,
                msg=f"frequency off for {signal.freq_mhz:.6f} MHz",
            )
            # the sign of the drift is the whole point: a sign error still
            # produces a strong peak at the wrong rate. Tolerance is about two
            # drift steps.
            self.assertAlmostEqual(match.drift_hz_s, signal.drift_hz_s, delta=0.05)
            self.assertGreater(match.snr, 10.0)

        # the 40x stationary carrier must have been blanked, not reported
        self.assertGreaterEqual(meta["blanked_channels"], 1)
        self.assertGreaterEqual(meta["blanked_rows"], 1)

    def test_weak_stationary_carrier_is_reported(self):
        """A moderate stationary tone is a candidate, not grounds for deletion."""
        path = os.path.join(self.tmp.name, "weak_carrier.fil")
        simulate(
            path,
            nchans=NCHANS,
            nints=NINTS,
            tsamp=TSAMP,
            fch1=FCH1,
            foff=FOFF,
            signals=(InjectedSignal(chan_freq(4000), drift_hz_s=0.0, amplitude=3.0),),
            rfi=(),
            glitch=False,
            seed=2,
        )
        hits, meta = self.run_search(path)
        self.assertEqual(meta["blanked_channels"], 0)
        self.assertEqual(
            len(hits), 1, msg=[(h.freq_mhz, h.drift_hz_s, h.snr) for h in hits]
        )
        self.assertAlmostEqual(hits[0].drift_hz_s, 0.0, delta=0.03)
        self.assertLess(abs(hits[0].freq_mhz - chan_freq(4000)), 2 * abs(FOFF))

    def test_overwhelming_carrier_is_excised(self):
        """A tone hundreds of sigma above the floor is terrestrial and is blanked."""
        path = os.path.join(self.tmp.name, "loud_carrier.fil")
        simulate(
            path,
            nchans=NCHANS,
            nints=NINTS,
            tsamp=TSAMP,
            fch1=FCH1,
            foff=FOFF,
            signals=(),
            rfi=(RFILine(channel=4000, amplitude=200.0),),
            glitch=False,
            seed=5,
        )
        hits, meta = self.run_search(path)
        self.assertEqual(hits, [])
        self.assertGreaterEqual(meta["blanked_channels"], 1)

    def test_noise_only_produces_no_candidates(self):
        path = os.path.join(self.tmp.name, "noise.fil")
        simulate(
            path,
            nchans=NCHANS,
            nints=NINTS,
            tsamp=TSAMP,
            fch1=FCH1,
            foff=FOFF,
            signals=(),
            rfi=(),
            glitch=False,
            seed=3,
        )
        hits, _ = self.run_search(path, snr_threshold=8.0)
        self.assertEqual(
            hits, [],
            msg="noise-only data produced candidates: "
            + ", ".join(f"{h.snr:.1f} sigma" for h in hits),
        )

    def test_wide_band_with_few_integrations_keeps_its_rows(self):
        """A wide band with few integrations must not have every row excised.

        The row baseline has to be the median of the row medians: taking the
        median of the whole array instead leaves every row low by the same few
        per cent, which exceeds the row threshold once n_f makes the row noise
        small, so a real fine-resolution file would be blanked end to end.
        """
        path = os.path.join(self.tmp.name, "wide_few_ints.fil")
        simulate(
            path,
            nchans=32768,
            nints=16,
            tsamp=TSAMP,
            fch1=FCH1,
            foff=FOFF,
            signals=(InjectedSignal(chan_freq(12000), drift_hz_s=0.0, amplitude=10.0),),
            rfi=(),
            glitch=True,
            seed=6,
        )
        # one chunk, so blanked_rows counts rows of this band rather than a sum
        # over chunk-sized pieces of it
        hits, meta = search_file(
            path,
            SearchConfig(snr_threshold=8.0, max_drift_hz_s=1.4, max_hits=20),
            chan_chunk=32768,
            max_ints=None,
        )
        self.assertEqual(
            meta["blanked_rows"], 2,
            msg=f"expected only the two glitch rows blanked, got {meta['blanked_rows']} of 16",
        )
        self.assertEqual(
            len(hits), 1, msg=[(h.freq_mhz, h.snr) for h in hits]
        )

    def test_drift_beyond_range_is_not_reported_at_wrong_rate(self):
        """A fast signal outside the searched range must not appear at a bogus drift."""
        path = os.path.join(self.tmp.name, "fast.fil")
        simulate(
            path,
            nchans=NCHANS,
            nints=NINTS,
            tsamp=TSAMP,
            fch1=FCH1,
            foff=FOFF,
            signals=(InjectedSignal(chan_freq(4000), drift_hz_s=0.9, amplitude=1.5),),
            rfi=(),
            glitch=False,
            seed=4,
        )
        # search a range that excludes 0.9 Hz/s; the signal smears away entirely
        hits, _ = self.run_search(path, max_drift_hz_s=0.2)
        for hit in hits:
            self.assertLess(abs(hit.drift_hz_s), 0.2 + 1e-6)


class MapTailTests(unittest.TestCase):
    """The bulk correction fixes the width; the tail still has to be counted."""

    def test_gaussian_map_tail_matches_the_model(self):
        rng = np.random.default_rng(7)
        snr_map = rng.standard_normal((200, 2000))
        tail = map_tail(snr_map)
        self.assertEqual(tail[4.0]["cells"], 400_000)
        # Poisson scatter on the ~13 cells a normal puts above 4 sigma, so this
        # is a loose band: it only rules out the ratio being wildly miscomputed
        self.assertGreater(tail[4.0]["ratio"], 0.35)
        self.assertLess(tail[4.0]["ratio"], 2.0)

    def test_heavy_tail_is_reported_as_excess(self):
        rng = np.random.default_rng(8)
        snr_map = rng.standard_normal((100, 2000))
        snr_map.flat[:2000] = 6.0
        tail = map_tail(snr_map)
        # the injected cells, plus the handful the untouched bulk contributes
        self.assertGreaterEqual(tail[4.0]["observed"], 2000)
        self.assertGreater(tail[4.0]["ratio"], 20.0)

    def test_blanked_channels_are_not_counted_as_tail_cells(self):
        snr_map = np.zeros((4, 1000), dtype=np.float32)
        live = np.ones(1000, dtype=bool)
        live[:600] = False
        snr_map[:, :600] = 9.0
        tail = map_tail(snr_map, live)
        self.assertEqual(tail[4.0]["cells"], 1600)
        self.assertEqual(tail[4.0]["observed"], 0)

    def test_search_carries_the_tail_in_its_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "tail.fil")
            simulate(
                path,
                nchans=4096,
                nints=16,
                tsamp=TSAMP,
                fch1=FCH1,
                foff=FOFF,
                signals=(),
                rfi=(),
                glitch=False,
                seed=9,
            )
            _hits, meta = search_file(
                path, SearchConfig(snr_threshold=5.0, max_drift_hz_s=1.4),
                chan_chunk=4096, max_ints=None,
            )
        self.assertIn("snr_tail", meta)
        self.assertEqual(sorted(meta["snr_tail"]), [4.0, 5.0])
        # a map is one cell per drift trial per live channel, not per integration
        n_drifts = drift_trials(
            n_ints=meta["n_ints"],
            tsamp_s=TSAMP,
            df_mhz=abs(FOFF),
            boxcar=3,
            max_drift_hz_s=1.4,
        )
        self.assertEqual(
            meta["snr_tail"][4.0]["cells"],
            n_drifts * (4096 - meta["blanked_channels"]),
        )


if __name__ == "__main__":
    unittest.main()
