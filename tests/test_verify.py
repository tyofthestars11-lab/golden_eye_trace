"""End-to-end ON/OFF tests.

Two files, same band, one hour apart. The ON scan points at the target and
contains a source at a fixed sky position plus a stationary carrier; the OFF
scan points away, so it contains the carrier and nothing else. A verifier that
works must keep the source and reject the carrier, and it must have got there by
predicting where the source would have been rather than by luck.
"""

import os
import tempfile
import unittest

from seti_search.cadence import (
    SITES,
    SkyPosition,
    predict_observation,
    topocentric_drift,
)
from seti_search.filterbank import InjectedSignal, RFILine, simulate
from seti_search.search import SearchConfig, search_file
from seti_search.verify import verify_cadence

FCH1 = 1420.4
FOFF = -2.86e-6
NCHANS = 4096
NINTS = 200
TSAMP = 1.0

MJD_ON = 58967.23
MJD_OFF = MJD_ON + 1.0 / 24.0  # one hour later
OBSERVER = SITES["ovro"]
SKY = SkyPosition(75.0, 20.0)
SIGNAL_CHANNEL = 1200
CARRIER_CHANNEL = 3000


def chan_freq(channel: int) -> float:
    return FCH1 + FOFF * channel


def run_search(path: str):
    config = SearchConfig(snr_threshold=10.0, max_drift_hz_s=1.4, max_hits=20)
    return search_file(path, config, chan_chunk=2048)


class CadenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def build_pair(self):
        """A sky source in ON only, and a carrier in both scans."""
        freq = chan_freq(SIGNAL_CHANNEL)
        drift_on = topocentric_drift(0.0, freq * 1e6, OBSERVER, MJD_ON, SKY)
        freq_off_hz, drift_off = predict_observation(
            freq * 1e6, drift_on, OBSERVER, MJD_ON, MJD_OFF, SKY
        )
        carrier = (RFILine(channel=CARRIER_CHANNEL, amplitude=3.0),)

        on_path = os.path.join(self.tmp.name, "on.fil")
        off_path = os.path.join(self.tmp.name, "off.fil")
        simulate(
            on_path,
            nchans=NCHANS,
            nints=NINTS,
            tsamp=TSAMP,
            fch1=FCH1,
            foff=FOFF,
            signals=(InjectedSignal(freq, drift_hz_s=drift_on, amplitude=2.0),),
            rfi=carrier,
            glitch=False,
            seed=11,
        )
        simulate(
            off_path,
            nchans=NCHANS,
            nints=NINTS,
            tsamp=TSAMP,
            fch1=FCH1,
            foff=FOFF,
            signals=(),
            rfi=carrier,
            glitch=False,
            seed=12,
        )
        return on_path, off_path, freq, drift_on, freq_off_hz / 1e6, drift_off

    def verify_default(self, on_hits, off_hits, on_meta, off_meta, **overrides):
        return verify_cadence(
            on_hits,
            off_hits,
            on_meta=on_meta,
            off_meta=off_meta,
            observer=OBSERVER,
            sky=SKY,
            on_mjd=overrides.pop("on_mjd", MJD_ON),
            off_mjd=overrides.pop("off_mjd", MJD_OFF),
            **overrides,
        )

    def test_sky_signal_survives_and_interference_does_not(self):
        on_path, off_path, freq, drift_on, freq_off, drift_off = self.build_pair()
        on_hits, on_meta = run_search(on_path)
        off_hits, off_meta = run_search(off_path)

        self.assertEqual(len(on_hits), 2, msg=[(h.freq_mhz, h.drift_hz_s) for h in on_hits])
        self.assertEqual(len(off_hits), 1, msg=[(h.freq_mhz, h.drift_hz_s) for h in off_hits])

        result = self.verify_default(on_hits, off_hits, on_meta, off_meta)
        self.assertTrue(result.predicted)

        confirmed, rejected = result.confirmed, result.rejected
        self.assertEqual(
            len(confirmed), 1, msg=[(v.on.freq_mhz, v.reason) for v in result.verdicts]
        )
        self.assertEqual(
            len(rejected), 1, msg=[(v.on.freq_mhz, v.reason) for v in result.verdicts]
        )

        survivor = confirmed[0]
        self.assertLess(abs(survivor.on.freq_mhz - freq), 2 * abs(FOFF))

        # The physics, not just the bookkeeping. Carrying the measured ON
        # observation forward an hour has to land where the model puts the
        # source, in both frequency and drift.
        self.assertIsNotNone(survivor.predicted_drift_hz_s)
        self.assertAlmostEqual(survivor.predicted_drift_hz_s, drift_off, delta=0.05)
        self.assertAlmostEqual(survivor.predicted_freq_mhz, freq_off, delta=10e-6)
        self.assertAlmostEqual(
            survivor.drift_discrimination_hz_s, abs(drift_off - drift_on), delta=0.05
        )
        # the source drifted the other way from the carrier, so a verifier that
        # ignored drift could not have got here by accident
        self.assertGreater(abs(drift_on), 0.02)

        carrier = rejected[0]
        self.assertLess(abs(carrier.on.freq_mhz - chan_freq(CARRIER_CHANNEL)), 2 * abs(FOFF))
        self.assertIsNotNone(carrier.off_hit)
        self.assertLess(carrier.on.snr, survivor.on.snr)

    def test_without_epoch_or_coordinates_only_frequency_is_compared(self):
        on_path, off_path, _, _, _, _ = self.build_pair()
        on_hits, on_meta = run_search(on_path)
        off_hits, off_meta = run_search(off_path)

        result = verify_cadence(on_hits, off_hits, on_meta=on_meta, off_meta=off_meta)
        self.assertFalse(result.predicted)
        self.assertTrue(any("drift prediction unavailable" in note for note in result.notes))
        for verdict in result.verdicts:
            self.assertIsNone(verdict.predicted_drift_hz_s)

        # The source is far enough from the carrier that the frequency-only
        # comparison still separates them; the carrier has not moved at all.
        self.assertEqual(len(result.confirmed), 1)
        self.assertEqual(len(result.rejected), 1)
        self.assertIsNotNone(result.rejected[0].off_hit)

    def test_short_cadence_is_reported_as_weak(self):
        on_path, off_path, _, _, _, _ = self.build_pair()
        on_hits, on_meta = run_search(on_path)
        off_hits, off_meta = run_search(off_path)

        result = self.verify_default(
            on_hits, off_hits, on_meta, off_meta, on_mjd=MJD_ON, off_mjd=MJD_ON
        )
        self.assertTrue(any("longer cadence" in note for note in result.notes))
        for verdict in result.verdicts:
            self.assertLess(verdict.drift_discrimination_hz_s, 1e-6)

    def test_frequency_window_is_what_decides(self):
        on_path, off_path, _, _, _, _ = self.build_pair()
        on_hits, on_meta = run_search(on_path)
        off_hits, off_meta = run_search(off_path)

        tight = self.verify_default(
            on_hits, off_hits, on_meta, off_meta, freq_tolerance_hz=1.0
        )
        wide = self.verify_default(
            on_hits, off_hits, on_meta, off_meta, freq_tolerance_hz=1e6
        )
        self.assertEqual(len(wide.confirmed), 0, msg="a huge window must match everything")
        self.assertGreater(len(tight.confirmed), len(wide.confirmed))


if __name__ == "__main__":
    unittest.main()
