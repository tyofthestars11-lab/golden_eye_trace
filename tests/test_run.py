"""Linking across a run.

The point of a run rather than a single cadence is that the drift a telescope
measures wanders as Earth turns, while the drift that belongs to the source does
not. These tests build a source at a fixed barycentric drift whose topocentric
drift differs between scans, and check that the run links it - and that linking
without the barycentric conversion, on the same files, does not.
"""

import os
import tempfile
import unittest

from seti_search.cadence import SITES, SkyPosition, doppler_factor, topocentric_drift
from seti_search.filterbank import InjectedSignal, RFILine, simulate
from seti_search.run import SEEN_IN_OFF, Scan, link_run
from seti_search.search import SearchConfig, search_file

FCH1 = 1420.4
FOFF = -2.86e-6
NCHANS = 4096
NINTS = 200
TSAMP = 1.0

MJD0 = 58967.23
OBSERVER = SITES["ovro"]
SKY = SkyPosition(75.0, 20.0)
CARRIER_CHANNEL = 3000
SIGNAL_CHANNEL = 1200

# the source's frequency as seen at MJD0; its barycentric frequency is fixed
TOPO_HZ_AT_MJD0 = (FCH1 + FOFF * SIGNAL_CHANNEL) * 1e6


def barycentric_hz() -> float:
    return TOPO_HZ_AT_MJD0 / doppler_factor(OBSERVER, MJD0, SKY)


def run_search(path: str):
    config = SearchConfig(snr_threshold=10.0, max_drift_hz_s=1.4, max_hits=20)
    return search_file(path, config, chan_chunk=2048)


def freq_tolerance_hz(span_days: float) -> float:
    """How well a referred frequency can be known across a run.

    Referring every hit to one epoch multiplies the drift error by the elapsed
    time, so a 12 h run with drifts good to 0.005 Hz/s knows the frequency to a
    few hundred Hz, not a few channels. Asserting tighter would be asserting
    something the method cannot deliver.
    """
    return max(4 * abs(FOFF) * 1e6, 0.02 * span_days * 86400.0)


class LinkTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def make_scan(
        self,
        name: str,
        kind: str,
        mjd: float,
        *,
        source: bool,
        seed: int,
        nints: int = NINTS,
        threshold: float = 10.0,
    ):
        path = os.path.join(self.tmp.name, name)
        signals = ()
        if source:
            topocentric = barycentric_hz() * doppler_factor(OBSERVER, mjd, SKY)
            drift = topocentric_drift(0.0, topocentric, OBSERVER, mjd, SKY)
            signals = (InjectedSignal(topocentric / 1e6, drift_hz_s=drift, amplitude=2.0),)
        simulate(
            path,
            nchans=NCHANS,
            nints=nints,
            tsamp=TSAMP,
            fch1=FCH1,
            foff=FOFF,
            signals=signals,
            rfi=(RFILine(channel=CARRIER_CHANNEL, amplitude=3.0),),
            glitch=False,
            seed=seed,
        )
        config = SearchConfig(
            snr_threshold=threshold, max_drift_hz_s=1.4, max_hits=20
        )
        hits, meta = search_file(path, config, chan_chunk=2048)
        return Scan(path=path, kind=kind, mjd=mjd, hits=hits, meta=meta)

    def source_event(self, result, span_days: float = 0.05):
        """The event at the source's barycentric frequency.

        Every scan also contains the stationary carrier, which links in its own
        right, so a bare count of events is not the thing to assert on.
        """
        matches = [
            e
            for e in result.events
            if abs(e.freq_hz - barycentric_hz()) <= freq_tolerance_hz(span_days)
        ]
        self.assertEqual(
            len(matches), 1, msg=[(e.freq_mhz, e.n_scans) for e in result.events]
        )
        return matches[0]

    def test_links_a_source_whose_topocentric_drift_wanders(self):
        # six hours apart: the drift a telescope measures moves by ~0.2 Hz/s
        scans = [
            self.make_scan("on1.fil", "on", MJD0, source=True, seed=21),
            self.make_scan("on2.fil", "on", MJD0 + 0.25, source=True, seed=22),
            self.make_scan("on3.fil", "on", MJD0 + 0.5, source=True, seed=23),
        ]
        result = link_run(scans, observer=OBSERVER, sky=SKY)
        self.assertTrue(result.predicted)

        span = 0.5  # days, first to last ON
        sources = [
            e
            for e in result.events
            if abs(e.freq_hz - barycentric_hz()) <= freq_tolerance_hz(span)
        ]
        self.assertEqual(
            len(sources), 1, msg=[(e.freq_mhz, e.n_scans) for e in result.events]
        )
        event = sources[0]
        self.assertEqual(event.n_scans, 3)
        self.assertAlmostEqual(event.drift_hz_s, 0.0, delta=0.05)

        # The other event is the stationary carrier, linked in the frame where it
        # *is* consistent - the telescope's own - so that it can be given a
        # verdict rather than quietly dropped.
        others = [e for e in result.events if e is not event]
        self.assertEqual(len(others), 1, msg=[(e.freq_mhz, e.n_scans) for e in result.events])
        carrier_mhz = FCH1 + FOFF * CARRIER_CHANNEL
        for member in others[0].members:
            self.assertLess(
                abs(member.candidate.freq_mhz - carrier_mhz),
                4 * abs(FOFF),
                msg="an event appeared that is neither the source nor the carrier",
            )

        # and the link really did need the conversion: the raw drifts disagree
        raw = [m.candidate.drift_hz_s for m in event.members]
        self.assertGreater(max(raw) - min(raw), 0.1)

        # Without a site and sky position the same files link only what does not
        # move: the stationary carrier. The source, whose measured drift wanders
        # by more than the tolerance, is never linked - which is precisely why a
        # run needs the barycentric conversion.
        blind = link_run(scans)
        self.assertFalse(blind.predicted)
        self.assertTrue(any("barycentric conversion unavailable" in n for n in blind.notes))
        carrier_mhz = FCH1 + FOFF * CARRIER_CHANNEL
        self.assertEqual(len(blind.events), 1, msg=[(e.freq_mhz, e.n_scans) for e in blind.events])
        for member in blind.events[0].members:
            self.assertLess(
                abs(member.candidate.freq_mhz - carrier_mhz), 4 * abs(FOFF),
                msg="a blind run linked something other than the stationary carrier",
            )

    def test_an_event_that_survives_pointing_away_is_rejected(self):
        scans = [
            self.make_scan("on1.fil", "on", MJD0, source=True, seed=31),
            self.make_scan("off1.fil", "off", MJD0 + 0.02, source=False, seed=32),
            self.make_scan("on2.fil", "on", MJD0 + 0.04, source=True, seed=33),
        ]
        result = link_run(scans, observer=OBSERVER, sky=SKY)
        statuses = sorted(e.status for e in result.events)
        self.assertIn(SEEN_IN_OFF, statuses)
        self.assertEqual(len(result.survivors), 1)

        survivor = result.survivors[0]
        self.assertLess(
            abs(survivor.freq_hz - barycentric_hz()), freq_tolerance_hz(0.04)
        )
        self.assertEqual(survivor.off_matches, [])

        rejected = [e for e in result.events if not e.survived][0]
        self.assertTrue(rejected.off_matches)
        self.assertEqual(rejected.status, SEEN_IN_OFF)

    def test_a_stationary_carrier_links_in_the_telescope_frame_and_is_rejected(self):
        # Six hours between ON scans breaks the carrier's *implied* barycentric
        # drift, so it has to be linked in the topocentric frame to get a verdict
        # at all. Silence would be the wrong answer: it is interference, and the
        # report should say so.
        scans = [
            self.make_scan("on1.fil", "on", MJD0, source=False, seed=61),
            self.make_scan("on2.fil", "on", MJD0 + 0.25, source=False, seed=62),
            self.make_scan("off1.fil", "off", MJD0 + 0.125, source=False, seed=63),
        ]
        result = link_run(scans, observer=OBSERVER, sky=SKY)
        self.assertEqual(
            len(result.events), 1, msg=[(e.freq_mhz, e.status) for e in result.events]
        )
        event = result.events[0]
        self.assertEqual(event.n_scans, 2)
        self.assertEqual(event.status, SEEN_IN_OFF)
        self.assertEqual(len(result.survivors), 0)

        carrier_mhz = FCH1 + FOFF * CARRIER_CHANNEL
        for member in event.members:
            self.assertLess(abs(member.candidate.freq_mhz - carrier_mhz), 4 * abs(FOFF))
        # the members disagree as soon as they are referred to the barycentric
        # frame, which is why the topocentric link was the one that found them
        implied = [m.drift_hz_s for m in event.members]
        self.assertGreater(abs(implied[0] - implied[1]), 0.05)

    def test_unequal_integration_times_are_weighted_by_sensitivity(self):
        # four times the integrations is twice the SNR for the same source
        scans = [
            self.make_scan("on1.fil", "on", MJD0, source=True, seed=71, nints=200),
            self.make_scan("on2.fil", "on", MJD0 + 0.02, source=True, seed=72, nints=800),
        ]
        result = link_run(scans, observer=OBSERVER, sky=SKY)
        event = self.source_event(result, 0.02)
        self.assertEqual(event.n_scans, 2)

        snrs = sorted(m.candidate.snr for m in event.members)
        self.assertGreater(snrs[1] / snrs[0], 1.6, msg="the deeper scan must score higher")

        # divided by its own sensitivity, the two scans agree on the source
        amplitudes = [m.amplitude for m in event.members]
        self.assertLess(
            abs(amplitudes[0] - amplitudes[1]) / max(amplitudes),
            0.25,
            msg=amplitudes,
        )

        # combining scans of unequal depth must beat the deeper one but not
        # simply add them
        self.assertGreater(event.score, max(snrs))
        self.assertLess(event.score, snrs[0] + snrs[1])

    def test_a_deep_scan_that_should_have_seen_the_event_is_flagged(self):
        scans = [
            self.make_scan("on1.fil", "on", MJD0, source=True, seed=81, nints=800),
            self.make_scan("on2.fil", "on", MJD0 + 0.02, source=True, seed=82, nints=800),
            self.make_scan("on3.fil", "on", MJD0 + 0.04, source=False, seed=83, nints=800),
        ]
        result = link_run(scans, observer=OBSERVER, sky=SKY)
        event = self.source_event(result, 0.04)
        self.assertEqual(event.n_scans, 2)
        self.assertIn(scans[2].path, event.missed_scans)
        self.assertGreater(event.expected_snr[scans[2].path], 10.0)

    def test_a_threshold_above_the_expected_snr_does_not_blame_a_scan(self):
        scans = [
            self.make_scan("on1.fil", "on", MJD0, source=True, seed=84, nints=800),
            self.make_scan("on2.fil", "on", MJD0 + 0.02, source=True, seed=85, nints=800),
            self.make_scan(
                "on3.fil", "on", MJD0 + 0.04, source=False, seed=86, nints=800,
                threshold=400.0,
            ),
        ]
        result = link_run(scans, observer=OBSERVER, sky=SKY)
        event = self.source_event(result, 0.04)
        # the scan was cut at 400 sigma, so its silence says nothing at the
        # amplitude this event implies
        self.assertEqual(event.missed_scans, [])
        self.assertIn(scans[2].path, event.expected_snr)
        self.assertLess(event.expected_snr[scans[2].path], 400.0)

    def test_a_single_on_scan_produces_no_event(self):
        scans = [self.make_scan("on1.fil", "on", MJD0, source=True, seed=41)]
        result = link_run(scans, observer=OBSERVER, sky=SKY)
        self.assertEqual(result.events, [], msg="one scan cannot make an event")


if __name__ == "__main__":
    unittest.main()
