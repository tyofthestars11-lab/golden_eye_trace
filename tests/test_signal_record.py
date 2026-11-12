"""Tests for the typed candidate record.

What matters here is that the record cannot be built into an unsupported claim:
a closed disposition vocabulary, a real date, a positive resolution, and a
positive significance somewhere. A record that breaks those is a bug in the
caller, and the constructor has to say so rather than store it.
"""

import unittest

from seti_search.signal_record import (
    DISPOSITIONS,
    SignalRecord,
    channel_rate,
)


def make(**overrides):
    base = dict(
        source_file="data/20240325.fits",
        observation_date="20240325",
        instrument="OVRO-LWA",
        frequency_hz=45.984284e6,
        drift_hz_per_s=41190.2,
        snr_global=953.98,
        snr_at_cell=30.60,
        local_scale=31.18,
        resolution_hz_per_s=392.76,
        channel_rate=0.1049,
        duration_s=256.0,
        on_source_status="on-source",
        transform_version="ovro-local-noise-v1",
        disposition="ASTROPHYSICAL",
    )
    base.update(overrides)
    return SignalRecord(**base)


class RecordConstructionTests(unittest.TestCase):
    def test_a_valid_record_keeps_both_significances(self):
        record = make()
        self.assertEqual(record.snr_global, 953.98)
        self.assertEqual(record.snr_at_cell, 30.60)
        # the ratio is the whole point of storing both
        self.assertLess(record.significance_ratio(), 0.05)

    def test_round_trip_through_a_dict(self):
        record = make()
        again = SignalRecord.from_dict(record.as_dict())
        self.assertEqual(again, record)

    def test_unknown_fields_are_rejected(self):
        data = make().as_dict()
        data["confidence"] = 0.99
        with self.assertRaises(ValueError):
            SignalRecord.from_dict(data)

    def test_disposition_vocabulary_is_closed(self):
        with self.assertRaises(ValueError):
            make(disposition="PROBABLY ALIENS")
        self.assertIn("NULL RESULT", DISPOSITIONS)
        self.assertIn("PIPELINE ARTIFACT", DISPOSITIONS)

    def test_observation_date_must_be_a_real_day(self):
        with self.assertRaises(ValueError):
            make(observation_date="20241332")
        with self.assertRaises(ValueError):
            make(observation_date="2024-03-25")

    def test_positive_scales_are_required(self):
        for bad in (0.0, -1.0, float("nan")):
            with self.assertRaises(ValueError):
                make(resolution_hz_per_s=bad)
            with self.assertRaises(ValueError):
                make(local_scale=bad)

    def test_a_record_with_no_significance_is_not_a_candidate(self):
        with self.assertRaises(ValueError):
            make(snr_global=0.0, snr_at_cell=-2.0)

    def test_on_source_status_is_constrained(self):
        with self.assertRaises(ValueError):
            make(on_source_status="somewhere")


class ChannelRateTests(unittest.TestCase):
    def test_rate_matches_the_survey_lane(self):
        rate = channel_rate(41190.2, 100547.19562243501, 0.256)
        self.assertAlmostEqual(rate, 0.10487, places=4)

    def test_a_rate_past_nyquist_folds_and_keeps_its_sign(self):
        # 0.4965 is inside; 0.5035 is its mirror and comes back negative
        inside = channel_rate(0.4965 * 100547.19562243501 / 0.256,
                              100547.19562243501, 0.256)
        outside = channel_rate(0.5035 * 100547.19562243501 / 0.256,
                               100547.19562243501, 0.256)
        self.assertAlmostEqual(inside, 0.4965, places=6)
        self.assertAlmostEqual(outside, -0.4965, places=6)

    def test_degenerate_scaling_is_rejected(self):
        with self.assertRaises(ValueError):
            channel_rate(100.0, 0.0, 0.256)
        with self.assertRaises(ValueError):
            channel_rate(100.0, 1e5, 0.0)


if __name__ == "__main__":
    unittest.main()
