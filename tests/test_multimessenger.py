"""Tests for the strict radio-time gate."""

from __future__ import annotations

import datetime as dt
import unittest
from types import SimpleNamespace

from seti_search.multimessenger import (
    Assessment,
    MultimessengerNotTestable,
    assess,
    enforce,
    require_testable,
)


def record(**metadata) -> SimpleNamespace:
    return SimpleNamespace(observation_date="20240322", duration_s=256.0,
                           frequency_hz=47.2416e6, metadata=metadata)


def coincidence(status: str) -> SimpleNamespace:
    return SimpleNamespace(status=status)


class MissingTimeTests(unittest.TestCase):
    def test_a_date_is_not_a_clock(self):
        result = assess(record())
        self.assertEqual(result.radio_time_status, "MISSING")
        self.assertEqual(result.multimessenger_status, "NOT TESTABLE")
        self.assertIsNone(result.window_s)
        self.assertTrue(any("a date is not a clock" in r for r in result.reasons))

    def test_the_missing_timestamp_is_priced(self):
        # a 24 h window is not free: the reason quotes what it would cost
        result = assess(record())
        joined = " | ".join(result.reasons)
        self.assertIn("chance probability", joined)
        self.assertIn("256 s", joined)
        self.assertIn("4.68", joined)

    def test_an_unusable_timestamp_says_why_it_is_unusable(self):
        result = assess(record(start_time_utc="2024-03-22T08:15:00"))
        self.assertEqual(result.radio_time_status, "MISSING")
        self.assertTrue(any("unusable" in r or "not a clock" in r
                            for r in result.reasons))

    def test_missing_time_cannot_be_exempted(self):
        with self.assertRaises(MultimessengerNotTestable):
            require_testable(assess(record()))


class GradedTests(unittest.TestCase):
    STAMP = "2024-03-22T08:15:00+00:00"

    def test_a_start_without_a_position_is_still_not_testable(self):
        result = assess(record(start_time_utc=self.STAMP))
        self.assertEqual(result.radio_time_status, "UTC_START")
        self.assertEqual(result.multimessenger_status, "NOT TESTABLE")
        self.assertEqual(result.window_s, 6 * 256.0)
        self.assertTrue(any("whole sky" in r for r in result.reasons))

    def test_a_start_and_a_position_are_testable(self):
        result = assess(record(start_time_utc=self.STAMP, ra_deg=100.0,
                               dec_deg=-20.0, radio_localisation_deg=0.5))
        self.assertEqual(result.radio_time_status, "FULL")
        self.assertEqual(result.multimessenger_status, "TESTABLE")
        self.assertTrue(result.testable)
        self.assertEqual(result.pointing_status, "PRESENT")

    def test_the_window_is_reported_in_seconds(self):
        result = assess(record(start_time_utc=self.STAMP, ra_deg=1.0,
                               dec_deg=2.0), sigma=3.0)
        self.assertAlmostEqual(result.window_s, 3.0 * 2.0 * 256.0)

    def test_a_missing_localisation_is_recorded_not_assumed(self):
        result = assess(record(start_time_utc=self.STAMP, ra_deg=1.0, dec_deg=2.0))
        self.assertTrue(any("no uncertainty" in r for r in result.reasons))

    def test_a_zero_duration_has_no_window(self):
        rec = record(start_time_utc=self.STAMP, ra_deg=1.0, dec_deg=2.0)
        rec.duration_s = 0.0
        result = assess(rec)
        self.assertEqual(result.multimessenger_status, "NOT TESTABLE")
        self.assertIsNone(result.window_s)


class ClosedVocabularyTests(unittest.TestCase):
    def test_bad_statuses_are_refused(self):
        with self.assertRaises(ValueError):
            Assessment(radio_time_status="SOMETIME",
                       multimessenger_status="NOT TESTABLE", window_s=None,
                       pointing_status="MISSING", reasons=("x",), start_utc=None)
        with self.assertRaises(ValueError):
            Assessment(radio_time_status="FULL",
                       multimessenger_status="PROBABLY", window_s=1.0,
                       pointing_status="PRESENT", reasons=("x",), start_utc=None)

    def test_a_testable_status_needs_its_reasons(self):
        with self.assertRaises(ValueError):
            Assessment(radio_time_status="FULL",
                       multimessenger_status="TESTABLE", window_s=1.0,
                       pointing_status="PRESENT", reasons=(), start_utc=None)


class EnforcementTests(unittest.TestCase):
    def test_a_coincidence_on_an_untestable_record_stops_the_run(self):
        gate = assess(record())
        with self.assertRaises(MultimessengerNotTestable) as caught:
            enforce(gate, coincidence("COINCIDENT"))
        self.assertIn("NOT TESTABLE", str(caught.exception))

    def test_a_missed_query_is_fine(self):
        gate = assess(record())
        enforce(gate, coincidence("UNINFORMATIVE"))
        enforce(gate, coincidence("NO EVENT IN WINDOW"))

    def test_a_testable_record_may_report_a_coincidence(self):
        gate = assess(record(start_time_utc="2024-03-22T08:15:00+00:00",
                             ra_deg=100.0, dec_deg=-20.0))
        enforce(gate, coincidence("COINCIDENT"))


if __name__ == "__main__":
    unittest.main()
