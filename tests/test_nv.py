"""Tests for the NV-centre lane: the Zeeman reach, the schema, and the detraction rule."""

from __future__ import annotations

import datetime as dt
import unittest

from seti_search.nv import (
    NVAsEvidence,
    NVSensorRecord,
    band_covered,
    field_from_magnetic,
    guard_support,
    magnetic_from_field,
    significance,
    verdict,
    zeeman_bias_t,
)

UTC = dt.timezone.utc


def sensor(**overrides) -> NVSensorRecord:
    base = dict(
        sensor_id="NV-1",
        instrument="ensemble diamond",
        sensing_mode="rf-dressed",
        band_low_hz=13.4e6,
        band_high_hz=86.9e6,
        sensitivity_t_per_sqrt_hz=1.0e-9,
        bandwidth_hz=7.5e6,
        bias_field_t=0.1,
        start_time_utc=dt.datetime(2024, 3, 22, 8, 15, 0, tzinfo=UTC),
        time_uncertainty_s=1.0,
        axis_tilt_deg=54.7,
        data_release="test-release",
        public_source_url="https://example.invalid/nv",
    )
    base.update(overrides)
    return NVSensorRecord(**base)


class ZeemanTests(unittest.TestCase):
    def test_fifty_megahertz_needs_about_a_tenth_of_a_tesla(self):
        self.assertAlmostEqual(zeeman_bias_t(50e6), 0.1007, places=3)

    def test_the_band_edges(self):
        self.assertGreater(zeeman_bias_t(13.4e6), zeeman_bias_t(86.9e6))
        self.assertAlmostEqual(zeeman_bias_t(13.4e6), 0.102, places=3)
        self.assertAlmostEqual(zeeman_bias_t(86.9e6), 0.0994, places=3)

    def test_a_frequency_above_the_splitting_is_unreachable(self):
        with self.assertRaises(ValueError):
            zeeman_bias_t(2.87e9)
        with self.assertRaises(ValueError):
            zeeman_bias_t(0.0)


class FieldTests(unittest.TestCase):
    def test_one_picotesla_is_three_tenths_of_a_millivolt_per_metre(self):
        # 1 pT * c = 0.2998 mV/m; the factor is c, not 3e8 rounded
        self.assertAlmostEqual(field_from_magnetic(1e-12) * 1e3, 0.29979,
                               places=4)

    def test_the_conversion_round_trips(self):
        self.assertAlmostEqual(magnetic_from_field(1.0) * 299_792_458.0, 1.0,
                               places=12)

    def test_band_coverage(self):
        self.assertAlmostEqual(band_covered(13.4e6, 86.9e6, 0.0, 1e9), 1.0)
        self.assertAlmostEqual(band_covered(13.4e6, 86.9e6, 0.0, 1e6), 0.0)
        half = band_covered(0.0, 100.0, 0.0, 50.0)
        self.assertAlmostEqual(half, 0.5)

    def test_an_empty_band_is_refused(self):
        with self.assertRaises(ValueError):
            band_covered(10.0, 10.0, 0.0, 100.0)


class SignificanceTests(unittest.TestCase):
    def test_the_bandwidth_is_in_the_denominator(self):
        # eta sqrt(df) = 1e-13 * 1e3 = 1e-10 T, so 2e-10 T is 2 sigma and
        # 2e-13 T is 0.002 sigma - the sqrt really does bite
        self.assertAlmostEqual(significance(2e-10, 1e-13, 1e6), 2.0)
        self.assertAlmostEqual(significance(2e-13, 1e-13, 1e6), 0.002)

    def test_a_wider_bandwidth_costs_signal_to_noise(self):
        self.assertLess(significance(1e-11, 1e-13, 1e8),
                        significance(1e-11, 1e-13, 1e6))

    def test_bad_inputs_are_refused(self):
        with self.assertRaises(ValueError):
            significance(-1.0, 1e-13, 1e6)
        with self.assertRaises(ValueError):
            significance(1.0, 0.0, 1e6)


class SchemaTests(unittest.TestCase):
    def test_a_valid_record_builds(self):
        self.assertEqual(sensor().sensor_id, "NV-1")

    def test_a_naive_start_is_refused(self):
        with self.assertRaises(ValueError):
            sensor(start_time_utc=dt.datetime(2024, 3, 22, 8, 15))

    def test_the_mode_vocabulary_is_closed(self):
        with self.assertRaises(ValueError):
            sensor(sensing_mode="ouija")

    def test_a_reversed_band_is_refused(self):
        with self.assertRaises(ValueError):
            sensor(band_low_hz=86.9e6, band_high_hz=13.4e6)

    def test_sensitivity_and_bandwidth_are_both_required(self):
        with self.assertRaises(ValueError):
            sensor(sensitivity_t_per_sqrt_hz=0.0)
        with self.assertRaises(ValueError):
            sensor(bandwidth_hz=0.0)

    def test_provenance_is_required(self):
        with self.assertRaises(ValueError):
            sensor(data_release="")
        with self.assertRaises(ValueError):
            sensor(public_source_url="")

    def test_the_e_field_sensitivity_is_reported(self):
        record = sensor(sensitivity_t_per_sqrt_hz=1e-12)
        self.assertAlmostEqual(record.field_sensitivity_v_m * 1e3, 0.29979,
                               places=4)

    def test_covers(self):
        record = sensor()
        self.assertTrue(record.covers(47.2416e6))
        self.assertFalse(record.covers(200e6))


class VerdictTests(unittest.TestCase):
    STAMP = dt.datetime(2024, 3, 22, 8, 15, 0, tzinfo=UTC)

    def test_no_sensor_is_not_testable(self):
        result = verdict(None, frequency_hz=47e6, radio_time_utc=self.STAMP)
        self.assertEqual(result["verdict"], "NOT TESTABLE")

    def test_no_radio_timestamp_is_not_testable(self):
        result = verdict(sensor(), frequency_hz=47e6, radio_time_utc=None)
        self.assertEqual(result["verdict"], "NOT TESTABLE")
        self.assertIn("timestamp", result["reasons"][0])

    def test_out_of_band_is_not_testable(self):
        result = verdict(sensor(band_low_hz=1e6, band_high_hz=10e6),
                         frequency_hz=47e6, radio_time_utc=self.STAMP)
        self.assertEqual(result["verdict"], "NOT TESTABLE")
        self.assertIn("MHz", result["reasons"][0])

    def test_a_sensor_acquisition_far_in_time_is_not_testable(self):
        later = dt.datetime(2024, 3, 23, 8, 15, 0, tzinfo=UTC)
        result = verdict(sensor(start_time_utc=later), frequency_hz=47e6,
                         radio_time_utc=self.STAMP)
        self.assertEqual(result["verdict"], "NOT TESTABLE")

    def test_a_null_in_window_is_no_information_not_a_negative(self):
        result = verdict(sensor(), frequency_hz=47e6, radio_time_utc=self.STAMP,
                         b_amplitude_t=None)
        self.assertEqual(result["verdict"], "NO LOCAL TRANSIENT (NO INFORMATION)")
        self.assertIn("not evidence", " ".join(result["reasons"]))

    def test_a_weak_amplitude_is_below_threshold(self):
        result = verdict(sensor(), frequency_hz=47e6, radio_time_utc=self.STAMP,
                         b_amplitude_t=1e-9)
        self.assertEqual(result["verdict"], "NO LOCAL TRANSIENT (NO INFORMATION)")
        self.assertLess(result["snr"], 3.0)

    def test_a_strong_local_transient_detracts(self):
        # threshold is 1e-9 * sqrt(7.5e6) = 2.7e-6 T, so a nT is far below it
        # and 10 uT is well above
        result = verdict(sensor(), frequency_hz=47e6, radio_time_utc=self.STAMP,
                         b_amplitude_t=1e-5)
        self.assertEqual(result["verdict"], "LOCAL TRANSIENT (DETRACTS)")
        self.assertGreater(result["snr"], 3.0)
        self.assertGreater(result["e_v_m"], 0.0)


class GuardTests(unittest.TestCase):
    def test_a_local_field_record_can_never_support_anything(self):
        with self.assertRaises(NVAsEvidence):
            guard_support("20240322 47.2416 MHz")

    def test_it_is_refused_with_or_without_a_claim(self):
        with self.assertRaises(NVAsEvidence):
            guard_support("20240322 47.2416 MHz",
                          claimed="CONFIRMED TECHNOSIGNATURE CANDIDATE")


if __name__ == "__main__":
    unittest.main()
