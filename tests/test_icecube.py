"""Tests for the IceCube ingest: units, refusals and the containment convention."""

from __future__ import annotations

import json
import math
import os
import unittest

from seti_search.icecube import (
    IceCubeSchemaError,
    from_alert,
    sigma_from_containment_deg,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAMPLE = os.path.join(ROOT, "data", "icecube_alerts_sample.json")

URL = "https://example.invalid/alert"
RELEASE = "test-release"


def gcn_row(**overrides) -> dict:
    row = {
        "event_name": ["IceCube-230416A"],
        "id": ["137840_57034692_0"],
        "alert_type": "initial",
        "alert_tense": "current",
        "analysis_pipeline": "IceCube Bronze Track alert",
        "alert_topology": "Track",
        "ra": 345.82,
        "dec": 9.01,
        "ra_dec_error": 0.5,
        "containment_probability": 0.9,
        "trigger_time": "2023-04-16T05:22:26.150574Z",
        "nu_energy": 127.29,
        "p_astro": 0.34064,
        "far": 8.029e-8,
        "far_unit": "hz",
    }
    row.update(overrides)
    return row


class ContainmentTests(unittest.TestCase):
    def test_ninety_percent_radius_is_not_a_sigma(self):
        # 0.5 deg at 90% is 0.233 sigma; using 0.5 directly would widen the cone
        # by 2.1x and every chance probability with it
        self.assertAlmostEqual(sigma_from_containment_deg(0.5, 0.9), 0.5 / 2.146,
                               places=3)
        self.assertAlmostEqual(sigma_from_containment_deg(1.0, 0.9), 0.466, places=3)

    def test_fifty_percent_radius(self):
        self.assertAlmostEqual(sigma_from_containment_deg(1.0, 0.5), 0.849, places=3)

    def test_the_conversion_is_the_gaussian_one(self):
        # solve the Rayleigh CDF back: 1 - exp(-R^2 / 2 sigma^2)
        sigma = sigma_from_containment_deg(2.0, 0.68)
        self.assertAlmostEqual(1 - math.exp(-(2.0 ** 2) / (2 * sigma ** 2)), 0.68,
                               places=9)

    def test_impossible_containment_is_refused(self):
        for bad in (0.0, 1.0, 1.5):
            with self.assertRaises(IceCubeSchemaError):
                sigma_from_containment_deg(1.0, bad)
        with self.assertRaises(IceCubeSchemaError):
            sigma_from_containment_deg(0.0, 0.9)


class ConversionTests(unittest.TestCase):
    def test_tev_becomes_gev(self):
        alert = from_alert(gcn_row(), public_source_url=URL, data_release=RELEASE)
        self.assertAlmostEqual(alert.event.energy_estimate_gev, 127290.0)

    def test_hz_becomes_per_year(self):
        alert = from_alert(gcn_row(), public_source_url=URL, data_release=RELEASE)
        self.assertAlmostEqual(alert.far_per_year, 8.029e-8 * 3.1557e7)
        self.assertGreater(alert.far_per_year, 2.0)
        self.assertLess(alert.far_per_year, 3.0)

    def test_a_per_year_rate_is_read_as_written(self):
        row = gcn_row()
        row.pop("far_unit")
        row["far"] = 2.1
        alert = from_alert(row, public_source_url=URL, data_release=RELEASE)
        self.assertAlmostEqual(alert.far_per_year, 2.1)

    def test_the_grade_comes_from_the_pipeline(self):
        bronze = from_alert(gcn_row(), public_source_url=URL, data_release=RELEASE)
        gold = from_alert(gcn_row(analysis_pipeline="IceCube Gold Track alert"),
                          public_source_url=URL, data_release=RELEASE)
        self.assertEqual(bronze.event.alert_status, "BRONZE")
        self.assertEqual(gold.event.alert_status, "GOLD")

    def test_an_ungraded_row_is_refused_rather_than_guessed(self):
        row = gcn_row(analysis_pipeline="IceCube Track alert")
        with self.assertRaises(IceCubeSchemaError):
            from_alert(row, public_source_url=URL, data_release=RELEASE)

    def test_the_conversions_are_recorded(self):
        alert = from_alert(gcn_row(), public_source_url=URL, data_release=RELEASE)
        joined = " | ".join(alert.conversions)
        self.assertIn("1 sigma", joined)
        self.assertIn("GeV", joined)
        self.assertIn("per year", joined)


class RefusalTests(unittest.TestCase):
    def test_a_retraction_is_not_an_event(self):
        with self.assertRaises(IceCubeSchemaError) as caught:
            from_alert(gcn_row(alert_type="retraction"),
                       public_source_url=URL, data_release=RELEASE)
        self.assertIn("retracted", str(caught.exception))

    def test_the_notice_time_is_not_the_event_time(self):
        row = gcn_row()
        row.pop("trigger_time")
        with self.assertRaises(IceCubeSchemaError) as caught:
            from_alert(row, public_source_url=URL, data_release=RELEASE)
        self.assertIn("alert_datetime", str(caught.exception))

    def test_a_naive_trigger_time_is_refused(self):
        with self.assertRaises(ValueError):
            from_alert(gcn_row(trigger_time="2023-04-16T05:22:26"),
                       public_source_url=URL, data_release=RELEASE)

    def test_a_hardware_injection_needs_permission(self):
        with self.assertRaises(IceCubeSchemaError):
            from_alert(gcn_row(alert_tense="injections"),
                       public_source_url=URL, data_release=RELEASE)
        alert = from_alert(gcn_row(alert_tense="injections"), allow_test=True,
                           public_source_url=URL, data_release=RELEASE)
        self.assertEqual(alert.event.alert_status, "test")
        self.assertTrue(alert.event.synthetic)

    def test_an_error_without_a_containment_is_refused(self):
        row = gcn_row()
        row.pop("ra_dec_error")
        with self.assertRaises(IceCubeSchemaError):
            from_alert(row, public_source_url=URL, data_release=RELEASE)


class LegacyLayoutTests(unittest.TestCase):
    def test_the_amon_columns_are_accepted(self):
        row = {
            "id": "143317_9424183",
            "record_number": 1,
            "notice_type": "BRONZE",
            "trigger_time": "2026-10-20T03:58:53.39Z",
            "ra": 224.12,
            "dec": 20.62,
            "error90_arcmin": 33.92,
            "error50_arcmin": 18.41,
            "energy_tev": 115.56,
            "p_astro": 0.37984,
            "far_per_year": 2.1081,
        }
        alert = from_alert(row, public_source_url=URL, data_release=RELEASE)
        self.assertEqual(alert.kind, "amon-table")
        self.assertEqual(alert.event.alert_status, "BRONZE")
        self.assertAlmostEqual(alert.event.angular_error_deg, 33.92 / 60 / 2.146,
                               places=3)
        self.assertAlmostEqual(alert.event.energy_estimate_gev, 115560.0)
        self.assertEqual(alert.event.ra_deg, 224.12)
        self.assertTrue(any("record_number 1" in c for c in alert.conversions))

    def test_a_row_without_a_name_uses_its_id(self):
        row = {
            "id": "143317_9424183",
            "notice_type": "BRONZE",
            "trigger_time": "2026-10-20T03:58:53.39Z",
            "ra": 1.0, "dec": 2.0, "error90_arcmin": 10.0,
            "energy_tev": 1.0,
        }
        alert = from_alert(row, public_source_url=URL, data_release=RELEASE)
        self.assertEqual(alert.event.event_id, "IC-143317_9424183")


class SampleArchiveTests(unittest.TestCase):
    def test_every_bundled_row_ingests(self):
        with open(SAMPLE, encoding="utf-8") as fh:
            rows = json.load(fh)
        self.assertGreaterEqual(len(rows), 2)
        for row in rows:
            alert = from_alert(row, public_source_url=str(row.get("_source")),
                               data_release=RELEASE)
            self.assertGreater(alert.event.angular_error_deg, 0.0)
            self.assertTrue(alert.event.public_source_url)
            # a real alert must never be labelled as a test
            self.assertNotEqual(alert.event.alert_status, "test")


if __name__ == "__main__":
    unittest.main()
