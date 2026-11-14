"""Tests for the neutrino lane: the schema, the null, and the refusals."""

from __future__ import annotations

import datetime as dt
import math
import unittest
from types import SimpleNamespace

from seti_search.neutrino import (
    Coincidence,
    DetectorExposure,
    NeutrinoAsEvidence,
    NeutrinoEvent,
    angular_separation_deg,
    chance_probability,
    guard_support,
    match_catalog,
    solid_angle_sr,
)

UTC = dt.timezone.utc


def event(**overrides) -> NeutrinoEvent:
    base = dict(
        event_id="EV-1",
        detector="TEST",
        event_time_utc=dt.datetime(2024, 3, 22, 8, 15, 0, tzinfo=UTC),
        time_uncertainty_s=0.5,
        ra_deg=100.0,
        dec_deg=-20.0,
        angular_error_deg=1.0,
        energy_estimate_gev=1e14,
        energy_class="track",
        event_type="through-going",
        alert_status="archival",
        public_source_url="https://example.invalid/EV-1",
        data_release="test-release",
    )
    base.update(overrides)
    return NeutrinoEvent(**base)


def exposure(n: int = 100) -> DetectorExposure:
    return DetectorExposure(livetime_s=3.156e7, n_events=n)


class SchemaTests(unittest.TestCase):
    def test_a_valid_event_builds(self):
        self.assertEqual(event().event_id, "EV-1")

    def test_a_naive_timestamp_is_refused(self):
        # the alternative is guessing the zone, which is how a 4-hour error gets
        # into a coincidence window
        with self.assertRaises(ValueError):
            event(event_time_utc=dt.datetime(2024, 3, 22, 8, 15))

    def test_a_zero_angular_error_is_refused(self):
        with self.assertRaises(ValueError):
            event(angular_error_deg=0.0)

    def test_ra_must_be_inside_the_circle(self):
        with self.assertRaises(ValueError):
            event(ra_deg=360.0)

    def test_vocabularies_are_closed(self):
        with self.assertRaises(ValueError):
            event(energy_class="muon")
        with self.assertRaises(ValueError):
            event(event_type="sideways")
        with self.assertRaises(ValueError):
            event(alert_status="PLATINUM")

    def test_an_event_without_a_release_is_not_evidence(self):
        with self.assertRaises(ValueError):
            event(data_release="")
        with self.assertRaises(ValueError):
            event(public_source_url="")

    def test_round_trip_keeps_the_timestamp_aware(self):
        back = NeutrinoEvent.from_dict(event().to_dict())
        self.assertEqual(back.event_time_utc, event().event_time_utc)
        self.assertEqual(back, event())


class GeometryTests(unittest.TestCase):
    def test_identical_positions_are_zero_apart(self):
        self.assertAlmostEqual(angular_separation_deg(10, 20, 10, 20), 0.0)

    def test_ra_wrap_is_one_degree_not_359(self):
        sep = angular_separation_deg(359.5, 0.0, 0.5, 0.0)
        self.assertAlmostEqual(sep, 1.0, places=6)

    def test_the_pole_separates_by_colatitude(self):
        sep = angular_separation_deg(0.0, 90.0, 180.0, 80.0)
        self.assertAlmostEqual(sep, 10.0, places=6)

    def test_a_declination_outside_the_sphere_is_refused(self):
        with self.assertRaises(ValueError):
            angular_separation_deg(0, 0, 0, 91)

    def test_solid_angle_grows_with_the_cone(self):
        self.assertLess(solid_angle_sr(0.5), solid_angle_sr(2.0))
        self.assertAlmostEqual(solid_angle_sr(180.0), 4 * math.pi, places=6)


class ChanceTests(unittest.TestCase):
    def test_no_exposure_means_no_null(self):
        # None, never zero: the absence of a null is not a small null
        self.assertIsNone(chance_probability(1.0, 86400.0, None))

    def test_the_density_is_per_steradian_per_second(self):
        self.assertAlmostEqual(exposure(100).density_per_sr_s,
                               100 / (3.156e7 * 4 * math.pi))

    def test_the_box_grows_with_angle_and_time(self):
        exp = exposure()
        self.assertLess(chance_probability(0.5, 3600, exp),
                        chance_probability(2.0, 3600, exp))
        self.assertLess(chance_probability(1.0, 3600, exp),
                        chance_probability(1.0, 86400, exp))

    def test_an_alert_grade_rate_gives_a_small_probability(self):
        # 100 events a year over the sky, a 3 sigma 1 deg cone, one day
        self.assertLess(chance_probability(3.0, 86400.0, exposure(100)), 1e-3)


class CoincidenceTests(unittest.TestCase):
    def test_absence_can_never_be_informative(self):
        with self.assertRaises(ValueError):
            Coincidence(status="NO EVENT IN WINDOW", absence_is_informative=True,
                        limits=("x",))

    def test_an_empty_query_must_name_its_limits(self):
        with self.assertRaises(ValueError):
            Coincidence(status="NO EVENT IN WINDOW")

    def test_a_coincidence_needs_its_event_and_null(self):
        with self.assertRaises(ValueError):
            Coincidence(status="COINCIDENT")

    def test_no_pointing_is_uninformative_not_absent(self):
        result = match_catalog([event()], ra_deg=None, dec_deg=None,
                               time_utc=dt.datetime(2024, 3, 22, tzinfo=UTC))
        self.assertEqual(result.status, "UNINFORMATIVE")
        self.assertIn("no pointing", result.limits[0])

    def test_no_timestamp_is_uninformative(self):
        result = match_catalog([event()], ra_deg=100.0, dec_deg=-20.0,
                               time_utc=None)
        self.assertEqual(result.status, "UNINFORMATIVE")
        self.assertIn("UTC", result.limits[0])

    def test_coincident_without_exposure_is_uninformative(self):
        # inside the box, but with no exposure there is no chance probability to
        # quote, so it is not reported as a coincidence
        result = match_catalog([event()], ra_deg=100.0, dec_deg=-20.0,
                               time_utc=dt.datetime(2024, 3, 22, 8, 15, 0,
                                                    tzinfo=UTC),
                               radio_time_uncertainty_s=10.0,
                               radio_localisation_deg=0.5)
        self.assertEqual(result.status, "UNINFORMATIVE")
        self.assertTrue(any("no detector exposure" in lim for lim in result.limits))

    def test_a_far_event_is_outside_the_box(self):
        far = event(ra_deg=200.0, dec_deg=40.0)
        result = match_catalog([far], ra_deg=100.0, dec_deg=-20.0,
                               time_utc=dt.datetime(2024, 3, 22, 8, 15, 0,
                                                    tzinfo=UTC),
                               radio_time_uncertainty_s=10.0,
                               exposure=exposure())
        self.assertEqual(result.status, "NO EVENT IN WINDOW")
        self.assertTrue(result.limits)

    def test_a_near_event_gives_a_coincidence_with_a_null(self):
        result = match_catalog([event()], ra_deg=100.2, dec_deg=-20.1,
                               time_utc=dt.datetime(2024, 3, 22, 8, 15, 0,
                                                    tzinfo=UTC),
                               radio_time_uncertainty_s=10.0,
                               radio_localisation_deg=0.3,
                               exposure=exposure())
        self.assertEqual(result.status, "COINCIDENT")
        self.assertIsNotNone(result.chance_probability)
        self.assertLess(result.separation_deg, 1.0)

    def test_the_best_event_is_the_least_chance_likely_not_the_closest(self):
        near = event(event_id="NEAR", angular_error_deg=8.0)
        better = event(event_id="BETTER", angular_error_deg=0.4)
        result = match_catalog([near, better], ra_deg=100.0, dec_deg=-20.0,
                               time_utc=dt.datetime(2024, 3, 22, 8, 15, 0,
                                                    tzinfo=UTC),
                               radio_time_uncertainty_s=10.0,
                               exposure=exposure())
        self.assertEqual(result.status, "COINCIDENT")
        self.assertEqual(result.event.event_id, "BETTER")


class GuardTests(unittest.TestCase):
    def record(self, disposition: str):
        return SimpleNamespace(disposition=disposition,
                               observation_date="20240325",
                               frequency_hz=45.9843e6)

    def coincidence(self) -> Coincidence:
        return match_catalog([event()], ra_deg=100.0, dec_deg=-20.0,
                             time_utc=dt.datetime(2024, 3, 22, 8, 15, 0,
                                                  tzinfo=UTC),
                             radio_time_uncertainty_s=10.0,
                             exposure=exposure())

    def test_an_invalid_candidate_cannot_be_restored(self):
        with self.assertRaises(NeutrinoAsEvidence) as caught:
            guard_support(self.record("MASK-MANUFACTURED SIGNAL"),
                          self.coincidence())
        self.assertIn("cannot restore", str(caught.exception))

    def test_an_invalid_candidate_is_refused_even_when_asking_for_more(self):
        with self.assertRaises(NeutrinoAsEvidence):
            guard_support(self.record("MASK-MANUFACTURED SIGNAL"),
                          self.coincidence(),
                          claimed="CONFIRMED TECHNOSIGNATURE CANDIDATE")

    def test_a_promotion_is_refused(self):
        with self.assertRaises(NeutrinoAsEvidence):
            guard_support(self.record("UNRESOLVED"), self.coincidence(),
                          claimed="CONFIRMED TECHNOSIGNATURE CANDIDATE")

    def test_a_coincidence_attaches_to_a_valid_record(self):
        attached = guard_support(self.record("UNRESOLVED"), self.coincidence())
        self.assertEqual(attached["neutrino_event_id"], "EV-1")
        self.assertIsNotNone(attached["neutrino_chance_probability"])

    def test_an_absence_cannot_be_attached_at_all(self):
        absent = match_catalog([event(ra_deg=200.0, dec_deg=40.0)],
                               ra_deg=100.0, dec_deg=-20.0,
                               time_utc=dt.datetime(2024, 3, 22, 8, 15, 0,
                                                    tzinfo=UTC),
                               radio_time_uncertainty_s=10.0,
                               exposure=exposure())
        with self.assertRaises(NeutrinoAsEvidence):
            guard_support(self.record("UNRESOLVED"), absent)


if __name__ == "__main__":
    unittest.main()
