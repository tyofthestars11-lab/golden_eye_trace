"""Tests for the power calculations: monotonicity, units and inverses."""

from __future__ import annotations

import math
import unittest

from seti_search.sensitivity import (
    background_events,
    chance_probability,
    energy_fluence_limit_erg_cm2,
    fluence_limit_cm2,
    horizon_mpc,
    power,
    required_events,
    signal_events,
    window_for_chance_probability,
)


class CountingTests(unittest.TestCase):
    def test_ninety_percent_power_needs_two_point_three_events(self):
        self.assertAlmostEqual(required_events(0.9), 2.3026, places=3)
        self.assertAlmostEqual(power(required_events(0.9)), 0.9, places=9)

    def test_power_is_monotone_and_bounded(self):
        self.assertLess(power(0.5), power(2.0))
        self.assertLess(power(0.0), 0.5)
        self.assertTrue(0.0 <= power(0.0) < 1e-12)

    def test_impossible_probabilities_are_refused(self):
        for bad in (0.0, 1.0, -0.1, 1.5):
            with self.assertRaises(ValueError):
                required_events(bad)
        with self.assertRaises(ValueError):
            chance_probability(-1.0)


class BackgroundTests(unittest.TestCase):
    def test_the_background_scales_with_the_window(self):
        one = background_events(100, 1.0, 1.0)
        two = background_events(100, 1.0, 256.0)
        self.assertAlmostEqual(two / one, 256.0, places=6)

    def test_the_background_scales_with_the_rate(self):
        self.assertAlmostEqual(background_events(200, 1.0, 100.0),
                               2 * background_events(100, 1.0, 100.0))

    def test_a_bigger_cone_is_a_bigger_background(self):
        self.assertLess(background_events(100, 1.0, 256.0),
                        background_events(100, 3.0, 256.0))

    def test_the_number_is_small_enough_to_be_believable(self):
        # 100 alerts a year, a one-degree cone, a 256 s window
        mu = background_events(100, 1.0, 256.0)
        self.assertLess(mu, 1e-5)
        self.assertGreater(mu, 0.0)

    def test_the_window_inverts(self):
        for rate in (30, 100, 1000):
            for cone in (0.5, 1.0, 4.68):
                window = window_for_chance_probability(rate, cone, 0.01)
                mu = background_events(rate, cone, window)
                self.assertAlmostEqual(chance_probability(mu), 0.01, places=9)

    def test_an_inverse_needs_a_rate(self):
        with self.assertRaises(ValueError):
            window_for_chance_probability(0, 1.0, 0.01)


class SensitivityTests(unittest.TestCase):
    def test_fluence_falls_with_area_and_duration(self):
        base = fluence_limit_cm2(0.9, 1e4, 100.0)
        self.assertAlmostEqual(fluence_limit_cm2(0.9, 2e4, 100.0), base / 2)
        self.assertAlmostEqual(fluence_limit_cm2(0.9, 1e4, 200.0), base / 2)

    def test_the_unit_conversion_is_in_the_number(self):
        # 1e-10 per cm^2 per s through 1e4 m^2 (1e8 cm^2) for 256 s is 2.56 events
        self.assertAlmostEqual(signal_events(1e-10, 1e4, 256.0), 2.56, places=12)

    def test_energy_fluence_is_the_number_fluence_times_the_pivot(self):
        number = fluence_limit_cm2(0.9, 1e4, 256.0)
        energy = energy_fluence_limit_erg_cm2(0.9, 1e4, 256.0, 1e6)
        self.assertAlmostEqual(energy, number * 1e6 * 1.602176634e-3, places=20)

    def test_the_horizon_grows_as_the_square_root_of_energy(self):
        limit = energy_fluence_limit_erg_cm2(0.9, 1e4, 256.0, 1e6)
        near = horizon_mpc(1e47, limit)
        far = horizon_mpc(4e47, limit)
        self.assertAlmostEqual(far / near, 2.0, places=9)

    def test_a_bigger_detector_reaches_further(self):
        small = energy_fluence_limit_erg_cm2(0.9, 1e3, 256.0, 1e6)
        large = energy_fluence_limit_erg_cm2(0.9, 1e4, 256.0, 1e6)
        self.assertGreater(horizon_mpc(1e47, large), horizon_mpc(1e47, small))

    def test_negatives_are_refused(self):
        with self.assertRaises(ValueError):
            signal_events(-1.0, 1e4, 10.0)
        with self.assertRaises(ValueError):
            fluence_limit_cm2(0.9, 0.0, 10.0)
        with self.assertRaises(ValueError):
            horizon_mpc(0.0, 1.0)

    def test_the_flux_limit_is_physically_sensible(self):
        # a PeV-scale detector with 1e4 m^2 over a day should reach a limit of
        # order 1e-11 erg/cm^2, not 1e-3 or 1e-17
        limit = energy_fluence_limit_erg_cm2(0.9, 1e4, 86400.0, 1e6)
        self.assertLess(limit, 1e-9)
        self.assertGreater(limit, 1e-13)


if __name__ == "__main__":
    unittest.main()
