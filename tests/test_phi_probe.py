"""The rung law must round-trip, and the phase test must have a null."""

from __future__ import annotations

import math
import unittest

from scripts.phi_probe import (
    PHI,
    angular_distance,
    null_phase,
    phase_deg,
    probe_expansion,
    probe_phase,
    round_trip,
    to_rung,
)


class RoundTripTests(unittest.TestCase):
    def test_the_map_inverts_to_machine_precision(self):
        for value in (1.0, 47.5358e6, 1e14, 7.3e-3):
            self.assertLess(round_trip(value), 1e-12)

    def test_phi_squared_is_two_rungs(self):
        self.assertAlmostEqual(to_rung(PHI ** 2), 2.0, places=12)

    def test_phi_itself_is_one_rung(self):
        self.assertAlmostEqual(to_rung(PHI), 1.0, places=12)


class PhaseTests(unittest.TestCase):
    def test_the_phase_is_the_fractional_rung_in_degrees(self):
        self.assertAlmostEqual(phase_deg(36.734329351), 264.359, places=2)

    def test_angular_distance_wraps(self):
        self.assertAlmostEqual(angular_distance(359.0, 1.0), 2.0)
        self.assertAlmostEqual(angular_distance(10.0, 200.0), 170.0)

    def test_the_null_depends_on_the_tolerance(self):
        self.assertLess(null_phase(0.5), null_phase(2.0))
        self.assertAlmostEqual(null_phase(1.0), 4.0 / 360.0)


class RecordTests(unittest.TestCase):
    def test_a_record_that_is_on_target_would_be_flagged(self):
        # construct the frequency whose phase is exactly 137.5 deg and check it
        rung = 36.0 + 137.5 / 360.0
        freq_mhz = (PHI ** rung) / 1e6
        result = probe_phase(freq_mhz)
        self.assertTrue(result["hit"]["0.5"])
        self.assertLess(result["distance_to_137.5"], 0.5)

    def test_an_arithmetic_comb_is_not_a_fibonacci_expansion(self):
        block = probe_expansion([45.9843, 46.9828, 47.9380, 48.9937])
        self.assertGreater(block["closest_ratio_gap"], 0.01)
        self.assertLess(block["mean_rung_step"], 0.5)


if __name__ == "__main__":
    unittest.main()
