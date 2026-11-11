"""Golden-angle rung tests.

The arithmetic here is exact and checkable: rung is a logarithm to base phi, so
rung(phi**k) = k with no remainder. What the tests also pin down is the limit:
fractional rungs of a smooth population are uniform, so a single small drift is
not a finding, and only repetition or non-uniformity says anything - and what
they say is that the *pipeline* quantized something.
"""

import math
import unittest

import numpy as np

from seti_search.rungs import (
    GOLDEN_ANGLE_DEG,
    PHI,
    address,
    delta_rung,
    ladder_step,
    repeated_scales,
    rung,
    rung_drift,
    rung_uniformity,
)


class RungArithmeticTests(unittest.TestCase):
    def test_exact_ladder_points(self):
        self.assertAlmostEqual(rung(1.0), 0.0, places=12)
        self.assertAlmostEqual(rung(PHI), 1.0, places=12)
        self.assertAlmostEqual(rung(PHI**2), 2.0, places=12)
        self.assertAlmostEqual(rung(1.0 / PHI), -1.0, places=12)
        for k in range(-8, 9):
            self.assertAlmostEqual(rung(PHI**k), float(k), places=10)

    def test_golden_angle_is_one_over_phi_squared_turn(self):
        # a full turn is exactly phi^2 golden angles, because phi^2 = phi + 1
        self.assertAlmostEqual(360.0 / GOLDEN_ANGLE_DEG, PHI**2, places=10)
        self.assertAlmostEqual(GOLDEN_ANGLE_DEG, 137.5077640, places=5)
        self.assertAlmostEqual(PHI**2, PHI + 1.0, places=12)

    def test_delta_rung_is_a_log_ratio(self):
        a, b = 1420.4e6, 1421.0e6
        self.assertAlmostEqual(delta_rung(a, b), math.log(a / b) / math.log(PHI), places=10)

    def test_drift_is_bounded_and_not_special(self):
        rng = np.random.default_rng(7)
        values = rng.lognormal(mean=10.0, sigma=3.0, size=2000)
        drifts = np.array([rung_drift(v) for v in values])
        self.assertLessEqual(drifts.max(), 0.5)
        # roughly uniform on [0, 0.5]: about half below 0.25, about a tenth
        # below 0.05. A framework that treats small drift as membership has to
        # live with this.
        self.assertGreater((drifts < 0.25).mean(), 0.40)
        self.assertLess((drifts < 0.25).mean(), 0.60)
        self.assertGreater((drifts < 0.05).mean(), 0.05)
        self.assertLess((drifts < 0.05).mean(), 0.16)

    def test_address_reports_nearest_rung_and_drift(self):
        spot = address(PHI**5 * 1.01, "Hz")
        self.assertEqual(spot.nearest, 5)
        self.assertAlmostEqual(spot.drift, abs(spot.rung - 5.0), places=12)
        self.assertIn("nearest 5", str(spot))

    def test_zero_and_negative_are_rejected(self):
        for bad in (0.0, -1.0):
            with self.assertRaises(ValueError):
                rung(bad)


class LadderTests(unittest.TestCase):
    def test_integer_steps_are_recognised(self):
        step = ladder_step(PHI**3 * 7.0, 7.0)
        self.assertTrue(step.integer_step)
        self.assertEqual(step.nearest_integer, 3)
        self.assertAlmostEqual(step.factor, PHI**3, places=10)

    def test_a_plain_scale_factor_is_not_a_ladder_step(self):
        step = ladder_step(3000.0, 2000.0)  # 3/2, not a power of phi
        self.assertFalse(step.integer_step)
        self.assertAlmostEqual(step.factor, 1.5, places=12)
        # the residual is exact arithmetic, so it is round-off, not a near miss
        self.assertAlmostEqual(step.residual, 0.0, places=10)
        self.assertIn("not a ladder step", str(step))


class ArtifactTests(unittest.TestCase):
    def test_a_repeated_scale_factor_is_caught(self):
        # one object measured at two scales, seven times over: ratios all 1.5
        values = [100.0 * 1.5**k for k in range(8)]
        found = repeated_scales(values, min_pairs=3)
        self.assertTrue(found)
        # disjoint pairs are required, so a short chain cannot reach the count
        self.assertTrue(
            any(abs(r.factor - 1.5) < 1e-3 and r.pairs >= 3 for r in found),
            msg=[str(r) for r in found],
        )
        # the factor is reported as a magnitude, never below 1
        self.assertTrue(all(r.factor >= 1.0 for r in found))

    def test_unrelated_values_show_no_repetition(self):
        rng = np.random.default_rng(11)
        values = list(rng.uniform(1e3, 1e5, size=40))
        self.assertEqual(repeated_scales(values, min_pairs=3), [])

    def test_quantized_population_fails_the_uniformity_test(self):
        values = [PHI**k * 3.0 for k in range(4, 20)]
        result = rung_uniformity(values)
        self.assertLess(result["p_value"], 1e-3)

    def test_smooth_population_passes(self):
        rng = np.random.default_rng(3)
        values = rng.lognormal(mean=14.0, sigma=2.0, size=200)
        result = rung_uniformity(values)
        self.assertGreater(result["p_value"], 0.01, msg=result)
        self.assertEqual(result["n"], 200)

    def test_too_few_values_is_not_a_test(self):
        result = rung_uniformity([1.0, 2.0, 3.0])
        self.assertNotEqual(result["p_value"], result["p_value"])  # NaN


if __name__ == "__main__":
    unittest.main()
