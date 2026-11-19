"""phi^2 = phi + 1 as a recurrence: the rate, the depth, and the rejected window."""

from __future__ import annotations

import math
import unittest
from fractions import Fraction

from seti_search.phi_fixedpoint import (
    CONTRACTION,
    PHI,
    STEPS_PER_DECADE,
    agreement_depth,
    continued_fraction,
    convergents,
    derivative,
    effective_depth,
    excess,
    fibonacci,
    is_phi_shaped,
    orbit,
    step,
)

LN_PHI = math.log(PHI)


def rung_width_hz(tolerance_rungs: float, freq_hz: float) -> float:
    """Physical width spanned by a rung tolerance, used by the fixture below."""
    return (PHI ** tolerance_rungs - PHI ** (-tolerance_rungs)) * freq_hz


class RecurrenceTests(unittest.TestCase):
    def test_phi_is_the_fixed_point(self):
        self.assertAlmostEqual(step(PHI), PHI, places=15)
        self.assertAlmostEqual(PHI ** 2, PHI + 1.0, places=15)

    def test_the_derivative_at_the_fixed_point_is_the_contraction(self):
        self.assertAlmostEqual(derivative(PHI), -CONTRACTION, places=15)
        self.assertAlmostEqual(CONTRACTION, 1.0 / PHI ** 2, places=16)

    def test_every_start_contracts_at_the_same_rate(self):
        for start in (0.5, 1.0, 2.0, 6.485, 100.0):
            ob = orbit(start, turns=30)
            self.assertAlmostEqual(ob.contraction, CONTRACTION, places=4)

    def test_the_residual_shrinks_two_turns_at_a_time(self):
        # T has a negative derivative at phi, so the orbit straddles the fixed
        # point: residuals alternate but every second one is smaller
        res = orbit(0.5, turns=20).residuals
        for i in range(1, len(res) - 2):
            self.assertLess(res[i + 2], res[i])

    def test_turns_per_decade_match_the_rate(self):
        self.assertAlmostEqual(STEPS_PER_DECADE,
                               1.0 / abs(math.log10(CONTRACTION)), places=12)
        self.assertAlmostEqual(STEPS_PER_DECADE, 2.3925, places=4)

    def test_non_positive_inputs_are_refused(self):
        with self.assertRaises(ValueError):
            step(0.0)
        with self.assertRaises(ValueError):
            derivative(-1.0)
        with self.assertRaises(ValueError):
            orbit(0.0)


class LadderTests(unittest.TestCase):
    def test_phi_reads_add_one_at_every_level(self):
        part = continued_fraction(PHI, depth=12)
        self.assertEqual(part[:2], [1, 1])
        self.assertTrue(is_phi_shaped(PHI))
        self.assertTrue(is_phi_shaped(1.0 / PHI))
        self.assertTrue(is_phi_shaped(PHI ** 2))

    def test_the_channel_grid_says_add_one_exactly_once(self):
        # (k+1)/k = 1 + 1/k, so [1; k]: the +1 is carried after one level
        part = continued_fraction(Fraction(318, 317))
        self.assertEqual(part, [1, 317])
        self.assertEqual(agreement_depth(part), 0)
        self.assertEqual(excess(part), 316)
        self.assertFalse(is_phi_shaped(318 / 317))

    def test_a_convergent_matches_until_its_terminating_two(self):
        part = continued_fraction(Fraction(13, 8))
        self.assertEqual(part, [1, 1, 1, 1, 2])
        self.assertEqual(agreement_depth(part), 3)
        self.assertEqual(excess(part), 1)

    def test_depth_is_monotone_in_closeness_to_phi(self):
        depths = [agreement_depth(continued_fraction(PHI * (1 + j)))
                  for j in (1e-1, 1e-2, 1e-3, 1e-4, 1e-6)]
        self.assertEqual(depths, sorted(depths))
        self.assertGreater(depths[-1], depths[0] + 5)

    def test_effective_depth_is_continuous_and_can_be_negative(self):
        self.assertTrue(math.isinf(effective_depth(PHI)))
        self.assertAlmostEqual(effective_depth(PHI * (1 + 1e-4)), 9.57, places=1)
        self.assertAlmostEqual(effective_depth(PHI * (1 + 1e-2)), 4.78, places=1)
        self.assertLess(effective_depth(869 / 134), 0.0)

    def test_bad_depth_is_refused(self):
        with self.assertRaises(ValueError):
            continued_fraction(0.0)
        with self.assertRaises(ValueError):
            continued_fraction(2.0, depth=0)


class ConvergentDecayTests(unittest.TestCase):
    def test_the_ladder_is_the_fibonacci_ladder(self):
        conv = convergents(PHI, depth=8)
        self.assertEqual([(p, q) for p, q, _ in conv[:5]],
                         [(1, 1), (2, 1), (3, 2), (5, 3), (8, 5)])
        fibs = {fibonacci(n) for n in range(1, 12)}
        for p, q, _ in conv:
            self.assertIn(p, fibs)
            self.assertIn(q, fibs)

    def test_each_level_costs_one_over_phi_squared(self):
        conv = convergents(PHI, depth=14)
        ratios = [conv[i][2] / conv[i - 1][2] for i in range(1, len(conv))]
        for r in ratios[-6:]:
            self.assertAlmostEqual(r, CONTRACTION, places=3)


class RejectedToleranceTest(unittest.TestCase):
    """The ±0.01 rung window, kept as a fixture for why it was withdrawn."""

    def test_a_rung_window_is_not_scale_invariant(self):
        low = rung_width_hz(0.01, 13.4e6)
        high = rung_width_hz(0.01, 86.9e6)
        self.assertAlmostEqual(high / low, 86.9 / 13.4, places=9)
        self.assertGreater(high / low, 6.0)

    def test_the_tolerated_width_is_not_a_constant_of_the_instrument(self):
        widths = [rung_width_hz(0.01, f) for f in (13.4e6, 47e6, 86.9e6)]
        self.assertAlmostEqual(widths[0] / 1e6, 0.1290, places=3)
        self.assertAlmostEqual(widths[-1] / 1e6, 0.8363, places=3)
        # a fixed rung tolerance therefore tests a different physical width
        # at every frequency, which is why the threshold was withdrawn
        self.assertGreater(max(widths) / min(widths), 6.0)

    def test_the_rate_itself_is_scale_free(self):
        # the rate, unlike a window, has no unit: it is the same number at
        # every frequency because it is the derivative at the fixed point
        for f in (13.4e6, 47e6, 86.9e6):
            self.assertAlmostEqual(math.exp(math.log(PHI ** 2) * 1.0),
                                   PHI ** 2, places=12)
            self.assertAlmostEqual(derivative(PHI), -1.0 / PHI ** 2, places=16)


if __name__ == "__main__":
    unittest.main()
