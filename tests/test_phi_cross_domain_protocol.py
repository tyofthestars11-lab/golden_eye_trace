"""The protocol must detect a planted interval, and must not invent one."""

from __future__ import annotations

import math
import unittest

from seti_search.phi_cross_domain import (
    BASES,
    PHI,
    TARGET,
    Series,
    bh,
    deltas,
    fit_ratio,
    hit_count,
    holdout_prediction,
    injection_power,
    permutation_p,
    rung_error,
)

LAW = PHI ** TARGET  # 1.2017833...


class TargetTests(unittest.TestCase):
    def test_the_target_is_one_over_phi_squared(self):
        self.assertAlmostEqual(TARGET, 1.0 / PHI ** 2, places=15)
        self.assertAlmostEqual(TARGET, 0.38196601125010515, places=15)

    def test_the_implied_ratio_is_1_2018(self):
        self.assertAlmostEqual(LAW, 1.2017833, places=6)

    def test_a_planted_series_has_deltas_on_target(self):
        planted = [1.0 * LAW ** i for i in range(5)]
        for delta in deltas(planted, PHI):
            self.assertAlmostEqual(delta, TARGET, places=12)


class NullTests(unittest.TestCase):
    def test_a_clean_series_gives_no_hits_and_p_one(self):
        comb = [45.984284 + i for i in range(8)]
        result = permutation_p(comb, PHI, 0.01, 200, seed=1)
        self.assertEqual(result["observed"], 0)
        self.assertAlmostEqual(result["p_value"], 1.0)

    def test_a_planted_series_is_detected(self):
        planted = [10.0 * LAW ** i for i in range(8)]
        result = permutation_p(planted, PHI, 0.01, 200, seed=1)
        self.assertEqual(result["observed"], 7)
        self.assertLess(result["p_value"], 0.01)

    def test_the_p_value_never_reaches_zero(self):
        planted = [10.0 * LAW ** i for i in range(6)]
        result = permutation_p(planted, PHI, 0.01, 100, seed=1)
        self.assertGreaterEqual(result["p_value"], 1.0 / 101.0)

    def test_a_single_permutation_count_is_refused(self):
        with self.assertRaises(ValueError):
            permutation_p([1.0, 2.0, 3.0], PHI, 0.01, 50)


class HoldoutTests(unittest.TestCase):
    def test_a_planted_law_predicts_perfectly(self):
        planted = [5.0 * LAW ** i for i in range(8)]
        result = holdout_prediction(planted, PHI, 0.5)
        self.assertAlmostEqual(result["hit_rate"], 1.0)
        self.assertLess(result["median_rung_error"], 1e-9)

    def test_a_comb_is_predicted_badly_by_phi(self):
        comb = [45.984284 + i for i in range(8)]
        result = holdout_prediction(comb, PHI, 0.5)
        self.assertGreater(result["median_rung_error"], 0.2)
        self.assertEqual(result["hit_rate"], 0.0)

    def test_the_fitted_constant_beats_every_fixed_base_on_a_comb(self):
        comb = [45.984284, 46.982763, 47.937962, 48.993707]
        fitted = fit_ratio(comb)
        test = fitted["test_ratios"]
        fitted_error = sorted(
            rung_error(fitted["train_ratio"], r, PHI) for r in test
        )[len(test) // 2]
        for base in BASES.values():
            fixed = holdout_prediction(comb, base, 0.5)
            self.assertLess(fitted_error, fixed["median_rung_error"])


class PowerTests(unittest.TestCase):
    def test_the_false_positive_arm_is_clean(self):
        comb = [45.984284 + i for i in range(8)]
        rate = injection_power(comb, PHI, 0.0, trials=20, permutations=200,
                               seed=3)
        self.assertEqual(rate, 0.0)

    def test_a_blanket_law_is_detected_on_a_long_series(self):
        comb = [10.0 + i for i in range(10)]
        rate = injection_power(comb, PHI, 1.0, trials=20, permutations=200,
                               seed=3)
        self.assertGreaterEqual(rate, 0.9)


class SeriesTests(unittest.TestCase):
    def test_a_unit_must_be_declared(self):
        with self.assertRaises(ValueError):
            Series("s", "radio", "  ", (1.0, 2.0, 3.0), "p")

    def test_values_must_be_positive_and_finite(self):
        with self.assertRaises(ValueError):
            Series("s", "radio", "MHz", (1.0, 0.0, 3.0), "p")
        with self.assertRaises(ValueError):
            Series("s", "radio", "MHz", (1.0, math.inf, 3.0), "p")

    def test_two_values_are_not_a_series(self):
        with self.assertRaises(ValueError):
            Series("s", "radio", "MHz", (1.0, 2.0), "p")


class MultiplicityTests(unittest.TestCase):
    def test_bh_is_monotone_and_bounded(self):
        adjusted = bh({"a": 0.001, "b": 0.02, "c": 0.5})
        self.assertLessEqual(adjusted["a"], adjusted["b"])
        self.assertLessEqual(adjusted["b"], adjusted["c"])
        self.assertTrue(all(0 <= v <= 1 for v in adjusted.values()))


if __name__ == "__main__":
    unittest.main()
