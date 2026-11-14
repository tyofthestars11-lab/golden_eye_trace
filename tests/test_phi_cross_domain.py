"""The preregistered test must fire on a planted law, or a null is meaningless."""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from phi_cross_domain_test import (  # noqa: E402
    PHI,
    PRIMARY_BASE,
    PRIMARY_TARGET,
    bh_adjust,
    fisher_combine,
    load_dataset,
    observed_hits,
    rung_delta,
    InputError,
)

STEP = PHI ** PRIMARY_TARGET  # 1.201783: one rung interval per step


class PowerTests(unittest.TestCase):
    def test_a_planted_law_is_detected_every_time(self):
        planted = [10.0 * STEP ** i for i in range(6)]
        hits, deltas = observed_hits(planted, PRIMARY_BASE, 0.01)
        self.assertEqual(hits, len(deltas))
        for delta in deltas:
            self.assertAlmostEqual(delta, PRIMARY_TARGET, places=9)

    def test_an_arithmetic_comb_is_not_detected(self):
        comb = [45.984284 + i for i in range(5)]
        hits, deltas = observed_hits(comb, PRIMARY_BASE, 0.01)
        self.assertEqual(hits, 0)
        for delta in deltas:
            self.assertLess(abs(delta - PRIMARY_TARGET), 10.0)
            self.assertLess(abs(delta), 0.1)

    def test_the_tolerance_is_two_sided_and_tight(self):
        # a ratio one rung interval away, plus 0.005, is still inside; plus
        # 0.05 is not
        near = [10.0, 10.0 * PHI ** (PRIMARY_TARGET + 0.005)]
        far = [10.0, 10.0 * PHI ** (PRIMARY_TARGET + 0.05)]
        self.assertEqual(observed_hits(near, PRIMARY_BASE, 0.01)[0], 1)
        self.assertEqual(observed_hits(far, PRIMARY_BASE, 0.01)[0], 0)

    def test_one_channel_of_a_megahertz_band_is_wider_than_the_tolerance(self):
        # d(delta_r)/df = 1/(f ln phi), so a 1 MHz step at 47 MHz moves the
        # rung interval by 0.044 - four times the tolerance. A survey at
        # channel resolution cannot resolve the prediction at all.
        per_mhz = 1.0 / (47.0e6 * math.log(PHI))
        self.assertGreater(per_mhz * 1e6, 0.01 * 4)


class StatisticsTests(unittest.TestCase):
    def test_rung_delta_is_zero_for_equal_values(self):
        self.assertAlmostEqual(rung_delta(1.0, PRIMARY_BASE), 0.0)

    def test_fisher_combines_to_one_for_null_p_values(self):
        _, p = fisher_combine([1.0, 1.0, 1.0])
        self.assertAlmostEqual(p, 1.0)

    def test_fisher_is_small_for_small_p_values(self):
        _, p = fisher_combine([0.001, 0.001, 0.001])
        self.assertLess(p, 1e-6)

    def test_bh_is_monotone_and_bounded(self):
        adjusted = bh_adjust({"a": 0.01, "b": 0.04, "c": 0.03})
        self.assertLessEqual(adjusted["a"], adjusted["b"])
        self.assertLessEqual(adjusted["a"], adjusted["c"])
        self.assertTrue(all(0.0 <= v <= 1.0 for v in adjusted.values()))


class LoaderTests(unittest.TestCase):
    def test_two_rows_are_refused(self):
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as fh:
            fh.write("order,value\n1,1.0\n2,2.0\n")
            path = Path(fh.name)
        with self.assertRaises(InputError):
            load_dataset(path)

    def test_a_non_positive_value_is_refused(self):
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as fh:
            fh.write("order,value\n1,1.0\n2,0.0\n3,3.0\n")
            path = Path(fh.name)
        with self.assertRaises(InputError):
            load_dataset(path)


if __name__ == "__main__":
    unittest.main()
