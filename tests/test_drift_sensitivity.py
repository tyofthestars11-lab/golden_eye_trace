"""Completeness falls with drift, and the guard's prediction is conservative.

The fixture is the measured artifact rather than a re-run: the sweep takes
ninety seconds and is not something to repeat in a unit test. What is asserted
here is the shape of the result - flat while the window holds the whole trail,
falling after, and never falling below the guard's prediction - so a change to
the dedoppler's shift convention or to the stabilizer's law shows up as a
failure here rather than as a quietly different number in a report.
"""

import json
import math
import os
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARTIFACT = os.path.join(ROOT, "prereg", "drift_sensitivity.json")


def load():
    with open(ARTIFACT, encoding="utf-8") as fh:
        return json.load(fh)


class DriftSensitivityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = load()
        cls.rows = {row["rate_hz_s"]: row for row in cls.data["rates"]}

    def test_the_ladder_reaches_the_fold(self):
        self.assertAlmostEqual(max(self.rows) / 1e3, 190.0, delta=1.0)
        self.assertAlmostEqual(self.data["nyquist_hz_s"] / 1e3, 196.4, delta=0.5)

    def test_the_window_limit_is_the_stabilizers_number(self):
        self.assertAlmostEqual(self.data["window_limit_hz_s"] / 1e3, 58.9, delta=0.3)

    def test_a_stationary_lane_is_its_own_reference(self):
        for bracket in self.rows[0.0]["brackets"]:
            self.assertAlmostEqual(bracket["ratio_observed"], 1.0, places=9)

    def test_the_ratio_is_flat_while_the_whole_trail_fits(self):
        # the median across brackets, because one cell is a median of five
        # injections and carries about +/-0.15 of its own scatter
        for rate, row in self.rows.items():
            if row["fraction_held"] >= 1.0 - 1e-9 and rate <= 20000.0:
                ratios = [b["ratio_observed"] for b in row["brackets"]]
                self.assertGreater(sorted(ratios)[len(ratios) // 2], 0.75,
                                   msg=f"{rate} kHz/s kept too little")

    def test_the_ratio_falls_once_the_trail_outgrows_the_window(self):
        fold = self.rows[190000.0]
        self.assertLess(fold["fraction_held"], 0.25)
        for bracket in fold["brackets"]:
            self.assertLess(bracket["ratio_observed"], 0.75)
            self.assertGreater(bracket["ratio_observed"], 0.4)

    def test_the_guard_overstates_the_loss_rather_than_understating_it(self):
        # Where the truncation is deep enough to stand out from the five-injection
        # scatter, the measured significance must not fall below sqrt(held): a
        # window the guard passes must not be hiding a shortfall.
        for rate, row in self.rows.items():
            if row["fraction_held"] > 0.65:
                continue
            for bracket in row["brackets"]:
                self.assertGreaterEqual(
                    bracket["ratio_observed"] + 0.02, bracket["predicted_ratio"],
                    msg=f"{rate} kHz/s fell below the guard's prediction")

    def test_completeness_falls_with_drift_at_the_threshold_strength(self):
        dim = self.rows[0.0]["brackets"][0]
        fold = self.rows[190000.0]["brackets"][0]
        self.assertEqual(dim["same-entry"]["fraction"], 1.0)
        self.assertLess(fold["same-entry"]["fraction"], dim["same-entry"]["fraction"])
        self.assertLessEqual(fold["same-entry"]["fraction"], 0.6)

    def test_a_bright_line_survives_what_a_dim_one_does_not(self):
        fold = self.rows[190000.0]["brackets"]
        self.assertGreaterEqual(fold[-1]["same-entry"]["fraction"],
                                fold[0]["same-entry"]["fraction"])
        self.assertLessEqual(fold[0]["uniform"]["fraction"], 0.2)

    def test_the_rows_agree_with_their_own_geometry(self):
        for row in self.data["rates"]:
            if row["trail_channels"] > 0:
                self.assertLessEqual(row["fraction_held"],
                                     min(1.0, 96.0 / row["trail_channels"]) + 1e-9)
            for bracket in row["brackets"]:
                self.assertAlmostEqual(
                    bracket["predicted_ratio"],
                    math.sqrt(bracket["same-entry"]["mean_held"]), places=9)


if __name__ == "__main__":
    unittest.main()
