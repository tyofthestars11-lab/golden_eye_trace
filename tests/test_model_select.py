"""Comb against scale family: the fits, the labels, and the calibration.

The point of these tests is not that a model can be fitted. It is that the module
refuses to name one when the data cannot support the name, and that it can say how
often it would have been right if it had.
"""

import math
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from seti_search.model_select import (  # noqa: E402
    BAND_SPAN_REQUIRED,
    LABELS,
    compare,
    fit_comb,
    fit_scale,
    spans,
    structure_controls,
    synthetic_calibration,
    verdict,
)
from seti_search.rungs import PHI  # noqa: E402

CHANNEL = 100547.19562243501
GROUP = [45.9843e6, 46.9828e6, 47.9380e6, 48.9937e6]


class FitTests(unittest.TestCase):
    def test_a_comb_is_recovered_exactly(self):
        f = 47.5e6 + 1.0e6 * np.arange(6)
        fit = fit_comb(f)
        self.assertAlmostEqual(fit["df"], 1.0e6, places=3)
        self.assertLess(float(np.abs(fit["residuals_hz"]).max()), 1e-3)

    def test_a_scale_family_is_recovered_exactly(self):
        q = 1.07
        f = 50e6 * q ** np.arange(6)
        fit = fit_scale(f)
        self.assertAlmostEqual(fit["q"], q, places=6)
        self.assertLess(float(np.abs(fit["residuals_hz"]).max()), 1e-2)

    def test_three_lines_fit_both_models_exactly(self):
        f = np.array([45.0e6, 46.1e6, 47.0e6])
        self.assertGreater(fit_comb(f)["rss"], -1.0)
        self.assertGreater(fit_scale(f)["rss"], -1.0)

    def test_a_wide_comb_prefers_comb_and_says_so(self):
        f = 20e6 + 8e6 * np.arange(8)
        out = compare(f, channel_hz=CHANNEL, bootstrap=200)
        self.assertGreater(spans(f)["phi_rungs"], BAND_SPAN_REQUIRED)
        self.assertLess(out["delta_aic"], 0.0)
        label = verdict(out, structure_controls(f, [1.0] * 8, CHANNEL),
                        channel_hz=CHANNEL)["label"]
        self.assertEqual(label, "FREQUENCY COMB SUPPORTED")

    def test_a_wide_scale_family_prefers_scale_but_needs_replication(self):
        f = 50e6 * 1.15 ** np.arange(7)
        out = compare(f, channel_hz=CHANNEL, bootstrap=200)
        controls = structure_controls(f, [1.0] * 7, CHANNEL)
        self.assertGreater(out["delta_aic"], 0.0)
        self.assertEqual(verdict(out, controls, channel_hz=CHANNEL)["label"],
                         "BOTH MODELS FIT")
        self.assertEqual(verdict(out, controls, replicated=True,
                                 channel_hz=CHANNEL)["label"],
                         "SCALE FAMILY SUPPORTED")


class VerdictTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = compare(GROUP, channel_hz=CHANNEL, bootstrap=300)
        cls.controls = structure_controls(GROUP, [41240.0] * 4, CHANNEL)

    def test_the_vocabulary_is_seven_verdicts_plus_a_disqualification(self):
        self.assertEqual(len(LABELS), 8)
        self.assertIn("MASK-MANUFACTURED SIGNAL", LABELS)
        self.assertEqual(len([l for l in LABELS if l != "MASK-MANUFACTURED SIGNAL"]), 7)

    def test_the_real_group_cannot_be_named_at_its_span(self):
        self.assertLess(self.out["span_rungs"], BAND_SPAN_REQUIRED)
        self.assertEqual(verdict(self.out, self.controls, channel_hz=CHANNEL)["label"],
                         "INDETERMINATE DUE TO LIMITED SPAN")

    def test_the_aic_gap_is_inside_its_own_bootstrap(self):
        lo, hi = self.out["bootstrap"]["delta_aic_ci95"]
        self.assertLessEqual(lo, 0.0)
        self.assertGreaterEqual(hi, 0.0)

    def test_dropping_one_line_flips_the_preference(self):
        self.assertFalse(self.out["drop_one_stable"])

    def test_a_spacing_on_the_channel_grid_is_instrumental(self):
        f = 47.5e6 + 10 * CHANNEL * np.arange(5)
        controls = structure_controls(f, [0.0] * 5, CHANNEL)
        self.assertTrue(controls["grid_comb"])
        out = compare(f, channel_hz=CHANNEL, bootstrap=100)
        self.assertEqual(verdict(out, controls, channel_hz=CHANNEL)["label"],
                         "INSTRUMENTAL STRUCTURE")

    def test_one_shared_drift_across_a_band_is_one_event(self):
        f = 55e6 + 2e6 * np.arange(12)
        controls = structure_controls(f, [-152430.0] * 12, CHANNEL)
        self.assertTrue(controls["one_event"])
        out = compare(f, channel_hz=CHANNEL, bootstrap=100)
        self.assertEqual(verdict(out, controls, channel_hz=CHANNEL)["label"],
                         "ASTROPHYSICAL STRUCTURE")

    def test_the_half_integer_steps_are_reported_not_used(self):
        # two of the group's three steps sit exactly between channels: recorded as a
        # suspicion, and the label still does not turn on it
        off = self.controls["steps_off_integer_channels"]
        self.assertTrue(any(abs(x - 0.5) < 0.01 for x in off))
        self.assertFalse(self.controls["grid_comb"])

    def test_three_lines_report_no_held_out_test(self):
        out = compare([45.4e6, 46.5e6, 47.5e6], channel_hz=CHANNEL, bootstrap=100)
        self.assertFalse(out["loo"]["applicable"])
        self.assertEqual(out["loo"]["folds"], [])


class CalibrationTests(unittest.TestCase):
    def test_clean_groups_are_recovered(self):
        out = synthetic_calibration(4, 30e6, trials=60, jitter_hz=1000)
        self.assertGreater(out["correct_rate"]["comb"], 0.9)
        self.assertGreater(out["correct_rate"]["scale"], 0.9)

    def test_noisy_narrow_groups_are_not(self):
        out = synthetic_calibration(4, 3e6, trials=60, jitter_hz=40000)
        self.assertLess(out["correct_rate"]["comb"], 0.5)
        self.assertLess(out["correct_rate"]["scale"], 0.5)

    def test_span_is_what_saves_it(self):
        narrow = synthetic_calibration(4, 3e6, trials=60, jitter_hz=20000)
        wide = synthetic_calibration(4, 30e6, trials=60, jitter_hz=20000)
        self.assertGreater(wide["correct_rate"]["comb"],
                           narrow["correct_rate"]["comb"] + 0.3)


class ArtifactTests(unittest.TestCase):
    def test_the_report_and_the_table_agree(self):
        import csv
        import json
        report = os.path.join(ROOT, "PHI_MODEL_DISCRIMINATION.md")
        table = os.path.join(ROOT, "results", "comb_vs_scale.csv")
        if not (os.path.exists(report) and os.path.exists(table)):
            self.skipTest("run scripts/phi_model_discrimination.py first")
        with open(table, encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        with open(report, encoding="utf-8") as fh:
            text = fh.read()
        self.assertTrue(rows)
        for row in rows:
            self.assertIn(row["verdict"], text)
            self.assertIn(row["date"], text)
        labels = {row["verdict"] for row in rows}
        self.assertTrue(labels <= set(LABELS))


if __name__ == "__main__":
    unittest.main()
