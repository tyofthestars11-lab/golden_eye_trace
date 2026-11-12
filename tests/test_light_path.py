"""The row test must not throw away most of the data, and the ledger must say so.

The synthetic case is built the way a solar dynamic spectrum is built: one
correlated gain across all channels, varying in time. Independent-channel
arithmetic cannot describe that, and the test pins the consequence.
"""

import csv
import json
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from seti_search.search import excise_rfi  # noqa: E402


class RowTestTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(11)
        n_t, n_f = 400, 64
        base = 1.0 + 0.05 * rng.standard_normal((n_t, n_f))
        common = 1.0 + 0.25 * np.sin(np.linspace(0, 6 * np.pi, n_t))[:, None]
        self.data = (base * common).astype(np.float32)

    def test_theoretical_scale_blanks_most_integrations(self):
        clean, rfi = excise_rfi(self.data, sigma=6.0, carrier_sigma=20.0,
                                row_noise="theory")
        self.assertGreater(rfi["blanked_rows"], 0.25 * self.data.shape[0])

    def test_empirical_scale_keeps_them(self):
        clean, rfi = excise_rfi(self.data, sigma=6.0, carrier_sigma=20.0,
                                row_noise="empirical")
        self.assertLessEqual(rfi["blanked_rows"], 2)
        self.assertEqual(rfi["row_noise"], "empirical")

    def test_an_actual_glitch_is_still_caught(self):
        glitchy = self.data.copy()
        glitchy[7, :] *= 12.0
        _, rfi = excise_rfi(glitchy, sigma=6.0, carrier_sigma=20.0,
                            row_noise="empirical")
        self.assertGreaterEqual(rfi["blanked_rows"], 1)

    def test_the_default_is_unchanged(self):
        _, rfi = excise_rfi(self.data)
        self.assertEqual(rfi["row_noise"], "theory")


class LedgerTests(unittest.TestCase):
    def test_the_ledger_carries_the_three_masks(self):
        path = os.path.join(ROOT, "results", "light_path_ledger.csv")
        if not os.path.exists(path):
            self.skipTest("run scripts/light_path.py first")
        with open(path, encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        masks = {row["mask"] for row in rows}
        self.assertEqual(masks, {"none", "theory", "empirical"})

    def test_the_group_is_not_a_candidate_under_an_honest_mask(self):
        path = os.path.join(ROOT, "results", "light_path_ledger.csv")
        if not os.path.exists(path):
            self.skipTest("run scripts/light_path.py first")
        with open(path, encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        def sigmas(mask, window, day="20240325"):
            return [float(r["E_sigma_calibrated"]) for r in rows
                    if r["mask"] == mask and r["window"] == window and r["date"] == day]
        theory = sigmas("theory", "96ch")
        empirical = sigmas("empirical", "96ch")
        none = sigmas("none", "96ch")
        self.assertGreater(min(theory), 100.0)
        self.assertLess(max(empirical), 2.0)
        self.assertLess(max(none), 2.0)


if __name__ == "__main__":
    unittest.main()
