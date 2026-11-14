"""Proportions must behave like proportions: monotone, bounded, scale-aware."""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from proportion_not_binary import (  # noqa: E402
    BAND_HIGH_HZ,
    BAND_LOW_HZ,
    CHANNEL_HZ,
    TARGET,
    channel_check,
    crossover_frequency,
    delta_r,
    posterior_table,
)

PHI = (1.0 + math.sqrt(5.0)) / 2.0


class PosteriorTests(unittest.TestCase):
    def test_zero_hits_is_not_a_verdict(self):
        post = posterior_table(0, 14)
        self.assertEqual(post["posterior"], "Beta(1, 15)")
        self.assertGreater(post["median"], 0.0)
        self.assertLess(post["median"], 0.1)
        self.assertGreater(post["upper_95"], post["median"])

    def test_the_survival_is_monotone_in_the_rate(self):
        post = posterior_table(0, 14)
        self.assertGreater(post["p_rate_above_0.02"],
                           post["p_rate_above_0.10"])
        self.assertGreater(post["p_rate_above_0.10"],
                           post["p_rate_above_0.381966"])

    def test_the_closed_form_matches_for_hits_out_of_pairs(self):
        # Beta(1, b) has CDF 1 - (1-x)^b, so the 95% bound is 1 - 0.05^(1/b)
        post = posterior_table(0, 9)
        self.assertAlmostEqual(post["upper_95"], 1 - 0.05 ** (1 / 10), places=6)

    def test_more_hits_raise_the_rate(self):
        low = posterior_table(1, 20)["median"]
        high = posterior_table(5, 20)["median"]
        self.assertGreater(high, low)

    def test_impossible_counts_are_refused(self):
        with self.assertRaises(ValueError):
            posterior_table(5, 3)


class ScaleTests(unittest.TestCase):
    def test_a_multiplicative_step_is_scale_invariant(self):
        ratio = PHI ** TARGET
        for freq in (BAND_LOW_HZ, 47e6, BAND_HIGH_HZ):
            step = (ratio - 1.0) * freq
            self.assertAlmostEqual(delta_r(step, freq), TARGET, places=12)

    def test_an_additive_step_is_not(self):
        low = delta_r(CHANNEL_HZ, BAND_LOW_HZ)
        high = delta_r(CHANNEL_HZ, BAND_HIGH_HZ)
        self.assertGreater(low / high, 6.0)

    def test_the_channel_check_spans_only_the_search_window(self):
        # 96 channels of 100.5 kHz is 9.65 MHz, so the rung interval moves by
        # the ratio of the window edges, about 1.3 - not the 6.5 of the full band
        block = channel_check()
        self.assertGreater(block["spread_factor"], 1.2)
        self.assertLess(block["spread_factor"], 2.0)
        whole_band = (delta_r(CHANNEL_HZ, BAND_LOW_HZ)
                      / delta_r(CHANNEL_HZ, BAND_HIGH_HZ))
        self.assertGreater(whole_band, 6.0)

    def test_the_crossover_frequency_solves_the_equation(self):
        block = crossover_frequency()
        self.assertAlmostEqual(
            delta_r(block["step_hz"], block["frequency_hz"]), TARGET, places=9)
        # and it is below the observed band, which is why nothing hit
        self.assertLess(block["frequency_hz"], BAND_LOW_HZ)


if __name__ == "__main__":
    unittest.main()
