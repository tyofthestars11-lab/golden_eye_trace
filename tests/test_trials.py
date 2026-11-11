import math
import unittest
from types import SimpleNamespace

import numpy as np

from seti_search.cli import _event_budget, _rejection_budget
from seti_search.search import dedoppler
from seti_search.trials import (
    drift_trials,
    expected_false_alarms,
    log10_tail,
    sigma_for_band_false_alarm,
    tail_probability,
)


class TailTests(unittest.TestCase):
    def test_known_probabilities(self):
        self.assertAlmostEqual(tail_probability(0.0), 0.5)
        self.assertAlmostEqual(tail_probability(1.0), 0.15865525393145707, places=12)
        self.assertAlmostEqual(tail_probability(5.0), 2.866515718791939e-07, places=18)
        self.assertAlmostEqual(tail_probability(-1.0), 1.0 - 0.15865525393145707, places=12)

    def test_stricter_threshold_means_fewer_false_alarms(self):
        self.assertLess(tail_probability(7.0), tail_probability(6.0))
        self.assertLess(
            expected_false_alarms(5.0, 10**8), expected_false_alarms(4.0, 10**8)
        )

    def test_log_form_survives_underflow(self):
        self.assertAlmostEqual(log10_tail(5.0), math.log10(tail_probability(5.0)))
        # a 60 sigma threshold is around 1e-783: no float holds it, the log does
        self.assertEqual(tail_probability(60.0), 0.0)
        self.assertLess(log10_tail(60.0), -770.0)
        self.assertGreater(log10_tail(60.0), -790.0)


class ThresholdTests(unittest.TestCase):
    def test_inverse_hits_the_requested_band_rate(self):
        for probability in (0.5, 0.01, 1e-4):
            cells = 1_000_000
            sigma = sigma_for_band_false_alarm(probability, cells)
            realised = 1.0 - (1.0 - tail_probability(sigma)) ** cells
            self.assertAlmostEqual(realised / probability, 1.0, places=6)

    def test_more_trials_need_a_higher_threshold(self):
        few = sigma_for_band_false_alarm(0.01, 10**4)
        many = sigma_for_band_false_alarm(0.01, 10**9)
        self.assertGreater(many, few)
        self.assertGreater(many - few, 0.4)

    def test_rejects_impossible_requests(self):
        for bad in (0.0, 1.0, -0.1, 1.5):
            with self.assertRaises(ValueError):
                sigma_for_band_false_alarm(bad, 100)
        with self.assertRaises(ValueError):
            sigma_for_band_false_alarm(0.01, 0)


class DriftTrialTests(unittest.TestCase):
    def test_matches_the_search_it_describes(self):
        # the same grid dedoppler builds, counted rather than assumed
        data = np.zeros((200, 64), dtype=np.float32)
        for boxcar, max_drift, step in ((3, 4.0, None), (1, 4.0, 0.01)):
            _, drifts = dedoppler(
                data,
                1.0,
                1e-6,
                max_drift_hz_s=max_drift,
                drift_step_hz_s=step,
                boxcar=boxcar,
            )
            self.assertEqual(
                drift_trials(
                    n_ints=200,
                    tsamp_s=1.0,
                    df_mhz=1e-6,
                    boxcar=boxcar,
                    max_drift_hz_s=max_drift,
                    drift_step_hz_s=step,
                ),
                drifts.size,
            )
            self.assertGreater(drifts.size, 10)

    def test_nyquist_clips_the_range(self):
        # half a channel per integration is as fast as shift-and-add can follow,
        # so asking for more trials cannot buy any
        self.assertEqual(
            drift_trials(n_ints=16, tsamp_s=1.0, df_mhz=1e-3, boxcar=1, max_drift_hz_s=1e9),
            drift_trials(n_ints=16, tsamp_s=1.0, df_mhz=1e-3, boxcar=1, max_drift_hz_s=500.0),
        )

    def test_longer_observations_buy_stricter_steps(self):
        short = drift_trials(n_ints=100, tsamp_s=1.0, df_mhz=1e-6, boxcar=3, max_drift_hz_s=0.5)
        long = drift_trials(n_ints=800, tsamp_s=1.0, df_mhz=1e-6, boxcar=3, max_drift_hz_s=0.5)
        self.assertGreater(long, short)

    def test_rejects_degenerate_inputs(self):
        with self.assertRaises(ValueError):
            drift_trials(n_ints=1, tsamp_s=1.0, df_mhz=1e-6)
        with self.assertRaises(ValueError):
            drift_trials(n_ints=10, tsamp_s=0.0, df_mhz=1e-6)


class RejectionTests(unittest.TestCase):
    def test_a_loud_hit_is_called_interference(self):
        text = _rejection_budget(SimpleNamespace(snr=30.5), 10.0, 1_089_536)
        self.assertIn("interference", text)
        self.assertNotIn("not decisive", text)

    def test_a_borderline_hit_is_called_undecided(self):
        # 4 sigma over ten million cells: noise does reach this strength, so a
        # matching OFF peak does not prove a transmitter
        text = _rejection_budget(SimpleNamespace(snr=4.0), 4.0, 10_000_000)
        self.assertIn("not decisive", text)
        self.assertNotIn("interference", text)

    def test_it_quotes_both_budgets(self):
        text = _rejection_budget(SimpleNamespace(snr=12.0), 10.0, 5000)
        self.assertIn("5,000 cells", text)
        self.assertIn("10.00 sigma threshold", text)
        self.assertIn("12.0 sigma hit", text)


META = {
    "n_ints": 200,
    "tsamp": 1.0,
    "df_mhz": 1e-6,
    "boxcar": 3,
    "nchans": 1024,
}
# 200 integrations at 1 s and 1 Hz channels, boxcar 3, Nyquist 0.5 Hz/s: 133
# drift trials per channel, so 133 * 1024 cells in a band this size
CELLS = 133 * 1024


def _scan(path: str):
    return SimpleNamespace(path=path, meta=dict(META))


def _member(path: str, snr: float):
    return SimpleNamespace(scan=path, candidate=SimpleNamespace(snr=snr))


def _event(*members):
    return SimpleNamespace(members=list(members), off_matches=[])


class EventBudgetTests(unittest.TestCase):
    args = SimpleNamespace(max_drift=4.0, drift_step=None, snr=10.0)
    by_path = {"a.fil": _scan("a.fil"), "b.fil": _scan("b.fil")}

    def test_one_scan_is_the_single_cell_chance(self):
        budget = _event_budget(_event(_member("a.fil", 10.0)), self.by_path, self.args)
        self.assertEqual(budget["cells_total"], CELLS)
        self.assertEqual(budget["per_scan_cells"], CELLS)
        self.assertAlmostEqual(
            budget["chance_log10"], math.log10(tail_probability(10.0) * CELLS), places=4
        )
        # 10 sigma over this band is about one chance in 10^18
        self.assertAlmostEqual(budget["chance_log10"], -17.98, places=2)

    def test_a_second_scan_multiplies_the_chance(self):
        one = _event_budget(_event(_member("a.fil", 10.0)), self.by_path, self.args)
        two = _event_budget(
            _event(_member("a.fil", 10.0), _member("b.fil", 10.0)),
            self.by_path,
            self.args,
        )
        self.assertAlmostEqual(two["chance_log10"], 2 * one["chance_log10"], places=4)
        self.assertEqual(two["cells_total"], 2 * CELLS)

    def test_an_off_match_is_quoted_per_window_not_per_band(self):
        event = _event(_member("a.fil", 10.0))
        event.off_matches = [
            SimpleNamespace(
                scan="b.fil",
                hit=SimpleNamespace(snr=8.0),
                # ten channels wide at 1 Hz spacing, at any drift
                window_hz=10.0,
            )
        ]
        budget = _event_budget(event, self.by_path, self.args)
        off = budget["off"]
        self.assertIsNotNone(off)
        self.assertEqual(off["channels"], 10)
        self.assertEqual(off["cells"], 10 * 133)
        self.assertLess(off["cells"], budget["cells_total"])
        self.assertAlmostEqual(
            off["chance_log10"], math.log10(tail_probability(8.0) * 10 * 133), places=4
        )
        self.assertLess(off["chance_log10"], 0.0)

    def test_a_survivor_has_no_off_budget(self):
        budget = _event_budget(_event(_member("a.fil", 10.0)), self.by_path, self.args)
        self.assertIsNone(budget["off"])

    def test_the_threshold_budget_is_reported_per_scan(self):
        budget = _event_budget(
            _event(_member("a.fil", 10.0), _member("b.fil", 10.0)),
            self.by_path,
            SimpleNamespace(max_drift=4.0, drift_step=None, snr=5.6),
        )
        self.assertAlmostEqual(budget["sigma"], 5.6)
        self.assertAlmostEqual(
            budget["expected_per_scan"], tail_probability(5.6) * CELLS, places=12
        )
        self.assertGreater(budget["expected_per_scan"], 0.0)


if __name__ == "__main__":
    unittest.main()
