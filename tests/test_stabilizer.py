"""The Lucas stabilizer refuses a window that cannot hold the drift trail."""

import unittest

from seti_search.stabilizer import (
    MAX_AMPLITUDE_LOSS,
    WindowTooNarrow,
    apparent_rate,
    check_window,
    fastest_rate_hz_s,
    fraction_held,
    narrowest_window,
    require_window,
    stabilize,
    trail_channels,
)

CHANNEL_HZ = 100547.19562243501
TSAMP_S = 0.256
SAMPLES = 1000
LANE = 41190.2        # Hz/s, the 20240325 burst group
ALIAS = 195015.7      # Hz/s, the 20240920 trio, 0.4965 channels per sample
SLOW = 500.0          # Hz/s, an ordinary slow line


class TrailGeometryTests(unittest.TestCase):
    def test_the_burst_group_crosses_105_channels(self):
        trail = trail_channels(LANE, CHANNEL_HZ, TSAMP_S, SAMPLES)
        self.assertAlmostEqual(trail, 104.9, delta=0.2)

    def test_the_fold_makes_a_negative_drift_the_same_trail(self):
        self.assertAlmostEqual(
            trail_channels(-LANE, CHANNEL_HZ, TSAMP_S, SAMPLES),
            trail_channels(LANE, CHANNEL_HZ, TSAMP_S, SAMPLES),
            places=9,
        )

    def test_the_alias_sits_just_under_half_a_channel_per_sample(self):
        trail = trail_channels(ALIAS, CHANNEL_HZ, TSAMP_S, SAMPLES)
        self.assertAlmostEqual(trail / SAMPLES, 0.4965, delta=0.001)

    def test_no_rate_can_cross_more_than_half_a_channel_per_sample(self):
        # the fold bounds the *apparent* trail, so a 502-channel window holds any
        # rate at all: that is the width at which this guard can never fire
        for rate in (LANE, ALIAS, 1e6, 1e9, -7.7e8):
            self.assertLessEqual(
                trail_channels(rate, CHANNEL_HZ, TSAMP_S, SAMPLES), 500.0)
        for rate in (1e6, 1e9, -7.7e8):
            self.assertEqual(
                check_window(rate, 502, SAMPLES, CHANNEL_HZ, TSAMP_S).verdict,
                "HOLDS")

    def test_the_fold_lands_inside_one_zone(self):
        self.assertAlmostEqual(apparent_rate(ALIAS, CHANNEL_HZ, TSAMP_S),
                               0.4965, delta=0.001)
        # one whole channel per sample is 392,762 Hz/s here: adding cycles does
        # not change which lane the dedoppler sees
        one_cycle_hz_s = CHANNEL_HZ / TSAMP_S
        self.assertAlmostEqual(
            apparent_rate(ALIAS + 1234 * one_cycle_hz_s, CHANNEL_HZ, TSAMP_S),
            apparent_rate(ALIAS, CHANNEL_HZ, TSAMP_S), delta=1e-6)
        self.assertAlmostEqual(
            apparent_rate(-ALIAS - 99 * one_cycle_hz_s, CHANNEL_HZ, TSAMP_S),
            apparent_rate(-ALIAS, CHANNEL_HZ, TSAMP_S), delta=1e-6)

    def test_a_stationary_line_needs_no_window_at_all(self):
        self.assertEqual(trail_channels(0.0, CHANNEL_HZ, TSAMP_S, SAMPLES), 0.0)
        self.assertEqual(fraction_held(1, 0.0), 1.0)

    def test_zero_samples_is_refused(self):
        with self.assertRaises(ValueError):
            trail_channels(LANE, CHANNEL_HZ, TSAMP_S, 0)


class WindowVerdictTests(unittest.TestCase):
    def test_the_survey_window_holds_the_burst_group(self):
        check = check_window(LANE, 96, SAMPLES, CHANNEL_HZ, TSAMP_S)
        self.assertTrue(check.fits)
        self.assertEqual(check.verdict, "HOLDS")
        self.assertAlmostEqual(check.fraction_held, 0.915, delta=0.005)
        self.assertAlmostEqual(check.amplitude_loss, 0.043, delta=0.005)
        self.assertIn("4.3%", check.note)

    def test_the_narrow_window_is_refused_for_the_burst_group(self):
        check = check_window(LANE, 32, SAMPLES, CHANNEL_HZ, TSAMP_S)
        self.assertFalse(check.fits)
        self.assertEqual(check.verdict, "TRUNCATES")
        self.assertAlmostEqual(check.fraction_held, 0.305, delta=0.005)
        with self.assertRaises(WindowTooNarrow) as caught:
            require_window(LANE, 32, SAMPLES, CHANNEL_HZ, TSAMP_S)
        message = str(caught.exception)
        self.assertIn("105-channel trail", message)
        self.assertIn("32-channel window", message)
        self.assertIs(caught.exception.check.verdict, check.verdict)

    def test_the_same_narrow_window_is_fine_for_a_slow_line(self):
        # the failure is the width relative to the drift, not the width
        self.assertTrue(check_window(SLOW, 32, SAMPLES, CHANNEL_HZ, TSAMP_S).fits)

    def test_the_alias_trio_is_window_limited_in_the_survey_window(self):
        check = check_window(ALIAS, 96, SAMPLES, CHANNEL_HZ, TSAMP_S)
        self.assertFalse(check.fits)
        self.assertLess(check.fraction_held, 0.20)

    def test_the_tolerance_is_the_boundary_it_claims_to_be(self):
        for window in (32, 48, 96, 128):
            check = check_window(LANE, window, SAMPLES, CHANNEL_HZ, TSAMP_S)
            self.assertEqual(check.fits,
                             check.amplitude_loss <= MAX_AMPLITUDE_LOSS + 1e-12)

    def test_a_stricter_tolerance_refuses_more(self):
        self.assertTrue(check_window(LANE, 96, SAMPLES, CHANNEL_HZ, TSAMP_S).fits)
        strict = check_window(LANE, 96, SAMPLES, CHANNEL_HZ, TSAMP_S,
                              max_amplitude_loss=0.01)
        self.assertFalse(strict.fits)

    def test_a_zero_window_is_a_bug_not_a_verdict(self):
        with self.assertRaises(ValueError):
            check_window(LANE, 0, SAMPLES, CHANNEL_HZ, TSAMP_S)


class RateLimitTests(unittest.TestCase):
    def test_the_limit_is_where_the_tolerance_is_satisfied_exactly(self):
        limit = fastest_rate_hz_s(96, SAMPLES, CHANNEL_HZ, TSAMP_S)
        self.assertAlmostEqual(limit / 1000.0, 58.9, delta=0.2)
        at = check_window(limit, 96, SAMPLES, CHANNEL_HZ, TSAMP_S)
        self.assertAlmostEqual(at.amplitude_loss, MAX_AMPLITUDE_LOSS, places=9)
        self.assertTrue(at.fits)
        just_over = check_window(limit * 1.05, 96, SAMPLES, CHANNEL_HZ, TSAMP_S)
        self.assertFalse(just_over.fits)

    def test_the_searched_drift_range_is_wider_than_the_window_can_hold(self):
        # the survey ran to the fold, 196.4 kHz/s, with a window that holds
        # 58.9 kHz/s: rates between those are found and under-measured
        nyquist_hz_s = 0.5 * CHANNEL_HZ / TSAMP_S
        limit = fastest_rate_hz_s(96, SAMPLES, CHANNEL_HZ, TSAMP_S)
        self.assertGreater(nyquist_hz_s, 3.0 * limit)
        self.assertTrue(check_window(LANE, 96, SAMPLES, CHANNEL_HZ,
                                     TSAMP_S).fits)
        self.assertFalse(check_window(0.9 * nyquist_hz_s, 96, SAMPLES,
                                      CHANNEL_HZ, TSAMP_S).fits)

    def test_a_narrow_window_has_a_narrower_drift_range(self):
        self.assertAlmostEqual(
            fastest_rate_hz_s(32, SAMPLES, CHANNEL_HZ, TSAMP_S) / 1000.0,
            19.6, delta=0.2)

    def test_a_bad_tolerance_is_refused(self):
        for bad in (-0.1, 1.0, 2.0):
            with self.assertRaises(ValueError):
                fastest_rate_hz_s(96, SAMPLES, CHANNEL_HZ, TSAMP_S, bad)


class StabilizeTests(unittest.TestCase):
    def test_the_narrowest_window_is_the_one_the_tolerance_needs(self):
        for rate in (SLOW, LANE, ALIAS):
            width = narrowest_window(SAMPLES, CHANNEL_HZ, TSAMP_S, rate)
            self.assertTrue(check_window(rate, width, SAMPLES, CHANNEL_HZ,
                                         TSAMP_S).fits)
            if width > 1:
                self.assertFalse(check_window(rate, width - 1, SAMPLES,
                                              CHANNEL_HZ, TSAMP_S).fits)

    def test_stabilize_reports_rather_than_repairs(self):
        check = stabilize(LANE, SAMPLES, CHANNEL_HZ, TSAMP_S)
        self.assertTrue(check.fits)
        self.assertEqual(check.window_channels, 68)
        self.assertIn("fits in the window", check.note)

    def test_the_record_round_trips(self):
        payload = check_window(LANE, 96, SAMPLES, CHANNEL_HZ, TSAMP_S).as_dict()
        self.assertEqual(payload["verdict"], "HOLDS")
        self.assertEqual(payload["window_channels"], 96)
        self.assertAlmostEqual(payload["trail_channels"], 104.9, delta=0.2)


if __name__ == "__main__":
    unittest.main()
