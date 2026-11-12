"""Permanent regression: a narrow search window truncates a fast drift trail.

This is the fixture for the invalid 0/10 run. The 20240325 lanes drift at about
0.105 channels per sample, which is 105 channels across a 1000-sample
integration. A 32-channel window cannot hold that trail, so the dedoppler sum is
truncated mid-flight and a 954 sigma line comes back at nothing.

The test asserts the geometry and the consequence, on synthetic data, so it runs
in a second and never needs the 700 MB archive:

  * the lane traverses about 105 channels, so at most about 30% of it is inside a
    32-channel window;
  * the wide window recovers the lane strongly;
  * the narrow window recovers it far more weakly.

If a future change makes the narrow window look fine, either the geometry moved
or the dedoppler started wrapping around, and both are worth failing on.
"""

import unittest

import numpy as np

from seti_search.search import SearchConfig, search_band

CHANNEL_HZ = 100547.19562243501
TSAMP_S = 0.256
RATE = 0.10487              # cycles per sample, the 20240325 lanes
SAMPLES = 1000
CHANNELS = 96


def synthetic_lane(channels: int = CHANNELS, samples: int = SAMPLES,
                   rate: float = RATE, amplitude: float = 60.0,
                   seed: int = 5) -> np.ndarray:
    """Exponential noise with a Gaussian lane drifting at `rate` channels/sample."""
    rng = np.random.default_rng(seed)
    t = np.arange(samples)
    centre = channels / 2.0 + rate * t
    j = np.arange(channels)[None, :]
    data = rng.exponential(1.0, (samples, channels))
    data += amplitude * np.exp(-0.5 * ((j - centre[:, None]) / 1.5) ** 2)
    return data.astype(np.float32)


def search(window: np.ndarray, threshold: float = 5.0):
    freqs = 40.0 + np.arange(window.shape[1]) * CHANNEL_HZ / 1e6
    cfg = SearchConfig(boxcar=3, max_drift_hz_s=0.5 * CHANNEL_HZ / TSAMP_S,
                       snr_threshold=threshold, max_hits=20)
    hits, _meta = search_band(window, freqs, TSAMP_S, cfg)
    return hits


class LaneGeometryTests(unittest.TestCase):
    def test_the_lane_crosses_about_a_hundred_channels(self):
        crossed = RATE * SAMPLES
        self.assertGreater(crossed, 100.0)
        self.assertLess(crossed, 110.0)
        # and therefore only about a third of it can sit in a 32-channel window
        self.assertLess(32.0 / crossed, 0.35)

    def test_survey_bandwidth_recovers_the_lane(self):
        hits = search(synthetic_lane())
        self.assertTrue(hits, msg="the 96-channel window must find this lane")
        self.assertGreater(max(h.snr for h in hits), 20.0)

    def test_narrow_window_loses_most_of_the_lane(self):
        wide = max(h.snr for h in search(synthetic_lane()))
        narrow_window = synthetic_lane()[:, CHANNELS // 2 - 16:CHANNELS // 2 + 16]
        narrow_hits = search(narrow_window)
        narrow = max((h.snr for h in narrow_hits), default=0.0)
        # the trail leaves the window, so the narrow search cannot be as strong
        self.assertLess(narrow, 0.8 * wide,
                        msg=f"narrow {narrow:.1f} vs wide {wide:.1f}")

    def test_the_failure_is_geometry_and_not_the_noise_reference(self):
        # the same 32-channel window searched with a slower lane is fine, so the
        # window itself is not broken, only its width relative to the drift
        slow = synthetic_lane(rate=0.005)
        narrow = slow[:, CHANNELS // 2 - 16:CHANNELS // 2 + 16]
        hits = search(narrow)
        self.assertTrue(hits, msg="a slow lane must survive the narrow window")


if __name__ == "__main__":
    unittest.main()
