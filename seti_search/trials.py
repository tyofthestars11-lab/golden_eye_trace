"""How many cells a search actually tests, and what that costs in false alarms.

A dedoppler search does not test one number per channel. It tests one number per
channel *per drift trial*, so a 1% chance of a noise peak anywhere in the band
corresponds to a far stricter per-cell threshold than 1% does. This module holds
the arithmetic: the trial count (mirroring `search.dedoppler` exactly), the
Gaussian tail and its inverse.

The Gaussian model is an idealisation - the search normalises with a robust
estimate taken from the same data - and the trials are correlated, since
neighbouring drift steps share integrations and the boxcar smears channels. Both
effects are real and they pull in opposite directions, so a rate quoted here is a
calibrated estimate of the false-alarm budget, not a guarantee.
"""

from __future__ import annotations

import math

__all__ = [
    "drift_trials",
    "expected_false_alarms",
    "log10_tail",
    "sigma_for_band_false_alarm",
    "tail_probability",
]

_SQRT2 = math.sqrt(2.0)
_SQRT_2PI = math.sqrt(2.0 * math.pi)
_LN10 = math.log(10.0)


def drift_trials(
    *,
    n_ints: int,
    tsamp_s: float,
    df_mhz: float,
    boxcar: int = 3,
    max_drift_hz_s: float = 4.0,
    drift_step_hz_s: float | None = None,
) -> int:
    """Drift trials per channel, using the same grid as `search.dedoppler`.

    Keeping this in step with the search matters: if the count drifts, every
    false-alarm rate quoted downstream is quietly wrong.
    """
    n_ints = int(n_ints)
    tsamp_s = float(tsamp_s)
    df_mhz = abs(float(df_mhz))
    if n_ints < 2:
        raise ValueError("need at least two integrations to dedoppler")
    if tsamp_s <= 0 or df_mhz <= 0:
        raise ValueError("integration time and channel width must be positive")

    total_time = n_ints * tsamp_s
    # beyond half a channel per integration the shift-and-add cannot align a
    # signal at all, so the search clips the range here too
    nyquist = 0.5 * df_mhz / tsamp_s * 1e6
    max_drift = min(float(max_drift_hz_s), nyquist)
    if drift_step_hz_s is None:
        step = 0.5 * int(boxcar) * df_mhz * 1e6 / total_time
    else:
        step = float(drift_step_hz_s)
    step = max(step, 1e-5)
    return 2 * int(math.floor(max_drift / step)) + 1


def tail_probability(sigma: float) -> float:
    """P(z > sigma) for a standard normal, or 0.0 once that underflows."""
    return 0.5 * math.erfc(float(sigma) / _SQRT2)


def log10_tail(sigma: float) -> float:
    """log10 of the same tail, finite even where the probability underflows.

    A 60 sigma threshold sits around 1e-783, which no float can hold, but the
    number is still the honest one to print.
    """
    p = tail_probability(sigma)
    if p > 0.0:
        return math.log10(p)
    s = float(sigma)
    # p ~ exp(-s^2/2) / (s * sqrt(2*pi)) in the far tail
    return (-0.5 * s * s - math.log(s * _SQRT_2PI)) / _LN10


def sigma_for_band_false_alarm(probability: float, n_trials: int) -> float:
    """Sigma that gives a band-wide chance ``probability`` of any noise peak.

    The probability is for the whole search, not per cell: with N cells the
    per-cell tail has to be ``1 - (1 - p)**(1/N)`` so that seeing *nothing* is
    the likely outcome.
    """
    probability = float(probability)
    if not 0.0 < probability < 1.0:
        raise ValueError("false-alarm probability must be between 0 and 1")
    n_trials = int(n_trials)
    if n_trials < 1:
        raise ValueError("need at least one trial cell")

    target = -math.expm1(math.log1p(-probability) / n_trials)
    if target <= 0.0:  # too small to represent; the tail is long past caring
        return 60.0
    # the tail falls by about forty decades over 0..40 sigma, so bisection on a
    # bracket that wide is both safe and more than accurate enough
    lo, hi = 0.0, 60.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if tail_probability(mid) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def expected_false_alarms(sigma: float, n_trials: int) -> float:
    """Noise peaks expected above ``sigma`` across ``n_trials`` cells."""
    return tail_probability(sigma) * int(n_trials)
