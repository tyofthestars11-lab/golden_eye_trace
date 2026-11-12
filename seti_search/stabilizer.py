"""The Lucas stabilizer: a search window has to hold the trail it is chasing.

A dedoppler integrates along a straight line. A lane drifting at ``nu`` channels
per sample sweeps ``|nu| * samples`` channels across one integration, so a window
of ``W`` channels can contain at most ``W / trail`` of that lane's coherent sum.
That is the *best* case, because it assumes the lane enters the window at the
right moment; in general it exits partway through and the rest of the
integration contributes noise on top of a partial sum.

The surviving significance goes as the **square root** of that fraction, and the
reason is an implementation detail worth stating plainly. ``dedoppler`` shifts
each integration by its accumulated offset and drops any integration whose shift
runs past the band rather than wrapping it, so a lane that leaves the window
takes its *noise* with it. With ``n_eff`` of ``n_t`` integrations surviving, the
signal sums to ``A * n_eff`` and the noise to ``sigma * sqrt(n_eff)``, so

    significance = A * n_eff / (sigma * sqrt(n_eff))
                 = sqrt(n_eff) * (A / sigma)
                 = sqrt(fraction) * significance_with_the_whole_trail

A dedoppler that instead kept integrating the noise from integrations the signal
had left would lose significance linearly in the fraction. This one does not, so
the law here is the square root, and ``scripts/drift_sensitivity.py`` measures it
rather than asserting it. If the shift convention in ``dedoppler`` ever changes
to wrap or to pad, that measurement is the thing that will catch the change.

The survey's numbers make the failure concrete. The 20240325 lanes drift
0.10487 channels per sample, which is a 104.9-channel trail across 256 s of
0.256 s integrations. A 96-channel window holds 91.5% of it and loses 4.3% of
the significance - a footnote. A 32-channel window holds 30.5% and loses 44.7%
of it, which is the invalid 0/10 run: the sigma it reported was a
statement about the window rather than about the sky.

So the guard computes four numbers before any search runs - chunk rate, trail,
fraction held, amplitude lost - and refuses the search when the window cannot
hold the trail to within a declared tolerance. The tolerance is on *amplitude*,
not on significance, because the loss is coherent: ``sqrt(fraction)`` of the
signal survives and the noise does not shrink with it.

Why 20%: at that loss the window keeps 64% of the trail (0.8^2), which is a
factor of 0.89 in the trail length the search can hold - the point past which the
stabilizer's own geometry, rather than the noise, is shaping the threshold. It is
also well inside the spread of the local noise scale these candidates are
measured against (0.30 to 31.18 of the band-wide estimate), so a smaller loss is
not resolvable and a larger one has to be carried in the record.

This is a pre-flight check, not a repair. It does not widen the window for you:
silently searching something other than what was asked for is the same class of
error as silently truncating. ``stabilize`` returns the narrowest window that
would hold the trail, so a caller can decide to use it and record that it did.

    check  = check_window(drift_hz_s=41190.2, window_channels=96, samples=1000)
    check.fits, check.fraction_held, check.amplitude_loss

    require_window(...)  # raises WindowTooNarrow instead of returning a verdict
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .signal_record import channel_rate

__all__ = [
    "MAX_SIGNIFICANCE_LOSS",
    "WindowCheck",
    "WindowTooNarrow",
    "apparent_rate",
    "check_window",
    "fastest_rate_hz_s",
    "fraction_held",
    "require_window",
    "significance_loss",
    "stabilize",
    "trail_channels",
]

#: Largest share of a lane's significance a window may remove before the search
#: is refused. See the module docstring for why 20%.
MAX_SIGNIFICANCE_LOSS = 0.20

#: A dedoppler cannot resolve a rate faster than half a channel per sample; past
#: that the drift folds and the apparent trail is shorter, not longer.
NYQUIST_CHANNELS_PER_SAMPLE = 0.5


def apparent_rate(drift_hz_per_s: float, channel_hz: float, tsamp_s: float) -> float:
    """Drift in cycles per sample, folded into the one zone a dedoppler can see.

    ``signal_record.channel_rate`` folds a rate once, which is right for the
    rates the survey searched (all inside one Nyquist zone) but not for an
    arbitrary input: 2546 cycles per sample would come back as 2545. Taking the
    modulo of the once-folded value puts any input in [-0.5, 0.5), where the
    apparent trail cannot exceed half the integration.
    """
    nu = channel_rate(drift_hz_per_s, channel_hz, tsamp_s)
    return (nu + 0.5) % 1.0 - 0.5


def trail_channels(drift_hz_per_s: float, channel_hz: float, tsamp_s: float,
                   samples: int) -> float:
    """How many channels a lane at ``drift_hz_per_s`` crosses in one integration.

    The rate is folded first, so a 0.4965 cycle/sample lane - which *is* a
    -0.5035 cycle/sample lane seen from the other side - gets the same trail as
    any other rate near Nyquist.
    """
    if samples <= 0:
        raise ValueError(f"samples must be positive, got {samples}")
    return abs(apparent_rate(drift_hz_per_s, channel_hz, tsamp_s)) * float(samples)


def fraction_held(window_channels: float, trail: float) -> float:
    """Best-case share of the integration during which the lane is inside.

    A lane with no drift is inside for the whole integration, whatever the
    window. A lane longer than the window leaves it, and cannot be present for
    more than ``window / trail`` of the time.
    """
    if window_channels <= 0:
        raise ValueError(f"window must be positive, got {window_channels}")
    if trail <= 0.0:
        return 1.0
    return min(1.0, float(window_channels) / trail)


def significance_loss(fraction: float) -> float:
    """Significance removed by holding only ``fraction`` of the trail.

    Square root, because ``dedoppler`` drops the integrations whose shift leaves
    the band, so the signal and the noise are truncated together and only the
    surviving ``sqrt(n_eff)`` of each raises the significance. See the module
    docstring for the derivation and for the linear law that would apply if the
    noise were not truncated with the signal.
    """
    if not 0.0 <= fraction <= 1.0:
        raise ValueError(f"fraction must be in [0, 1], got {fraction}")
    return 1.0 - math.sqrt(fraction)


def fastest_rate_hz_s(window_channels: float, samples: int, channel_hz: float,
                      tsamp_s: float,
                      max_significance_loss: float = MAX_SIGNIFICANCE_LOSS) -> float:
    """The fastest drift this window can search without exceeding the tolerance.

    This is the honest drift range of a search, as opposed to the Nyquist range
    the dedoppler will happily evaluate. For the survey's 96-channel, 1000-sample
    window it is 58.9 kHz/s, while the dedoppler runs to the fold at 196.4 kHz/s:
    a rate between those two is searched, found and under-measured.
    """
    if not 0.0 <= max_significance_loss < 1.0:
        raise ValueError("max_significance_loss must be in [0, 1)")
    if samples <= 0 or channel_hz <= 0 or tsamp_s <= 0 or window_channels <= 0:
        raise ValueError("window, samples, channel width and sample time must be positive")
    hold = (1.0 - max_significance_loss) ** 2
    nu = min(NYQUIST_CHANNELS_PER_SAMPLE,
             float(window_channels) / (float(samples) * hold))
    return nu * channel_hz / tsamp_s


def narrowest_window(samples: int, channel_hz: float, tsamp_s: float,
                     drift_hz_per_s: float,
                     max_significance_loss: float = MAX_SIGNIFICANCE_LOSS) -> int:
    """The narrowest window that holds this lane to the tolerance, in channels."""
    hold = (1.0 - max_significance_loss) ** 2
    trail = trail_channels(drift_hz_per_s, channel_hz, tsamp_s, samples)
    return max(1, int(math.ceil(trail * hold - 1e-9)))


@dataclass(frozen=True)
class WindowCheck:
    """The geometry of one search window against one drift rate."""

    window_channels: int
    samples: int
    channel_rate: float
    trail_channels: float
    fraction_held: float
    significance_loss: float
    fits: bool
    verdict: str
    note: str

    def as_dict(self) -> dict:
        return {
            "window_channels": self.window_channels,
            "samples": self.samples,
            "channel_rate": self.channel_rate,
            "trail_channels": self.trail_channels,
            "fraction_held": self.fraction_held,
            "significance_loss": self.significance_loss,
            "fits": self.fits,
            "verdict": self.verdict,
            "note": self.note,
        }


def check_window(drift_hz_per_s: float, window_channels: int, samples: int,
                 channel_hz: float, tsamp_s: float,
                 max_significance_loss: float = MAX_SIGNIFICANCE_LOSS) -> WindowCheck:
    """Decide whether this window can hold a lane at this drift."""
    if window_channels <= 0:
        raise ValueError(f"window must be positive, got {window_channels}")
    nu = apparent_rate(drift_hz_per_s, channel_hz, tsamp_s)
    trail = abs(nu) * float(samples)
    held = fraction_held(window_channels, trail)
    lost = significance_loss(held)
    fits = lost <= max_significance_loss + 1e-12
    if fits:
        note = (f"{held:.1%} of the trail fits in the window, so the dedoppler "
                f"loses at most {lost:.1%} of its significance to the window edge")
        verdict = "HOLDS"
    else:
        note = (f"a {trail:.0f}-channel trail in a {window_channels}-channel "
                f"window: at most {held:.1%} of the integration can contain the "
                f"lane, so a sigma from this window is a statement about the "
                f"window (limit {max_significance_loss:.0%} of the significance)")
        verdict = "TRUNCATES"
    return WindowCheck(
        window_channels=int(window_channels),
        samples=int(samples),
        channel_rate=nu,
        trail_channels=trail,
        fraction_held=held,
        significance_loss=lost,
        fits=fits,
        verdict=verdict,
        note=note,
    )


class WindowTooNarrow(ValueError):
    """Raised instead of searching a window that cannot hold the trail."""

    def __init__(self, check: WindowCheck) -> None:
        self.check = check
        super().__init__(f"{check.verdict}: {check.note}")


def require_window(drift_hz_per_s: float, window_channels: int, samples: int,
                   channel_hz: float, tsamp_s: float,
                   max_significance_loss: float = MAX_SIGNIFICANCE_LOSS) -> WindowCheck:
    """``check_window``, but refusing to continue on a truncated search."""
    check = check_window(drift_hz_per_s, window_channels, samples, channel_hz,
                         tsamp_s, max_significance_loss)
    if not check.fits:
        raise WindowTooNarrow(check)
    return check


def stabilize(drift_hz_per_s: float, samples: int, channel_hz: float,
              tsamp_s: float,
              max_significance_loss: float = MAX_SIGNIFICANCE_LOSS) -> WindowCheck:
    """The narrowest window that would hold this lane, checked.

    Returned rather than applied: a caller that widens its own search has to say
    so, because the wider window tests more cells and reports a different sigma.
    """
    width = narrowest_window(samples, channel_hz, tsamp_s, drift_hz_per_s,
                             max_significance_loss)
    return check_window(drift_hz_per_s, width, samples, channel_hz, tsamp_s,
                        max_significance_loss)
