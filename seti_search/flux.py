"""Turning a measured amplitude into a flux density in Jy.

Everything upstream is in units of the file's own noise: a candidate's SNR, and
an event's amplitude, say "this many noise standard deviations". That is enough
to compare scans of one receiver and useless for comparing two. Flux density is
the common currency, and it costs three numbers per scan - channel width,
integration time and system equivalent flux density (SEFD) - plus one assumption
about how wide the line is.

The radiometer equation
-----------------------
For one channel of width ``df`` integrated for ``tsamp`` seconds, the noise
fluctuation expressed as a flux density is

    sigma_S = SEFD / sqrt(df * tsamp)        [Jy]

and a narrowband signal of flux density ``S`` filling that channel therefore
appears at a per-sample signal-to-noise of

    a = S * sqrt(df * tsamp) / SEFD.

So ``S = a * SEFD / sqrt(df * tsamp)``. SEFD is the bridge to hardware:

    SEFD = Tsys / gain = 2 k Tsys / A_eff        [Jy]

with ``2k/Jy = 2761.3 K m^2 / Jy`` for an aperture area in square metres. It is
the flux density that doubles the system temperature, so it folds in Tsys, the
aperture, the receiver and the spillover in one measured number - which is why
it is the quantity to calibrate rather than any of them separately.

The boxcar penalty
------------------
An event's amplitude is defined by ``SNR = A * sqrt(N * W)`` for ``N``
integrations and a boxcar of ``W`` channels, which is the response to a signal
*matched* to the boxcar. A line narrower than that scores lower: the boxcar adds
noise across all ``W`` channels but signal across only ``w`` of them, so

    A = a * w / W,        hence        a = A * W / w.

``w`` is the one thing the pipeline does not measure, so it has to be assumed.
The default is one channel (an unresolved line, the usual technosignature case),
which is the conservative choice: assuming a wide line would inflate the flux by
up to the boxcar width.

What this is not
----------------
A flux density computed this way is only as good as the SEFD handed to it, and
SEFD varies with receiver, elevation, weather and the band. Nothing here
calibrates; it converts with a number you supply. It also assumes the file's
noise is radiometer noise, which a band full of interference is not.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = [
    "COMPLETENESS_MARGIN",
    "K_M2_PER_JY",
    "FluxEstimate",
    "boxcar_scale",
    "combine_fluxes",
    "flux_density_jy",
    "limit_jy",
    "noise_per_sample_jy",
    "run_limit",
    "sefd_from_gain",
    "sefd_from_tsys",
]

K_BOLTZMANN = 1.380649e-23
JY = 1e-26
# 2k/Jy, in K m^2 per Jy: SEFD[Jy] = K_M2_PER_JY * Tsys[K] / A_eff[m^2]
K_M2_PER_JY = 2.0 * K_BOLTZMANN / JY
# 95% one-sided Gaussian: a source this far above the threshold is missed 5% of
# the time, which is what makes an upper limit a promise rather than a rumour.
COMPLETENESS_MARGIN = 1.6448536269514722


def sefd_from_tsys(tsys_k: float, effective_area_m2: float) -> float:
    """SEFD in Jy from system temperature and effective collecting area."""
    if effective_area_m2 <= 0:
        raise ValueError("effective area must be positive")
    return K_M2_PER_JY * float(tsys_k) / float(effective_area_m2)


def sefd_from_gain(tsys_k: float, gain_k_per_jy: float) -> float:
    """SEFD in Jy from system temperature and gain in K/Jy."""
    if gain_k_per_jy <= 0:
        raise ValueError("gain must be positive")
    return float(tsys_k) / float(gain_k_per_jy)


def noise_per_sample_jy(sefd_jy: float, df_hz: float, tsamp_s: float) -> float:
    """Radiometer noise in one channel of one integration, as a flux density."""
    if df_hz <= 0 or tsamp_s <= 0:
        raise ValueError("channel width and integration time must be positive")
    return float(sefd_jy) / math.sqrt(float(df_hz) * float(tsamp_s))


def boxcar_scale(boxcar: int, line_width_channels: float = 1.0) -> float:
    """W/w: how much a line narrower than the boxcar is penalised.

    A line wider than the boxcar gains nothing further, so the ratio is capped
    at 1 for w >= W rather than being allowed below it.
    """
    width = max(float(boxcar), 1.0)
    line = min(max(float(line_width_channels), 1.0), width)
    return width / line


def flux_density_jy(
    *,
    amplitude: float,
    sensitivity: float,
    sefd_jy: float,
    df_hz: float,
    tsamp_s: float,
    boxcar: int,
    line_width_channels: float = 1.0,
) -> tuple[float, float]:
    """Convert one scan's amplitude to (flux density, uncertainty) in Jy.

    ``sensitivity`` is the scan's sqrt(N*W) factor; the amplitude estimate has
    variance 1/sensitivity^2, so its uncertainty in Jy follows the same scale.
    """
    if sensitivity <= 0:
        raise ValueError("sensitivity must be positive")
    scale = (
        boxcar_scale(boxcar, line_width_channels)
        * noise_per_sample_jy(sefd_jy, df_hz, tsamp_s)
    )
    return amplitude * scale, scale / float(sensitivity)


def limit_jy(
    *,
    threshold: float,
    sensitivity: float,
    sefd_jy: float,
    df_hz: float,
    tsamp_s: float,
    boxcar: int,
    line_width_channels: float = 1.0,
    margin_sigma: float = 0.0,
) -> float:
    """Flux density at which a source would reach ``threshold + margin_sigma``.

    This is the other half of a non-detection. If nothing was found above the
    threshold, then nothing brighter than this is present in the searched band -
    at ``margin_sigma = 0`` a source sitting exactly here is caught about half
    the time, so quote ``COMPLETENESS_MARGIN`` alongside it for a limit that
    holds nearly always.
    """
    if sensitivity <= 0:
        raise ValueError("sensitivity must be positive")
    scale = (
        boxcar_scale(boxcar, line_width_channels)
        * noise_per_sample_jy(sefd_jy, df_hz, tsamp_s)
    )
    return (float(threshold) + float(margin_sigma)) * scale / float(sensitivity)


def run_limit(limits: dict[str, float], min_scans: int) -> tuple[float, str] | None:
    """The flux a whole run rules out, given each ON scan's own limit.

    A run only claims a source that shows up in ``min_scans`` scans, so the
    depth that matters is the min_scans-th shallowest scan, not the deepest: a
    source brighter than that is seen enough times to be linked, while a
    fainter one might reach only the best scans and be discarded.
    """
    if min_scans <= 0 or len(limits) < min_scans:
        return None
    ordered = sorted((value, name) for name, value in limits.items())
    value, name = ordered[min_scans - 1]
    return float(value), name


@dataclass
class FluxEstimate:
    """A flux density combined across scans, in Jy."""

    jy: float
    uncertainty_jy: float
    n_scans: int
    per_scan: dict[str, float]

    def __str__(self) -> str:
        return (
            f"{self.jy:.3g} +/- {self.uncertainty_jy:.2g} Jy"
            f" from {self.n_scans} scan(s)"
        )


def combine_fluxes(
    estimates: list[tuple[float, float, str]],
    *,
    line_width_channels: float = 1.0,
) -> FluxEstimate | None:
    """Inverse-variance combination of per-scan flux densities.

    Taking ``estimates`` as (flux, sigma, scan) triples. This is the whole point
    of calibrating: a 10 Jy receiver and a 100 Jy receiver contribute in
    proportion to their variances, so the shallower scan nudges the answer
    instead of dragging it.
    """
    usable = [(f, s, name) for f, s, name in estimates if s > 0 and np.isfinite(f)]
    if not usable:
        return None
    weights = np.array([1.0 / sigma**2 for _, sigma, _ in usable], dtype=np.float64)
    values = np.array([flux for flux, _, _ in usable], dtype=np.float64)
    total = float(weights.sum())
    combined = float(np.sum(values * weights) / total)
    return FluxEstimate(
        jy=combined,
        uncertainty_jy=float(1.0 / math.sqrt(total)),
        n_scans=len(usable),
        per_scan={name: float(flux) for flux, _, name in usable},
    )
