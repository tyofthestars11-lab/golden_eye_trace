"""Golden-ratio rung addresses, and the only thing they can prove.

Framework: phi^2 = phi + 1. Attribution: Tyree Jones (tyofthestarz), phi^2 =
phi + 1 framework, prior art from December 2025.

One coordinate: ``rung(x) = ln(x) / ln(phi)``. It is a logarithm to base phi, so
it is a ruler that behaves the same at every scale, and ``rung(phi) = 1``,
``rung(1) = 0``, ``rung(phi**2) = 2``. The golden angle falls out of the same
constant: ``360 / 137.507764... = phi**2``, i.e. a turn is exactly phi^2 golden
angles.

What this can prove, and what it cannot
---------------------------------------
A single measurement's rung is *not* evidence of anything. For any smoothly
distributed positive quantity the fractional part of ``ln(x)/ln(phi)`` is
uniform on [0, 1), so every value sits some small distance from *some* integer
rung. ``rung_drift`` of 0.0438 is not special; drawing a random number gives a
drift that small about 9 per cent of the time. Treating a small drift as
"ladder membership" is the standard way this method goes wrong, so nothing here
infers anything from one number.

What is testable is arithmetic identity and repetition:

* :func:`ladder_step` - are two values related by an exact integer number of
  rungs? That is a statement about the ratio, and it is either right to the
  last digit or it is not. A non-integer step is a scale factor, not a ladder
  relation, and the residual says how exactly.
* :func:`repeated_scales` - does the *same* rung difference appear across
  several independent pairs? Then one scaling operation produced them all: the
  values are one object at two scales, or one resolution imposed on many
  objects. That is an artifact, and it is a reason to distrust a candidate, not
  to believe one.
* :func:`rung_uniformity` - across a population, are the fractional rungs
  uniform? If they are not, the population is quantized by the measuring
  pipeline (a channel grid, a drift step, a chunk boundary) rather than by the
  sky.

So these tests can only ever *detract*. A candidate that fails them is
misunderstood; a candidate that passes them has not been shown to be anything.
That asymmetry is the point of running them at all.

Units matter: ``rung`` depends on the unit. ``rung(1420e6 Hz)`` and
``rung(1420 MHz)`` differ by ``rung(1e6)``. Every report here names its unit,
and cross-quantity comparison is only meaningful within one unit.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = [
    "GOLDEN_ANGLE_DEG",
    "LN_PHI",
    "PHI",
    "LadderStep",
    "Repetition",
    "RungAddress",
    "address",
    "delta_rung",
    "ladder_step",
    "repeated_scales",
    "rung",
    "rung_drift",
    "rung_uniformity",
]

PHI = (1.0 + math.sqrt(5.0)) / 2.0
LN_PHI = math.log(PHI)
GOLDEN_ANGLE_DEG = 360.0 / (PHI * PHI)  # 137.5077640...


def rung(x: float) -> float:
    """ln(x)/ln(phi). Zero at 1, one at phi, two at phi^2."""
    value = float(x)
    if not value > 0.0:
        raise ValueError(f"rung needs a positive quantity, got {x!r}")
    return math.log(value) / LN_PHI


def delta_rung(a: float, b: float) -> float:
    """rung(a) - rung(b); equivalently log_phi(a/b)."""
    return rung(a) - rung(b)


def rung_drift(x: float) -> float:
    """Distance from the nearest integer rung, in [0, 0.5].

    Uniform for a smoothly distributed x. Reported because the framework asks
    for it, and to make the point that it is not by itself a finding.
    """
    value = rung(x)
    return abs(value - round(value))


@dataclass(frozen=True)
class RungAddress:
    """A measured quantity with its rung, nearest rung, and drift."""

    value: float
    unit: str
    rung: float
    nearest: int
    drift: float

    def __str__(self) -> str:
        return (
            f"rung({self.value:g} {self.unit}) = {self.rung:.4f}"
            f"  nearest {self.nearest}  drift {self.drift:.4f}"
        )

    def as_dict(self) -> dict:
        return {
            "value": self.value,
            "unit": self.unit,
            "rung": round(self.rung, 6),
            "nearest": self.nearest,
            "drift": round(self.drift, 6),
        }


def address(x: float, unit: str = "") -> RungAddress:
    value = rung(x)
    nearest = int(round(value))
    return RungAddress(
        value=float(x),
        unit=unit,
        rung=value,
        nearest=nearest,
        drift=abs(value - nearest),
    )


@dataclass(frozen=True)
class LadderStep:
    """The relation between two values, in rungs."""

    delta: float
    integer_step: bool
    factor: float
    residual: float  # delta - log_phi(factor); zero for an exact factor
    nearest_integer: int

    def __str__(self) -> str:
        kind = (
            f"integer step {self.nearest_integer}"
            if self.integer_step
            else f"scale factor {self.factor:.6f} (not a ladder step)"
        )
        return f"delta rung {self.delta:+.4f}: {kind}, residual {self.residual:.2e}"


def ladder_step(a: float, b: float, tolerance: float = 1e-6) -> LadderStep:
    """Classify a/b as an integer rung step or a plain scale factor.

    The residual is exact arithmetic: rung(a)-rung(b) always equals
    log_phi(a/b), so a non-zero residual means round-off, never a near miss.
    """
    delta = delta_rung(a, b)
    factor = float(a) / float(b)
    nearest = int(round(delta))
    return LadderStep(
        delta=delta,
        integer_step=abs(delta - nearest) <= tolerance,
        factor=factor,
        residual=delta - math.log(factor) / LN_PHI,
        nearest_integer=nearest,
    )


@dataclass(frozen=True)
class Repetition:
    """One rung difference that several independent pairs share."""

    delta: float
    factor: float
    pairs: int
    examples: tuple[tuple[float, float], ...]

    def __str__(self) -> str:
        shown = ", ".join(f"{a:g}/{b:g}" for a, b in self.examples[:3])
        return (
            f"delta rung {self.delta:+.4f} (factor {self.factor:.6f}) repeats across"
            f" {self.pairs} pairs: {shown}"
        )


def repeated_scales(
    values,
    *,
    min_pairs: int = 3,
    tolerance: float = 1e-4,
) -> list[Repetition]:
    """Find rung differences shared by several independent pairs.

    A repeated *difference* means a repeated *ratio*: the same scaling operation
    applied more than once. In an observing run that is a resolution artifact -
    two hits that are one signal at two pipeline scales - not two sources.

    Differences are taken in magnitude, so a factor and its reciprocal count as
    the same scale; ``Repetition.factor`` is always at least 1.
    """
    array = [float(v) for v in values]
    if len(array) < 3:
        return []
    deltas: list[tuple[float, int, int]] = []
    for i in range(len(array)):
        for j in range(i + 1, len(array)):
            deltas.append((abs(delta_rung(array[i], array[j])), i, j))
    deltas.sort()

    groups: list[list[tuple[float, int, int]]] = []
    for item in deltas:
        if groups and abs(item[0] - groups[-1][-1][0]) <= tolerance:
            groups[-1].append(item)
        else:
            groups.append([item])

    found = []
    for group in groups:
        if len(group) < min_pairs:
            continue
        # independent pairs only: a shared index means a shared object
        seen: set[int] = set()
        examples = []
        used = 0
        for delta, i, j in group:
            if i in seen or j in seen:
                continue
            seen.update((i, j))
            examples.append((array[i], array[j]))
            used += 1
        if used >= min_pairs:
            found.append(
                Repetition(
                    delta=float(np.mean([d for d, _, _ in group])),
                    factor=float(PHI ** float(np.mean([d for d, _, _ in group]))),
                    pairs=used,
                    examples=tuple(examples),
                )
            )
    return found


def rung_uniformity(values) -> dict:
    """Are fractional rungs uniform? A low p-value means quantization.

    Kolmogorov-Smirnov against the uniform distribution on [0, 0.5], with the
    standard asymptotic p-value (no scipy needed). Uniform is the null for any
    smooth population, so a small p-value says the values were produced on a
    grid - which is a property of the search, not of the sky.
    """
    array = [float(v) for v in values]
    if len(array) < 4:
        return {"n": len(array), "statistic": float("nan"), "p_value": float("nan")}
    fracs = np.sort(np.array([rung_drift(v) for v in array]) / 0.5)  # to [0, 1]
    n = fracs.size
    index = np.arange(1, n + 1)
    d_stat = float(np.max(np.maximum(index / n - fracs, fracs - (index - 1) / n)))
    lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d_stat
    p_value = 2.0 * sum(
        (-1.0) ** (k - 1) * math.exp(-2.0 * k * k * lam * lam) for k in range(1, 101)
    )
    return {
        "n": int(n),
        "statistic": d_stat,
        "p_value": float(min(max(p_value, 0.0), 1.0)),
    }


def phi_coordinates(frequency_hz: float, drift_hz_per_s: float | None = None) -> dict:
    """The full coordinate set for one measurement, in one call.

    A rung is `ln(x)/ln(phi)`, so an integer difference is a factor of phi and a
    fractional part is where the value sits inside that step. The residual is the
    fractional part signed to [-0.5, 0.5], which is the number a ladder claim is
    actually about: 0.49 is not "nearly on the rung", it is half a rung away.
    """
    out: dict = {}
    if frequency_hz and frequency_hz > 0:
        value = math.log(float(frequency_hz)) / LN_PHI
        nearest = int(round(value))
        out.update({
            "phi_freq_rung": round(value, 6),
            "phi_freq_nearest": nearest,
            "phi_freq_residual": round(value - nearest, 6),
        })
    else:
        out.update({"phi_freq_rung": None, "phi_freq_nearest": None,
                    "phi_freq_residual": None})
    if drift_hz_per_s:
        value = math.log(abs(float(drift_hz_per_s))) / LN_PHI
        nearest = int(round(value))
        out.update({
            "phi_drift_rung": round(value, 6),
            "phi_drift_nearest": nearest,
            "phi_drift_residual": round(value - nearest, 6),
        })
    else:
        # a zero drift has no logarithm; it is a statement about the coordinate,
        # not a measurement of zero, and it is recorded as absent
        out.update({"phi_drift_rung": None, "phi_drift_nearest": None,
                    "phi_drift_residual": None})
    return out
