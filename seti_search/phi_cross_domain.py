"""The cross-domain protocol, with the parts that can fail left in.

The claim under test is narrow and falsifiable:

    adjacent ratios in an ordered series sit at Δr = 1/φ² = 0.381966 rungs,
    measured in base φ, within ±0.01.

Three things a test of it must contain, each of which can only make the claim
harder to support:

1. **A null.** The permutation null shuffles the observed values across the
   fixed order, which preserves the value multiset and destroys only the
   sequence. A hit rate that the shuffle reproduces is not structure.
2. **Held-out prediction.** A fixed interval predicts the *next* value from the
   current one. If that prediction is no better on data the interval was not
   derived from - and here the interval is fixed a priori, so every pair is
   held-out by construction - then the interval is not describing the series.
   The baseline that matters is not "no structure": it is a *fitted constant
   ratio*, which is what a series with any regular spacing will obey. φ must
   beat the fitted constant, not merely beat zero.
3. **Power.** With n pairs the test can only see a law that applies to nearly
   every step. Injecting the law at a known fraction and measuring the
   detection rate turns "0 hits" into a statement with a scope.

Every function here reports a number beside the null that gives it meaning, and
nothing in this module can promote a disposition: φ may only detract. A positive
result here would support the claim that the tested series share a spacing, not
that the spacing is intelligent, alive, or universal.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field

__all__ = [
    "BASES",
    "PHI",
    "Series",
    "TARGET",
    "TOLERANCE",
    "bh",
    "deltas",
    "fit_ratio",
    "hit_count",
    "holdout_prediction",
    "injection_power",
    "permutation_p",
    "rung_error",
]

PHI = (1.0 + math.sqrt(5.0)) / 2.0
TARGET = 1.0 / PHI ** 2
TOLERANCE = 0.01

#: The controls. Each is given the *same* fixed rung interval in its own
#: coordinates, so the comparison asks whether the interval is special in base φ
#: rather than whether it is special at all.
BASES = {
    "phi": PHI,
    "base_2": 2.0,
    "base_e": math.e,
    "base_10": 10.0,
    "base_sqrt2": math.sqrt(2.0),
}


@dataclass(frozen=True)
class Series:
    """One ordered series, with its units and its provenance stated.

    ``unit`` and ``normalisation`` are required because a ratio of two values in
    different units is not a ratio: the protocol only ever divides a value by
    another value from the same series, and saying which unit that is costs
    nothing and prevents a silent mix.
    """

    name: str
    domain: str
    unit: str
    values: tuple[float, ...]
    provenance: str
    normalisation: str = "none needed: only within-series ratios are taken"
    synthetic: bool = False
    #: Name of a series this one is a subset of, if any. A nested series cannot
    #: be counted as an independent domain, and the report says so rather than
    #: combining it into a Fisher statistic.
    nested_with: str | None = None
    labels: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if len(self.values) < 3:
            raise ValueError(f"{self.name}: at least three values are needed")
        for value in self.values:
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{self.name}: values must be finite and positive")
        if not self.unit.strip():
            raise ValueError(f"{self.name}: a unit must be declared")


def deltas(values, base: float) -> list[float]:
    """Adjacent rung intervals in the given base."""
    if base <= 1.0:
        raise ValueError("a logarithmic base must exceed 1")
    return [math.log(b / a) / math.log(base)
            for a, b in zip(values, values[1:])]


def hit_count(values, base: float, tolerance: float = TOLERANCE) -> int:
    return sum(abs(d - TARGET) <= tolerance for d in deltas(values, base))


def rung_error(predicted_ratio: float, actual_ratio: float,
               base: float) -> float:
    """Distance between two ratios, in rungs of the base under test."""
    if predicted_ratio <= 0 or actual_ratio <= 0:
        raise ValueError("ratios must be positive")
    return abs(math.log(actual_ratio / predicted_ratio)) / math.log(base)


def permutation_p(values, base: float, tolerance: float = TOLERANCE,
                  permutations: int = 10_000,
                  seed: int = 20261114) -> dict:
    """P(>= observed hits) under a shuffle of the values across the fixed order."""
    if permutations < 100:
        raise ValueError("at least 100 permutations are needed to say anything")
    observed = hit_count(values, base, tolerance)
    rng = random.Random(seed)
    pool = list(values)
    at_least = 0
    counts = []
    for _ in range(permutations):
        rng.shuffle(pool)
        count = hit_count(pool, base, tolerance)
        counts.append(count)
        if count >= observed:
            at_least += 1
    p = (at_least + 1) / (permutations + 1)  # add-one: never report p = 0
    return {"observed": observed, "p_value": p,
            "null_mean": sum(counts) / len(counts),
            "null_max": max(counts), "permutations": permutations}


def fit_ratio(values, held_out_fraction: float = 0.5) -> dict:
    """The strong baseline: a constant ratio fitted on the training half.

    The median is used rather than the mean so that a single cross-day jump does
    not decide the baseline; the point of this baseline is to be the strongest
    honest predictor a regular series admits.
    """
    ratios = [b / a for a, b in zip(values, values[1:])]
    k = max(1, int(round(len(ratios) * (1.0 - held_out_fraction))))
    train, test = ratios[:k], ratios[k:]
    if not test:
        test = ratios[-1:]
        train = ratios[:-1]
    # The median, not the geometric mean: one large cross-day jump otherwise
    # sets the baseline below 1 and makes the comparison about that jump.
    ordered = sorted(train)
    middle = len(ordered) // 2
    ratio = (ordered[middle] if len(ordered) % 2
             else 0.5 * (ordered[middle - 1] + ordered[middle]))
    return {"train_ratio": ratio, "n_train": len(train),
            "n_test": len(test), "test_ratios": test}


def holdout_prediction(values, base: float,
                       held_out_fraction: float = 0.5) -> dict:
    """Score a fixed rung interval against the held-out adjacent ratios.

    The interval is fixed a priori, so "training" only sets the offset: each
    prediction is ``previous value x base**TARGET``. What is measured is whether
    that interval lands on the next value, which is exactly the claim.
    """
    ratios = [b / a for a, b in zip(values, values[1:])]
    k = max(1, int(round(len(ratios) * (1.0 - held_out_fraction))))
    test = ratios[k:] or ratios[-1:]
    predicted = base ** TARGET
    errors = [rung_error(predicted, actual, base) for actual in test]
    hits = sum(1 for actual in test if abs(math.log(actual / predicted))
               / math.log(base) <= TOLERANCE)
    return {"predicted_ratio": predicted, "n_test": len(test),
            "hit_rate": hits / len(test),
            "median_rung_error": sorted(errors)[len(errors) // 2],
            "mean_rung_error": sum(errors) / len(errors)}


def injection_power(values, base: float, fraction: float, *,
                    trials: int = 200, permutations: int = 100,
                    tolerance: float = TOLERANCE,
                    threshold: float = 0.05,
                    seed: int = 20261114) -> float:
    """Detection rate when the law holds on ``fraction`` of the steps.

    The remaining steps keep the series' own geometric drift, so the injection
    adds the law without destroying the series it is added to. A detection is
    ``p <= threshold`` from the permutation null of the injected series.
    """
    if not (0.0 <= fraction <= 1.0):
        raise ValueError("fraction is a proportion")
    ratios = [b / a for a, b in zip(values, values[1:])]
    drift = math.exp(sum(math.log(r) for r in ratios) / len(ratios))
    law = base ** TARGET
    rng = random.Random(seed)
    detected = 0
    n = len(ratios)
    for trial in range(trials):
        n_law = int(round(fraction * n))
        which = set(rng.sample(range(n), n_law))
        cs = [1.0]
        for j in range(n):
            cs.append(cs[-1] * (law if j in which else drift))
        series = [values[0] * c for c in cs]
        result = permutation_p(series, base, tolerance, permutations,
                               seed=seed + trial)
        if result["p_value"] <= threshold:
            detected += 1
    return detected / trials


def bh(p_values: dict[str, float]) -> dict[str, float]:
    """Benjamini-Hochberg adjustment over every test actually performed."""
    ordered = sorted(p_values.items(), key=lambda item: item[1])
    m = len(ordered)
    out: dict[str, float] = {}
    running = 1.0
    for rank_from_end, (key, p) in enumerate(reversed(ordered), start=1):
        rank = m - rank_from_end + 1
        running = min(running, p * m / rank)
        out[key] = min(running, 1.0)
    return out
