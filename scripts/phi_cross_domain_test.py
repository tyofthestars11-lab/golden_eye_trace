#!/usr/bin/env python3
"""
phi_cross_domain_test.py

Preregistered cross-domain test of a fixed phi-rung interval.

Prediction
----------
For an ordered sequence of positive measurements x_i in each independent domain,
adjacent multiplicative ratios should show an excess near:

    delta_phi = 1 / phi**2 = 0.38196601125010515

where:

    delta_r = log(x[i+1] / x[i]) / log(base)

The primary test uses base=phi and tolerance=0.01. Control bases are 2, e, 10,
and sqrt(2), with the same tolerance in their respective log-base rung
coordinates. The script never compares dimensional quantities across domains; it
compares within-domain ratios only.

Input CSV format
----------------
Each CSV must contain one row per ordered measurement and these columns:

    order,value

Optional columns are retained in output:

    id,label,uncertainty

Rows are sorted by numeric `order` if possible, otherwise lexicographically.
The order is part of the preregistration: do not reorder rows after seeing the
results. Missing, non-finite, or non-positive values are rejected.

Example
-------
python3 phi_cross_domain_test.py \\
  --dataset genomics=genomics.csv \\
  --dataset astronomy=astronomy.csv \\
  --dataset physics=physics.csv \\
  --out results/phi_cross_domain \\
  --permutations 100000 \\
  --seed 20261114

The script uses only Python's standard library.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

PHI = (1.0 + math.sqrt(5.0)) / 2.0
PRIMARY_BASE = PHI
PRIMARY_TARGET = 1.0 / PHI**2
PRIMARY_TOLERANCE = 0.01
CONTROL_BASES = {
    "base_2": 2.0,
    "base_e": math.e,
    "base_10": 10.0,
    "base_sqrt2": math.sqrt(2.0),
}


@dataclass(frozen=True)
class Observation:
    order: str
    value: float
    row_number: int
    record_id: str = ""
    label: str = ""
    uncertainty: float | None = None


@dataclass(frozen=True)
class Pair:
    left_index: int
    right_index: int
    left_value: float
    right_value: float
    ratio: float


class InputError(ValueError):
    pass


def finite_positive(value: float) -> bool:
    return math.isfinite(value) and value > 0.0


def parse_optional_float(raw: str, field: str, row_number: int) -> float | None:
    raw = raw.strip()
    if raw == "":
        return None
    try:
        value = float(raw)
    except ValueError as exc:
        raise InputError(
            f"row {row_number}: {field!r} is not numeric: {raw!r}") from exc
    if not math.isfinite(value) or value < 0:
        raise InputError(f"row {row_number}: {field!r} must be finite and >= 0")
    return value


def load_dataset(path: Path) -> list[Observation]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if (reader.fieldnames is None or "order" not in reader.fieldnames
                or "value" not in reader.fieldnames):
            raise InputError(
                f"{path}: CSV must contain columns 'order' and 'value'")
        observations: list[Observation] = []
        for row_number, row in enumerate(reader, start=2):
            raw_order = (row.get("order") or "").strip()
            raw_value = (row.get("value") or "").strip()
            if not raw_order:
                raise InputError(f"{path}: row {row_number}: empty order")
            try:
                value = float(raw_value)
            except ValueError as exc:
                raise InputError(
                    f"{path}: row {row_number}: value is not numeric: "
                    f"{raw_value!r}") from exc
            if not finite_positive(value):
                raise InputError(
                    f"{path}: row {row_number}: value must be finite and > 0")
            observations.append(
                Observation(
                    order=raw_order,
                    value=value,
                    row_number=row_number,
                    record_id=(row.get("id") or "").strip(),
                    label=(row.get("label") or "").strip(),
                    uncertainty=parse_optional_float(
                        row.get("uncertainty", ""), "uncertainty", row_number),
                )
            )
    if len(observations) < 3:
        raise InputError(f"{path}: at least 3 observations are required")
    # Numeric order is preferred; otherwise preserve a deterministic lexical order.
    try:
        sorted_observations = sorted(observations, key=lambda x: float(x.order))
    except ValueError:
        sorted_observations = sorted(observations, key=lambda x: x.order)
    return sorted_observations


def adjacent_pairs(values: list[float]) -> list[Pair]:
    return [
        Pair(i, i + 1, values[i], values[i + 1], values[i + 1] / values[i])
        for i in range(len(values) - 1)
    ]


def rung_delta(ratio: float, base: float) -> float:
    return math.log(ratio) / math.log(base)


def target_for_base(base: float) -> float:
    # The fixed prediction is expressed as 1/phi^2 in rung units for every
    # logarithmic base, so controls test whether the effect is base-specific.
    return PRIMARY_TARGET


def observed_hits(values: list[float], base: float,
                  tolerance: float) -> tuple[int, list[float]]:
    deltas = [rung_delta(pair.ratio, base) for pair in adjacent_pairs(values)]
    target = target_for_base(base)
    hits = sum(abs(delta - target) <= tolerance for delta in deltas)
    return hits, deltas


def permutation_null(
    values: list[float],
    base: float,
    tolerance: float,
    permutations: int,
    rng: random.Random,
) -> tuple[int, int, float, list[int]]:
    """Shuffle values across the fixed order while preserving the value multiset."""
    observed, _ = observed_hits(values, base, tolerance)
    shuffled = list(values)
    null_counts: list[int] = []
    at_least = 0
    for _ in range(permutations):
        rng.shuffle(shuffled)
        count, _ = observed_hits(shuffled, base, tolerance)
        null_counts.append(count)
        if count >= observed:
            at_least += 1
    # Add-one correction avoids reporting p=0 for finite Monte Carlo samples.
    p = (at_least + 1) / (permutations + 1)
    return observed, at_least, p, null_counts


def odds_ratio(observed_hits_count: int, n_pairs: int,
               null_counts: list[int]) -> float | None:
    if not null_counts or n_pairs == 0:
        return None
    null_mean = statistics.mean(null_counts)
    # Haldane-style continuity correction for zero cells.
    a = observed_hits_count + 0.5
    b = (n_pairs - observed_hits_count) + 0.5
    c = null_mean + 0.5
    d = (n_pairs - null_mean) + 0.5
    return (a * d) / (b * c)


def fisher_combine(p_values: Iterable[float]) -> tuple[float, float]:
    ps = [min(max(float(p), 1e-300), 1.0) for p in p_values]
    statistic = -2.0 * sum(math.log(p) for p in ps)
    # For an even number of degrees of freedom 2k, chi-square survival is
    # exp(-x/2) * sum_{j=0}^{k-1} (x/2)^j / j!.
    k = len(ps)
    if k == 0:
        return 0.0, 1.0
    x = statistic / 2.0
    survival = math.exp(-x) * sum(x**j / math.factorial(j) for j in range(k))
    return statistic, min(max(survival, 0.0), 1.0)


def bh_adjust(p_values: dict[str, float]) -> dict[str, float]:
    ordered = sorted(p_values.items(), key=lambda item: item[1])
    m = len(ordered)
    adjusted: dict[str, float] = {}
    running = 1.0
    for rank_from_end, (key, p) in enumerate(reversed(ordered), start=1):
        rank = m - rank_from_end + 1
        value = min(running, p * m / rank)
        running = value
        adjusted[key] = min(value, 1.0)
    return adjusted


def dataset_result(
    name: str,
    path: Path,
    observations: list[Observation],
    permutations: int,
    seed: int,
) -> dict[str, Any]:
    values = [obs.value for obs in observations]
    n_pairs = len(values) - 1
    bases = {"phi": PRIMARY_BASE, **CONTROL_BASES}
    results: dict[str, Any] = {
        "name": name,
        "path": str(path),
        "n_observations": len(values),
        "n_adjacent_pairs": n_pairs,
        "order": [obs.order for obs in observations],
        "values": values,
        "tests": {},
    }
    for offset, (base_name, base) in enumerate(bases.items()):
        rng = random.Random(seed + 1009 * offset)
        hits, above, p, null_counts = permutation_null(
            values, base, PRIMARY_TOLERANCE, permutations, rng)
        _, deltas = observed_hits(values, base, PRIMARY_TOLERANCE)
        results["tests"][base_name] = {
            "base": base,
            "target_rung_interval": PRIMARY_TARGET,
            "tolerance": PRIMARY_TOLERANCE,
            "hit_count": hits,
            "pair_count": n_pairs,
            "hit_fraction": hits / n_pairs,
            "observed_deltas": deltas,
            "null_mean_hit_count": statistics.mean(null_counts),
            "null_sd_hit_count": statistics.pstdev(null_counts),
            "null_at_least_as_extreme": above,
            "permutation_count": permutations,
            "p_value": p,
            "odds_ratio_vs_null_mean": odds_ratio(hits, n_pairs, null_counts),
        }
    return results


def flatten_results(dataset_results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result in dataset_results:
        for base_name, test in result["tests"].items():
            rows.append(
                {
                    "dataset": result["name"],
                    "base": base_name,
                    "hit_count": test["hit_count"],
                    "pair_count": test["pair_count"],
                    "hit_fraction": test["hit_fraction"],
                    "null_mean_hit_count": test["null_mean_hit_count"],
                    "p_value": test["p_value"],
                    "odds_ratio_vs_null_mean": test["odds_ratio_vs_null_mean"],
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_report(
    output_dir: Path,
    dataset_results: list[dict[str, Any]],
    combined: dict[str, Any],
    config: dict[str, Any],
) -> None:
    lines = [
        "# Preregistered Cross-Domain φ-Structure Test",
        "",
        "## Frozen prediction",
        "",
        f"The primary target is Δr = 1/φ² = {PRIMARY_TARGET:.15f}, with tolerance "
        f"±{PRIMARY_TOLERANCE:.6f}. Each Δr is computed from an adjacent "
        "within-domain ratio: log(x[i+1]/x[i]) / log(base).",
        "",
        "Dimensional quantities are never compared directly across domains. Each "
        "dataset is tested only through ratios of its own positive measurements "
        "in its declared order.",
        "",
        "## Dataset results",
        "",
        "| dataset | φ hits | pairs | φ p-value | φ odds ratio |",
        "|---|---:|---:|---:|---:|",
    ]
    for result in dataset_results:
        test = result["tests"]["phi"]
        lines.append(
            f"| {result['name']} | {test['hit_count']} | "
            f"{test['pair_count']} | {test['p_value']:.6g} | "
            f"{test['odds_ratio_vs_null_mean']:.6g} |")
    lines += [
        "",
        f"Fisher combined statistic: {combined['phi']['fisher_statistic']:.6g}",
        f"Fisher combined p-value: {combined['phi']['fisher_p_value']:.6g}",
        "",
        "## Control bases",
        "",
        "The control bases use the same fixed φ target in their own logarithmic "
        "coordinates and the same tolerance. This is a preregistered sensitivity "
        "comparison, not a refit of the target.",
        "",
        "| base | combined Fisher p-value |",
        "|---|---:|",
    ]
    for base_name, summary in combined.items():
        lines.append(f"| {base_name} | {summary['fisher_p_value']:.6g} |")
    lines += [
        "",
        "## Decision rule",
        "",
        "The prediction is supported only if every domain has an odds ratio ≥ 2, "
        "every domain has p < 0.01, the combined preregistered p-value is "
        "< 1e-4, and the φ result is stronger than the control-base results "
        "after the declared multiplicity adjustment. This script reports those "
        "quantities; it does not convert them into a claim about life or "
        "intelligence.",
        "",
        "## Reproducibility",
        "",
        f"Configuration: `{json.dumps(config, sort_keys=True)}`",
        "",
        "A significant result would support only the narrower claim that the "
        "tested datasets share a reproducible multiplicative spacing near the "
        "preregistered φ interval. It would not by itself establish a universal "
        "law, consciousness, or a technosignature.",
    ]
    (output_dir / "REPORT.md").write_text("\n".join(lines) + "\n",
                                          encoding="utf-8")


def parse_dataset_spec(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise InputError(f"--dataset must be NAME=PATH, got {spec!r}")
    name, raw_path = spec.split("=", 1)
    name = name.strip()
    if not name or not raw_path.strip():
        raise InputError(f"--dataset must be NAME=PATH, got {spec!r}")
    return name, Path(raw_path).expanduser().resolve()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        required=True,
        metavar="NAME=CSV",
        help="Repeat exactly three times, e.g. genomics=data/genomics.csv",
    )
    parser.add_argument("--out", default="results/phi_cross_domain",
                        help="Output directory")
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20261114)
    args = parser.parse_args(argv)

    if len(args.dataset) != 3:
        parser.error("provide exactly three --dataset arguments")
    if args.permutations < 1000:
        parser.error("--permutations must be at least 1000")

    try:
        specs = [parse_dataset_spec(spec) for spec in args.dataset]
        names = [name for name, _ in specs]
        if len(set(names)) != 3:
            raise InputError("dataset names must be unique")
        loaded = [(name, path, load_dataset(path)) for name, path in specs]
        results = [
            dataset_result(name, path, observations, args.permutations,
                           args.seed + i * 100_003)
            for i, (name, path, observations) in enumerate(loaded)
        ]
    except (OSError, InputError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    out_dir = Path(args.out).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    p_values_by_base: dict[str, list[float]] = {
        "phi": [], **{key: [] for key in CONTROL_BASES}}
    for result in results:
        for base_name, test in result["tests"].items():
            p_values_by_base[base_name].append(test["p_value"])

    combined: dict[str, Any] = {}
    for base_name, p_values in p_values_by_base.items():
        statistic, p_value = fisher_combine(p_values)
        combined[base_name] = {
            "per_dataset_p_values": p_values,
            "fisher_statistic": statistic,
            "fisher_p_value": p_value,
        }

    all_p = {
        f"{result['name']}::{base_name}": test["p_value"]
        for result in results
        for base_name, test in result["tests"].items()
    }
    adjusted = bh_adjust(all_p)
    for result in results:
        for base_name, test in result["tests"].items():
            test["bh_adjusted_p_value_across_all_domain_base_tests"] = adjusted[
                f"{result['name']}::{base_name}"
            ]

    config = {
        "phi": PHI,
        "primary_target_rung": PRIMARY_TARGET,
        "tolerance": PRIMARY_TOLERANCE,
        "control_bases": CONTROL_BASES,
        "permutations": args.permutations,
        "seed": args.seed,
        "input_order_rule": "numeric order if possible, else lexical order",
        "null": "shuffle observed values across fixed declared order",
    }
    payload = {
        "config": config,
        "datasets": results,
        "combined": combined,
    }
    (out_dir / "results.json").write_text(json.dumps(payload, indent=2) + "\n",
                                          encoding="utf-8")
    write_csv(out_dir / "summary.csv", flatten_results(results))
    make_report(out_dir, results, combined, config)

    primary = combined["phi"]
    print(f"Wrote {out_dir}")
    print(f"Primary target Δr = {PRIMARY_TARGET:.15f} ± {PRIMARY_TOLERANCE}")
    print(f"Combined φ Fisher p = {primary['fisher_p_value']:.6g}")
    for result in results:
        test = result["tests"]["phi"]
        print(
            f"{result['name']}: {test['hit_count']}/{test['pair_count']} hits, "
            f"p={test['p_value']:.6g}, OR={test['odds_ratio_vs_null_mean']:.6g}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
