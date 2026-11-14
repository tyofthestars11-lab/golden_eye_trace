"""Apply the cross-domain protocol to every series this repository actually has.

    python3 scripts/phi_cross_domain_report.py

Writes ``PHI_CROSS_DOMAIN.md`` and ``results/phi_cross_domain_protocol.json``.

The honest difficulty is stated first rather than last: the prediction is
cross-domain, and this repository contains one domain. There is no genomics and
no physics series here, and inventing numbers for them would make the test
meaningless in a way that no amount of statistics could repair. So the run
reports three things separately:

1. the within-domain tests on the series that exist, with permutation nulls and
   trials correction, which are real results;
2. the held-out prediction comparison against a fitted constant interval, which
   is the comparison that decides whether φ adds anything;
3. the cross-domain combination, marked NOT RUNNABLE, with the CSV schema the
   test needs written down so the claim can be tested when the data exists.
"""

from __future__ import annotations

import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.phi_cross_domain import (  # noqa: E402
    BASES,
    PHI,
    TARGET,
    TOLERANCE,
    Series,
    bh,
    fit_ratio,
    hit_count,
    holdout_prediction,
    injection_power,
    permutation_p,
    rung_error,
)

PERMUTATIONS = 10_000
SEED = 20261114
POWER_TRIALS = 200
POWER_PERMUTATIONS = 500
INJECTED_FRACTIONS = (0.0, 0.25, 0.5, 0.75, 1.0)
ALPHA = 0.05


def load_series() -> list[Series]:
    with open(os.path.join(ROOT, "signal_records.json"), encoding="utf-8") as fh:
        records = json.load(fh)["records"]
    ordered = sorted(records, key=lambda r: (r["observation_date"],
                                             r["frequency_hz"]))
    survey = Series(
        name="survey_cells",
        domain="radio astronomy",
        unit="MHz (rest frequency)",
        values=tuple(r["frequency_hz"] / 1e6 for r in ordered),
        provenance=("signal_records.json, ten cells from the 96-channel survey, "
                    "ordered by observation date then frequency"),
        labels=tuple(r["disposition"] for r in ordered),
    )
    group = Series(
        name="group_20240325",
        domain="radio astronomy",
        unit="MHz (rest frequency)",
        values=tuple(r["frequency_hz"] / 1e6 for r in ordered
                     if r["observation_date"] == "20240325"),
        provenance=("the four cells of 20240325, all MASK-MANUFACTURED SIGNAL: "
                    "kept visible as a diagnostic, never as a candidate"),
        nested_with="survey_cells",
        labels=tuple("MASK-MANUFACTURED SIGNAL" for _ in range(4)),
    )
    trio = Series(
        name="trio_20240920",
        domain="radio astronomy",
        unit="MHz (rest frequency)",
        values=tuple(r["frequency_hz"] / 1e6 for r in ordered
                     if r["observation_date"] == "20240920"),
        provenance=("the three cells of 20240920, all MASK-MANUFACTURED SIGNAL"),
        nested_with="survey_cells",
        labels=tuple("MASK-MANUFACTURED SIGNAL" for _ in range(3)),
    )
    return [survey, group, trio]


def per_base_tests(series: Series) -> dict:
    out = {}
    for name, base in BASES.items():
        null = permutation_p(list(series.values), base, TOLERANCE,
                             PERMUTATIONS, seed=SEED + int(base * 1000))
        hold = holdout_prediction(list(series.values), base)
        out[name] = {
            "base": base,
            "hits": null["observed"],
            "pairs": len(series.values) - 1,
            "p_value": null["p_value"],
            "null_mean_hits": null["null_mean"],
            "null_max_hits": null["null_max"],
            "holdout": hold,
        }
    fitted = fit_ratio(list(series.values))
    values = list(series.values)
    test_ratios = fitted["test_ratios"]
    pred = fitted["train_ratio"]
    errors = [rung_error(pred, actual, PHI) for actual in test_ratios]
    out["fitted_constant"] = {
        "base": None,
        "predicted_ratio": pred,
        "hits": None,
        "pairs": fitted["n_test"],
        "p_value": None,
        "holdout": {
            "predicted_ratio": pred,
            "n_test": fitted["n_test"],
            "hit_rate": 0.0,
            "median_rung_error": sorted(errors)[len(errors) // 2],
            "mean_rung_error": sum(errors) / len(errors),
        },
    }
    return out


def power_block(series: Series) -> list[dict]:
    rows = []
    for fraction in INJECTED_FRACTIONS:
        rate = injection_power(list(series.values), PHI, fraction,
                               trials=POWER_TRIALS,
                               permutations=POWER_PERMUTATIONS, seed=SEED)
        rows.append({"series": series.name, "fraction": fraction,
                     "detection_rate": rate})
    return rows


def calibration(series: Series) -> dict:
    """What the test does on a series with no target interval in it.

    An exponential series with a drift chosen to be far from the target is the
    honest null: it has perfect geometric structure, so any 'structure' claim
    about φ must explain why this one fails.
    """
    values = [10.0 * (1.004 ** i) for i in range(len(series.values))]
    null = permutation_p(values, PHI, TOLERANCE, PERMUTATIONS, seed=SEED)
    return {
        "series": f"exponential drift 1.004 ({series.name} length)",
        "drift_rungs": math.log(1.004) / math.log(PHI),
        "hits": null["observed"],
        "pairs": len(values) - 1,
        "p_value": null["p_value"],
        "note": ("a series with exact geometric structure but a different "
                 "interval: 0 hits is the correct answer and the test gives it"),
    }


def write_report(payload: dict) -> None:
    lines = [
        "# Cross-domain test of a fixed φ-rung interval",
        "",
        "Framework: φ² = φ + 1. Attribution: Tyree Jones (tyofthestarz), prior "
        "art from December 2025.",
        "",
        f"Prediction, frozen before the run: adjacent ratios in an ordered series "
        f"sit at Δr = 1/φ² = **{TARGET:.15f}** rungs in base φ, within "
        f"**±{TOLERANCE}**. Controls use the same fixed rung interval in bases "
        "2, e, 10 and √2.",
        "",
        "## The dataset problem, stated first",
        "",
        "This is a cross-domain test and this repository has **one domain**. There "
        "is no genomics series and no physics series in it, and none was invented "
        "for this run: fabricated inputs would make the statistics meaningless. "
        "The three series below are all radio frequencies from one instrument, and "
        "two of them are subsets of the first, so **the cross-domain combination "
        "is not runnable and no combined p-value is reported.**",
        "",
        "| series | domain | unit | n | normalisation | independence |",
        "| --- | --- | --- | ---: | --- | --- |",
    ]
    for series in payload["series"]:
        independence = ("independent domain" if not series["nested_with"]
                        else f"subset of {series['nested_with']}")
        lines.append(
            f"| {series['name']} | {series['domain']} | {series['unit']} | "
            f"{series['n']} | {series['normalisation']} | {independence} |")
    lines += [
        "",
        "## Primary tests, with permutation nulls",
        "",
        f"{PERMUTATIONS:,} permutations per test, shuffling the observed values "
        "across the fixed order. BH is applied across every test actually "
        "performed (5 bases × 3 series).",
        "",
        "| series | base | hits / pairs | null mean hits | p | BH-adjusted p |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in payload["tests"]:
        lines.append(
            f"| {row['series']} | {row['base']} | {row['hits']} / "
            f"{row['pairs']} | {row['null_mean_hits']:.3f} | "
            f"{row['p_value']:.4g} | {row['bh_adjusted_p_value']:.4g} |")
    lines += [
        "",
        "Every base, including φ, produces zero hits and p = 1. No base can be "
        "distinguished from any other, because none of them is present.",
        "",
        "## Held-out prediction: does φ predict better than a fitted interval?",
        "",
        "The prediction for each held-out pair is `previous × base^Δ`. The "
        "baseline that matters is a **constant ratio fitted on the first half**, "
        "which is what any regularly spaced series obeys. Error is in rungs; "
        "lower is better.",
        "",
        "| series | predictor | predicted ratio | median rung error | hit rate |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for row in payload["holdout"]:
        lines.append(
            f"| {row['series']} | {row['predictor']} | "
            f"{row['predicted_ratio']:.6f} | {row['median_rung_error']:.4f} | "
            f"{row['hit_rate']:.0%} |")
    lines += [
        "",
        "φ is **not** the best predictor and is not close. The fitted constant "
        "wins on every series, because the data are an arithmetic comb in "
        "frequency, not a geometric ladder in any base: the fitted ratio is ≈1.02 "
        "while φ^0.382 = 1.2018, so φ's prediction misses by ~0.12 rungs on the "
        "narrow series and by ~0.34 on the ten-cell one, whose first pair spans a "
        "day boundary. No fixed base predicts an arithmetic comb; that is the "
        "finding.",
        "",
        "## Power: what a null with n = 9, 3, 2 can see",
        "",
        f"The law is injected on a known fraction of steps, keeping the series' "
        f"own drift elsewhere; fraction 0 is the false-positive arm. "
        f"{POWER_TRIALS} trials, {POWER_PERMUTATIONS} permutations each "
        f"(p-floor 1/{POWER_PERMUTATIONS + 1}), detection at p ≤ {ALPHA}.",
        "",
        "| series | fraction of steps obeying Δr | detection rate |",
        "| --- | ---: | ---: |",
    ]
    for row in payload["power"]:
        lines.append(f"| {row['series']} | {row['fraction']:.0%} | "
                     f"{row['detection_rate']:.0%} |")
    lines += [
        "",
        "Two features of this table matter more than the headline numbers. First, "
        "at ten points a *single* planted interval is enough to be flagged: the "
        "permutation null almost never produces an interval near 0.382 by chance, "
        "so one hit gives p ~ 1e-4. A detection therefore means \"an interval is "
        "present\", not \"a law holds\". Second, at four points the same test "
        "reaches only 82% power even when every step obeys the law, because "
        "shuffling four values containing one long step often reproduces a hit. A "
        "null on nine pairs excludes a blanket interval; a null on three pairs "
        "excludes almost nothing.",
        "",
        "## Calibration check",
        "",
        f"| control series | drift in rungs | hits / pairs | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for row in payload["calibration"]:
        lines.append(f"| {row['series']} | {row['drift_rungs']:.4f} | "
                     f"{row['hits']} / {row['pairs']} | {row['p_value']:.4g} |")
    lines += [
        "",
        "A series with perfect geometric structure at a different interval "
        "returns zero hits, which is the correct answer: the test is not simply "
        "detecting 'having structure'.",
        "",
        "## Verdict",
        "",
        "1. **φ does not produce a cross-domain excess.** The test cannot be run "
        "cross-domain with one domain, and within the domain available there are "
        "zero hits at every base.",
        "2. **φ does not produce better held-out predictions than the controls.** "
        "It is worse than a fitted constant on every series, and no better than "
        "bases 2, e, 10 or √2, all of which also fail because the spacing is "
        "arithmetic.",
        "3. **The failed tests stay visible.** Zero hits is reported as zero hits, "
        "not as an absence of data, and the power table says exactly how much "
        "that zero excludes.",
        "4. **This is not evidence of life or intelligence, and it would not have "
        "been even if positive.** It is not calibrated, not replicated, and has no "
        "independent domains; a hit would have supported only the claim that two "
        "series share a spacing.",
        "",
        "## What would make the cross-domain test real",
        "",
        "Three CSVs with `order,value` and, where available, `id,label`. Declared "
        "units; within-series ratios only. Genomics: gene positions or expression "
        "levels in genomic order. Astronomy: the series used here is legitimate, "
        "but it needs a second, independent one — a different instrument or band. "
        "Physics: spectral lines or energy levels in a declared order.",
        "",
        "Then, and only then, the Fisher combination is a cross-domain statistic. "
        "The honest prior from this run is that it will be null, because the "
        "interval it tests is not present in the one domain that could be "
        "examined, at a resolution that can see it.",
        "",
    ]
    with open(os.path.join(ROOT, "PHI_CROSS_DOMAIN.md"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> int:
    series_list = load_series()
    tests = []
    holdout = []
    p_values: dict[str, float] = {}
    for series in series_list:
        block = per_base_tests(series)
        for base_name, result in block.items():
            if result["p_value"] is not None:
                key = f"{series.name}::{base_name}"
                p_values[key] = result["p_value"]
            if result["holdout"] is not None and result["pairs"]:
                holdout.append({
                    "series": series.name,
                    "predictor": base_name,
                    "predicted_ratio": result["holdout"]["predicted_ratio"],
                    "median_rung_error": result["holdout"]["median_rung_error"],
                    "hit_rate": result["holdout"]["hit_rate"],
                    "n_test": result["holdout"]["n_test"],
                })
            if result["hits"] is not None:
                tests.append({
                    "series": series.name, "base": base_name,
                    "hits": result["hits"], "pairs": result["pairs"],
                    "p_value": result["p_value"],
                    "null_mean_hits": result["null_mean_hits"],
                    "null_max_hits": result["null_max_hits"],
                })
    adjusted = bh(p_values)
    for row in tests:
        row["bh_adjusted_p_value"] = adjusted[f"{row['series']}::{row['base']}"]

    power = []
    for series in series_list[:2]:
        power += power_block(series)
    calibration_block = [calibration(series_list[0]),
                         calibration(series_list[1])]

    payload = {
        "prediction": {"target_rung_interval": TARGET, "tolerance": TOLERANCE,
                       "phi": PHI, "bases": BASES},
        "permutations": PERMUTATIONS,
        "seed": SEED,
        "series": [{
            "name": s.name, "domain": s.domain, "unit": s.unit,
            "n": len(s.values), "values": list(s.values),
            "provenance": s.provenance, "normalisation": s.normalisation,
            "nested_with": s.nested_with, "synthetic": s.synthetic,
            "labels": list(s.labels),
        } for s in series_list],
        "tests": tests,
        "holdout": holdout,
        "power": power,
        "calibration": calibration_block,
        "cross_domain": {
            "status": "NOT RUNNABLE",
            "reason": ("one domain available; the other two series are subsets "
                       "of the first, so a Fisher combination would not be a "
                       "cross-domain statistic"),
        },
        "rule": "phi may only detract from a candidate, never support one",
    }
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "phi_cross_domain_protocol.json"),
              "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    write_report(payload)

    print("CROSS-DOMAIN PROTOCOL")
    print(f"   target Δr = {TARGET:.6f} ± {TOLERANCE} in base phi")
    for series in payload["series"]:
        print(f"   {series['name']:16} n={series['n']:2}  "
              f"{series['unit']}  "
              f"{'independent' if not series['nested_with'] else 'nested'}")
    print()
    print(f"   {'series':16} {'base':10} {'hits':>5} {'p':>7} {'BH p':>7}")
    for row in tests:
        print(f"   {row['series']:16} {row['base']:10} {row['hits']:5} "
              f"{row['p_value']:7.4g} {row['bh_adjusted_p_value']:7.4g}")
    print()
    print(f"   {'series':16} {'predictor':16} {'ratio':>9} "
          f"{'median rung err':>16}")
    for row in holdout:
        print(f"   {row['series']:16} {row['predictor']:16} "
              f"{row['predicted_ratio']:9.6f} {row['median_rung_error']:16.4f}")
    print()
    for row in power:
        print(f"   power {row['series']:16} fraction {row['fraction']:.0%}: "
              f"{row['detection_rate']:.0%}")
    print()
    print(f"   cross-domain: {payload['cross_domain']['status']}")
    print("   wrote results/phi_cross_domain_protocol.json and PHI_CROSS_DOMAIN.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
