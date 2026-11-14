# Cross-domain test of a fixed φ-rung interval

Framework: φ² = φ + 1. Attribution: Tyree Jones (tyofthestarz), prior art from December 2025.

Prediction, frozen before the run: adjacent ratios in an ordered series sit at Δr = 1/φ² = **0.381966011250105** rungs in base φ, within **±0.01**. Controls use the same fixed rung interval in bases 2, e, 10 and √2.

## The dataset problem, stated first

This is a cross-domain test and this repository has **one domain**. There is no genomics series and no physics series in it, and none was invented for this run: fabricated inputs would make the statistics meaningless. The three series below are all radio frequencies from one instrument, and two of them are subsets of the first, so **the cross-domain combination is not runnable and no combined p-value is reported.**

| series | domain | unit | n | normalisation | independence |
| --- | --- | --- | ---: | --- | --- |
| survey_cells | radio astronomy | MHz (rest frequency) | 10 | none needed: only within-series ratios are taken | independent domain |
| group_20240325 | radio astronomy | MHz (rest frequency) | 4 | none needed: only within-series ratios are taken | subset of survey_cells |
| trio_20240920 | radio astronomy | MHz (rest frequency) | 3 | none needed: only within-series ratios are taken | subset of survey_cells |

## Primary tests, with permutation nulls

10,000 permutations per test, shuffling the observed values across the fixed order. BH is applied across every test actually performed (5 bases × 3 series).

| series | base | hits / pairs | null mean hits | p | BH-adjusted p |
| --- | --- | ---: | ---: | ---: | ---: |
| survey_cells | phi | 0 / 9 | 0.000 | 1 | 1 |
| survey_cells | base_2 | 0 / 9 | 0.000 | 1 | 1 |
| survey_cells | base_e | 0 / 9 | 0.000 | 1 | 1 |
| survey_cells | base_10 | 0 / 9 | 0.000 | 1 | 1 |
| survey_cells | base_sqrt2 | 0 / 9 | 0.200 | 1 | 1 |
| group_20240325 | phi | 0 / 3 | 0.000 | 1 | 1 |
| group_20240325 | base_2 | 0 / 3 | 0.000 | 1 | 1 |
| group_20240325 | base_e | 0 / 3 | 0.000 | 1 | 1 |
| group_20240325 | base_10 | 0 / 3 | 0.000 | 1 | 1 |
| group_20240325 | base_sqrt2 | 0 / 3 | 0.000 | 1 | 1 |
| trio_20240920 | phi | 0 / 2 | 0.000 | 1 | 1 |
| trio_20240920 | base_2 | 0 / 2 | 0.000 | 1 | 1 |
| trio_20240920 | base_e | 0 / 2 | 0.000 | 1 | 1 |
| trio_20240920 | base_10 | 0 / 2 | 0.000 | 1 | 1 |
| trio_20240920 | base_sqrt2 | 0 / 2 | 0.000 | 1 | 1 |

Every base, including φ, produces zero hits and p = 1. No base can be distinguished from any other, because none of them is present.

## Held-out prediction: does φ predict better than a fitted interval?

The prediction for each held-out pair is `previous × base^Δ`. The baseline that matters is a **constant ratio fitted on the first half**, which is what any regularly spaced series obeys. Error is in rungs; lower is better.

| series | predictor | predicted ratio | median rung error | hit rate |
| --- | --- | ---: | ---: | ---: |
| survey_cells | phi | 1.201783 | 0.3401 | 0% |
| survey_cells | base_2 | 1.303116 | 0.3529 | 0% |
| survey_cells | base_e | 1.465162 | 0.3618 | 0% |
| survey_cells | base_10 | 2.409717 | 0.3732 | 0% |
| survey_cells | base_sqrt2 | 1.141541 | 0.3239 | 0% |
| survey_cells | fitted_constant | 0.986714 | 0.0731 | 0% |
| group_20240325 | phi | 1.201783 | 0.3367 | 0% |
| group_20240325 | base_2 | 1.303116 | 0.3505 | 0% |
| group_20240325 | base_e | 1.465162 | 0.3602 | 0% |
| group_20240325 | base_10 | 2.409717 | 0.3725 | 0% |
| group_20240325 | base_sqrt2 | 1.141541 | 0.3191 | 0% |
| group_20240325 | fitted_constant | 1.021022 | 0.0020 | 0% |
| trio_20240920 | phi | 1.201783 | 0.3402 | 0% |
| trio_20240920 | base_2 | 1.303116 | 0.3529 | 0% |
| trio_20240920 | base_e | 1.465162 | 0.3618 | 0% |
| trio_20240920 | base_10 | 2.409717 | 0.3732 | 0% |
| trio_20240920 | base_sqrt2 | 1.141541 | 0.3239 | 0% |
| trio_20240920 | fitted_constant | 1.025722 | 0.0110 | 0% |

φ is **not** the best predictor and is not close. The fitted constant wins on every series, because the data are an arithmetic comb in frequency, not a geometric ladder in any base: the fitted ratio is ≈1.02 while φ^0.382 = 1.2018, so φ's prediction misses by ~0.12 rungs on the narrow series and by ~0.34 on the ten-cell one, whose first pair spans a day boundary. No fixed base predicts an arithmetic comb; that is the finding.

## Power: what a null with n = 9, 3, 2 can see

The law is injected on a known fraction of steps, keeping the series' own drift elsewhere; fraction 0 is the false-positive arm. 200 trials, 500 permutations each (p-floor 1/501), detection at p ≤ 0.05.

| series | fraction of steps obeying Δr | detection rate |
| --- | ---: | ---: |
| survey_cells | 0% | 0% |
| survey_cells | 25% | 100% |
| survey_cells | 50% | 100% |
| survey_cells | 75% | 100% |
| survey_cells | 100% | 100% |
| group_20240325 | 0% | 0% |
| group_20240325 | 25% | 0% |
| group_20240325 | 50% | 0% |
| group_20240325 | 75% | 0% |
| group_20240325 | 100% | 82% |

Two features of this table matter more than the headline numbers. First, at ten points a *single* planted interval is enough to be flagged: the permutation null almost never produces an interval near 0.382 by chance, so one hit gives p ~ 1e-4. A detection therefore means "an interval is present", not "a law holds". Second, at four points the same test reaches only 82% power even when every step obeys the law, because shuffling four values containing one long step often reproduces a hit. A null on nine pairs excludes a blanket interval; a null on three pairs excludes almost nothing.

## Calibration check

| control series | drift in rungs | hits / pairs | p |
| --- | ---: | ---: | ---: |
| exponential drift 1.004 (survey_cells length) | 0.0083 | 0 / 9 | 1 |
| exponential drift 1.004 (group_20240325 length) | 0.0083 | 0 / 3 | 1 |

A series with perfect geometric structure at a different interval returns zero hits, which is the correct answer: the test is not simply detecting 'having structure'.

## Verdict

1. **φ does not produce a cross-domain excess.** The test cannot be run cross-domain with one domain, and within the domain available there are zero hits at every base.
2. **φ does not produce better held-out predictions than the controls.** It is worse than a fitted constant on every series, and no better than bases 2, e, 10 or √2, all of which also fail because the spacing is arithmetic.
3. **The failed tests stay visible.** Zero hits is reported as zero hits, not as an absence of data, and the power table says exactly how much that zero excludes.
4. **This is not evidence of life or intelligence, and it would not have been even if positive.** It is not calibrated, not replicated, and has no independent domains; a hit would have supported only the claim that two series share a spacing.

## What would make the cross-domain test real

Three CSVs with `order,value` and, where available, `id,label`. Declared units; within-series ratios only. Genomics: gene positions or expression levels in genomic order. Astronomy: the series used here is legitimate, but it needs a second, independent one — a different instrument or band. Physics: spectral lines or energy levels in a declared order.

Then, and only then, the Fisher combination is a cross-domain statistic. The honest prior from this run is that it will be null, because the interval it tests is not present in the one domain that could be examined, at a resolution that can see it.

