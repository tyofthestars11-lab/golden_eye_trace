# Preregistered Cross-Domain φ-Structure Test

## Frozen prediction

The primary target is Δr = 1/φ² = 0.381966011250105, with tolerance ±0.010000. Each Δr is computed from an adjacent within-domain ratio: log(x[i+1]/x[i]) / log(base).

Dimensional quantities are never compared directly across domains. Each dataset is tested only through ratios of its own positive measurements in its declared order.

## Dataset results

| dataset | φ hits | pairs | φ p-value | φ odds ratio |
|---|---:|---:|---:|---:|
| survey_cells | 0 | 9 | 1 | 1 |
| group_20240325 | 0 | 3 | 1 | 1 |
| trio_20240920 | 0 | 2 | 1 | 1 |

Fisher combined statistic: -0
Fisher combined p-value: 1

## Control bases

The control bases use the same fixed φ target in their own logarithmic coordinates and the same tolerance. This is a preregistered sensitivity comparison, not a refit of the target.

| base | combined Fisher p-value |
|---|---:|
| phi | 1 |
| base_2 | 1 |
| base_e | 1 |
| base_10 | 1 |
| base_sqrt2 | 1 |

## Decision rule

The prediction is supported only if every domain has an odds ratio ≥ 2, every domain has p < 0.01, the combined preregistered p-value is < 1e-4, and the φ result is stronger than the control-base results after the declared multiplicity adjustment. This script reports those quantities; it does not convert them into a claim about life or intelligence.

## Reproducibility

Configuration: `{"control_bases": {"base_10": 10.0, "base_2": 2.0, "base_e": 2.718281828459045, "base_sqrt2": 1.4142135623730951}, "input_order_rule": "numeric order if possible, else lexical order", "null": "shuffle observed values across fixed declared order", "permutations": 100000, "phi": 1.618033988749895, "primary_target_rung": 0.38196601125010515, "seed": 20261114, "tolerance": 0.01}`

A significant result would support only the narrower claim that the tested datasets share a reproducible multiplicative spacing near the preregistered φ interval. It would not by itself establish a universal law, consciousness, or a technosignature.
