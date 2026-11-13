> **Disqualified premise.** The two multi-line groups below were selected by
> the excision mask, not by a valid detection: their cells read 1.44, 0.54,
> -0.15 and -0.87 sigma without the theory-scaled row test, and 954, 604, 447
> and 171 sigma with it (`LIGHT_PATH.md`). Their spacing arithmetic is kept as
> the worked example of selection bias and is **not** a regularity claim.
> `prereg/phi_address.json` carries `status` per group: none of the two is a
> candidate.

# Rung addresses for every measured quantity

Framework: φ² = φ + 1. Attribution: Tyree Jones (tyofthestarz), prior art from December 2025. Coordinate: `rung(x) = ln x / ln φ`.

Every record's quantities are addressed in `prereg/phi_address.json`. What follows is what the addressing *buys*: the four tests, on the real candidates, with a null beside each statistic.

## What the tests can and cannot do

A single rung drift is not a finding. For any smooth population the fractional part of `ln x / ln φ` is uniform, so a value lands within 0.05 rungs about a tenth of the time, and one number close to an integer is not evidence of anything. The table below shows that this is not a property of φ: the same median fractional distance comes out of base 2, e, 10 and √2 within the noise, because they are all logarithms. Where the coordinate earns its place is in *relations* — a rung difference is a ratio, and a repeated ratio is a repeated scaling.

| base | median fractional distance | band span (rungs) | uniform expectation |
| --- | --- | --- | --- |
| phi | 0.2770 | 0.351 | 0.25 |
| 2 | 0.4794 | 0.244 | 0.25 |
| e | 0.3284 | 0.169 | 0.25 |
| 10 | 0.3253 | 0.073 | 0.25 |
| sqrt2 | 0.0412 | 0.488 | 0.25 |

Across bases the median moves by only 0.4382, so no address is φ-specific at the level of one value.

## The clarification that does work: comb or scale family

A group of lines can be spaced by a fixed **frequency step** — a clock, a channel grid, an intermodulation product, all of which are comb-like — or by a fixed **factor**, which is one thing measured at several scales. On a linear axis inside a narrow window the two are nearly indistinguishable, and they mean opposite things. Rung differences separate them: a comb has constant Δf and drifting Δrung, a scale family the reverse.

| date | lines | step spread | rung spread | p(comb) | p(scale) | read |
| --- | --- | --- | --- | --- | --- | --- |
| 20240325 | 4 | 0.1002 | 0.0784 | 0.0025 | 0.0018 | regular; comb and scale family fit equally |
| 20240920 | 3 | 0.2094 | 0.2319 | 0.0973 | 0.1082 | not regular at 5% |

Each p compares the observed spread against 4000 groups drawn uniformly across the same band, so a regular-looking arrangement of four lines in 9.55 MHz has to beat luck before it is read as structure. `p(comb)` is the chance of a spread that tight in frequency, `p(scale)` the same in logs.

## Ladder steps between neighbours

Per test A, a pair is a ladder step only if the rung difference is an integer; otherwise it is a plain scale factor and the residual says how exactly.

| date | pair (MHz) | Δrung | integer step? | factor |
| --- | --- | --- | --- | --- |
| 20240325 | 45.9843 → 46.9828 | +0.0446 | no | 1.021713 |
| 20240325 | 46.9828 → 47.9380 | +0.0418 | no | 1.020331 |
| 20240325 | 47.9380 → 48.9937 | +0.0453 | no | 1.022023 |
| 20240920 | 45.4206 → 46.5890 | +0.0528 | no | 1.025722 |
| 20240920 | 46.5890 → 47.5358 | +0.0418 | no | 1.020323 |

No pair is an integer step, so no pair is a ladder relation. That is the expected outcome and it is the honest one: the lenses here are scale factors, and a scale factor is a description of how two lines differ, not a claim about what produced them.

## Repeated scales

None. `repeated_scales` needs three *independent* pairs sharing a rung difference, and it applies a shared-index rule: pairs that reuse a line are one object, not two. A chain of four lines gives at most two independent pairs, so a same-factor chain is not by itself a repeated-scaling artifact — the test says so rather than manufacturing a hit.

## The drift population

Across the 10 non-zero drifts, fractional rungs are uniform with p = 0.165. A low p here would mean the drift values were quantized by the search grid rather than by the sky, and there would be nothing left to interpret.

## Reading

The framework is a ruler, and a ruler's job is to make two things comparable rather than to decide between them. What it adds to a binary method is precisely the comb-versus-scale question: an accept/reject on each line cannot see it, and it is the difference between a group of sources and one source counted four times. What it cannot do is turn a small drift into a detection — that is arithmetic, not caution, and the base table above is the demonstration.

