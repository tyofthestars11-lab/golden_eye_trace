# Proportions, not thresholds

Framework: φ² = φ + 1. Attribution: Tyree Jones (tyofthestarz), prior art from December 2025.

A ratio does not guess and does not terminate. Two consequences are measured here.

## 1. Zero out of fourteen is a distribution, not a verdict

Counting every adjacent pair tested for the frozen interval: **0 hits in 14 pairs**. With a uniform prior that gives a posterior over the per-pair rate:

**Rate ~ Beta(1, 15)**, median **0.045**.

| question | answer |
| --- | ---: |
| P(rate ≤ 0.019) — the largest rate still consistent with 1 hit in 50 | 1 − 0.739 = 0.261 |
| P(rate ≤ 0.10) | 0.794 |
| P(rate ≤ 1/φ² = 0.382) | 0.999 |
| 95% upper bound on the rate | 0.181 |

So the statement supported by the data is a constraint, not a refutation: **a per-pair occurrence rate below 0.18 at 95%**, with a 21% posterior probability that the rate exceeds 0.02 and a 0% probability that it exceeds the frozen 0.382. A blanket law - every step on the interval - is excluded far harder than that: if each pair obeyed with 90% probability, zero hits in fourteen has likelihood 10⁻¹⁴.

## 2. Scale invariance: the property a φ description must have

A multiplicative step has the same rung interval at every frequency. An additive step does not: the same 1 MHz is a different fraction of the band at 13 MHz and at 87 MHz. Measured:

| scale | step | at | ratio | Δr (base φ) |
| --- | ---: | ---: | ---: | ---: |
| one channel | 0.100547 MHz | 13.40 MHz | 1.007504 | 0.0155 |
| one channel | 0.100547 MHz | 30.00 MHz | 1.003352 | 0.0070 |
| one channel | 0.100547 MHz | 47.00 MHz | 1.002139 | 0.0044 |
| one channel | 0.100547 MHz | 60.00 MHz | 1.001676 | 0.0035 |
| one channel | 0.100547 MHz | 86.90 MHz | 1.001157 | 0.0024 |
| comb step (1 MHz) | 1 MHz | 13.40 MHz | 1.074627 | 0.1496 |
| comb step (1 MHz) | 1 MHz | 45.98 MHz | 1.021747 | 0.0447 |
| comb step (1 MHz) | 1 MHz | 47.24 MHz | 1.021168 | 0.0435 |
| comb step (1 MHz) | 1 MHz | 53.78 MHz | 1.018593 | 0.0383 |
| comb step (1 MHz) | 1 MHz | 86.90 MHz | 1.011507 | 0.0238 |
| whole band | 13.4 MHz | 13.40 MHz | 6.485075 | 3.8850 |

At one channel step the interval moves from 0.0024 to 0.0155 rungs across the band - a factor of 6.5. At the comb step it moves from 0.0238 to 0.1496, a factor of 6.3.

The frozen target, by contrast, has one value everywhere it is evaluated:

| interval | at 13.4 MHz | at 47 MHz | at 86.9 MHz |
| --- | ---: | ---: | ---: |
| 1/φ² (multiplicative) | 0.3820 | 0.3820 | 0.3820 |
| 1 MHz (additive) | 0.1496 | 0.0438 | 0.0238 |

**This is the whole difference.** A φ description predicts that the measured interval is the same number at every scale; an arithmetic comb predicts that it falls as 1/f. The instrument's step is 100.5 kHz, fixed in Hz, so the measured intervals *must* fall as 1/f, and they do. The data cannot be scale-invariant in this coordinate, and no choice of base repairs that: it is a property of the spacing, not of the logarithm.

## What that means for the earlier zero

The zero from the cross-domain run is therefore not a failed guess. It is the proportion of pairs landing on a fixed multiplicative interval in a series whose spacing is additive, and that proportion is 0 by construction at all but one frequency - the one where this comb's step happens to equal the interval exactly, **4.96 MHz**. The survey searched no cell there.

## The honest form of every result in this repository

| result | as a proportion |
| --- | --- |
| mask dependence | 954σ becomes 1.44σ when the mask changes: the ratio of the two is 663, and it is a property of the mask |
| comb vs scale | ΔAIC +2.14 with a bootstrap interval spanning zero: the evidence is a distribution, and its 95% mass includes 0 |
| φ phase | 0 of 10 phases within tolerance, posterior rate < 0.18 |
| multimessenger | 0 of 10 testable: a ratio of 0/10, not a negative detection |
| rung law | exact to 3.1×10⁻¹⁵: the one statement here that terminates |

Nothing in this document raises a disposition. φ remains a coordinate that can only detract, and the proportions above are the form that detraction takes.

