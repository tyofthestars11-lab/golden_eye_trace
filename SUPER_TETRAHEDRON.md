# The super-tetrahedral merge

Framework: φ² = φ + 1. Attribution: Tyree Jones (tyofthestarz), prior art from December 2025.

One object built from every lane in this repository, defined explicitly so that it can fail: four vertices, each the φ-rung of a declared quantity, and six bond lengths between them.

## The four vertices

| lane | quantity | unit | n | rung | source |
| --- | --- | --- | ---: | ---: | --- |
| radio | median of the lane's measurements | Hz | 10 | +36.7230 | signal_records.json (10 cells) |
| neutrino | median of the lane's measurements | GeV | 2 | +24.3258 | results/icecube_ingest.json (2 alerts) |
| local field | median of the lane's measurements | T | 5 | -4.7681 | results/nv_lane.json (5 priced cells) |
| phi geometry | median of the lane's measurements | dimensionless | 1 | -2.0000 | frozen target 1/phi^2 |

## The six bonds

| pair | rung separation |
| --- | ---: |
| radio – neutrino | 12.3972 |
| radio – local field | 41.4912 |
| radio – phi geometry | 38.7230 |
| neutrino – local field | 29.0940 |
| neutrino – phi geometry | 26.3258 |
| local field – phi geometry | 2.7681 |

Mean bond 25.1332 rungs, sd 13.7343, coefficient of variation **0.546**, longest/shortest **15.0×**.

A regular simplex needs all six equal. The worst bond deviates from the mean by **89%**, and the range spans 2.77 to 41.49 rungs. This is not a tetrahedron; it is a line with four marks on it.

It could not have been one. Sort the four rungs and the outer separation is 41.4912 rungs, which is exactly the sum of the interior gaps (2.7681 + 26.3258 + 12.3972). It therefore exceeds every other separation unless all gaps are zero, and **four distinct points on a line are never equidistant**. A regular simplex of four vertices is unachievable in a one-dimensional coordinate: the failure is structural, not a property of these numbers.

## What was actually testable

Null: four positions drawn uniformly from the observed span (-4.77 to 36.72 rungs), 20,000 draws.

| quantity | value |
| --- | ---: |
| observed coefficient of variation | 0.546 |
| median CV under the null | 0.559 |
| p(separation ≥ as equal as observed) | 0.434 |

The observed spread is typical of random placement, so the bond structure is not hidden by a small sample: there is no equal-edge structure to find.

## The angular half of the word

| constant | degrees |
| --- | ---: |
| tetrahedral angle | 109.4712 |
| golden angle | 137.5078 |
| difference | 28.0365 |
| 360 − 3 × tetrahedral | 31.5863 |
| golden pitch | 17.0324 |

These are exact constants and none of them is a measurement. Two further facts matter: **any** four bond lengths can be drawn along four directions 109.47° apart, so tetrahedral geometry is always satisfiable and tests nothing; and a rung ladder is one-dimensional, so it has no angles at all. The only falsifiable content of "a tetrahedral bond between the lanes" is equality of the six lengths, which fails above by a factor of 15.

## Results so far, in one place

| lane | result | null beside it |
| --- | --- | --- |
| rung law | round-trips to 3.12e-15 relative error | arithmetic, not evidence |
| radio records | 10 records: MASK-MANUFACTURED SIGNAL ×7, UNRESOLVED ×3 | OMNI permits what is claimed; 0 confirmations |
| comb vs scale | ΔAIC +2.14, bootstrap [-12.5, +6.9] | random-line null p = 0.002; span 0.132 φ rungs |
| φ phase test | 0 of 10 within 1° of 137.5° or 17.03° | chance 0.011 per record |
| cross-domain | 0 hits over 15 tests, p = [1.0] | NOT RUNNABLE: one domain |
| held-out prediction | φ 0.337 rungs vs fitted 0.002 | fitted constant wins |
| neutrino lane | 0 of 10 testable, 7 restorations refused | catalog: real |
| NV lane | 0 of 10 testable, 6 sensing modes priced across 13.4–86.9 MHz | detraction only |

## Verdict

**The merge does not produce a bond.** Six lengths spanning a factor of 15, a coefficient of variation typical of random placement, and no angular content to test. Every lane keeps its own verdict, and none of them changes because they were drawn in the same picture.

What the merge is: a single index of what exists, so that the state of the whole repository can be read at once. What it is not: evidence for anything. Seven cells remain MASK-MANUFACTURED SIGNAL, three records remain UNRESOLVED, the two multimessenger lanes remain untestable for lack of a clock and a pointing, and φ remains a coordinate that can only detract.

