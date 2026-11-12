# The light path

The measured field, followed from the file to the disposition, one stage at
a time. The question was where a candidate's significance comes from. The
answer is not the sky.

## Headline

**The 954 sigma group and the 8.65 sigma trio are both manufactured by the
excision mask.** At the candidates' own drift, in the survey's own window, the
same cells read:

| day | no mask | theory row test | empirical row test |
| --- | --- | --- | --- |
| 20240325 | 1.44 0.54 -0.15 -0.87 | **954.18 604.16 447.18 171.27** | 1.44 0.54 -0.15 -0.87 |
| 20240920 | -1.20 -0.42 0.11 | **8.65 7.90 6.49** | -1.20 -0.42 0.11 |

with integrations blanked:

| day | no mask | theory | empirical |
| --- | --- | --- | --- |
| 20240325 | 0 | 505 | 0 |
| 20240920 | 0 | 228 | 1 |

The theory test flags half the integrations; the empirical one flags none.

## Why

`excise_rfi`'s row test asks whether an integration's band median is further
from the median of medians than `sigma * 1.25/sqrt(n_f) * typical_std`. The
`1.25/sqrt(n_f)` term is the standard error of a median *if the channels were
independent*. On a real OVRO-LWA solar product the channel gains are
correlated - the Sun is one source across the band - so the row medians
scatter far more than that, and a two-sided 6 sigma threshold on the
theoretical scale flags integrations that are perfectly ordinary. It rejects
both bright *and* faint rows, keeping a subset selected on the same quantity
the search then measures, and the dedoppler sum of that subset is a
different - much quieter - statistic.

The empirical test keeps the structure and replaces the theoretical scale
with `1.4826 * MAD(row medians)`, which is what "one integration far
brighter than the rest" has to be measured against. It is opt-in
(`row_noise="empirical"`); the default is unchanged, because flipping it
would silently rewrite every sigma this repository has published.

## The 2x2

Declared design: frequency span x integration duration, with the drift window
guard printed beside each cell.

| window | duration | integrations | drift window held | candidates |
| --- | --- | --- | --- | --- |
| narrow-96 | short-1000 | 1000 | - | 10 |
| narrow-96 | full | 2802-7009 | 3%-7% | 6 |
| wide-731 | short-1000 | 1000-1000 | 100%-100% | 32 |
| wide-731 | full | 2802-7009 | 21%-52% | 32 |

Candidate counts in this table are counts produced by the theory-scaled mask,
so they are an upper bound on what the uncorrected pipeline reports, not a
measurement of the sky. The wide-band cells truncate the drift window because
a Nyquist drift needs more channels than the band has at these durations; that
is the stabilizer's warning, not a null result.

## Every stage, both windows

| day | window | mask | A raw | C sigma (map) | D rows/ch blanked | E sigma | F factor | G local sigma | H recovered |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20240325 | narrow-96 | none | 1.680e-05 | 1.26 | 0/0 | 1.44 | 0.870 | 1.68 | False |
| 20240325 | narrow-96 | theory | 1.680e-05 | 359.23 | 505/0 | 954.18 | 0.377 | - | True |
| 20240325 | narrow-96 | empirical | 1.680e-05 | 1.26 | 0/0 | 1.44 | 0.870 | 1.68 | False |
| 20240325 | wide-731 | none | 1.992e-05 | 0.04 | 0/156 | 0.03 | 1.527 | - | False |
| 20240325 | wide-731 | theory | 1.992e-05 | -0.57 | 530/208 | -0.29 | 1.965 | - | False |
| 20240325 | wide-731 | empirical | 1.992e-05 | 0.32 | 0/208 | 0.20 | 1.587 | - | False |
| 20240920 | narrow-96 | none | 7.263e-01 | -1.02 | 0/0 | -1.20 | 0.854 | -4.09 | False |
| 20240920 | narrow-96 | theory | 7.263e-01 | 7.73 | 228/0 | 8.65 | 0.893 | 18.61 | True |
| 20240920 | narrow-96 | empirical | 7.263e-01 | -1.02 | 1/0 | -1.20 | 0.854 | -4.08 | False |
| 20240920 | wide-731 | none | 2.619e+00 | -0.31 | 0/0 | -0.18 | 1.702 | -1.68 | False |
| 20240920 | wide-731 | theory | 2.619e+00 | -0.16 | 504/79 | -0.09 | 1.675 | - | False |
| 20240920 | wide-731 | empirical | 2.619e+00 | -0.27 | 0/79 | -0.16 | 1.667 | -0.77 | False |

One line per configuration shown at its first cell; every cell is in
`results/light_path_ledger.csv`, 42 rows.

## What the stages say

- **A to B** (bandpass) is reversible and moves nothing much: the cells sit at
  the band's own level.
- **C** (dedoppler, nothing removed) puts every cell inside 1.5 sigma. This is
  the number that does not depend on a mask.
- **D** is where the story is. In the narrow window the theory test removes
  half the integrations and no channels; in the wide window it also removes
  a quarter to half the channels. Which set survives is decided by the band's
  brightness at that moment, which correlates with what the Sun was doing.
- **E to F** (calibration) is a factor of 0.38 to 0.87 - real, but not the
  effect. The effect is in D's selection, which C never sees.
- **G** (local scale) is undefined where a neighbourhood has been blanked to
  nothing, which is itself a diagnostic: on 20240325 the local reference has
  no live data to measure.
- **H** (peak extraction) returns the cells under the theory mask and not
  otherwise.

## Verdict, corrected

The four frequencies at 45.98-48.99 MHz are **not validated candidates**. The
earlier `INDETERMINATE DUE TO LIMITED SPAN` was a verdict about whether a
regularity is a comb or a scale family; it presupposed that the four
frequencies were selected by a valid detection. They were selected by a mask
that manufactures 954 sigma out of 1.44. Under an honestly scaled mask they
are at 1.44, 0.54, -0.15 and -0.87 sigma, which is not a detection at all.

The mathematical result stands on its own terms and keeps no evidential
weight: four frequencies 0.988-1.056 MHz apart are regular at p = 0.002
against random lines, and a comb and a scale family fit them equally. Four
numbers chosen by an invalid selection can be regular for the same reason any
four chosen numbers can.

## Polarization

The product carries two Stokes planes and the reader kept one. Both are now
carried through the same aligned sum, and the ratio is in the ledger. It is
reported, not interpreted: without the header's CTYPE cards there is no
guarantee which two Stokes parameters these are, so a ratio is a measurement
whose meaning is not established - which is the standard this record holds
everything else to.

## Reproduce

```
python3 scripts/light_path.py          # ~11 s, both days, 2x2 x 3 masks
```
