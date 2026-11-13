# Comb or scale family

> **Superseded premise (see [LIGHT_PATH.md](LIGHT_PATH.md)).** The four
> frequencies this document fits are not validated candidates. Under an
> honestly scaled excision mask they sit at 1.44, 0.54, -0.15 and -0.87 sigma,
> against 954 sigma under the theory-scaled mask that produced them. Everything
> below is a correct statement about four chosen frequencies and carries no
> evidential weight about the sky.

Does the regularly spaced group at 45.98-48.99 MHz keep its regularity when
the frequency coverage is widened, and does the wider coverage say whether it
is a frequency comb or a scale family?

**Answer: the group is regular and undecided.** It is well described by both
models, and at its span and residual size this pipeline recovers the true model
in roughly a quarter to a third of synthetic trials - worse than a coin toss -
so naming either one would be choosing, not measuring. It is not evidence of a
scale law and it is not evidence of anything artificial.

## 1. What was expanded, and what was declared first

The survey dedopplered a 96-channel window in the middle of each file. Every
file is 731 channels wide (13.4-86.9 MHz), so 87% of the frequency coverage had
never been searched. `prereg/expanded_prereg.json` declares the whole-band,
all-integration search before it was run; `prereg/expanded_prereg_bc.json`
declares two further variants and records plainly that it was written *after*
the first run, because that run changed frequency and duration together and
could not attribute a lost candidate to either.

| variant | channels | integrations | drift window held | candidates |
| --- | --- | --- | --- | --- |
| A | 731 | all integrations | 21%-52% | 32 |
| B | 731 | first 1000 | 100%-100% | 32 |
| C | 96 | all integrations | 3%-7% | 6 |

Variant A is the honest whole-band search, and the window guard says it cannot
be read as a null result: a Nyquist drift window needs 0.5 channel per sample,
which over 2802-7009 integrations is 1400-3500 channels of trail in a
731-channel band, so between half and four fifths of the drift range is
truncated. Variant B holds the survey's 1000 integrations and keeps the window
at 100%; variant C holds the survey's 96 channels and expands duration only.

## 2. Why the group cannot yet be tested in the wider band

Two window-dependent gates, both measured with the same data, same threshold
and same integration count:

| gate | 96-channel window | 731-channel band |
| --- | --- | --- |
| noise scale (`calibrate_map`) | 0.377 | 1.965 |
| channels blanked by the carrier test | 0/96 | 208/731 |
| time samples blanked | 505/1000 | 530/1000 |

The scale differs by a factor
5.2, so the same cell is scored 954 sigma in one window and 14 in the other.
And the channels holding the four lines are among those the wide-band carrier
test blanks: the quarter of the band that the carrier test removes is not the
same quarter in a 96-channel window. The group is therefore **masked**, not
refuted, by the wider band - and masked is not a result either.

The near-Nyquist trio behaves differently: its channels survive in both
windows, but its significance goes from 8.65 sigma to -0.04.
504 of 1000 time samples are blanked in the wide band against 228 in the narrow one, so the integrations the peak was built from are not in the
wide-band clean data. It was already classified an instrumental alias, and
this does not change that.

## 3. Every group, one label each

| source | date | lines | span (phi rungs) | dAIC (+ = scale) | boot CI | held-out | null p | label |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| survey-96ch-1000int | 20240325 | 4 | 0.132 | 2.14 | [-12.5, 6.9] | scale | 0.002 | MASK-MANUFACTURED SIGNAL |
| survey-96ch-1000int | 20240920 | 3 | 0.095 | -0.62 | [-1.7, 2.4] | undefined | - | MASK-MANUFACTURED SIGNAL |
| expanded-A-731ch-allint | 20240326 | 14 | 0.707 | 5.87 | [1.2, 14.1] | tie | - | ASTROPHYSICAL STRUCTURE |
| expanded-A-731ch-allint | 20240326 | 4 | 0.095 | 0.30 | [-0.6, 0.9] | scale | 0.038 | NEITHER MODEL FIT |

Seven labels were available; four verdicts were returned. No group returned
SCALE FAMILY SUPPORTED or FREQUENCY COMB SUPPORTED, because none of them spans
the single phi rung the rule requires, and a scale claim additionally needs the
same q on two separated groups - and there is only this one group.

## 4. The 20240325 group in full

- 4 lines, span 3.009 MHz = 0.132 phi rungs
- comb: step 0.998347 MHz, AIC 95.17, BIC 93.94, rms 0.214 channels
- scale: q = 1.02125582, AIC 93.03, BIC 91.80, rms 0.164 channels
- delta AIC 2.14 favours the scale family, and its 95%
  sign-flip bootstrap interval is [-12.5, 6.9], which contains zero
- held out one line at a time: scale wins 3 of 4 folds, comb 1, on mean errors of 27.4 kHz and 39.3 kHz
- dropping a single line flips the preference in 2 of 4 cases, so the ranking is not robust to the set
- random lines drawn inside the group's own span exceed this AIC gap in 0.2% of 1500 draws (median gap 0.09)

So the group is *more* regular than random lines (p = 0.002) while the choice between the two kinds of
regularity is unresolvable. Both statements are true at once and the record
keeps them apart.

## 5. What the pipeline can actually resolve

Synthetic groups with a known generator, fitted both ways, correct rate:

| residual jitter | comb correct | scale correct |
| --- | --- | --- |
| 1 kHz | 100% | 100% |
| 5 kHz | 83% | 81% |
| 10 kHz | 56% | 54% |
| 20 kHz | 33% | 26% |
| 40 kHz | 18% | 17% |
| 80 kHz | 7% | 8% |

| span (4 lines, 20 kHz jitter) | comb correct | scale correct |
| --- | --- | --- |
| 3 MHz | 33% | 26% |
| 10 MHz | 100% | 100% |
| 30 MHz | 100% | 100% |
| 60 MHz | 100% | 100% |

The real group's residuals are 0.21 and 0.16 channels, about 16-22 kHz. On the first table that
sits between the 10 and 20 kHz rows: a correct call about half to a quarter
of the time, and a wrong one the rest. The span table shows the other way
out - at the same 20 kHz jitter, a 10 MHz span recovers the model every
time, and 30 MHz does too. The binding limit is span, not line count.

What would decide it, in order of cost:

1. a continuation of this group over 10 MHz or more - the same regularity
   carried by more lines, at the residual this data already has
2. residuals below about 5 kHz at the present 3 MHz span, which needs a
   narrower channel or a longer integration on a signal that holds it
3. the same q on a second, separated group, which the scale claim requires
   before it can be made at all

## 6. Physical controls

| date | steps off an integer channel | one event? | drift spread | span / spacing |
| --- | --- | --- | --- | --- |
| 20240325 | 0.070, 0.500, 0.500 | no | 0.03% | 3.0 |
| 20240920 | 0.380, 0.417 | no | 0.06% | 2.0 |
| 20240326 | 0.490, 0.313, 0.425, 0.225, 0.154, 0.357, 0.099, 0.237, 0.152, 0.230, 0.390, 0.446, 0.031 | yes | 0.02% | 13.4 |
| 20240326 | 0.046, 0.330, 0.323 | no | 0.00% | 2.9 |

No group is a comb on the channel grid: no spacing lands within
0.02 channel of an integer. The 20240325 steps are 0.069, 0.500 and 0.500
channels off integer - two of the three sit exactly between channels, which is
recorded as a suspicion rather than used as a verdict, since the tolerance
that would make it fire was chosen after seeing it.

The 20240326 event is the useful negative control: 14 lines over 22 MHz, all
with the same drift to within 0.02%, which is one broadband event crossing the
band, not fourteen lines at a scale. The pipeline calls it ASTROPHYSICAL
STRUCTURE rather than manufacturing a model for it, which is the behaviour the
test was built to have.

## 7. Reproduce

```
python3 scripts/expanded_search.py --prereg        # refuses a second time
python3 scripts/expanded_search.py --variant A     # ~90 s, all six days
python3 scripts/expanded_search.py --variant B     # ~8 s
python3 scripts/expanded_search.py --variant C     # ~45 s
python3 scripts/band_calibration_check.py          # window controls
python3 scripts/phi_model_discrimination.py        # this document and the tables
```

Standing rule unchanged: phi is a coordinate, not a criterion. A rung address
can make two things comparable and it can put a group on a ladder; it cannot
promote a candidate, and nothing here does.
