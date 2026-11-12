# Methods

Every method below states its statistic, the assumption it rests on, the null it
is compared against, and the limit that would make it lie. The limits are not
footnotes; several of them have already produced wrong answers in this project.

## 1. Ingest and channelisation

`dynamic()` reads a FITS primary HDU, infers the axis order from sizes (the
product carries no CTYPE/CRVAL/CDELT), and checks the result against the
documented span and cadence. If the layout check fails it raises instead of
guessing.

- channel width `100547.19562243501` Hz, sample time `0.256` s
- 731 channels x 73.5 MHz; 7009 samples x 0.50 h in the largest file
- Nyquist drift `0.5 * df / tsamp` = 196,381 Hz/s

**Limit.** Without WCS there is no source position, so no barycentric correction.
Only a bound can be quoted: Earth's rotation and orbit induce a few Hz/s, three
orders of magnitude below every candidate here.

## 2. Incoherent dedoppler

`dedoppler()` shift-and-adds rows at trial drifts on a grid whose step keeps the
smear below half a boxcar, then `_smoothed_snr()` reports
`(smoothed - median) / (1.4826 * MAD)` and `calibrate_map()` rescales the map so
its bulk is one sigma wide.

- drift step `0.5 * boxcar * df / T`; at 1000 samples and boxcar 3 that is
  **589 Hz/s**, while the intrinsic width of a line in drift space is
  `df/T` = **392.8 Hz/s**. Both are quoted because they answer different
  questions: the step limits how finely two drifts can be separated, the width
  limits how precisely one drift is known. Confusing the two is how a 393 Hz/s
  "resolution" gets attached to a search that cannot separate drifts that close
- a rate beyond 0.5 cycles/sample folds: 0.6 is recorded as **-0.4**, sign
  included

**Limit.** The statistic is incoherent: it sums power, so it cannot distinguish a
narrowband carrier from a broad feature that happens to drift.

## 3. Significance: band-wide versus local

The survey's denominator is one MAD across the whole searched window. The local
run replaces it with a per-channel spread from a ±16 channel neighbourhood across
every drift trial, with one round of clipping.

| reference | what it assumes | failure |
| --- | --- | --- |
| band-wide | the band is uniformly quiet | a bright sub-band inflates significance |
| local | the neighbourhood is source-free | a bright source contaminates its own scale |

Measured range of the ratio across ten candidates: **0.30 to 31.18**. For the
20240325 group the ratio is 12-31 because the neighbourhood *is* the source, so
those values are lower bounds, not estimates.

## 4. Nonlinear event time

Five coordinates, each fitted on the first half of a track and scored on the
second: `linear`, `phase` (cumulative |dφ|), `chirp` (∫(f-f0)/B dt),
`accel` (chirp minus a survey-grid drift model), `energy` (cumulative excess
power). `tau_phi` is a feature, not a coordinate.

**Statistic.** Matched-filter excess, `sqrt(N)(mean - median)/sd`, defined on
intensity. A Hilbert phase-locking value is also reported but is an **envelope
proxy**: the product carries no carrier phase, and the PLV is not used to decide
anything.

**Null.** Best-of-five against best-of-five randomised monotone warps. A
shuffled-time null is degenerate for an amplitude statistic — a permutation keeps
the values, and the measured deviation is 1.8e-15 — so it cannot be the control.

**Limit.** `chirp` is fitted to the drift that defines the candidate, so a
constant-acceleration line is linearised by construction. Only the split-half
value and the warp comparison mean anything.

## 5. φ lock (feature only)

Lock is a phase advance: `nu = drift * tsamp / df` and
`theta = 1/phi^2 = 0.381966`. Distance is `|nu/theta - k|` with **k >= 1**;
`k = 0` is excluded because a stationary line is not an advance and including it
lets a near-zero rate score arbitrarily close to the lattice.

**Null.** 2,000 randomized-phase draws, best distance per draw.

**Result.** No candidate within 0.05 at k >= 1; control p = 1.0000. The φ
columns exist to be checkable, not to rank.

## 6. Gates

| gate | test | status in this archive |
| --- | --- | --- |
| repeat | same frequency and drift on another day | testable |
| persistence | same line in a longer integration, and in other windows | testable |
| off-source | a pointing away from the source | **UNTESTABLE** |
| barycentric | corrected drift needs a sky position | **UNTESTABLE** |
| independent telescope | a second instrument | **UNTESTABLE** |

Three gates are untestable, not failed. A candidate cannot be promoted past them
and the headline stays NULL RESULT.

## 7. Injection recovery

Synthetic unit-mean exponential noise, amplitude fixed at 3, rates
`0, 0.02, theta, 0.45, 0.6`. `0.6` is tracked at the folded `-0.4`.

**Why synthetic.** An earlier version injected into the real solar spectra and
recovered 0/8: the bandpass and excision stages act on that structure, so the
test measured them rather than the search.

## 8. Disposition

Closed vocabulary in `seti_search/signal_record.py`. Order of precedence follows
evidence, not preference: a near-Nyquist rate is an alias regardless of how
strong it is, four parallel lanes at one drift are astrophysical regardless of
how significant they are, and "confirmed" is never available here.

## Comb against scale family (`seti_search/model_select.py`)

**Input.** A group of three or more line frequencies in Hz, and the channel width.

**Transform.** Fit `f_n = f0 + n*df` and `f_n = f0*q**n` to the same quantity by
least squares - the scale fit is profiled, so for each `q` the best `f0` is closed
form and the single nonlinear parameter is searched exhaustively. Both models take
two parameters, so AIC and BIC are comparable; leave-one-out refits from scratch
per fold and predicts the omitted line in Hz and in φ rungs; a sign-flip bootstrap
of the better model's residuals puts a 95% interval on the AIC gap; dropping each
line in turn tests whether the preference is an artefact of one measurement. Two
physical readings outrank the fit: spacings on the channel grid are instrumental,
and one drift shared across many spacings is one broadband event.

**Output.** `results/comb_vs_scale.csv`, `results/heldout_predictions.csv`,
`results/expanded_phi_groups.json`, two figures, and exactly one of seven labels
per group.

**Claim.** Nothing about a source. A named model requires a group spanning at
least one φ rung, a decisive AIC gap, agreement from held-out prediction, a
bootstrap interval excluding zero, and - for a scale family - the same `q` on a
second separated group. `synthetic_calibration` reports the rate at which the
pipeline recovers a known model at the data's own residual size and span; when
that rate is near half, the label is INDETERMINATE DUE TO LIMITED SPAN and the
number is printed beside it.

## The light path (`scripts/light_path.py`)

**Input.** A day file, a set of frequencies, and the drift each peaks at.

**Transform.** Eight stages, each recording the number it produces: A raw power,
B bandpass ratio, C dedoppler sigma with nothing removed, D the excision mask and
its counts, E sigma after blanking, F the calibration factor, G the local-scale
sigma, H whether peak extraction returns the cell. Each stage is repeated for
three masks - none, the theory-scaled row test, and an empirical one scaled by the
measured scatter of row medians - in both a 96-channel and a 731-channel window.
Both Stokes planes go through the same aligned sum.

**Output.** `LIGHT_PATH.md`, `results/light_path_ledger.csv` (42 rows),
`results/light_path_stages.json`.

**Claim.** Only about the processing. It attributes a candidate's significance to
a stage; it does not confirm or refute a source. What it found is that the
excision mask, not the sky, supplies the sigma, and that is a statement about this
pipeline's arithmetic.
