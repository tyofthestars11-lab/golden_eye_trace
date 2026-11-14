# Failure modes

Rejected methods are kept as fixtures. Each entry names the wrong result, the
tell that exposed it, and what replaced it. The machine-readable form is
`seti_search/regressions.py`; the executable form is in `tests/`.

## 1. Narrow window truncates a drift trail — **the 0/10 run**

| | |
| --- | --- |
| method | 32-channel search windows with a local noise reference |
| result | **0 of 10 survive**, including a 954σ lane |
| tell | a 954σ line reporting 0σ is a statement about the window |
| cause | the 20240325 lanes drift 0.105 channels/sample = **105 channels** across 1000 samples; a 32-channel window holds 30% of the trail |
| fix | keep the survey's 96-channel window, change only the noise scale |
| after | **9 of 10 survive** |
| fixture | `tests/test_window_truncation.py`, `prereg/regression_narrow_window.json` |

The test asserts the geometry (trail ≈ 105 channels, ≤ 30% inside the window) and
the consequence on synthetic data, plus a slow-lane control showing the 32-channel
window is not broken *per se*.

## 2. k = 0 in the φ lock

| | |
| --- | --- |
| result | control p = 0.0265 from a 172 Hz/s line |
| cause | k = 0 is the stationary lattice point; distance to it shrinks with the rate |
| fix | `k = max(1, round(nu/theta))`, written into `prereg/phi_lock.json` |
| after | p = 1.0000 |

## 3. Alias fold sign

| | |
| --- | --- |
| result | 0/2 aliased injections recovered |
| cause | 0.6 folds to **−0.4**, and the test compared against +0.4 |
| tell | a recovered alias at −157 kHz/s scored against +157 |
| fix | signed fold, verified exactly at six rates |
| after | fold exact at every probed rate |

## 4. Injection into real data

| | |
| --- | --- |
| result | 0/8 recovered, every rate failing equally |
| cause | the solar spectra carry bandpass structure the excision stages act on, so the test measured those stages |
| fix | synthetic unit-variance exponential noise |
| after | 6/8 recovered, 0/3 false positives |

## 5. Shuffled null as the primary control

| | |
| --- | --- |
| result | null percentile 1.000 for every candidate and coordinate |
| cause | a permutation preserves the sample values, so an amplitude statistic cannot move; measured deviation 1.8e-15 |
| tell | every candidate passing every null at the 100th percentile |
| fix | best-of-five randomised monotone warps, equal degrees of freedom |
| after | 0 of 12 beat the warp-gain control |

## 6. Centroid in chunks

| | |
| --- | --- |
| result | the two days came out **swapped** |
| cause | ~3σ per chunk, so the centroid followed whichever channel was brightest |
| tell | each day recovering the other's frequency |
| fix | full-integration dedoppler plus a within-day window scan |
| after | neither line persists; **NOT A REPEAT** |

## 7. Band-wide noise reference — **superseded, not invalid**

| | |
| --- | --- |
| result | significances wrong by 0.03× to 3.0× |
| cause | one MAD across 9.66 MHz is not a constant false-alarm rate |
| fix | per-channel local scale |
| after | 9 of 10 survive; the set changes by one candidate |

This one is a scaling error, not a logical one. It is listed here because a fixed
threshold applied to a mis-scaled statistic is indistinguishable from a real
detection until you measure the scale.

## 8. φ treated as evidence

| | |
| --- | --- |
| result | no lock survives any control |
| evidence | population adjusted p = 0.455 and 0.935; echo harness nulls overlap; chance rate 0.264 vs 2 of 3 |
| rule | φ may only detract from a candidate, never promote one |
| fixture | `tests/test_rungs.py` |

## 9. External claims that are not measurements

The `phi-quantum-pulses` repository contains statements of the form "all 99
blocks verified with φ parameters, zero failures" and a NIST application
proposal. These are arithmetic consistency checks against a parameter set and a
proposal, respectively. They are not detections, are not used to produce any
number here, and are preserved unaltered in `phi-quantum-pulses/`.

The general lesson, which is why they are listed at all: a verification that can
only pass is not a test. Every fixture above failed at least once.

## 10. Reading a candidate out of one window width

Widening the band from the survey's 96 channels to all 731 moves the same four
lines from 954 sigma to 14 sigma with no change to the data, the threshold or the
integration count, and at full band they are gone entirely. Two gates do this and
neither is the sky:

- `calibrate_map` measures **one** spread for the whole map. The 73.5 MHz band's
  spread is `1.965` where the 96-channel window's is `0.377` - a factor 5.2, so
  every significance scales with the window.
- `excise_rfi`'s carrier test compares each channel against the band around it.
  In the full band it blanks **208 of 731** channels, and the channels holding the
  group are among them. In the 96-channel window it blanks none.

The 20240920 trio shows the third version: its channels survive, but the wide
band blanks **504** of 1000 time samples against **228** in the narrow window, so
the integrations the peak was built from are not in the wide-band data.

**Tell:** a candidate whose significance or existence depends on the window,
reported as if the window were the sky. **Fix:** print the window guard beside
every result, keep the window fixed within a comparison, and report a masked
candidate as masked - `prereg/band_calibration.json` records both gates per
candidate, and regression fixture `window-dependent-visibility` pins it.

## 11. An excision threshold that fabricates candidates

The row test in `excise_rfi` flags an integration when its band median departs
from the median of medians by more than `sigma * 1.25/sqrt(n_f) * typical_std`.
That denominator assumes independent channels; in a solar dynamic spectrum the
Sun spans the band, the gains are correlated, and the row medians scatter far
more. The threshold is also two-sided, so it deletes bright *and* faint
integrations and keeps a subset selected on the very quantity the search measures.

Measured at the candidates' own drift, same window, same calibration:

| day | no mask | theory row test | empirical row test |
| --- | --- | --- | --- |
| 20240325 | 1.44 0.54 -0.15 -0.87 | **954.18 604.16 447.18 171.27** | 1.44 0.54 -0.15 -0.87 |
| 20240920 | -1.20 -0.42 0.11 | **8.65 7.90 6.49** | -1.20 -0.42 0.11 |
| rows blanked | 0 / 0 | 505 / 228 | 0 / 1 |

**Tell:** a candidate whose sigma collapses when only the excision threshold
changes, and a mask that removes a large fraction of integrations rather than a
handful. **Fix:** `row_noise="empirical"` scales the same test by the measured
scatter of the row medians. It is opt-in, because switching the default rewrites
every sigma this repository has published - a decision to take deliberately, not
as a side effect.

## 12. Absence of a neutrino treated as evidence of artificiality

A neutrino detector has an energy threshold of order TeV, an angular resolution
of 0.5-15 degrees, a duty cycle below unity and a sky coverage that depends on
whether the event is up-going or down-going. A radio feature with no coincident
neutrino is therefore compatible with a real source that is neutrino-quiet, one
below threshold, one in a part of the sky seen through the Earth, and one that
sent its event while the detector was not taking data. "No neutrino, so it must
be RFI" is not an inference; it is a preference.

**Tell:** a claim that a signal is artificial because a catalog query returned
nothing, with no statement of the detector's exposure, threshold or coverage.
**Fix:** `Coincidence.absence_is_informative` is a field fixed at `False`, an
empty query must name the limit that made it empty, and `guard_support` refuses
to attach an absence to a record at all. When the radio record has no pointing
or no UTC timestamp the lane reports `UNINFORMATIVE`, not a null.

## 13. A neutrino coincidence used to restore an invalid candidate

The seven cells at 45-49 MHz were manufactured by the excision mask. A neutrino
that happens to point at that part of the sky is a real particle coincident with
an artifact, which is a statement about the catalog's density and nothing about
the radio data. Allowing it to restore the record would launder a selection
effect through an independent instrument, and the independence of that instrument
is exactly what would make the laundering hard to see.

**Tell:** a disposition that rises when a second messenger is added, with the
radio gates that failed still failing. **Fix:** `guard_support` raises
`NeutrinoAsEvidence` for any record in `INVALID_CANDIDATE`, before it considers
any claim; `scripts/neutrino_crossmatch.py` calls the guard on all seven and
records the refusals (7 of 7) in `results/neutrino_crossmatch.json`; and the
OMNI gate `neutrino_counterpart` is not required by any disposition, so a
coincidence can never be the reason for one.

## 14. A date treated as a clock

The survey product is named by day and carries no UTC start, so any "coincidence
with a neutrino" computed from it has a window of 24 hours or an unbounded one,
depending on how generous the reader is. Both are wrong in the same direction:
they manufacture a coincidence rate out of the filing convention. In a 4.68 deg
cone at 100 alerts/yr, a 24 h window carries a chance probability of 4.6e-4 -
about 300 times the survey's own 256 s integration - so substituting the date for
a timestamp buys a coincidence with a box nobody measured.

**Tell:** a crossmatch table whose window is "the observation" or "the day the
file is named after", or a `COINCIDENT` row next to a record that has no start
time. **Fix:** `seti_search/multimessenger.py` grades every record as `MISSING`,
`UTC_START` or `FULL`, never `PARTIAL` for a date; `assess` prices the 24 h
substitution in the reason string; `require_testable` refuses to proceed; and
`enforce` raises `MultimessengerNotTestable` if a coincidence is ever reported
for a record the gate called `NOT TESTABLE`. The OMNI gate
`neutrino_counterpart` reads the same verdict, so the table and the document
cannot disagree.

## Checklist before believing a result

1. Does any single number move by more than 2× when the noise reference changes?
2. Does a bright line vanish? If so, suspect the window before the sky.
3. Is the null degenerate (identical for every candidate)? Then it is not a null.
4. Is the transform fitted to the same quantity that defines the candidate? Then
   its in-sample agreement is circular.
5. Are the untestable gates marked UNTESTABLE rather than passed?
