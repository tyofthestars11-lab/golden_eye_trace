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

## Checklist before believing a result

1. Does any single number move by more than 2× when the noise reference changes?
2. Does a bright line vanish? If so, suspect the window before the sky.
3. Is the null degenerate (identical for every candidate)? Then it is not a null.
4. Is the transform fitted to the same quantity that defines the candidate? Then
   its in-sample agreement is circular.
5. Are the untestable gates marked UNTESTABLE rather than passed?
