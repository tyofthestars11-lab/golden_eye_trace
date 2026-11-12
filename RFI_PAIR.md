# Is the 47.242 / 47.314 MHz pair terrestrial RFI?

**Verdict: NOT A REPEAT (configuration-dependent transients; RFI unconfirmed).**

## Why the pairing needed re-testing

The survey paired these lines on 0.3 MHz in frequency and 800 Hz/s in drift, at a drift resolution of 392.8 Hz/s over a 256 s window. Those tolerances admit almost any two lines in the band. Three measurements replace them: a full-integration dedoppler search per day, a scan of successive 256 s windows through each day, and the Doppler bounds of the terrestrial mechanisms that could produce either line.

## Full integration

| date | integration s | drift resolution Hz/s | line found | frequency MHz | drift Hz/s | SNR | hits in window |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 20240322 | 717 | 140.2 | no |  |  |  | 2 |
| 20240324 | 1017 | 98.9 | no |  |  |  | 1 |

## Window dependence

The survey's 96-channel window and this 24-channel window use the same samples, the same dedoppler and the same drift range. The survey reports 5.5σ threshold detections here; the narrow window finds nothing at 4σ. A single 3-channel boxcar sums three channels, so a per-channel peak of p sigma becomes p*sqrt(3) after the boxcar:

| date | survey SNR | peak sigma, own band | peak x sqrt(3) | agrees |
| --- | --- | --- | --- | --- |
| 20240322 | 6.67 | 3.95 | 6.85 | yes |
| 20240324 | 11.67 | 0.83 | 1.43 | no |

For 20240322 the boxcar accounts for the survey number almost exactly; for 20240324 it does not, by a factor of eight. A line whose significance changes with the width of the window used to estimate the noise is not a robust property of the data, and that is the case here.

## Persistence through the day

| date | windows scanned | windows with a line near the candidate | best per window |
| --- | --- | --- | --- |
| 20240322 | 3 | 0 | 0 s: -; 256 s: -; 512 s: - |
| 20240324 | 4 | 0 | 0 s: -; 256 s: -; 512 s: -; 768 s: - |

## Doppler bounds at 47.28 MHz

| mechanism | speed m/s | accel m/s² | frequency offset Hz | drift Hz/s |
| --- | --- | --- | --- | --- |
| ionosphere, quiet, 0.01 TECU/min |  |  | 0.00 | 0.000 |
| ionosphere, active, 0.1 TECU/min |  |  | 0.05 | 0.000 |
| ionosphere, storm, 1 TECU/min |  |  | 0.47 | 0.003 |
| ionosphere, severe, 10 TECU/min |  |  | 4.74 | 0.026 |
| aircraft, 250 m/s, 1 g | 250 | 9.8 | 39.43 | 1.547 |
| aircraft, 250 m/s, 5 g | 250 | 49.0 | 39.43 | 7.735 |
| LEO satellite, 7.7 km/s, 9 m/s2 | 7700 | 9.0 | 1214.31 | 1.419 |
| balloon, 30 m/s | 30 | 0.5 | 4.73 | 0.079 |

The ionospheric row uses df = 40.3 dTEC/dt / (c f) with TEC in electrons per square metre. A quiet ionosphere moves a 47 MHz carrier by 0.05 Hz and a severe storm by 4.7 Hz: the whole ionospheric range is under 5 Hz and under 0.03 Hz/s. Aircraft, satellites and balloons reach 1,214 Hz of offset, because 7.7 km/s of orbital speed is real, but every one of them is bounded at a few Hz/s of drift by its own acceleration.

## What the measurements say

- **20240322** full integration: no line within 0.3 MHz at 4.0σ over 717 s, so the survey's 6.67σ detection does not persist. Its -172.8 Hz/s, taken at face value, would need 1,097 m/s² (112 g).
- **20240324** full integration: no line within 0.3 MHz at 4.0σ over 1017 s, so the survey's 11.67σ detection does not persist. Its -504.0 Hz/s, taken at face value, would need 3,193 m/s² (326 g).

## Direct answer

The premise fails before the classification can be made. These two lines cannot be tested for terrestrial RFI by frequency and drift stability, because neither survives re-measurement: each is present in one 256 s window and absent from a full-integration search and from every other window of its own day. Their survey significance is not reproducible at a different window width. So the pair is not a repeat, not a stable carrier, and no mechanism - terrestrial or otherwise - is needed to explain it. If it is terrestrial, then its quoted drifts are artifacts of the same noise estimate, since 112 g and 326 g are reached by no aircraft, satellite or balloon in the bounds table.

## Reading

Two things follow, and neither is the answer the tolerance-based pairing implied.

First, the pair is not a repeat. The two detections come from one 256 s window each, on consecutive-but-one days, and neither survives a full-integration search at 4σ even though the integration is three to four times longer, which for a coherent line scales the significance by 1.7 and would put it near 11σ and 19σ. A line that is present for four minutes and absent for the rest of the day is a transient, not a carrier that repeats.

Second, the drifts cannot be used to classify anything on their own. At the full-integration resolution the survey's -172.8 and -504.0 Hz/s are not significant; and at face value they would require 1,097 and 3,193 m/s², which is beyond every terrestrial bound in the table (the largest is 7.7 Hz/s). So the drifts neither support a terrestrial origin nor exclude one.

The honest classification is therefore not TERRESTRIAL RFI but **not a repeat**. Transient single-window detections in a crowded part of the band are most plausibly terrestrial or instrumental, but this data set cannot say which, and the previous disposition's `repeats_same_position: yes` was an artifact of the 0.3 MHz and 800 Hz/s tolerances.

## What would settle it

A re-observation at kilohertz frequency resolution on several days, with the array pointed away for an off-source check. The 100.5 kHz channelisation fixes any frequency comparison at 0.7 channel, and a single 0.5 h track per day cannot distinguish a four-minute terrestrial transient from a four-minute astrophysical one.
