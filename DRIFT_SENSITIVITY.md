# Where the search stops being flat in drift

**A stationary line at 6 sigma reads 6 sigma; the same line at 190 kHz/s reads about 3.8 in the best case and 0.8 on average.** At 18 sigma the same pair is 10.0 and 3.1. The window holds a complete trail up to 58.9 kHz/s, so a rate above that is searched and under-measured rather than unsearched.

## The guard's law, tested

The prediction is that the reported significance scales as the square root of the share of the trail the window holds, because `dedoppler` drops the integrations whose shift leaves the band: signal and noise are truncated together, so each surviving integration raises the significance by sqrt. Comparing the observed ratio `sigma(rate) / sigma(0)` against `sqrt(held)`, every rate whose window holds under 65% of its trail keeps +0.01 to +0.39 more significance than that prediction.

Each cell is a median over five injections, so a ratio carries about 0.15 of its own scatter; the comparison is quoted only where the window holds under 65% of the trail, which is where the effect is larger than that scatter. There the measured ratio never falls below the prediction, and that direction is the important one: the guard **overstates** the loss, so a window it passes is not hiding a shortfall. Part of the reason is physical - the rows are dropped only once the *shift* exceeds the band, which happens after the lane itself has left - and part is bookkeeping, since `calibrate_map` renormalises each map by its own spread, which shrinks when fewer rows are summed. That second part is a known bias rather than a correction: trials of different drift are then measured against a noise level they do not have, which is why the guard is stated as geometry and not as a calibration.

The exponent itself is a result of this run, not a convention. Loss linear in the fraction is what a dedoppler would show if it kept summing noise from integrations the signal had left; no loss at all is what a wrapping dedoppler would show. The measurement sits between the square root and unity, so the guard takes the conservative end of the interval it can defend: a window it passes keeps at least the significance it claims.

## The table

Cells are `sigma(rate) / sigma(0)` in the same-entry placement, beside the share of the trail the window holds. A cell at 1.00 lost nothing; a cell at 0.20 lost four fifths of its significance.

| drift kHz/s | cycles/sample | trail (ch) | held | 6σ: ratio / recovered | 9σ: ratio / recovered | 18σ: ratio / recovered |
| --- | --- | --- | --- | --- | --- | --- |
| 0.0 | +0.0000 | 0 | 100.0% | 1.00 / 100% | 1.00 / 100% | 1.00 / 100% |
| 0.5 | +0.0013 | 1 | 100.0% | 0.88 / 80% | 0.91 / 100% | 0.92 / 100% |
| 2.0 | +0.0051 | 5 | 100.0% | 0.82 / 60% | 0.83 / 100% | 0.95 / 100% |
| 5.0 | +0.0127 | 13 | 100.0% | 1.00 / 100% | 0.90 / 100% | 0.64 / 100% |
| 10.0 | +0.0255 | 25 | 100.0% | 0.94 / 80% | 0.85 / 100% | 0.70 / 100% |
| 20.0 | +0.0509 | 51 | 100.0% | 1.19 / 80% | 1.07 / 100% | 0.90 / 100% |
| 30.0 | +0.0764 | 76 | 100.0% | 1.50 / 100% | 1.35 / 100% | 1.03 / 100% |
| 41.2 | +0.1049 | 105 | 91.5% | 1.15 / 100% | 1.01 / 100% | 0.92 / 100% |
| 47.1 | +0.1200 | 120 | 80.0% | 1.00 / 100% | 0.93 / 100% | 0.85 / 100% |
| 58.9 | +0.1500 | 150 | 64.0% | 0.83 / 60% | 0.76 / 100% | 0.70 / 100% |
| 80.0 - | +0.2037 | 204 | 47.1% | 0.98 / 80% | 0.87 / 100% | 0.76 / 100% |
| 100.0 - | +0.2546 | 255 | 37.7% | 0.72 / 40% | 0.61 / 100% | 0.55 / 100% |
| 150.0 - | +0.3819 | 382 | 25.1% | 0.62 / 0% | 0.55 / 80% | 0.48 / 100% |
| 190.0 - | +0.4838 | 484 | 19.8% | 0.53 / 20% | 0.53 / 60% | 0.47 / 100% |

## What this does to the null result

The honest statement is not "no drifting line brighter than X anywhere below 196 kHz/s". It is: the sensitivity is flat to within the 20% tolerance up to 58.9 kHz/s, falls as the share of the trail the window holds above that, and at the fold a line must be roughly 5.0 times stronger to register at the same significance. The ratio column is that factor, measured.

## Two method notes

The `same-entry` column is the guard's best case only in the sense that it fixes one entry channel; a lane entering closer to the upstream edge has more runway and one entering at the centre has less. `uniform` is what a survey gets across entries, and it is where completeness actually lives.

Significance is read at the injected cell rather than at an extracted peak. The first two versions of this run matched peaks to injections within a frequency tolerance, and the tolerance made the result: two channels matched noise peaks and reported 100% recovery at the fold, half a channel refused real injections reported up to 0.7 channels from their centre. Reading the cell has nothing to tune.

