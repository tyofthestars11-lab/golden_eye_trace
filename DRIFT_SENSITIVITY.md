# Where the search stops being flat in drift

**A stationary line at 6 sigma reads 6 sigma; the same line at 190 kHz/s reads about 4.2 in the best case and 2.2 on average.** At 15 sigma the same pair is 12.1 and 4.1. The window holds a complete trail up to 58.9 kHz/s, so a rate above that is searched and under-measured rather than unsearched.

## The guard's law, tested

The prediction is that the reported significance scales as the square root of the share of the trail the window holds, because `dedoppler` drops the integrations whose shift leaves the band: signal and noise are truncated together, so each surviving integration raises the significance by sqrt. Comparing the observed ratio `sigma(rate) / sigma(0)` against `sqrt(held)`, every rate whose window holds under 65% of its trail keeps +0.09 to +0.33 more significance than that prediction.

Each cell is a median over five injections, so a ratio carries about 0.15 of its own scatter; the comparison is quoted only where the window holds under 65% of the trail, which is where the effect is larger than that scatter. There the measured ratio never falls below the prediction, and that direction is the important one: the guard **overstates** the loss, so a window it passes is not hiding a shortfall. Part of the reason is physical - the rows are dropped only once the *shift* exceeds the band, which happens after the lane itself has left - and part is bookkeeping, since `calibrate_map` renormalises each map by its own spread, which shrinks when fewer rows are summed. That second part is a known bias rather than a correction: trials of different drift are then measured against a noise level they do not have, which is why the guard is stated as geometry and not as a calibration.

The exponent itself is a result of this run, not a convention. Loss linear in the fraction is what a dedoppler would show if it kept summing noise from integrations the signal had left; no loss at all is what a wrapping dedoppler would show. The measurement sits between the square root and unity, so the guard takes the conservative end of the interval it can defend: a window it passes keeps at least the significance it claims.

## The table

Cells are `sigma(rate) / sigma(0)` in the same-entry placement, beside the share of the trail the window holds. A cell at 1.00 lost nothing; a cell at 0.20 lost four fifths of its significance.

| drift kHz/s | cycles/sample | trail (ch) | held | 6σ: ratio / recovered | 9σ: ratio / recovered | 15σ: ratio / recovered |
| --- | --- | --- | --- | --- | --- | --- |
| 0.0 | +0.0000 | 0 | 100.0% | 1.00 / 100% | 1.00 / 100% | 1.00 / 100% |
| 5.0 | +0.0127 | 13 | 100.0% | 0.99 / 60% | 1.00 / 100% | 0.93 / 100% |
| 20.0 | +0.0509 | 51 | 100.0% | 1.28 / 100% | 1.26 / 100% | 1.18 / 100% |
| 41.2 | +0.1049 | 105 | 91.5% | 1.27 / 100% | 1.21 / 100% | 1.01 / 100% |
| 58.9 | +0.1500 | 150 | 64.0% | 0.78 / 40% | 0.84 / 100% | 0.88 / 100% |
| 100.0 - | +0.2546 | 255 | 37.7% | 0.77 / 40% | 0.75 / 100% | 0.69 / 100% |
| 150.0 - | +0.3819 | 382 | 25.1% | 0.73 / 0% | 0.69 / 80% | 0.66 / 100% |
| 190.0 - | +0.4838 | 484 | 19.8% | 0.62 / 40% | 0.67 / 60% | 0.72 / 100% |

## What this does to the null result

The honest statement is not "no drifting line brighter than X anywhere below 196 kHz/s". It is: the sensitivity is flat to within the 20% tolerance up to 58.9 kHz/s, falls as the share of the trail the window holds above that, and at the fold a line must be roughly 5.0 times stronger to register at the same significance. The ratio column is that factor, measured.

## Two method notes

The `same-entry` column is the guard's best case only in the sense that it fixes one entry channel; a lane entering closer to the upstream edge has more runway and one entering at the centre has less. `uniform` is what a survey gets across entries, and it is where completeness actually lives.

Significance is read at the injected cell rather than at an extracted peak. The first two versions of this run matched peaks to injections within a frequency tolerance, and the tolerance made the result: two channels matched noise peaks and reported 100% recovery at the fold, half a channel refused real injections reported up to 0.7 channels from their centre. Reading the cell has nothing to tune.

