# Nonlinear event-time lane

**Verdict: NULL RESULT.** The linear result is unchanged and remains NULL RESULT. This lane runs beside it, on the same candidates, and overwrites nothing.

## What can be measured on this product

`spec_fits` carries intensity, not complex voltages, so there is no carrier phase to lock. The phase-locking value reported here is computed from a Hilbert analytic signal of the de-drifted track: it measures envelope structure and is **an envelope proxy, not carrier phase**, and it carries no decision weight. The decision statistic is a matched-filter excess score defined directly on intensity, `sqrt(N) * (mean - median) / sd` of the resampled track.

Shuffled-time nulls are reported in the CSV and are degenerate in the opposite direction from the usual one: a permutation keeps the sample values, so an amplitude statistic cannot move, and the null sits on the observation rather than below it. Measured worst case over all candidates and coordinates: 1.8e-15. That is, no coordinate can add coherence to the samples themselves. A coordinate changes the fitting basis; a warp that resamples changes the values by interpolation. Neither is evidence about the sky.

## Coordinates

| coordinate | definition |
| --- | --- |
| linear | `tau = t` |
| phase | `tau = cumulative |dphi| / (2 pi f_ref)` |
| chirp | `tau = integral (f - f0)/B dt = drift t^2 / (2B)` |
| accel | `tau with the survey-grid drift model removed` |
| energy | `tau = cumulative excess power / total excess` |
| phi | `tau_phi = frac(tau_event / PHI^2), feature only` |

Every transform is fitted on the **first half** of the track and evaluated on the **second half**. No coherence is measured on the samples that estimated the transform.

## The circularity this had to survive

`chirp` is fitted to the same linear drift that defines the candidate, so a constant-acceleration line is linearised *by construction* and its in-sample coherence carries no information. Two tests carry the weight. The split-half value applies the first half's drift to the second, so a burst whose drift changes fails even though the transform would have made it look coherent. And the gain over linear time is compared with the gain randomised monotone warps achieve, best-of-five against best-of-five, so transform and control have the same number of degrees of freedom.

## Per-candidate results

| date | MHz | raw drift Hz/s | linear score | best coordinate | best score | gain | warp gain p95 | beats control | split half | envelope PLV | repeats | disposition |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20240207 | 54.4530 | 8,194.7 | +6.76 | chirp | +7.33 | +0.57 | +1.41 | no | yes | 0.123 | 0 | UNRESOLVED |
| 20240207 | 53.7833 | 6,970.9 | +5.12 | phase | +9.34 | +4.22 | +9.51 | no | yes | 0.826 | 0 | UNRESOLVED |
| 20240207 | 55.9010 | 7,866.9 | +6.80 | accel | +7.39 | +0.59 | +5.06 | no | yes | 0.161 | 0 | UNRESOLVED |
| 20240322 | 47.2416 | -172.8 | +7.54 | chirp | +7.52 | -0.02 | +1.44 | no | yes | 0.767 | 1 | UNRESOLVED |
| 20240324 | 47.3143 | -504.0 | +5.09 | energy | +12.65 | +7.56 | +8.49 | no | yes | 0.904 | 1 | UNRESOLVED |
| 20240325 | 45.9843 | 41,190.2 | +6.93 | phase | +8.57 | +1.64 | +2.76 | no | yes | 0.134 | 0 | ASTROPHYSICAL |
| 20240325 | 46.9828 | 41,185.4 | +8.55 | chirp | +8.56 | +0.02 | +0.97 | no | yes | 0.319 | 0 | ASTROPHYSICAL |
| 20240325 | 47.9380 | 41,184.1 | +4.12 | chirp | +4.02 | -0.10 | +2.11 | no | yes | 0.085 | 0 | ASTROPHYSICAL |
| 20240325 | 48.9937 | 41,179.1 | -1.58 | phase | +7.04 | +8.62 | +10.36 | no | yes | 0.095 | 0 | ASTROPHYSICAL |
| 20240920 | 45.4206 | 195,015.7 | -11.83 | energy | -6.60 | +5.23 | +5.73 | no | no | 0.125 | 0 | INSTRUMENTAL ALIAS |
| 20240920 | 46.5890 | 194,896.8 | -13.03 | energy | -7.71 | +5.32 | +8.14 | no | no | 0.199 | 0 | INSTRUMENTAL ALIAS |
| 20240920 | 47.5358 | 194,905.6 | -13.89 | energy | -7.19 | +6.70 | +36.08 | no | no | 0.219 | 0 | INSTRUMENTAL ALIAS |

Scores are second-half only, in each coordinate. Envelope PLV is shown for completeness and is not used to rank anything.

## Controls

| date | MHz | warp p95 | warp envelope p95 | phase scrambled | reversed | shuffled channels | off band |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 20240207 | 54.4530 | +8.17 | 0.904 | +1.51 | +7.70 | +0.00 | +6.80 |
| 20240207 | 53.7833 | +14.63 | 0.968 | -0.73 | +8.72 | +7.82 | +5.15 |
| 20240207 | 55.9010 | +11.85 | 0.589 | +0.40 | +9.09 | +0.00 | +7.24 |
| 20240322 | 47.2416 | +8.98 | 0.833 | -1.14 | +10.63 | +6.94 | +8.94 |
| 20240324 | 47.3143 | +13.57 | 0.920 | -0.49 | +8.23 | +5.93 | +3.66 |
| 20240325 | 45.9843 | +9.69 | 0.304 | -1.73 | +3.87 | +0.00 | +7.42 |
| 20240325 | 46.9828 | +9.51 | 0.880 | -0.92 | +12.00 | +0.00 | +6.14 |
| 20240325 | 47.9380 | +6.23 | 0.438 | -2.06 | +4.69 | +0.00 | +6.69 |
| 20240325 | 48.9937 | +8.78 | 0.277 | -0.94 | +7.97 | +0.00 | +8.24 |
| 20240920 | 45.4206 | -6.10 | 0.963 | +3.24 | +5.65 | +0.00 | +20.93 |
| 20240920 | 46.5890 | -4.89 | 0.967 | +4.82 | +4.96 | +0.00 | +22.40 |
| 20240920 | 47.5358 | +22.19 | 0.969 | +4.14 | +4.66 | +0.00 | +20.89 |

## The 20240325 group, solar template first

4 lanes with a drift spread of 11.1 Hz/s, which is 0.028 of one resolution element, so the lanes are parallel at this resolution. A common drift across lanes separated by about 1 MHz is a burst crossing the band, not independent transmitters. A coronal exciter at 125 km/s gives 42,632 Hz/s at these frequencies against the measured 41,179-41,190 Hz/s. Parallel: yes.

## The 20240920 group, exact alias mapping

| injected cycles/sample | reported | fold |
| --- | --- | --- |
| 0.1000 | +0.1000 | +0.1000 |
| 0.3000 | +0.3000 | +0.3000 |
| 0.4965 | +0.4965 | +0.4965 |
| 0.5035 | -0.4965 | -0.4965 |
| 0.6000 | -0.4000 | -0.4000 |
| 0.9000 | -0.1000 | -0.1000 |

0.6 is recorded as -0.4: the fold changes the sign as well as the magnitude, which is why the earlier linear search could not report it at 0.6 and why a reported 0.4965 is equally consistent with a true -0.5035. The three 20240920 lines share one folded rate and one sign, the signature of a common sampling origin rather than of three sources.

## Injection recovery

| injected | reported | seed | linear score | best score | gain | shuffled p95 | detected |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.0000 | +0.0000 | 11 | +6.28 | +6.28 | +0.00 | +6.29 | yes |
| 0.0000 | +0.0000 | 12 | +5.99 | +5.99 | +0.00 | +5.17 | yes |
| 0.0200 | +0.0200 | 11 | +5.89 | +6.76 | +0.87 | +5.82 | yes |
| 0.0200 | +0.0200 | 12 | +5.15 | +5.80 | +0.65 | +5.85 | yes |
| 0.3820 | +0.3820 | 11 | +6.61 | +6.70 | +0.09 | +6.86 | yes |
| 0.3820 | +0.3820 | 12 | +6.33 | +6.33 | +0.00 | +6.20 | yes |
| 0.4500 | +0.4500 | 11 | +4.88 | +5.62 | +0.74 | +5.99 | yes |
| 0.4500 | +0.4500 | 12 | +4.77 | +5.47 | +0.70 | +5.43 | yes |
| 0.6000 | -0.4000 | 11 | +5.59 | +5.59 | +0.00 | +5.06 | yes |
| 0.6000 | -0.4000 | 12 | +5.47 | +6.24 | +0.77 | +5.63 | yes |

Synthetic, unit-mean exponential noise, amplitude 3. `detected` means the best score exceeded 5, the search threshold; the shuffled p95 is shown but, as above, is degenerate for an amplitude statistic. These are coherence numbers, not a sensitivity curve.

## Confirmation status

| date | MHz | nonlinear coherent | absent from controls | repeat | later observation | off source | barycentric | second telescope |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20240207 | 54.4530 | no | no | no | no | UNTESTABLE | UNTESTABLE (no WCS, no sky position) | UNTESTABLE |
| 20240207 | 53.7833 | no | yes | no | no | UNTESTABLE | UNTESTABLE (no WCS, no sky position) | UNTESTABLE |
| 20240207 | 55.9010 | no | no | no | no | UNTESTABLE | UNTESTABLE (no WCS, no sky position) | UNTESTABLE |
| 20240322 | 47.2416 | no | no | yes | yes | UNTESTABLE | UNTESTABLE (no WCS, no sky position) | UNTESTABLE |
| 20240324 | 47.3143 | no | yes | yes | yes | UNTESTABLE | UNTESTABLE (no WCS, no sky position) | UNTESTABLE |
| 20240325 | 45.9843 | no | yes | no | no | UNTESTABLE | UNTESTABLE (no WCS, no sky position) | UNTESTABLE |
| 20240325 | 46.9828 | no | no | no | no | UNTESTABLE | UNTESTABLE (no WCS, no sky position) | UNTESTABLE |
| 20240325 | 47.9380 | no | no | no | no | UNTESTABLE | UNTESTABLE (no WCS, no sky position) | UNTESTABLE |
| 20240325 | 48.9937 | no | no | no | no | UNTESTABLE | UNTESTABLE (no WCS, no sky position) | UNTESTABLE |
| 20240920 | 45.4206 | no | no | no | no | UNTESTABLE | UNTESTABLE (no WCS, no sky position) | UNTESTABLE |
| 20240920 | 46.5890 | no | no | no | no | UNTESTABLE | UNTESTABLE (no WCS, no sky position) | UNTESTABLE |
| 20240920 | 47.5358 | no | no | no | no | UNTESTABLE | UNTESTABLE (no WCS, no sky position) | UNTESTABLE |

Three gates are UNTESTABLE rather than failed: `spec_fits` is one disk-integrated beam per day with no off-source pointing, the product carries no WCS so a barycentric correction is undefined, and there is no second instrument in this archive.

## Dispositions

- **20240207 54.453 MHz** - UNRESOLVED; band: inside the 54-60 MHz TV channel 2 allocation (55.25 video, 59.75 audio); best coordinate `chirp` at score +7.33 against a linear score of +6.76 and a warp-gain control of +1.41; no astrophysical or instrumental explanation established.
- **20240207 53.783 MHz** - UNRESOLVED; band: inside the 50-54 MHz amateur 6 m allocation; best coordinate `phase` at score +9.34 against a linear score of +5.12 and a warp-gain control of +9.51; no astrophysical or instrumental explanation established.
- **20240207 55.901 MHz** - UNRESOLVED; band: inside the 54-60 MHz TV channel 2 allocation (55.25 video, 59.75 audio); best coordinate `accel` at score +7.39 against a linear score of +6.80 and a warp-gain control of +5.06; no astrophysical or instrumental explanation established.
- **20240322 47.242 MHz** - UNRESOLVED; band: below the 6 m band; no allocated broadcast service here; best coordinate `chirp` at score +7.52 against a linear score of +7.54 and a warp-gain control of +1.44; no astrophysical or instrumental explanation established.
- **20240324 47.314 MHz** - UNRESOLVED; band: below the 6 m band; no allocated broadcast service here; best coordinate `energy` at score +12.65 against a linear score of +5.09 and a warp-gain control of +8.49; no astrophysical or instrumental explanation established.
- **20240325 45.984 MHz** - ASTROPHYSICAL; band: below the 6 m band; no allocated broadcast service here; best coordinate `phase` at score +8.57 against a linear score of +6.93 and a warp-gain control of +2.76; solar burst: parallel lanes sharing one drift, and the drift matches a 125 km/s coronal exciter.
- **20240325 46.983 MHz** - ASTROPHYSICAL; band: below the 6 m band; no allocated broadcast service here; best coordinate `chirp` at score +8.56 against a linear score of +8.55 and a warp-gain control of +0.97; solar burst: parallel lanes sharing one drift, and the drift matches a 125 km/s coronal exciter.
- **20240325 47.938 MHz** - ASTROPHYSICAL; band: below the 6 m band; no allocated broadcast service here; best coordinate `chirp` at score +4.02 against a linear score of +4.12 and a warp-gain control of +2.11; solar burst: parallel lanes sharing one drift, and the drift matches a 125 km/s coronal exciter.
- **20240325 48.994 MHz** - ASTROPHYSICAL; band: below the 6 m band; no allocated broadcast service here; best coordinate `phase` at score +7.04 against a linear score of -1.58 and a warp-gain control of +10.36; solar burst: parallel lanes sharing one drift, and the drift matches a 125 km/s coronal exciter.
- **20240920 45.421 MHz** - INSTRUMENTAL ALIAS; band: below the 6 m band; no allocated broadcast service here; best coordinate `energy` at score -6.60 against a linear score of -11.83 and a warp-gain control of +5.73; rate 0.4965 cycles/sample, 0.0035 below Nyquist; the sampler cannot distinguish it from -0.5035.
- **20240920 46.589 MHz** - INSTRUMENTAL ALIAS; band: below the 6 m band; no allocated broadcast service here; best coordinate `energy` at score -7.71 against a linear score of -13.03 and a warp-gain control of +8.14; rate 0.4965 cycles/sample, 0.0035 below Nyquist; the sampler cannot distinguish it from -0.5035.
- **20240920 47.536 MHz** - INSTRUMENTAL ALIAS; band: below the 6 m band; no allocated broadcast service here; best coordinate `energy` at score -7.19 against a linear score of -13.89 and a warp-gain control of +36.08; rate 0.4965 cycles/sample, 0.0035 below Nyquist; the sampler cannot distinguish it from -0.5035.

## What the second lane adds

A coordinate in which a constant-acceleration line is linear by construction is a statement about the transform, not about the sky. What is a property of the sources is whether the drift persists between halves, and the split-half score measures exactly that. No candidate is both coherent beyond the flexibility-matched control and confirmed, because the gates that would carry a technosignature verdict are untestable with this product.
