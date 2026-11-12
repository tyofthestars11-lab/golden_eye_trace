# Confirmation report

**Headline: NULL RESULT** — no candidate passes the confirmation gates, so no signal is claimed as evidence of intelligent life.

## Gate matrix

| date | MHz | repeat_testable | repeats_same_position | off_source_available | independent_telescope | barycentric_reducible | physically_possible_orbital | explained_by_rfi |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20240207 | 54.4530 | pass | fail | fail | fail | fail | fail | pass |
| 20240207 | 53.7833 | pass | fail | fail | fail | fail | fail | fail |
| 20240207 | 55.9010 | fail | fail | fail | fail | fail | fail | pass |
| 20240322 | 47.2416 | pass | pass | fail | fail | fail | fail | fail |
| 20240324 | 47.3143 | pass | pass | fail | fail | fail | fail | fail |
| 20240325 | 45.9843 | pass | fail | fail | fail | fail | fail | fail |
| 20240325 | 46.9828 | pass | fail | fail | fail | fail | fail | fail |
| 20240325 | 47.9380 | pass | fail | fail | fail | fail | fail | fail |
| 20240325 | 48.9937 | pass | fail | fail | fail | fail | fail | fail |
| 20240920 | 45.4206 | pass | fail | fail | fail | fail | fail | fail |
| 20240920 | 46.5890 | pass | fail | fail | fail | fail | fail | fail |
| 20240920 | 47.5358 | pass | fail | fail | fail | fail | fail | fail |

Three gates cannot be attempted with this product at all, which is the central result: `off_source_available` and `independent_telescope` are structurally false for one disk-integrated beam from one array, and `barycentric_reducible` is false for every candidate whose drift exceeds the few Hz/s that Earth's motion can supply.

## Injection-recovery

| rate cycles/sample | recovered |
| --- | --- |
| 0.0000 | 1.00 |
| 0.0200 | 1.00 |
| 0.3820 | 1.00 |
| 0.6000 | 0.00 |

Rates 0 and 0.02 test stationary and slow lines, theta tests a lock-rate line, and 0.6 folds past Nyquist to test aliasing handling. Amplitude is fixed, so these are completeness at one signal strength, not a curve.

Noise-only runs: 3 trials at the same threshold, 0 with any hit.

## Dispositions

- **20240207 54.453 MHz, 8,195 Hz/s, SNR 7.7** — UNRESOLVED (in a TV band, but no repeat across days); inside the 54-60 MHz TV channel 2 allocation (55.25 video, 59.75 audio); nearest phi integer k = 1 at distance 0.945.
- **20240207 53.783 MHz, 6,971 Hz/s, SNR 7.2** — UNRESOLVED; inside the 50-54 MHz amateur 6 m allocation; nearest phi integer k = 1 at distance 0.954.
- **20240207 55.901 MHz, 7,867 Hz/s, SNR 7.0** — UNRESOLVED (in a TV band; repeat untestable, outside survey coverage); inside the 54-60 MHz TV channel 2 allocation (55.25 video, 59.75 audio); nearest phi integer k = 1 at distance 0.948.
- **20240322 47.242 MHz, -173 Hz/s, SNR 6.7** — UNRESOLVED; below the 6 m band; no allocated broadcast service here; nearest phi integer k = 1 at distance 0.999.
- **20240324 47.314 MHz, -504 Hz/s, SNR 11.7** — UNRESOLVED; below the 6 m band; no allocated broadcast service here; nearest phi integer k = 1 at distance 0.997.
- **20240325 45.984 MHz, 41,190 Hz/s, SNR 954.0** — ASTROPHYSICAL (bright burst group); below the 6 m band; no allocated broadcast service here; nearest phi integer k = 1 at distance 0.725.
- **20240325 46.983 MHz, 41,185 Hz/s, SNR 604.0** — ASTROPHYSICAL (bright burst group); below the 6 m band; no allocated broadcast service here; nearest phi integer k = 1 at distance 0.725.
- **20240325 47.938 MHz, 41,184 Hz/s, SNR 447.1** — ASTROPHYSICAL (bright burst group); below the 6 m band; no allocated broadcast service here; nearest phi integer k = 1 at distance 0.725.
- **20240325 48.994 MHz, 41,179 Hz/s, SNR 171.2** — ASTROPHYSICAL (bright burst group); below the 6 m band; no allocated broadcast service here; nearest phi integer k = 1 at distance 0.726.
- **20240920 45.421 MHz, 195,016 Hz/s, SNR 8.7** — INSTRUMENTAL ALIASING (rate at Nyquist); below the 6 m band; no allocated broadcast service here; nearest phi integer k = 1 at distance 0.300.
- **20240920 46.589 MHz, 194,897 Hz/s, SNR 7.9** — INSTRUMENTAL ALIASING (rate at Nyquist); below the 6 m band; no allocated broadcast service here; nearest phi integer k = 1 at distance 0.299.
- **20240920 47.536 MHz, 194,906 Hz/s, SNR 6.5** — INSTRUMENTAL ALIASING (rate at Nyquist); below the 6 m band; no allocated broadcast service here; nearest phi integer k = 1 at distance 0.299.

## What would change this

An off-source observation, a second telescope, and a candidate whose drift survives barycentric correction. None of the three exists in this data set, so the honest statement is that the archive contains drifting lines and none of them is a technosignature candidate.
