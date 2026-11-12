# The survey re-run with a local noise reference

**9 of 10 candidates survive.**

## What changed

The survey's statistic is (smoothed - median)/MAD, with one median and one MAD for the whole searched band, rescaled so the map's bulk is one sigma wide. That is a single noise scale for 9.66 MHz. This run changes only the scale: same window, same samples, same dedoppler, same boxcar, same calibration, same 5.5σ threshold, but the scale is measured per channel from a ±16 channel neighbourhood across every drift trial, with one round of clipping so a bright line cannot set the noise it is measured against.

Significance is read at the map cell the survey reported, not from a re-extracted peak. Dividing each column by its own scale moves the map's maximum, so a line can be strong and still not be the top hit; reading its own cell asks the question that matters.

| date | MHz | drift Hz/s | survey σ | at-cell σ | local scale | local σ | outcome |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 20240207 | 53.7833 | 6,970.9 | 7.16 | 7.16 | 0.74 | 9.66 | survives |
| 20240322 | 47.2416 | -172.8 | 6.67 | 6.67 | 1.29 | 5.16 | below threshold |
| 20240324 | 47.3143 | -504.0 | 11.67 | 11.67 | 1.30 | 8.97 | survives |
| 20240325 | 45.9843 | 41,190.2 | 953.98 | 954.00 | 31.18 | 30.60 | survives |
| 20240325 | 46.9828 | 41,185.4 | 604.03 | 604.05 | 27.69 | 21.81 | survives |
| 20240325 | 47.9380 | 41,184.1 | 447.08 | 447.09 | 19.54 | 22.88 | survives |
| 20240325 | 48.9937 | 41,179.1 | 171.23 | 171.23 | 12.26 | 13.97 | survives |
| 20240920 | 45.4206 | 195,015.7 | 8.66 | 8.66 | 0.59 | 14.74 | survives |
| 20240920 | 46.5890 | 194,896.8 | 7.91 | 7.91 | 0.34 | 23.35 | survives |
| 20240920 | 47.5358 | 194,905.6 | 6.49 | 6.67 | 0.30 | 22.05 | survives |

## The one that does not survive

**20240322 47.242 MHz** falls from 6.67σ to 5.16σ: the noise in its own ±16 channel neighbourhood is 1.29 of the band-wide estimate, so it does not clear 5.5σ. This is the line that RFI_PAIR.md also finds absent over the full 717 s integration, so the two tests agree that it is not a persistent feature.

## Where the local scale is not a noise estimate

For the 20240325 burst group the local scale comes out at 12 to 31 times the band-wide value, because the ±16 channel neighbourhood contains the bright lanes themselves and their smear cones. Their 14-31σ values are therefore lower bounds on a reference that is partly made of the source, not estimates of its significance. The astrophysical classification does not rest on the number: it rests on four lanes parallel to 0.028 of a resolution element with a drift matching a 125 km/s coronal exciter.

The 20240920 group moves the other way. Its neighbourhood is three times quieter than the band-wide average, so the significance rises to 14-23σ from 6.5-8.7σ. That does not make it astrophysical: those three lines sit at 0.4965 cycles per sample, within 0.0035 of Nyquist, and an injected 0.5035 appears at -0.4965. Their strength is a statement about the noise, not about their origin.

## Reading

The survey's significances were neither uniformly inflated nor uniformly conservative: the local scale ranges from 0.30 to 31.18 of the band-wide value across these ten candidates. What the re-run changes is the confidence attached to each line; what it does not change is a single classification. The bright group stays astrophysical, the near-Nyquist trio stays instrumental, the 53.783 MHz line stays an unexplained single line, and the 47.242 MHz line drops below threshold and joins the pair that fails its own repeat test.

One correction to an earlier attempt, recorded because it mattered: a first version of this test searched narrow 32-channel windows and returned 0 of 10, including a 954σ burst group. That was a broken test, not a result. The 20240325 lanes drift about 105 channels over the 1000-sample integration, so a 32-channel window truncates them mid-trail. Keeping the survey's 96-channel window and changing only the noise scale avoids that entirely.
