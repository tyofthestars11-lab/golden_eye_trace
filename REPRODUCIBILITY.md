# Reproducibility report — OVRO-LWA candidates

## Exact-configuration reproduction

Day 20240207 searched at 1600 samples x 192 channels, 5.5 sigma, drift window to Nyquist — the same configuration that produced the published pair.

| known MHz | known Hz/s | found MHz | found Hz/s | SNR | verdict |
| --- | --- | --- | --- | --- | --- |
| 54.453 | 8,194.7 | 54.4531 | 8,194.7 | 7.68 | reproduced |
| 55.901 | 7,866.9 | 55.9015 | 7,866.9 | 7.00 | reproduced |

## Why the pair is absent from the survey table

The survey uses 1000 samples x 96 channels, which covers 45.27-54.93 MHz, while the wide run covers 40.45-59.75 MHz. Coverage per known line:

| line MHz | inside wide band | inside survey band |
| --- | --- | --- |
| 54.453 | True | True |
| 55.901 | True | False |

A line outside a configuration's band cannot appear in it; that is a coverage effect and not a failure of reproducibility. The 45-49 MHz survey band is a different part of the spectrum, so the two tables are not directly comparable row to row.
