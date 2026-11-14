# What the neutrino lane could detect

Framework: φ² = φ + 1. Attribution: Tyree Jones (tyofthestarz), prior art from December 2025.

This document exists so that `UNINFORMATIVE` is a measured statement and
not an excuse. It does not change any disposition: no number here can
promote a record, and the seven mask-manufactured cells stay excluded
whatever the flux limit turns out to be.

## Assumptions

- whole-sky alert rate: 100 per year (range [30, 100, 1000])
- effective area: 1e+04 m^2 (range [100.0, 1000.0, 10000.0])
- pivot energy: 1e+06 GeV
- flare: 1e+47 erg over 10 s

order-of-magnitude values: the alert rate is a published-alert rate, the effective area is a PeV track value for a cubic-kilometre detector, and the flare energy is a bright transient. Each scales the answer linearly and is stated so it can be disputed.

## 1. How often a random alert lands in the box

| window | 1.0 deg (alert-grade) | 3.0 deg (3 sigma of 1 deg) | 4.68 deg (3 sigma of the demo geometry) |
| --- | --- | --- | --- |
| 1 s | mu 2.41e-10, p 2.41e-10 | mu 2.17e-09, p 2.17e-09 | mu 5.28e-09, p 5.28e-09 |
| 256 s (survey window) | mu 6.18e-08, p 6.18e-08 | mu 5.56e-07, p 5.56e-07 | mu 1.35e-06, p 1.35e-06 |
| 1 h | mu 8.69e-07, p 8.69e-07 | mu 7.82e-06, p 7.82e-06 | mu 1.90e-05, p 1.90e-05 |
| 1 day | mu 2.08e-05, p 2.08e-05 | mu 1.88e-04, p 1.88e-04 | mu 4.56e-04, p 4.56e-04 |
| 30 days | mu 6.25e-04, p 6.25e-04 | mu 5.63e-03, p 5.61e-03 | mu 1.37e-02, p 1.36e-02 |

The background is tiny even for a one-degree cone and a day, because the
alert rate is only ~100 per year: 100/yr over 41,253 square degrees is one
event per 413 square degrees per year, and a one-degree cone is 3.1 square
degrees. **A coincidence is therefore easy and means little.** The
limiting quantity is not the background, it is the signal, which is why
the next two tables matter more.

## 2. Timing precision needed to keep the null under 1%

| alert rate (/yr) | cone | longest window |
| --- | --- | --- |
| 30 | 1.00 deg | 1.39e+08 s (1606.8 days) |
| 30 | 3.00 deg | 1.54e+07 s (178.6 days) |
| 30 | 4.68 deg | 6.34e+06 s (73.4 days) |
| 100 | 1.00 deg | 4.16e+07 s (482.0 days) |
| 100 | 3.00 deg | 4.63e+06 s (53.6 days) |
| 100 | 4.68 deg | 1.9e+06 s (22.0 days) |
| 1000 | 1.00 deg | 4.16e+06 s (48.2 days) |
| 1000 | 3.00 deg | 4.63e+05 s (5.4 days) |
| 1000 | 4.68 deg | 1.9e+05 s (2.2 days) |

The radio product has **no timestamp at all** - the files are named by
date - so its window is not 256 s and not a day, it is undefined. That is
the whole of the `UNINFORMATIVE` verdict: not a null result, an
unbounded box.

## 3. Sensitivity: flux and reach

| A_eff (m^2) | duration | fluence (cm^-2) | energy fluence (erg/cm^2) | horizon for 1e47 erg |
| --- | --- | --- | --- | --- |
| 1e+02 | 1 s | 2.30e-06 | 3.69e-03 | 0 Mpc |
| 1e+02 | 10 s | 2.30e-07 | 3.69e-04 | 2 Mpc |
| 1e+02 | 256 s | 8.99e-09 | 1.44e-05 | 8 Mpc |
| 1e+02 | 1 day | 2.67e-11 | 4.27e-08 | 140 Mpc |
| 1e+03 | 1 s | 2.30e-07 | 3.69e-04 | 2 Mpc |
| 1e+03 | 10 s | 2.30e-08 | 3.69e-05 | 5 Mpc |
| 1e+03 | 256 s | 8.99e-10 | 1.44e-06 | 24 Mpc |
| 1e+03 | 1 day | 2.67e-12 | 4.27e-09 | 442 Mpc |
| 1e+04 | 1 s | 2.30e-08 | 3.69e-05 | 5 Mpc |
| 1e+04 | 10 s | 2.30e-09 | 3.69e-06 | 15 Mpc |
| 1e+04 | 256 s | 8.99e-11 | 1.44e-07 | 76 Mpc |
| 1e+04 | 1 day | 2.67e-13 | 4.27e-10 | 1399 Mpc |

At 1e+04 m^2, the same detector is **256 times more sensitive** over a 256 s window than over a 1 s window, and a
correspondingly smaller share of the sky is excluded per second. This is
the trade the lane cannot avoid: the window that maximises the chance of
catching a flare is the window that maximises the background.

## 4. The null a real alert already carries

| event | FAR (/yr) | p within 1 s | p within 256 s | p within 1 day |
| --- | --- | --- | --- | --- |
| IceCube-230416A | 2.53 | 8.03e-08 | 2.06e-05 | 6.91e-03 |
| IC-143317_9424183 | 2.11 | 6.68e-08 | 1.71e-05 | 5.76e-03 |

These come from the ingested alerts' own false alarm rates, which is a
better null than the published rate used above: it already accounts
for the alert threshold. Note what it buys - even a day-long window
gives a small probability - and what it does not: a small probability
for a *background event*, which says nothing about whether the radio
feature is real.

## What would make the lane informative

1. A radio pointing: RA/Dec per record, from a beam-formed or targeted
   follow-up. Without it the cone does not exist.
2. A UTC start time and an integration length. That turns the window from
   undefined into 256 s, at which point the background is ~1e-7 and the
   background stops being the problem.
3. Enough effective area or enough flare energy to put one event in the
   box - the table above says how much, and the answer scales linearly in
   both.

Until (1) and (2) exist, the honest report is that the lane cannot be
run, and this document is the measurement of that sentence.

