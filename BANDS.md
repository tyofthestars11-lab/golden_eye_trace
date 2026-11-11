# Pulling a second band, reading it, and registering the golden axis

The pre-registered test in `prereg/phi_2band.json` needs two bands at least a
decade apart, each with at least three resolved drifting lines. This note records
which open archives were reachable, what each actually delivers once read, and
how the registration was revised around the golden axis. Every number here was
produced by a script in this repository on a file fetched during the session.

```bash
python3 scripts/fetch_band.py --list
python3 scripts/fetch_band.py --get lwa-drx --with-payload
python3 scripts/fetch_band.py --gate
python3 scripts/presto_dat.py --dat <f>.dat --inf <f>.inf --inspect
python3 scripts/psrfits.py data/bands/lwa_drspec.fits --info
python3 scripts/psrfits.py data/bands/lwa_drspec.fits --to-fil data/lwa_drspec.fil
python3 scripts/phi_prereg_amend.py
python3 scripts/phi_twoband.py --band-a data/stare2_ovro.fil --band-b data/lwa_drspec.fil
```

## The survey

| source | band | centre | decades from 1405 MHz | status here |
| --- | --- | --- | --- | --- |
| STARE2 OVRO (local) | 1280–1530 MHz | 1405 | 0.00 | drift-blind: 122 kHz channels, 0.4 s |
| **LWA1 DRX, B0329+54** | **40.6–60.2 MHz** | **50.4** | **1.45** | **read: 1024 channels, 910 s** |
| LWA1 DRX, tuning 1 (time series) | 32.2–51.8 MHz | 42.0 | 1.52 | `.dat` is a power series, no axis |
| LOFAR HBA | 110–190 MHz | 150 | 0.97 | no connection (000); LTA login |
| MWA | 72–300 MHz | 186 | 0.88 | no connection (000); Pawsey login |
| BL Open Data (GBT/Parkes X band) | up to 11.2 GHz | 6150 | 0.64 | 200 frontend, `/api/*` 404 here |
| Zenodo 22848961 (Parkes) | 1360–1680 MHz | 1520 | 0.03 | FITS despite `.sf`; rejected by the gate |
| FAST | 1000–1500 MHz | 1250 | 0.05 | 200 landing page; login |

## What the two LWA products actually are

`presto_dat.py` decides from the numbers, not the extension. The `DM0.00` `.dat`
is 4,034,560 float32 samples, exactly the `.inf` bin count, with lag-1
autocorrelation +0.475 and channel means flat to 3 × 10⁻⁴ under the waterfall
reading: a **single dedispersed power series**, no frequency axis.

`psrfits.py` then reads the 36 MB `drspec` product, which is a genuine PSRFITS
`SUBINT` table: 1024 channels of **19.14 kHz** across 40.6–60.18 MHz, 709 rows,
each row `TSUBINT` = 1.284 s, so 910 s of integration and a drift resolution of
**21 Hz/s** — 15,000 times finer than STARE2. Two details cost a debugging pass
each and are worth recording: FITS tables are **big-endian**, and `TDIM` lists
dimensions fastest-varying first, so the numpy shape is that list reversed. The
row's 32 stored values are phase bins; a row-mean spectrum advances by `TSUBINT`,
not by one bin.

That gives the band axis and the time axis, at 1.45 decades. What it does not
give is drifting lines: the 3σ search finds **1 candidate and 0 with resolved
drift**. At 50 MHz the band is full of stationary terrestrial carriers, which is
a statement about the sky and the site rather than about the reader.

## Band A: why a longer cadence cannot fix it

Drift resolution is `Δν / T`, so with channels fixed at 122.07 kHz:

| drift | contiguous scan needed | channel width needed in a 300 s scan |
| --- | --- | --- |
| 0.1 Hz/s | 14.1 d | 30 Hz |
| 1 Hz/s | 1.41 d | 300 Hz |
| 10 Hz/s | 3.39 h | 3.0 kHz |
| 100 Hz/s | 20.3 min | 30 kHz |
| 1000 Hz/s | 2.03 min | 300 kHz |

STARE2 records 0.4 s, so resolving a 10 Hz/s drift would need 12,207 s at that
channel width. **The channel width is the limit, not the cadence**, which is why
revision 2 replaces the high band instead of extending it.

## The golden axis in the registration

The axis is used twice, both times as arithmetic that can fail. Its pitch,
17.0324°, read as a log-log slope gives tan μ = 0.30635, i.e. **1.46587 rungs per
decade** against **4.78497** for a flat rung rate. And it fixes where the bands
should sit: at an integer rung offset the predicted drift ratio is a pure power
of φ.

Between the reachable bands the offset is **6.915 rungs**. The axis would place
them 7 apart; the archive falls short by **0.085 rungs**. At that offset the two
statements predict, for drift in the high band over the low band:

| statement | exponent | predicted ratio |
| --- | --- | --- |
| flat rung rate | 1.0000 | **27.877** |
| golden axis | 0.30635 | **2.772** |

A slope known to better than **1.599 rungs/decade** separates them at 3σ, and
the two ratios above are what the measurement would have to choose between. The
deficit is recorded rather than absorbed: the prediction is evaluated at the
offset the data have, not the offset that would be convenient.

`scripts/phi_prereg_amend.py` writes `prereg/phi_2band_r2.json`, which keeps the
statistic, the two predictions, the 5% decision rule and the seed byte-for-byte
and records why the high band changed.

## Status

- Band axis: **satisfied and read** — 1.45 decades, 1024 × 19.14 kHz, 910 s.
- Drifting lines in that band: **0**, so the test remains **NOT MEASURED**, now
  for the only reason left: there is nothing drifting to measure.
- High band: needs ≤ 3 kHz channels over ≥ 300 s (≤ 10 Hz/s). The BL L band
  satisfies this at 2.9 Hz channels; its archive API returns 404 here.
