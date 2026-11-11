# Pulling a second band: what is open, what qualifies, and what band A cannot do

The pre-registered test in `prereg/phi_2band.json` needs two bands at least a
decade apart, each with at least three resolved drifting lines. This note records
which open archives were reachable from this workspace, what each actually
delivers, and the two structural obstacles that remain. Every status was checked
on 2026-11-11 with `scripts/fetch_band.py --check`; nothing is copied from a
search result.

```bash
python3 scripts/fetch_band.py --list                  # registry, decade column
python3 scripts/fetch_band.py --check lwa-drx mwa     # HTTP probe, no download
python3 scripts/fetch_band.py --get lwa-drx --with-payload   # 1.2 kB .inf + 16 MB .dat
python3 scripts/fetch_band.py --gate                  # does the ledger clear a decade?
python3 scripts/presto_dat.py --dat <f>.dat --inf <f>.inf --inspect
python3 scripts/presto_dat.py --dat <f>.dat --inf <f>.inf --out data/lwa_dm0.fil
python3 scripts/phi_twoband.py --band-a data/stare2_ovro.fil --band-b <second.fil>
```

## The survey

| source | band | centre | decades from 1405 MHz | status here |
| --- | --- | --- | --- | --- |
| STARE2 OVRO (local) | 1280–1530 MHz | 1405 | 0.00 | first band; 3 candidates at 3σ |
| **LWA1 DRX, B0329+54** | **32.2–51.8 MHz** | **42.0** | **1.52** | **200; `.inf` + 16.1 MB `.dat`** |
| LWA1 DRX, tuning 2 | 64.2–83.8 MHz | 74.0 | 1.28 | 200; same directory |
| LOFAR HBA | 110–190 MHz | 150 | 0.97 | no connection (000); LTA login |
| MWA | 72–300 MHz | 186 | 0.88 | no connection (000); Pawsey login |
| BL Open Data (GBT/Parkes X band) | up to 11.2 GHz | 6150 | 0.64 | 200 frontend, `/api/*` 404 here |
| Zenodo 22848961 (Parkes) | 1360–1680 MHz | 1520 | 0.03 | 200, 377 280 bytes, FITS despite `.sf` |
| FAST | 1000–1500 MHz | 1250 | 0.05 | 200 landing page; login |

Only the LWA clears a decade. Its `.inf` reports 4096 channels of 4.785 kHz
across 19.6 MHz and an 843.14 s scan, so its drift resolution is **5.68 Hz/s**
against STARE2's 3.17e+05 Hz/s — the LWA is 5.6e4 times finer. The ledger at
`bands/lwa-drx.json` records the centre, the channel count, the resolution and
the SHA-256 of both files; both hashes were recomputed from the downloads and
match.

## What the payload actually is

`scripts/presto_dat.py` reads the float32 payload, checks it against the `.inf`,
and decides the shape from the numbers:

- 4,034,560 float32 samples **exactly equals** the `.inf` bin count, and
  4,034,560 × 208.98 µs = 843.14 s, the stated observation length;
- lag-1 autocorrelation **+0.475**, the red noise of a 4.8 kHz-wide signal
  sampled at 209 µs;
- read as 985 rows × 4096 channels instead, the sample time becomes 0.856 s and
  the channel means are flat to **3 × 10⁻⁴** — there is no channel structure to
  be a frequency axis.

So the `DM0.00_topo.dat` product is a **single dedispersed power series** at
42 MHz, not a band. Converting it gives a one-channel SIGPROC file of 843.14 s
whose every sample round-trips exactly through this project's parser — a valid
conversion of a product that cannot carry a drift measurement.

The 2-D product in the same directory is `drx_56797_B0329+54_b3t1_subs_0001.fits`
(494 MB, search mode, sub-banded), with a 36 MB `drspec_*` dynamic spectrum
beside it. Those are PSRFITS and need a reader this project does not have yet;
the `.dat` route cannot substitute for it, because a time series has no
frequency axis at all.

## What the gate rejects

`--get zenodo-vela` downloads cleanly and looks like a gift: 377 kB named `.sf`,
i.e. SIGPROC by extension. It is FITS — `OBSFREQ 1520`, `OBSBW 320` — and **0.03
decades** from the first band. Name and extension were both misleading; the
header decided, and the ledger marks it `qualifies: false`.

## Band A: why a longer cadence cannot fix it

Drift resolution is channel width over integration time: `Δν / T`. STARE2 stores
122.07 kHz channels and 0.4 s integrations, so it resolves 3.17e+05 Hz/s — larger
than any plausible narrowband drift, which is why every candidate in it has zero
resolved drift. Holding the channel width fixed, resolving a drift `d` needs

| drift | contiguous scan needed at 122.07 kHz | channel width needed in a 300 s scan |
| --- | --- | --- |
| 0.1 Hz/s | 14.1 d | 30 Hz |
| 1 Hz/s | 1.41 d | 300 Hz |
| 10 Hz/s | 3.39 h | 3.0 kHz |
| 100 Hz/s | 20.3 min | 30 kHz |
| 1000 Hz/s | 2.03 min | 300 kHz |

At the drifts real beacons are expected to show, band A would need hours to days
of contiguous integration, against the 0.4 s it records; STARE2's design is
millisecond-cadence transients, not long narrowband stares. **The channel width
is the limit, not the cadence**: to resolve 10 Hz/s in a 300 s scan the
spectrometer would need 3 kHz channels, 40 times narrower than the 122 kHz it
writes. `scripts/phi_twoband.py` prints this line on every run that has too few
drifts, so the requirement cannot be forgotten.

That leaves the registered test needing a *second* fine-channel band on the high
side as well — a fine-channel 1.4 GHz dataset (BL L band writes 2.9 Hz channels,
0.01 Hz/s over 300 s) paired with the LWA at 1.52 decades. The LWA half now has a
verified band; the high half needs the BL archive's API, which returns 404 from
here.

## Why a decade, quantitatively

The two statements differ by 3.3191 rungs per decade (4.7849 flat against 1.4659
for the pitch), so a slope measured to better than `1.1064 × span` rungs/decade
separates them at 3σ:

| span (decades) | slope precision needed |
| --- | --- |
| 0.10 | 0.1106 |
| 0.30 | 0.3319 |
| 0.60 | 0.6638 |
| 0.90 | 0.9957 |
| 1.00 | 1.1064 |
| 1.52 | 1.6816 |

At the LWA's 5.68 Hz/s resolution, a line drifting at 50 Hz/s is measured to
about 11% (0.11 in ln), comfortably inside 1.68; a 10 Hz/s drift is only
marginally resolved (57%), so the lines that carry weight are the faster ones.

## Status

- Band axis: **satisfied**, verified, ledgered at 1.52 decades.
- Readable waterfall for that band: **missing** — needs a PSRFITS reader for the
  494 MB `_subs` file, because the `.dat` is a time series.
- Band A drift resolution: **structurally insufficient** at 122 kHz channels.
- Registered test: **NOT MEASURED**, and it now says which of these is blocking
  each time it runs.
