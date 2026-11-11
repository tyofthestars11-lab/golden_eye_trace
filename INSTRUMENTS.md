# Instruments that could carry a φ phase lock, and the ones that cannot

The lock is one ratio. A drifting line advances `ν = drift · t / Δf` turns per
sample, and the lock asks for `ν = θ = 1/φ² = 0.381966`. Inverting it gives the
only instrument property that matters:

```
Δf / t = drift / θ        (channel width over integration time)
```

So a **1 Hz/s** line needs `Δf/t = 2.62 Hz per second of sample`: about 2.6 Hz
channels with 1 s integrations, or 260 Hz channels with 100 s integrations. Every
instrument below is judged on that ratio alone, and the ratio is why most of them
fail.

Run `python3 scripts/phi_geometry.py` for the radio table and
`python3 scripts/ovro_lwa_test.py` for the open-archive test.

## Radio: what exists versus what is needed

| archive / product | Δf | t | Δf/t | θ-locked drift | open? |
| --- | --- | --- | --- | --- | --- |
| OVRO-LWA `spec_fits` (731 ch, 256 ms) | 100.5 kHz | 0.256 s | 3.9e5 | 150 kHz/s | **yes**, S3, CC BY 4.0 |
| OVRO-LWA `image_hdf` (144 ch, 10 s) | 384 kHz | 10 s | 3.8e4 | — | yes, but one cube per file |
| LWA1 drspec (read here) | 19.14 kHz | 1.284 s | 1.5e4 | 5,694 Hz/s | yes |
| STARE2 (read here) | 122.07 kHz | 0.4 s | 3.1e5 | 116,566 Hz/s | yes, one integration |
| CAMRAS Dwingeloo | 136.7 kHz | 468 µs | 2.9e8 | 1.1e8 Hz/s | yes, direct GET |
| Breakthrough Listen L band (nominal) | 2.9 Hz | 300 s | 9.7e-3 | 0.0037 Hz/s | metadata only; files 404 here |
| **needed for a 1 Hz/s chirp** | — | — | **2.62** | 1 Hz/s | — |

The open products sit four to eight orders of magnitude too coarse, which is what
pushes their locked drift onto chirps of 10⁴–10⁸ Hz/s. The OVRO-LWA spectrogram is
the best of them and is the one actually tested: a Type III solar burst drifts at
MHz/s, i.e. tens of channels per 256 ms sample, so even the fastest real drifting
lines in that archive alias past Nyquist rather than sitting on the lock.

## JWST: the same ratio, and why it is further away

JWST is not a radio telescope, but the lock is scale-free — a spectral channel and
an exposure time play the same role as a frequency channel and an integration. The
question is whether an astrophysical line moves 0.38 channels per exposure.

| instrument | bin spacing | shortest integration | drift per exposure ÷ spacing |
| --- | --- | --- | --- |
| NIRSpec G140H | 0.24 nm/px | ~1.1 s | ~3e-8 per second of exposure |
| NIRSpec PRISM | ~11 nm/px | ~1.1 s | ~7e-10 per second |
| NIRISS SOSS order 1 | 0.96 nm/px | ~5.5 s | ~2e-9 per second |
| MIRI LRS | ~23 nm/px | ~0.3–1.6 s | ~3e-9 per second |
| **a full hot-Jupiter transit** | — | 3–8 h | **~3e-4 channels total** |

The arithmetic behind the last row: a transiting hot Jupiter changes radial
velocity by about 0.1 m/s per minute, so a line at 1.1 µm moves
`Δλ = λ·Δv/c ≈ 6×10⁻⁹ nm/s`. Against a 0.24 nm G140H pixel that is `2.6×10⁻⁸` of
a channel per second, and even the whole ~20 m/s transit span is `7×10⁻⁵ nm`,
**three parts in ten thousand of one channel**.

To reach 0.382 channels per exposure you would need the line to cross a third of a
channel between reads — for G140H, a velocity change of about **90 km/s per
exposure**, three hundred times a planet's orbital speed and of order the escape
speed of a white dwarf. JWST's genuine time-variable systematic — trace drift at
the few-ppm level, `10⁻⁶` to `10⁻⁵` of a pixel — is closer to the mark than any
Doppler signal, which is the honest summary: **the lock would be measuring the
instrument, not the sky.**

That is not a defect of JWST. It follows from `Δf/t`: JWST trades channel width
for very long integrations, so its ratio is tiny, while the lock wants the
opposite. The mirror that costs ten billion dollars is further from this
particular test than a 19 kHz-channel radio array.

## What would actually carry the test

1. **A radio spectrometer at `Δf/t ≈ 2.6 Hz/s`.** Fine channels, second-scale
   integrations: 3 kHz channels at 1,150 s, or 700 Hz at 270 s. The Breakthrough
   Listen L band nominally qualifies, at 2.9 Hz channels over 300 s, giving a
   locked drift of 0.0037 Hz/s — an entirely plausible beacon drift. Its archive
   answers 404 from this workspace, so the geometry is right and the data are not
   reachable.
2. **≥987 samples** for a marginal echo and **≥3,948** for the ten-lag
   configuration, because lags are capped at `n/4`.
3. **`t` fixed by the instrument**, not chosen by dedoppler-summing blocks, or the
   lock constrains the binning instead of the source (see `GEOMETRY.md` §4).
4. **`|k| = 1`**, since aliasing otherwise lets any rate match some integer
   advance.

With those, the test is decisive and cheap: a few hours on an existing
fine-channel spectrometer, analysed by `scripts/phi_echo.py` and scored by the
rule in `prereg/phi_lock.json`.

## What the open-archive run found

`scripts/ovro_lwa_test.py` pulls a day of `spec_fits` from the public bucket,
reads the cube without astropy, searches a sub-band to Nyquist and scores every
candidate against the lock. Its numbers — candidates, resolved drifts, lock
distances, and whether the verdict is NOT MEASURED — are written to
`prereg/ovro_lwa.json` rather than quoted here, so this note cannot drift out of
step with the data.

Two properties of the product are worth recording because they cost time to find.
The axis order is `(time, freq, singleton, Stokes)`, not the `(Stokes, 1, freq,
time)` the reader code implies. And the files carry **no `CTYPE`, `CRVAL` or
`CDELT` cards**, so there is no world coordinate system to read: the frequency
axis has to be inferred from axis sizes and then validated against the documented
73.5 MHz span. The script prints that check (`731 x 100.5 kHz = 73.5 MHz`) and
refuses to run if it fails.

`scripts/ovro_lwa_population.py` then widens the search: it lists `spec_fits/`,
picks days inside a size window and searches each to Nyquist, scoring every
drift against the lock. Three 2024 days searched this way — the widest at 1,600
samples x 192 channels, 410 s of data at 246 Hz/s resolution — contain **two
drifting lines**, at **8,195 and 7,867 Hz/s**. Their rates are **0.0200 and 0.0209
cycles per sample, eighteen times below θ = 0.382**, and neither is within 0.05 θ
of an integer, against **0.15 expected by chance**. So the archive does hold
drifting lines, and they sit far from the lock; the earlier NOT MEASURED verdict
was a statement about coverage, not about the sky.
