# Narrowband drift-rate search

A small, dependency-light implementation of the standard technosignature search:
look for radio signals that are narrow in frequency and drift at a constant rate,
which is what a transmitter on a moving planet would look like and what nature
does not produce.

Pure Python + NumPy. It reads SIGPROC filterbank (`.fil`) files directly, so it
runs on Breakthrough Listen and other public radio data without conversion.

## Quick start

```bash
# build a file with three known injections buried in noise
python -m seti_search simulate data/demo.fil --truth data/demo.truth.json

# find them
python -m seti_search search data/demo.fil --out data/demo.hits.csv

# inspect any filterbank file
python -m seti_search info data/demo.fil

# test the candidates against an OFF scan: ON points at the target, OFF away
python -m seti_search verify data/on.fil data/off.fil \
  --site gbt --ra 75.0 --dec 20.0 -o data/verified.csv
```

`search` prints a table and writes a CSV:

```
  #   frequency (MHz)      drift (Hz/s)    SNR   channel
--------------------------------------------------------
  1      1420.397140          +0.5034     34.2      1001
  2      1420.379972          +0.0511     31.8      7001
  3      1420.385682          -0.9978     28.4      5000
```

## What the numbers mean

| column | meaning |
| --- | --- |
| `freq_mhz` | frequency at the first integration of the file, sky-referenced only to the receiver (no barycentric correction) |
| `drift_hz_s` | rate of change of frequency with time; positive means rising. A transmitter on a planet drifts at least from the observer's rotation, typically `~0.01-1 Hz/s` |
| `snr` | detection statistic: boxcar-smoothed dedoppler power above the band median, in units of the robust noise level (1.4826 × MAD). Not a Gaussian sigma in the strict sense, but comparable |
| `channel` | channel index in the file's own channel order (which descends in frequency when `foff < 0`) |

## Pipeline

1. **Bandpass control.** Estimate the bandpass with a block median across
   frequency and interpolate it back to every channel, then divide. Narrowband
   features do not bias a block median, so carriers survive this step and are
   handled next rather than being absorbed into the bandpass.
2. **RFI excision.** Blank only what cannot be a signal: dead or saturated
   channels (the noise variance is gone) and overwhelming carriers (median
   power tens of sigma above the noise floor), plus broadband glitches, which
   are caught along time. Stationary tones at ordinary amplitude are left
   alone and surface as zero-drift candidates — one file cannot tell a beacon
   from a carrier, so the pipeline reports rather than guesses, and an
   early "clean up the obvious RFI" pass does not quietly eat slow real
   signals.
3. **Dedoppler.** For each trial drift `d`, integration `i` is shifted by
   `round(-d · t_i / Δf)` channels and the rows are summed. A signal drifting at
   `d` adds coherently; noise does not. Shifts are sign-correct for bands that
   descend in frequency (`foff < 0`), which is most of them.
4. **Detection.** The summed spectrum is boxcar-smoothed (variance-normalised so
   the noise level is flat at the band edges), centred and scaled by a robust
   noise estimate, and peaks are taken above the threshold with suppression in
   both frequency and drift. Suppressing in drift matters: one real signal also
   scores across neighbouring drift trials.
5. **Noise calibration.** Step 4 divides by `1.4826 × MAD`, which is a Gaussian
   estimate. Total-power data is wider and heavier-tailed than that: on the
   public fast-transient file the dedoppler map's bulk is `1.53×` wider than the
   robust estimate, so a peak printed as `5.6σ` is a `3.7σ` peak, and the flux
   limit and false-alarm budget printed beside it inherit the optimism. Each
   band's map is therefore measured against its own bulk width - iterated
   sigma-clipping around the median, so a real signal cannot set the noise level
   it is judged against - and divided by it. `search` prints the factor whenever
   it is more than 2 per cent from 1. Because the correction happens before peak
   extraction, everything downstream (thresholds, Jy limits, chance
   probabilities, linking) is in the sigma the file actually has. Correcting the
   bulk does not correct the tail, and a false-alarm budget is a claim about the
   tail, so `search` also prints each map's measured upper tail - how many cells
   sit above 4 and 5 sigma against the Gaussian count - and the counts travel in
   the metadata. A ratio well above 1 there means the budget is a lower bound on
   the chance of a noise peak, not an estimate of it. Read it knowing that a
   bright candidate is included: its smear cone fills cells above every
   threshold across many drift trials, so a map holding a strong signal shows a
   large excess that is the signal, not the noise.

Cost is `O(n_drifts × n_ints × n_channels)` per band, vectorised in NumPy. On
one core, a 262144-channel band with 200 integrations comes to 133 drift trials
and takes about 10 seconds (measured).

## Options worth knowing

| flag | why |
| --- | --- |
| `--max-drift` | search range in Hz/s, default ±4, clamped to the channel-crossing limit `0.5·Δf/Δt` |
| `--boxcar` | frequency smoothing; wider is more sensitive to broadened signals, narrower resolves close pairs |
| `--snr` | threshold, default 10 |
| `--max-ints` | integrations to use, default 4096, `0` for all |
| `--tsamp` | override for headers that omit the integration time |

Lower drift rates are the most sensitive region (the signal stays aligned), which
is also where terrestrial interference lives.

## ON/OFF verification

One file gives candidates; `verify` decides which of them survive a cadence. For
each ON candidate it predicts where that signal should be in the OFF scan, using
the observer's motion, and then looks for it there:

- early in the cadence the source is at one frequency and drift; by the OFF scan
  Earth's rotation has changed the projected acceleration, so the same source
  sits at a slightly different drift (and, for a long gap, a shifted frequency);
- a terrestrial signal does not care where the telescope points, so it is still
  there, at the frequency it always had.

So each candidate is checked against two hypotheses with different windows: the
predicted sky position (window grows with the gap, because the prediction
inherits the uncertainty in the measured drift) and the unchanged terrestrial
position (a couple of channels, because a carrier does not move). A candidate
with an OFF counterpart at either is **rejected**; one with neither is
**confirmed**, which means it behaved like a source, not that it is one.

```
  2 candidate(s) in ON, 1 in OFF at 10 sigma
  site +37.2340, -118.2820, target RA 75.000 Dec +20.000
  epochs 58967.230000 / 58967.271667 MJD  (cadence 60.0 min)

  #   frequency (MHz)  drift (Hz/s)     SNR  verdict
----------------------------------------------------------------
  1       1420.389998       -0.0471    49.7  confirmed
  2       1420.391418       +0.0099    30.5  rejected
  #1 confirmed: absent from the OFF scan at its predicted frequency
      predicted in OFF: 1420.389776 MHz, -0.0759 Hz/s  (drift moves 0.0288 Hz/s between scans)
  #2 rejected: present in the OFF scan at the same frequency and drift: not
      switched off by pointing away, so not a sky signal
      predicted in OFF: 1420.391401 MHz, -0.0189 Hz/s  (drift moves 0.0288 Hz/s between scans)
      OFF hit: 1420.391418 MHz @ +0.0098 Hz/s, 28.4 sigma

1 confirmed, 1 rejected
note: the scans predict drifts differing by only 0.029 Hz/s, less than the 0.05
Hz/s tolerance: a longer cadence would make the drift test decisive
```

The physics is in `seti_search/cadence.py`: a low-precision solar ephemeris plus
Earth rotation gives the barycentric velocity and acceleration, which convert
between the drift a transmitter has and the drift the telescope sees. That
prediction is worth a test in itself, and `tests/test_cadence.py` checks it
against things that are not derived from it - Earth's orbital speed and its
perihelion excess, the size of the equatorial centripetal term, sidereal time at
J2000, and the fact that a fixed source cannot drift faster than the observer's
projected acceleration allows.

Supply `--ra/--dec` for the pointing (read from `src_raj`/`src_dej` if the header
has them), `--site`/`--lat`/`--lon` for the telescope, and the epochs. Epochs come
from the header when they are there (`tstart`, or `MJD_start` plus `MJD_hour`,
`MJD_minute`, `MJD_second`), and `--on-mjd`/`--off-mjd` override them. Without
them verification falls back to comparing frequencies only, says so in its
notes, and checks no drift at all.

The check assumes the OFF scan is about as sensitive as the ON one. If it is
shorter, has fewer usable channels, or used a higher threshold, it will fail to
show real signals and confirm everything - compare the integration times before
believing a confirmation.

The `drift_discrimination_hz_s` in each verdict is how far the drift moved
between the scans. If that is smaller than the drift tolerance, the cadence is
too short for the drift test to mean anything, and `verify` says so instead of
reporting a confirmation it cannot support.

## A whole run

`verify` tests one cadence. `link` tests a run — several ON scans of the same
target with OFF scans in between — by demanding that the signal hold the same
**barycentric** drift every time. That is a much stronger statement than
agreeing as measured, because the drift a telescope sees moves as Earth turns:
two scans six hours apart see the same source at drifts differing by ~0.2 Hz/s.

```bash
python -m seti_search link \
  --on data/on1.fil --on data/on2.fil --on data/on3.fil \
  --off data/off1.fil --site ovro --ra 75 --dec 20
```

Every hit is referred to two run-fixed quantities: its barycentric drift, and
its barycentric frequency at one reference epoch (a drifting source's own
frequency changes, so raw frequencies do not line up). Hits that agree on both,
coming from several distinct ON scans, form an event. The event is then tested
against every OFF scan at the position the prediction gives for that scan's
epoch *and* at the position its own hits had, since a carrier does not move. A
hit that does not move at all is linked in the topocentric frame instead — it
is consistent there even when its implied barycentric drift is not — so it gets
a verdict instead of vanishing from the report.

```
run: 3 ON, 1 OFF, 7 candidate(s)

  #   frequency (MHz)  drift (Hz/s)    score   ON  scans  status
------------------------------------------------------------------------
  1       1420.301605       -0.0014     87.0    3      3  survived
  2       1420.304160       +0.0695     50.4    3      3  seen in off
  #1 survived: the same barycentric drift in 3 ON scan(s)
      as measured: -0.0457  -0.1494  +0.0108   (spread 0.1603 Hz/s)
      frequency scatter 0.113 kHz, drift scatter 0.0089 Hz/s
  #2 seen in off: the same barycentric drift in 3 ON scan(s)
      OFF hit in data/run_off1.fil: 1420.391422 MHz @ -0.0096 Hz/s, 30.3 sigma
```

The first event is the one that matters: its drift moves by 0.16 Hz/s between
scans as measured, which no single tolerance would link, and reduces to
`-0.0014` once the rotation is removed.

Scans are rarely equally deep, so scoring is sensitivity-aware rather than a
plain sum. The detection statistic scales as `SNR = A·√(N·W)` for `N`
integrations and a boxcar of `W` channels, so each hit is divided by its own
sensitivity to recover the source amplitude and the scans are combined by
inverse variance — with equal scans that is exactly `√ΣSNR²`. The threshold is
deliberately *not* part of the weight: it decides which hits exist, not how
strong one is. It does decide whether an absence means anything, so from the
combined amplitude the SNR a missing scan should have shown is predicted, and if
that beats the scan's own threshold the event is flagged as missed, not quietly
scored lower.

Two things worth knowing. The frequency window grows with the run, because
referring to one epoch multiplies the drift error by the elapsed time: 0.005
Hz/s over 12 hours is 216 Hz, not a channel or two. And an epoch can be pinned
per file with `PATH@MJD` when a header does not carry one.

Without a site, a position and epochs the run still links, but on topocentric
drift — which only works while the drift has not moved much, so long runs stop
linking, and the output says which mode it used.

## Flux density

Every number above is in units of a scan's own noise: a 50 sigma hit in a deep
scan and a 50 sigma hit in a shallow one are not the same source. To put a run
on a physical scale - and to combine two telescopes into one answer - you need
three numbers per scan plus one assumption.

    S = a * SEFD / sqrt(df * tsamp)          [Jy]

with `df` the channel width, `tsamp` the integration time, and `SEFD` the system
equivalent flux density: the flux that doubles the system temperature,
`SEFD = Tsys / gain = 2761 * Tsys / A_eff` for an aperture area in square metres
and Tsys in kelvin. Give it either way:

    # one receiver: one value for every scan
    python3 -m seti_search link --on data/a.fil --on data/b.fil --sefd 10

    # two receivers in one run: override per file
    python3 -m seti_search link --on data/gbt.fil --on data/parkes.fil \
        --sefd data/gbt.fil=2 --sefd data/parkes.fil=35

    # or supply Tsys and gain and let it compute the SEFD
    python3 -m seti_search link --on data/a.fil --tsys 25 --gain 0.1

Scans are combined by inverse variance, so the quiet receiver carries the answer
and the noisy one only nudges it. A scan given no SEFD still counts towards the
score; it is just reported as absent from the flux, and the run says so.

The assumption is line width. An event's amplitude is defined for a signal
matched to the boxcar, and a line narrower than that loses by `W/w` - the boxcar
adds noise across all `W` channels but signal across only `w`. The default
`--line-width 1` treats every line as unresolved, which is what a technosignature
search looks for and the conservative choice. A genuine 8-channel line analysed
with the default is overstated by a factor of 8.

The same four flags work on `search` and on `verify`, so a single file and an
ON/OFF cadence can quote a flux density too - the formula is the same one, with
the file's own integration count fixing its sensitivity. `verify` calibrates the
ON scan (the OFF is a rejection test, not a measurement of the source).

A null result is still a result. Once a SEFD is supplied, `search` and `link`
report what the observation *ruled out*: with the detection threshold at
`T` sigma, a source of that strength has flux

    S_limit = T * SEFD * (W/w) / sqrt(df * tsamp * N * W)

so anything brighter in the searched band, at any drift inside the search range,
would have been found. Two numbers are quoted: the threshold itself (a source
sitting exactly there is caught about half the time) and the threshold plus
1.64 sigma (about 95%, which is the one to quote as a limit). For a run the
limit comes from the `--min-scans`-th shallowest ON scan, not the deepest - a
source must beat that one to be linked at all.

Those limits are also written to `<results>.limits.csv` beside whatever CSV the
run produced local (`run.linked.limits.csv`, `data/x.fil.hits.limits.csv`): one
row per calibrated ON scan, plus a `run` row naming the scan that set it. The
row carries the threshold, line width, drift range, SEFD, integration count,
boxcar, channel width, integration time and searched fraction, so a limit can be
recomputed months later instead of being taken on faith. On the `run` row the
trial columns describe the scan testing the most cells - the one the threshold
had to hold for - while `set_by` names the scan that set the flux limit, since
those are not always the same file.

What this is not: a calibration.

## False-alarm rates and trial counting

A search does not test one number per channel. It tests one per channel *per
drift trial*, so "10 sigma" and "10 sigma with nothing found in ten million
tries" are different claims. `--false-alarm P` starts from the second one:

    p_cell = 1 - (1 - P)^(1/N),    N = channels * drift trials

with `N` from the same grid the dedoppler builds, `2*floor(max_drift/step) + 1`
per channel (a test asserts the count and the search agree). For a run the
threshold comes from the file testing the most cells, because it has to hold for
every scan searched, and `--false-alarm` replaces `--snr` rather than competing
with it - the substitution is printed.

Every limit states what its threshold costs, whether or not the flag was used:
`expected_false_alarms` is the tail probability times the cell count, and
`false_alarm_log10` keeps the per-cell figure legible where it underflows (60
sigma is about 1e-783).

`verify` reports the same numbers for candidates rather than limits. Calling a
candidate rejected is a claim that the OFF scan saw the same transmitter, and the
claim is only as strong as the ON hit stood above noise: each rejected verdict
carries the noise count expected at the threshold and at the hit's own
significance, and says which reading applies. At 30 sigma over a million cells
that is interference; at a hair over a threshold that fires on noise anyway it is
two fluctuations meeting by chance. Every verdict row in the CSV carries
`threshold_sigma`, `trial_cells`, `snr_tail_log10` and `expected_noise_cells`,
confirmed or not.

Linked events carry the same budget, multiplied. Each member is a peak in one
scan's search, so the chance that noise alone puts a peak that strong somewhere
in that band is `cells * tail(snr)`, and independent scans multiply - which is
what a second ON scan actually buys. Because the product ignores the frequency
and drift agreement the link also demands, it overstates a chance coincidence
rather than flattering the event. The linked-event table prints it as
`log10 p` beside the threshold and cell count, a null run has no row to print,
and the events CSV keeps `threshold_sigma`, `trial_cells` and `chance_log10`.

An event that also appears in OFF is rejected, and that rejection gets the same
treatment. The OFF test accepts a hit anywhere in a frequency window, at any
drift, so the cells it searched are the window's channels times the OFF scan's
drift trials - far fewer than the whole file. The detail block quotes the noise
peaks expected in that window at the threshold and the probability that noise
alone reaches the OFF hit's own strength there: below 1, the OFF peak is not
chance noise and the rejection stands; at or above 1, the window makes peaks
that strong anyway and the rejection is only as good as the prior that
interference explains it. The events CSV keeps `off_window_hz`, `off_cells` and
`off_chance_log10`.

Two caveats that pull in opposite directions. Neighbouring drift trials reuse
integrations and the boxcar smears channels, so `N` overstates the number of
independent cells, which makes the quoted rate pessimistic. The tail itself is
assumed Gaussian while the search normalises with a robust estimate drawn from
the same data, which can make it optimistic. Treat it as a budget, not a
promise. SEFD moves with elevation, weather and the band;
the conversion assumes the file's noise is radiometer noise, which a band full of
interference is not; and the answer is only as good as the SEFD you supplied. It
converts, it does not measure.

## Golden-angle rungs

The φ framework in `README.md` is addressed on every event: `rung(x) =
ln(x)/ln(φ)`, a logarithm to base φ. The golden angle belongs to the same
constant — a full turn is exactly φ² of them (360 / 137.507764… = 2.618) — and
`rung(φ) = 1`, `rung(φ²) = 2`.

As a ruler it says nothing by itself, and the tool says so rather than hinting
otherwise. For any smoothly distributed quantity the fractional part of a rung
is uniform, so every value sits near *some* integer rung: a drift of 0.04 turns
up about one time in twelve by chance, and reporting that as ladder membership
would be reporting noise. What is actually testable is arithmetic and
repetition:

- a repeated rung difference across independent pairs means a repeated ratio —
  one scaling operation, i.e. two hits that are one signal at two pipeline
  scales, or one resolution imposed on many objects;
- a population whose fractional rungs are *not* uniform is quantized by the
  measuring pipeline — a channel grid, a drift step, a chunk boundary — rather
  than by the sky, and a Kolmogorov–Smirnov test against uniformity says how
  strongly.

Both can only ever *detract* from a candidate. Passing them shows nothing; that
asymmetry is the whole value of running them.

Rung addressing follows Tyree Jones (tyofthestarz), the φ² = φ + 1 framework,
prior art from December 2025.

## NV centres: the same read on an exactly known object

A nitrogen-vacancy centre is the best-behaved target this arithmetic can be
pointed at: its transition frequencies are tabulated, its linewidths are known
and its geometry is exact. `scripts/nv_check.py` runs the read on ten of its
constants - D = 2.87 GHz (plus/minus 1 MHz, drifting -74.2 kHz/K), 14N
A_parallel 2.16 MHz, A_perp 2.70 MHz, quadrupole 5.01 MHz, 15N 3.03 and
3.65 MHz, 13C nearest-neighbour 130 MHz, and zero-phonon lines at 637 nm
(1.945 eV, 15687 cm^-1), 575 nm (2.156 eV) and 1042 nm.

- **Address.** Expressed in six units - Hz, MHz, GHz, eV, cm^-1, nm - the same
  sixty tests give **3 hits against 6.0 expected** from uniform phase,
  binomial P(>= 3) = 0.947. Per unit: Hz 1, MHz 0, GHz 1, eV 1, cm^-1 0, nm 0.
  Rewriting the same frequency in another unit shifts every phase by
  log(10^6)/log(phi) = 28.7098, so the verdict is a property of the ruler. The
  14N quadrupole sits 0.006 rungs from an integer in GHz and 0.35 rungs away in
  MHz.
- **Geometry.** The NV axis is <111>, so axis-to-bond is arccos(-1/3) =
  109.4712 deg, exactly. Every angle a bond or axis pair can subtend is an
  arccosine of 0, +/-1/2, +/-1/sqrt(2), +/-1/sqrt(3), +/-1/3, +/-sqrt(2/3) or
  +/-2/3, and the nearest of those to 137.5078 deg is **135.0000 deg, off by
  2.5078 deg exactly**. Allowing high-index directions (|hkl| <= 4) produces 2654
  distinct angles, one of which lands 0.0420 deg away - a statement about how
  many pairs exist, not a near-hit. cos(137.5078 deg) = -0.737369 is the
  arccosine of a fifth-degree algebraic number, and a periodic lattice admits no
  five-fold rotation at all, phi being the five-fold divergence.
- **Pitch.** D/nu(HI) = 2.020549. The nearest simple forms are 2 (off 1.017%)
  and sqrt(phi^3) (off 1.862%), both chosen after the number was known; the two
  quantities are unrelated couplings and no baseline enters, so the pitch has
  nothing to act on.
- **Tolerance.** 0.05 rungs at 2.87 GHz is a 69.9 MHz cone, +/-35 MHz, against a
  2 MHz ensemble ODMR line: the criterion is ~35x coarser than the instrument,
  so it is falsifiable - a line 100 MHz off its rung is caught - but
  uninformative, because one frequency in ten satisfies it. At 1e-3 rungs the
  cone is 1.38 MHz, comparable to an ensemble line and a thousand times wider
  than a single centre resolves, which puts the tight version within reach of a
  bench; there the temperature term bites, since D moves 74.2 kHz/K and holding
  a 1.38 MHz cone needs ~19 K of control.

### The tightened claim

`scripts/phi_tight.py` runs the version that could be wrong: unit fixed to Hz,
tolerance 1e-3 rungs. Two numbers decide it. **Content** is 2t, so 0.2% per
test. **Threshold** is tau = ln(1 + delta/f) / ln(phi): while the cone is
narrower than the error bar no measurement can refute an on-rung value, so a
claim at tolerance t is checkable only when **t > tau**. Measured thresholds:
Cs 133 hyperfine 0 (exact, the SI second), the 21 cm line 1.46e-9, Rb 87
3.04e-10, NV D 7.24e-7, 14N A_parallel 1.92e-2, 13C 1.59e-2.

At 1e-3 rungs, **0 of 6 pass** (0.012 expected) and 4 of 6 could have been
checked - the two that could not are the 1% hyperfine terms, where one percent
of the frequency is twenty times the cone. At the framework's 0.05 all six can
be checked, and every pass is a 10% event, which is the problem. Content falls
as t falls while checkability needs t > tau, so the informative region is
tau < t << 1: inside it the claim is testable and empty, and outside it a pass
means nothing.

### Fitting instead of assuming

`scripts/phi_fit.py` takes the other route: let the scale float. On the three
candidates in the STARE2 file, with the offset free and the ratio fixed at phi
(one parameter), the best offset leaves **0.1438 rungs** where a joint pass
needs every gap within 0.10 - it fails, and uniform triples in the same band
pass 0.178 of the time (P(<= observed residual) = 0.334). Fitting the ratio as
well and searching the integers reproduces the candidates with rms 0.0000 rungs
- and reproduces uniform noise with median rms 0.0012, forty times inside the
0.05 tolerance. That is the trade: two continuous parameters plus free integers
absorb the spacing, so the fit cannot separate a candidate set from noise, which
is the discrimination the framework exists to provide. A model with a free ratio
per gap is an interpolant, not a law.

## Real data

Any SIGPROC `.fil` works. A small, public example is the STARE2 triggered burst
dataset (CaltechDATA record `jhn25-fsd29`, CC0), which is a real filterbank file
written by an instrument that is not Breakthrough Listen:

```bash
curl -sL -o data/stare2_ovro.fil \
  "https://data.caltech.edu/api/records/jhn25-fsd29/files/candidate_ovro_20200428.fil/content"
python -m seti_search info data/stare2_ovro.fil
```

That file is a fast radio burst, not a SETI target — it is useful for exercising
the reader against real header quirks (no `header_size` field, extra `MJD_*` and
`*_sample` keys, keys in a non-standard order), not for finding anything.

## Limits (read before trusting a candidate)

- **A candidate is not a detection.** A real signal has to reappear in an ON scan
  and vanish in an OFF scan, and its drift has to match the sky position. One
  file, one telescope, one time: that is a candidate.
- **Barycentric correction is used for decisions, not for reporting.** `search`
  reports topocentric drifts. `verify` converts them to barycentric, predicts the
  OFF scan and converts back, which needs the site, the pointing and the epochs;
  without all three it compares frequencies only and says so. The solar ephemeris
  is low-precision (0.01 degrees), which is ample for drift *differences* over a
  cadence but means an absolute barycentric frequency is not what you get here.
- **A confirmation is not a detection.** `verify` removes the objection that
  kills most candidates; it does not prove the survivor is artificial. More
  cadences, a second telescope and a check against known transmitters are still
  required.
- **Incoherent dedoppler.** Sensitivity is lower than a coherent de-doppler
  search; this trades a few percent for a large speed-up.
- **The RFI excision is a heuristic.** It removes dead channels, broadband
  glitches and carriers tens of sigma above the noise, but a strong stationary
  tone that is not quite that loud survives as a candidate, and so does
  terrestrial interference that drifts. Nothing here decides what is real:
  that is the ON/OFF cadence's job.
- **16-bit samples are read as unsigned** and polarisations are averaged, so
  Stokes parameters are not available.
