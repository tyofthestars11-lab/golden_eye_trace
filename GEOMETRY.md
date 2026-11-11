# The geometry a φ phase lock demands

The lock in `scripts/phi_echo.py` says a drifting line advances by
`θ = 1/φ² = 0.381966` turns per integration. In a waterfall that is a statement
about three measured quantities — drift `D` in Hz/s, channel width `Δf` in Hz, and
integration time `t` in seconds:

```
ν = D · t / Δf          cycles per sample
ν = θ   ⇔   t = θ · Δf / D
```

So the lock is not one drift, it is a family: every band has its own φ-locked
drift, set by its channel width and cadence. This note works out what that family
is, what it asks of an instrument, and what is left of the claim once the
integration time is allowed to float. All numbers come from
`python3 scripts/phi_geometry.py`, which writes `prereg/phi_geometry.json`.

## 1. The integration time the lock demands

`t = θΔf/D`, in seconds:

| channel | D=0.01 Hz/s | 0.1 | 1 | 10 | 100 | 1000 |
| --- | --- | --- | --- | --- | --- | --- |
| LWA DRX, 19.14 kHz | 7.31e5 (8.5 d) | 7.31e4 (20 h) | **7,311 (2.03 h)** | 731.1 | 73.11 | 7.311 |
| STARE2, 122.07 kHz | 4.66e6 | 4.66e5 | 4.66e4 (12.9 h) | 4,663 | 466.3 | 46.63 |
| 3 kHz | 1.15e5 | 1.15e4 | 1,146 | 114.6 | 11.46 | 1.146 |
| 700 Hz | 2.67e4 | 2,674 | 267.4 | 26.74 | 2.674 | 0.267 |
| 100 Hz | 3,820 | 382 | 38.2 | 3.82 | 0.382 | 0.038 |

Read down the 1 Hz/s column: to put a realistic beacon on the lock you need
integrations of minutes to days, not seconds. Read across a row: faster drifts
need shorter integrations, and the product `D·t = θΔf` is fixed by the band.

## 2. The lock is always at 76.4% of Nyquist

Since `θ / 0.5 = 0.7639`, the locked drift is always the same fraction of the
fastest drift the sampling can represent, whatever the instrument:

| band | channel | t | φ-locked drift | Nyquist | fraction | largest drift observed |
| --- | --- | --- | --- | --- | --- | --- |
| LWA DRX | 19.14 kHz | 1.284 s | **5,694 Hz/s** | 7,454 Hz/s | 76.4% | 19.1 Hz/s |
| STARE2 | 122.07 kHz | 0.4 s | 116,566 Hz/s | 152,588 Hz/s | 76.4% | 0.0 Hz/s |

That invariance is the useful answer to "which geometry": no choice of band or
channel width moves the lock off 76.4% of Nyquist, because θ is a fraction of a
channel per sample and 0.5 is the sampling limit. Only the drift in Hz/s changes,
scaling with `Δf/t`.

Two consequences for the data in hand. The LWA dynamic spectrum can host the lock
— 5,694 Hz/s is inside the ±7,454 Hz/s window that `phi_echo_recovery.py` searched,
and that search returned **11 candidates at 4σ, every drift between 10 and 19 Hz/s,
none above the 21 Hz/s resolution**. The observed drifts are 298× below the lock.
STARE2 cannot host it at all: the file spans 0.385 s with a 0.4 s cadence, one
integration, so no drift is defined and the rate does not exist.

So the answer is not "no instrument was found". The low band already had the right
geometry, was searched across the whole window that contains the lock, and holds
no resolved drifting line at any rate, let alone at θ.

## 3. What a slow beacon would need

For the population SETI actually searches — narrowband carriers drifting at
roughly 0.1 to 10 Hz/s — the lock demands:

| target | t per sample | dwell for ten Fibonacci lags | bandwidth crossed |
| --- | --- | --- | --- |
| 1 Hz/s, 19.14 kHz channels | 7,311 s (2.03 h) | **334 days** (3,948 samples) | 28.9 MHz |
| 1 Hz/s, real track length 5,872 | 7,311 s | 496.9 days | 42.9 MHz |
| 0.1 Hz/s, 19.14 kHz channels | 73,108 s (20.3 h) | 3,341 days | 28.9 MHz |

Ten lags is the count the echo is calibrated at, and it requires 3,948 samples
because lags are capped at `n/4`; at 987 samples only seven lags exist.

An 83-day continuous observation at 2-hour integrations is not a hard limit in
principle but it is far outside any survey cadence in `BANDS.md`, and the signal
has to stay coherent across every block of it. Two hours per integration is also
a demand on the transmitter's stability, not just the telescope's.

## 4. Where the claim loses its content

Two ways, both arithmetic rather than opinion.

**The integer advance must be bounded independently.** Within Nyquist the drift
lattice `ν = kθ` has exactly one physical member, `k = 1` at 0.382, because
`k = 2` at 0.764 is outside. But a discrete-time rate is defined modulo one cycle
per sample, and θ is irrational, so `kθ mod 1` is equidistributed in `[0, 0.5]`.
Allowing `k` up to **16** already leaves no gap wider than 5% of θ anywhere in the
window: every drift then matches some integer advance. The lock is a sharp test
only while `k` is restricted — and the restriction has to come from the sampling
geometry, not from the data, or it is a free parameter.

**The integration time is the analyst's choice.** `t = θΔf/D` can always be
satisfied by dedoppler-summing blocks of the right length: a 1 Hz/s carrier in a
19.14 kHz channel sits exactly on the lock when the blocks are 7,311 s long, and a
10 Hz/s carrier when they are 731 s. Unless `t` is fixed by the instrument before
the data are seen, the lock constrains the binning rather than the source, and no
observation can refute it.

## 5. Requirements for a decisive test

1. **Fix `t` from the raw dump** and record it before the search; no coherent
   re-binning, no block averaging chosen after the drift is known.
2. **Search to 76.4% of Nyquist.** A pipeline that caps drift at a few Hz/s never
   reaches the lock; the window is `Δf/(2t)`, which for 19 kHz channels at 1.284 s
   is 7,454 Hz/s. `scripts/dedoppler_audit.py` checks the search can live there:
   the grid reaches **99.0% of Nyquist**, a θ chirp comes back **0.0023 θ** from
   its injected rate at **27.9σ** against zero hits in noise, and the worst rate
   error over five injections is **0.0256 θ** — inside the 0.05 tolerance with a
   factor of two to spare, though two and a half times the half-step rule. The
   rate grid is `0.5·boxcar/n_ints` in cycles per sample, `0.0196 θ` here, so the
   grid is not what limits the test.
3. **Bound the integer `k`** to `|k| = 1` on the grounds that `k = 2` is outside
   Nyquist, and state that bound in the pre-registration.
4. **Collect at least 3,948 samples** for ten Fibonacci lags, the count the echo is
   calibrated at; 987 samples gives only seven. With `t` fixed, that sets the
   dwell, so the requirement is a real observing cost.
5. **State the chance.** A random drift within 5% of θ in the `k = 1` window is a
   **0.0764** event per candidate; multiplicity over candidate lists is applied on
   top of that.

With those five fixed, the lock is falsifiable and sharp. Without them it is
satisfiable on demand, which is the same as untestable.

## 6. What open archives can actually offer

The lock needs `channel/tsamp = drift/θ`, so a 1 Hz/s chirp needs **2.62 Hz of
channel per second of sample** — either about 1 s integrations with 2.6 Hz
channels, or 100 s integrations with 260 Hz channels. Nothing openly downloadable
from this workspace sits near that:

| archive | channel | tsamp | channel/tsamp | θ-locked drift | access |
| --- | --- | --- | --- | --- | --- |
| LWA1 drspec (read here) | 19.14 kHz | 1.284 s | 1.5e4 | 5,694 Hz/s | open, already searched |
| STARE2 (read here) | 122.07 kHz | 0.4 s | 3.1e5 | 116,566 Hz/s | open, one integration |
| CAMRAS Dwingeloo | 136.7 kHz | 468 µs | 2.9e8 | 1.1e8 Hz/s | open, direct GET |
| OVRO-LWA spec_fits (AWS) | 100.5 kHz | 0.256 s | 3.9e5 | 1.5e5 Hz/s | open S3, tested here |
| BL L band (nominal) | 2.9 Hz | 300 s | 9.7e-3 | 0.0037 Hz/s | metadata only, 404 here |

So the pre-registered test cannot run on any openly reachable data found here, and
the reason is now precise: `channel/tsamp` is four to eight orders of magnitude
too coarse, which pushes the lock onto chirps of 10⁴ to 10⁸ Hz/s that no
transmitter produces. The one geometry that would put θ at a plausible drift is
fine-channel beamforming — the Breakthrough Listen L band at 2.9 Hz channels over
300 s gives θ = 0.0037 Hz/s — and that archive's file API answers 404 from here.

The low-frequency options are worth naming because they are genuinely open. CAMRAS
publishes search-mode filterbanks with per-observation header sidecars and no
login, but at 136.7 kHz channels. The OVRO-LWA solar archive is a public S3 bucket
under CC BY 4.0 with 24 kHz resolution and 64 ms cadence — the right resolution,
but an imaging product whose θ-locked drift is 1.4e5 Hz/s. Neither meets the
ratio, and no archive reachable here meets both the ratio and the sample count.
