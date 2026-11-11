# A falsification campaign on the φ rung framework

This report measures the quantitative claims of the golden-ratio framework in
`README.md` — rung address `r(x) = ln(x)/ln(φ)`, golden pitch `μ = 17.0324°`,
golden-angle sky address `137.5078°` — and reports every statistic next to its
own null. Framework and prior art: Tyree Jones (tyofthestarz), φ² = φ + 1,
December 2025.

Every harness here is written so it can only *detract* from a candidate. A pass
proves nothing; the point of the exercise is that a coordinate which cannot
lose is not a criterion. Each number below is reproducible with the command in
its section, and each script exits 0 on the data named.

## Summary of verdicts

| claim | test | result |
| --- | --- | --- |
| coordinates separate signal from noise | 24 φ-free injections vs real candidates | same pass rates in both |
| rung address is a property of a frequency | six units | membership flips with the unit, 28.7098 rungs per factor 10⁶ |
| lattice rate flat in frequency (p = 1.0000) | coordinate identity | true by definition, cannot fail |
| pitch rate p = 0.3063 | same measurement | refuted at 3.26×, outside a pre-registered 5% |
| pitch names a baseline | drift → T* | every drift is the pitch at some T*; no baseline is selected |
| golden angle is a sky address | 3 frames, one source | offsets spread 93.5°, chance 0.081 |
| NV constants sit on the lattice | 60 tests | 3 hits, 6.0 expected, P = 0.947 |
| golden angle is a crystal angle | cubic geometry | exact gap 2.5078° to 135.0000°; no five-fold in a periodic lattice |
| tightened claim (Hz, 1e-3 rungs) | 6 precisely known frequencies | 0 pass, 4 of 6 checkable |
| fitting the scale rescues it | offset and ratio free | fits candidates to 0.0000 rungs and noise to 0.0012: no discrimination |
| drift exponent across bands | pre-registered | **NOT MEASURED** — no second band, no resolved drifts |

## 1. The read on real candidates

`python3 scripts/phi_read.py --fil data/stare2_ovro.fil` (STARE2 OVRO
filterbank, 1280.183–1530.061 MHz, 3 candidates at 3σ: 1287.622252 MHz at
3.7σ, 1379.895070 at 3.1σ, 1356.243742 at 3.0σ).

- **Rung.** 2 of 3 sit within 0.05 rungs. Uniform noise in the same band passes
  0.264 of the time. Expressing the same frequencies in Hz shifts every rung by
  28.7098, moves the chance rate to 0.000, and flips 2 of the 3 verdicts. The
  band spans 0.3705 rungs, so every channel shares one address.
- **Pitch.** All three candidates have zero drift, so `T* = ∞`: the pitch is
  satisfied by any baseline longer than 3.45 yr. The file's drift resolution is
  3.172e+05 Hz/s against a ±4 Hz/s search window, so the whole window is the
  pitch at some baseline. The pitch names a baseline, not a source.
- **Angle.** The header has no `src_raj`/`src_dej`, so the read falls back to
  Sgr A* as the reference. One fixed direction gives offsets 128.909°,
  −137.564° and 129.344° across the three candidates: spread 93.5°, with chance
  0.081 over three frames. A frame-dependent quantity is not an address.

## 2. The same read on signals containing no φ

`READ 4` injects 24 beacons with no φ anywhere in their design, in-band, and
recovers 24 of 24 at 24 of 24 injected frequencies. 3 of the 24 land within
0.05 rungs, against 2.4 ± 1.5 expected, and 24 of 24 are pitch-accepted at some
baseline.

That is the decisive comparison: the coordinate pass rates for injected noise
signals and for real candidates are the same, so the read is describing its own
arithmetic rather than the source.

## 3. The hunt, in Hz, at both tolerances

`python3 scripts/phi_hunt.py data/stare2_ovro.fil data/demo.fil data/fine_on.fil`

The null is the fraction of the band the cones actually cover, not `N × 2t`:
a band narrower than one rung never approaches a rung at all, and `N × 2t`
would invent opportunities that do not exist.

- **1280.183–1530.061 MHz, 0.3705 rungs.** The nearest rung is 0.051 rungs
  above the top edge. Coverage is 0.00% at 0.05 rungs and 0.00% at 1e-3, so no
  channel of this file is testable at either tolerance in Hz.
- **Narrow files** (`demo.fil` 0.000034 rungs, `fine_on.fil` 0.0001 rungs) sit
  inside a single gap end to end; their 1e-3 cone is 239 044 channels wide.

Two consequences. First, the "2 of 3 within 0.05" of §1 is a property of the
MHz ruler: in Hz the same band misses every rung. Second, a hunt for
technosignatures on the lattice is vacuous on all data here, at any tolerance
tight enough to be interesting — not because nothing was found, but because the
bands never reach a rung.

## 4. The one test that can fail, pre-registered

`python3 scripts/phi_prereg.py --fil data/stare2_ovro.fil` writes
`prereg/phi_2band.json` — statistic, predictions, 5% tolerance, decision rule,
evidence requirements, seed 1729 — and prints its SHA-256 (`5e4435f1ed0292fc…`
for the STARE2 band) before any measurement.

The statement: the drift exponent `p = d ln|drift| / d ln f` is 1.0000 (4.7849
rungs/decade) if the rung rate is flat, and `tan(17.0324°) = 0.3063` (1.4659
rungs/decade) if the pitch sets the scale.

- **Coordinate identity:** 4.7850 rungs/decade measured — the definition, not
  evidence.
- **Pitch:** predicted 1.4659 against 4.7850, a factor 3.2642, well outside 5%.
  On the framework's own competing statements, the pitch is refuted by one band
  before any signal is examined.
- **Signal side:** **NOT MEASURED.** Requirements: a band B at least a decade
  away, and ≥3 drifting lines resolved in each band. The real file supplies
  neither (drift resolution 3.17e+05 Hz/s, all drifts zero). A slope fitted from
  one band would describe that band's noise, so none is quoted.

## 5. NV centres: an exactly known object

`python3 scripts/nv_check.py`. Constants: D = 2.87 GHz (±1 MHz, −74.2 kHz/K),
¹⁴N A∥ 2.16 MHz, A⊥ 2.70 MHz, P 5.01 MHz, ¹⁵N 3.03/3.65 MHz, ¹³C nearest 130 MHz,
ZPLs 637 nm (1.945 eV, 15687 cm⁻¹), 575 nm (2.156 eV), IR 1042 nm.

- **Address.** Ten constants × six units = 60 tests: **3 hits against 6.0
  expected**, binomial P(≥3) = 0.947. Per unit: Hz 1, MHz 0, GHz 1, eV 1,
  cm⁻¹ 0, nm 0. The ¹⁴N quadrupole is 0.006 rungs from an integer in GHz and
  0.35 in MHz.
- **Geometry.** Axis-to-bond is arccos(−1/3) = 109.4712°, exact. Every angle a
  bond or axis pair subtends is an arccosine of 0, ±½, ±1/√2, ±1/√3, ±1/3,
  ±√(2/3) or ±2/3; the nearest to the golden angle is **135.0000°, off by
  2.5078° exactly**. Allowing |hkl| ≤ 4 produces 2654 exact angles and one lies
  0.0420° away — a multiplicity artefact, not a near-hit. cos(137.5078°) =
  −0.737369 is the arccosine of a fifth-degree algebraic number, and a periodic
  lattice admits no five-fold rotation, φ being the five-fold divergence.
- **Pitch.** D/ν(HI) = 2.020549; nearest simple forms 2 (off 1.017%) and √(φ³)
  (off 1.862%), both chosen after the number was known, with no baseline to act
  on.
- **Tolerance.** 0.05 rungs at 2.87 GHz is a 69.9 MHz cone against a 2 MHz
  ensemble line: falsifiable, but satisfied by one frequency in ten.

## 6. The tightened claim

`python3 scripts/phi_tight.py`. Unit fixed to Hz, tolerance 1e-3 rungs. Content
is `2t = 0.2%` per test. Threshold is `τ = ln(1 + δ/f)/ln φ`, the tightest a
measurement can test: a claim is checkable only when `t > τ`.

| quantity | Hz | τ | checkable at 1e-3 |
| --- | --- | --- | --- |
| Cs 133 hyperfine (SI second) | 9.192632e9 | 0 (exact) | yes |
| H 21 cm line | 1.420406e9 | 1.46e-9 | yes |
| Rb 87 hyperfine | 6.834683e9 | 3.04e-10 | yes |
| NV⁻ D | 2.870000e9 | 7.24e-7 | yes |
| ¹⁴N A∥ | 2.16e6 | 1.92e-2 | no |
| ¹³C nearest | 1.30e8 | 1.59e-2 | no |

At 1e-3 rungs: 0 of 6 pass, 0.012 expected, 4 of 6 checkable. At the
framework's 0.05: 0 of 6 pass, 6 of 6 checkable, and every pass is a 10% event.
Content falls as `t` falls while checkability needs `t > τ`, so the informative
window is `τ < t ≪ 1`: inside it the claim is testable and empty, outside it a
pass means nothing.

## 7. Fitting instead of assuming

`python3 scripts/phi_fit.py --fil data/stare2_ovro.fil`, on the three
candidates (gaps 0.1079 and 0.1438 rungs).

- **Offset free, ratio φ fixed** (one parameter): best offset leaves **0.1438
  rungs** where a joint pass needs every gap within 0.10. Fails. Uniform triples
  in the same band pass 0.178 of the time; P(≤ observed residual) = 0.334.
- **Offset and ratio free, integers searched:** integers (3,4) reproduce the
  candidates to **rms 0.0000 rungs**, and reproduce uniform noise to a median
  rms of **0.0012** — forty times inside the 0.05 tolerance. Random triples fit
  that exactly 1.3% of the time.

Two continuous parameters plus two free integers explain three frequencies, and
the ratio is fitted from the same gaps it explains. The fit is an interpolant,
not a law: it cannot be wrong about data it was fitted to, so it cannot
discriminate.

## 8. What would change the answer

1. **Two bands, ≥1 decade apart, ≥3 resolved drifting lines each.** The
   pre-registration in `prereg/phi_2band.json` fixes the statistic, the two
   predictions and the 5% rule. Either outcome is reportable; the script prints
   which statement survives.
2. **NV spectroscopy at 1e-3 rungs**, which needs ~19 K temperature control and
   a single-centre linewidth. There the claim is checkable and a pass is a
   1-in-500 event.
3. **A stated unit.** Any frequency can be placed on a rung by choosing Hz,
   MHz, eV or nm; the framework must name one before an address means anything.

## Reproduce

```bash
python3 scripts/phi_read.py --fil data/stare2_ovro.fil
python3 scripts/phi_prereg.py --fil data/stare2_ovro.fil
python3 scripts/phi_hunt.py data/stare2_ovro.fil data/demo.fil data/fine_on.fil
python3 scripts/nv_check.py
python3 scripts/phi_tight.py
python3 scripts/phi_fit.py --fil data/stare2_ovro.fil
python3 scripts/phi_check.py && python3 scripts/phi_native.py && python3 scripts/phi_superpose.py
```

Earlier harnesses, for completeness: the coordinate round trip reproduces to
1.1e-15 relative, one baseline is 435.2 MHz, ±4 Hz/s needs 4.71e-05° of angle,
the sky exponent is 1.0000 against the pitch's 0.306 (4.9×/decade), folding the
stack gives arc null R = 0.7972 against 0.9276 with circle R = 0.1232
(p = 0.955), and the superposition finds 3 false alarms against 2.8 expected.

## What this is not

Not a detection, and not a refutation of φ as mathematics. The rung coordinate
is a genuine logarithm; the golden angle is a genuine constant; the tetrahedral
angle really is arccos(−1/3). What these tests show is narrower and more
useful: the framework's *addresses* are unprivileged, its *tolerance* is either
too loose to fail or too tight to pass, and its *scale* is free in a way that
removes the test. The drift pipeline is untouched by all of this and stands on
its own measurements.
