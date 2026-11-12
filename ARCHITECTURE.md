# Architecture

Two repositories, one pipeline, and six layers that are never merged into a
single number.

| repository | what it actually contains | role here |
| --- | --- | --- |
| `golden_eye_trace` (this tree) | a numpy-and-stdlib SETI search: OVRO-LWA ingest, dedoppler, φ ledger, confirmation gates | the measurement system |
| `phi-quantum-pulses` | text corpora, a print-based φ-constants demonstration, a Torch I3D model, Vesuvius scroll outputs | source of φ vocabulary and worked examples; **not** a measurement pipeline |

The second repository is imported as *concepts and constants*, never as
detection code. Its own claims are audit findings, not results: see
[FAILURE_MODES.md](FAILURE_MODES.md) for which of them are testable and which are
not.

## Layers

Each layer answers one question, writes its own artifact, and is checkable
without the layers above it.

| # | layer | question | artifact |
| --- | --- | --- | --- |
| 1 | raw measurement | what did the instrument record, where, when | `prereg/ovro_survey.json`, `signal_records.json` |
| 2 | calibrated significance | how far above the *local* noise | `prereg/local_noise.json`, `LOCAL_NOISE.md` |
| 3 | nonlinear event time | is it coherent in a coordinate other than sample time | `nonlinear_candidates.json`, `NONLINEAR_TIME.md` |
| 4 | φ-coordinate ranking | where does it sit on the golden ladder | `prereg/phi_lock.json`, `PHI_FALSIFICATION.md` |
| 5 | physical interpretation | which mechanism explains it | `CANDIDATES.md`, `CONFIRMATION.md` |
| 6 | independent confirmation | did anyone else see it | gate matrix in `CONFIRMATION.md` |

A layer may restrict what a lower one is allowed to claim. Layer 4 may only ever
*detract* from a candidate; it can never promote one. Layer 6 is the only layer
allowed to use the word "confirmed", and in this archive it always returns
UNTESTABLE, because there is one beam, no off-source pointing, and one telescope.

## Data flow

```
spec_fits/YYYY/YYYYMMDD.fits        public bucket, CC BY 4.0, no WCS
        |
        |  dynamic()               (nfreq, ntime) Stokes-I float64
        v
  band window 96 ch x 1000 samples   survey configuration
        |
        |  control_bandpass -> excise_rfi -> dedoppler -> calibrate_map
        v
  dedoppler map (drift x channel)     snr_map, drifts
        |
        +--> extract_hits .......... candidates at 5.5 sigma      layer 1
        |
        +--> local_scale ........... per-channel noise, +/-16 ch  layer 2
        |
        +--> track_at + transforms . phase/chirp/accel/energy     layer 3
        |
        +--> phi_features .......... k >= 1 lock distance          layer 4
        |
        +--> gates ................. repeat/off-source/barycentric layer 5-6
        v
  signal_records.json                 ten typed records
```

## Script map

| script | input | transform | output | claim it is allowed to make |
| --- | --- | --- | --- | --- |
| `scripts/fetch_band.py` | bucket listing | download + cache | `data/*.fits` | none; transport only |
| `scripts/ovro_lwa_test.py` | one day file | FITS read, band window, dedoppler | `prereg/ovro_lwa.json` | a line exists in that window |
| `scripts/ovro_lwa_population.py` | day listing | multi-day dedoppler | `prereg/ovro_lwa_population.json` | population statistics |
| `scripts/ovro_survey.py` | 6 days | one configuration across days | `prereg/ovro_survey.json` | candidates at 5.5σ, band-wide noise |
| `scripts/ovro_confirm.py` | survey rows | reproduce, correct, gate, inject | `CONFIRMATION.md`, `CANDIDATES.md` | dispositions, never confirmation |
| `scripts/local_noise_survey.py` | survey rows | per-channel local scale | `LOCAL_NOISE.md` | significance under local noise |
| `scripts/nonlinear_time.py` | confirm rows | 5 time coordinates, split-half | `NONLINEAR_TIME.md` | coherence vs a warp null |
| `scripts/rfi_pair_test.py` | 2 day files | full-integration + window scan + Doppler bounds | `RFI_PAIR.md` | repeat or not-repeat |
| `scripts/build_signal_records.py` | the three JSONs | join, validate | `signal_records.json/.csv` | nothing new; it only joins |
| `scripts/phi_echo.py` etc. | band samples | φ lock, echo, geometry | `prereg/phi_*.json` | feature values and their nulls |
| `seti_search/search.py` | waterfall | dedoppler | SNR map, hits | significance of a cell |
| `seti_search/signal_record.py` | — | validation | typed record | none; it constrains claims |

## Invariants

1. **Two significances, never one.** `snr_global` is the band-wide reference the
   first survey used; `snr_at_cell` is against the candidate's own neighbourhood.
   Both are stored, with `local_scale` between them, so the calibration is visible.
2. **φ is a feature with a control beside it.** Any φ column is accompanied by the
   randomized-phase null from the same run.
3. **Every rejection is written down.** `seti_search/regressions.py` holds the
   rejected methods with the tell that exposed them.
4. **The verdict vocabulary is closed.** `DISPOSITIONS` in
   `seti_search/signal_record.py`; anything else fails validation.
5. **No layer may promote a candidate past an untestable gate.** In practice that
   means the headline stays NULL RESULT.

## Running it

```bash
python3 scripts/ovro_survey.py          # needs data/ and the bucket
python3 scripts/ovro_confirm.py
python3 scripts/local_noise_survey.py
python3 scripts/nonlinear_time.py
python3 scripts/rfi_pair_test.py
python3 scripts/build_signal_records.py
python3 -m unittest discover -s tests -q   # 120 tests, ~23 s, no network
```
