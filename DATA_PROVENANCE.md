# Data provenance

## Source

Public, anonymous, no credentials.

```
https://ovro-lwa-solar.s3.us-west-2.amazonaws.com/?list-type=2&prefix=spec_fits/
```

- **Instrument:** OVRO-LWA solar observing mode, Stokes I dynamic spectra.
- **Product:** `spec_fits/YYYY/YYYYMMDD.fits`, one primary HDU per day.
- **Licence:** CC BY 4.0 as published by the OVRO-LWA solar team.
- **Access:** `scripts/fetch_band.py` / `ensure_file()`, cached in `data/`.
- **Integrity:** SHA-256 of every cached file in
  [`prereg/data_checksums.txt`](prereg/data_checksums.txt). Verify with
  `cd data && sha256sum -c ../prereg/data_checksums.txt`.

## What was used

| file | size | days | samples | role |
| --- | --- | --- | --- | --- |
| `20240207.fits` | 82.1 MB | 2024-02-07 | 7009 | published pair reproduction |
| `20240322.fits` | 32.9 MB | 2024-03-22 | 2802 | 47.242 MHz line |
| `20240324.fits` | 46.5 MB | 2024-03-24 | 3971 | 47.314 MHz line |
| `20240325.fits` | 60.2 MB | 2024-03-25 | 5140 | four-lane solar group |
| `20240326.fits` | 60.2 MB | 2024-03-26 | 5140 | no candidates |
| `20240920.fits` | 41.1 MB | 2024-09-20 | 3504 | near-Nyquist trio |

Selection: every 2024 day inside a **20-140 MB** size window, which is 6 of 347
listed keys. That is a size-filtered subset, not a complete survey of the year,
and every claim about "the archive" means "these six days".

## Structure

- axes `(time, freq, 1, Stokes)`; layout inferred from axis sizes because the
  header carries **no CTYPE/CRVAL/CDELT**
- 731 channels, width `100547.19562243501` Hz, span 73.5 MHz
- the survey searched 96 of those 731 channels per day (a 9.7 MHz window at the
  band centre, channels 317-412); `scripts/expanded_search.py` searches all 731,
  and `prereg/expanded_search*.json` records which variant produced which rows
- cadence `0.256` s; the largest file is `0.50` h
- one **disk-integrated beam per file**: there is no off-source product, no
  second beam, and no second instrument in this archive

Barycentric correction is undefined without WCS and source coordinates. The build
quotes a bound instead: a few Hz/s, computed in `scripts/ovro_confirm.py`.

## Derived artifacts

| artifact | produced by | meaning |
| --- | --- | --- |
| `prereg/ovro_survey.json` | `scripts/ovro_survey.py` | 10 candidates, band-wide noise |
| `prereg/local_noise.json` | `scripts/local_noise_survey.py` | same 10, local noise |
| `prereg/ovro_confirm.json` | `scripts/ovro_confirm.py` | gates, corrections, injections |
| `prereg/rfi_pair.json` | `scripts/rfi_pair_test.py` | persistence and Doppler bounds |
| `nonlinear_candidates.json` | `scripts/nonlinear_time.py` | five coordinates, controls |
| `signal_records.json` | `scripts/build_signal_records.py` | typed join of the above |
| `prereg/regression_narrow_window.json` | manual, from a real run | the invalid 0/10 fixture |

## The φ corpus

`phi-quantum-pulses` contributes text corpora (`b64*.txt`, `chunks*/`), a
print-based constants demonstration (`phi_scroll_decode.py`), a Torch I3D model
(`i3dallnl.py`, requiring PyTorch, which this project does not use), and Vesuvius
scroll outputs. None of it is measurement data about the sky, and none of it is
used to produce a number in this pipeline. Its identity and licensing records
(`identity.json`, `LICENSING.md`) are preserved as published.

Its "all 99 blocks verified with φ parameters, zero failures" claim is a
consistency check of arithmetic constants against a parameter set. It is not a
detection and is not treated as one; see [FAILURE_MODES.md](FAILURE_MODES.md).

## Reproducing

```bash
python3 scripts/ovro_survey.py && python3 scripts/ovro_confirm.py
python3 scripts/local_noise_survey.py && python3 scripts/nonlinear_time.py
python3 scripts/rfi_pair_test.py && python3 scripts/build_signal_records.py
python3 scripts/ingest_icecube.py         # published alerts -> data/neutrino_events.json
python3 scripts/neutrino_crossmatch.py   # SYNTHETIC catalog unless data/neutrino_events.json exists
python3 scripts/neutrino_power.py && python3 scripts/nv_lane.py
python3 scripts/phi_cross_domain_test.py \
  --dataset survey_cells=data/phi_cross_domain/survey_cells.csv \
  --dataset group_20240325=data/phi_cross_domain/group_20240325.csv \
  --dataset trio_20240920=data/phi_cross_domain/trio_20240920.csv \
  --out results/phi_cross_domain --permutations 100000
python3 scripts/phi_cross_domain_report.py && python3 scripts/super_tetrahedron.py
python3 -m unittest discover -s tests -q
```

Network is needed only for the first fetch; everything after that reads `data/`.
