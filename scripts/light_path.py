"""Follow the measured field through the whole path, and find where it changes.

A candidate's sigma is a product of eight transformations, and this records the
number each one produces, so a change of significance can be attributed to a stage
instead of to "the pipeline".

The finding that made this necessary: `excise_rfi`'s row test compares each
integration's band median against `1.25/sqrt(n_f) * typical_std`, a noise figure
that assumes independent channels. On a real solar product the channel gains are
correlated, so row medians scatter far more than that, and a two-sided threshold at
6 sigma of that value flags half the integrations. Blanking half the integrations -
both those brighter *and* those fainter than the median - leaves a subset selected
on the very thing being searched for, and the dedoppler sum of that subset is a
different, much quieter statistic.

Measured here at the candidates' own drift, in the survey's own window:

  20240325   no mask  1.44  0.54 -0.15 -0.87      (0 rows blanked)
             theory   954.18 604.16 447.18 171.27 (505 of 1000 rows blanked)
             empirical  1.44  0.54 -0.15 -0.87      (0 rows blanked)

The empirical test is the same one scaled by the measured scatter of the row
medians, which is what "one integration far brighter than the rest" has to be
measured against. With it, nothing at these frequencies reaches 1.5 sigma.

Every stage is recorded with whether it can be inverted, and both Stokes planes are
carried because the product ships two and the earlier records kept only one.

    python3 scripts/light_path.py
"""

from __future__ import annotations

import csv
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import numpy as np  # noqa: E402

from ovro_lwa_population import read_primary  # noqa: E402
from ovro_lwa_test import TSAMP_S  # noqa: E402
from seti_search.rungs import rung, rung_drift  # noqa: E402
from seti_search.search import (  # noqa: E402
    SearchConfig,
    calibrate_map,
    control_bandpass,
    dedoppler,
    excise_rfi,
    extract_hits,
)

CHANNEL_HZ = 100547.19562243501
LOW_MHZ = 13.4
BOX = 3
NYQUIST = 0.5 * CHANNEL_HZ / TSAMP_S
INTS = 1000
WINDOWS = {"narrow-96": 96, "wide-731": 731}
MASKS = {
    "none": (1e9, 1e9, "theory"),
    "theory": (6.0, 20.0, "theory"),
    "empirical": (6.0, 20.0, "empirical"),
}
TARGETS = {
    "20240325": {"mhz": [45.9843, 46.9828, 47.9380, 48.9937], "drift": 41240.0},
    "20240920": {"mhz": [45.4206, 46.5890, 47.5358], "drift": 195007.0},
}


def mad_sigma(values) -> float:
    flat = np.asarray(values, dtype=np.float64).ravel()
    med = float(np.median(flat))
    return float(1.4826 * np.median(np.abs(flat - med)))


def read_planes(path):
    header, cube = read_primary(path)
    plane0 = np.nan_to_num(cube[:, :, 0, 0].astype(np.float32), nan=0.0)
    plane1 = np.nan_to_num(cube[:, :, 0, 1].astype(np.float32), nan=0.0)
    return plane0, plane1, {"nfreq": plane0.shape[1], "ntime": plane0.shape[0],
                            "stokes_planes": int(cube.shape[3])}


def trail_sum(block, channel, drift_hz_s):
    """Aligned sum along the cell's own trail - what the dedoppler forms."""
    shift = drift_hz_s * TSAMP_S / CHANNEL_HZ
    n_t = block.shape[0]
    total, count = 0.0, 0
    for t in range(n_t):
        j = int(round(channel - shift * (n_t - 1 - t)))
        if 0 <= j < block.shape[1]:
            total += float(block[t, j])
            count += 1
    return total / max(count, 1)


def trace(day, cells, drift, plane0, plane1, n_ch, c0):
    block = plane0[:INTS, c0:c0 + n_ch].astype(np.float32)
    p1 = plane1[:INTS, c0:c0 + n_ch].astype(np.float32)
    normalised, _ = control_bandpass(block)
    band_level = float(np.median(normalised))
    freqs = LOW_MHZ + (c0 + np.arange(n_ch)) * CHANNEL_HZ / 1e6
    out = {}
    for mask, (sig, car, rn) in MASKS.items():
        clean, rfi = excise_rfi(normalised, sigma=sig, carrier_sigma=car, row_noise=rn)
        snr_map, drifts = dedoppler(clean, TSAMP_S, CHANNEL_HZ / 1e6,
                                    max_drift_hz_s=NYQUIST, drift_step_hz_s=None,
                                    boxcar=BOX)
        live = ~(clean == 0).all(axis=0)
        cal, factor = calibrate_map(snr_map, live)
        row = int(np.argmin(np.abs(drifts - drift)))
        hits = extract_hits(snr_map, drifts, freqs, threshold=5.5, boxcar=BOX,
                            total_time=block.shape[0] * TSAMP_S, max_hits=400)
        cells_out = []
        for freq_mhz in cells:
            ch = int(round((freq_mhz - LOW_MHZ) * 1e6 / CHANNEL_HZ)) - c0
            neighbourhood = clean[:, max(0, ch - 8): ch + 9]
            local = mad_sigma(neighbourhood)
            global_sigma = mad_sigma(clean)
            cells_out.append({
                "date": day, "mask": mask, "window": f"{n_ch}ch",
                "freq_mhz": round(freq_mhz, 6),
                "freq_hz": round(freq_mhz * 1e6, 3),
                "channel_index": ch + c0,
                "channel_offset": round((freq_mhz - LOW_MHZ) * 1e6 / CHANNEL_HZ - c0 - ch, 4),
                "channel_hz": round(CHANNEL_HZ, 4),
                "drift_hz_s": round(drift, 1),
                "drift_trial": row,
                "n_drift_trials": int(len(drifts)),
                "A_raw_power": round(trail_sum(block, ch, drift), 10),
                "B_over_band_level": round(trail_sum(normalised, ch, drift)
                                           / band_level if band_level else 0.0, 6),
                "C_sigma_map_units": round(float(snr_map[row, ch]), 2),
                "D_channel_masked": bool(not live[ch]),
                "D_blanked_channels": int(rfi["blanked_channels"]),
                "D_blanked_samples": int(rfi["blanked_rows"]),
                "D_total_samples": int(clean.shape[0]),
                "D_row_test": rfi["row_noise"],
                "E_sigma_calibrated": round(float(cal[row, ch]), 2),
                "F_calibration_factor": round(float(factor), 4),
                "G_local_scale": round(global_sigma / local, 4) if local > 0 else None,
                "G_sigma_local": (round(float(cal[row, ch]) * global_sigma / local, 2)
                                  if local > 0 else None),
                "H_recovered": any(abs(h.freq_mhz - freq_mhz) < 0.01 for h in hits),
                "polarization_plane_ratio": (
                    round(trail_sum(p1, ch, drift) / max(trail_sum(block, ch, drift), 1e-30), 4)),
                "phi_freq_rung": round(rung(freq_mhz * 1e6), 6),
                "phi_drift_rung": round(rung_drift(drift), 6),
            })
        out[mask] = {"row_test": rn, "cells": cells_out,
                     "blanked_rows": rfi["blanked_rows"],
                     "blanked_channels": rfi["blanked_channels"],
                     "calibration_factor": float(factor)}
    return out


def two_by_two():
    """The declared design, four cells, read from the artifacts."""
    rows = []
    sources = {
        ("narrow-96", "short-1000"): "signal_records.json",
        ("narrow-96", "full"): "prereg/expanded_search_C.json",
        ("wide-731", "short-1000"): "prereg/expanded_search_B.json",
        ("wide-731", "full"): "prereg/expanded_search.json",
    }
    for (window, duration), path in sources.items():
        full = os.path.join(ROOT, path)
        if not os.path.exists(full):
            continue
        with open(full, encoding="utf-8") as fh:
            doc = json.load(fh)
        if "records" in doc:
            cands = len(doc["records"])
            held = None
            ints = 1000
        else:
            cands = sum(d["n_candidates"] for d in doc["days"])
            held = [d["window_check"]["fraction_held"] for d in doc["days"]]
            ints = [d["n_integrations"] for d in doc["days"]]
        rows.append({
            "window": window, "duration": duration, "artifact": path,
            "candidates": cands,
            "integrations": (ints if isinstance(ints, int) else f"{min(ints)}-{max(ints)}"),
            "drift_window_held": ("-" if held is None
                                  else f"{min(held):.0%}-{max(held):.0%}"),
        })
    return rows



def write_report(by_mask, rows, path):
    lines = []
    add = lines.append
    add("# The light path")
    add("")
    add("The measured field, followed from the file to the disposition, one stage at")
    add("a time. The question was where a candidate's significance comes from. The")
    add("answer is not the sky.")
    add("")
    add("## Headline")
    add("")
    add("**The 954 sigma group and the 8.65 sigma trio are both manufactured by the")
    add("excision mask.** At the candidates' own drift, in the survey's own window, the")
    add("same cells read:")
    add("")
    add("| day | no mask | theory row test | empirical row test |")
    add("| --- | --- | --- | --- |")
    for day, windows in by_mask.items():
        d = windows["narrow-96"]
        def fmt(m):
            return " ".join(f"{c['E_sigma_calibrated']:.2f}" for c in d[m]["cells"])
        add(f"| {day} | {fmt('none')} | **{fmt('theory')}** | {fmt('empirical')} |")
    add("")
    add("with integrations blanked:")
    add("")
    add("| day | no mask | theory | empirical |")
    add("| --- | --- | --- | --- |")
    for day, windows in by_mask.items():
        d = windows["narrow-96"]
        add(f"| {day} | {d['none']['blanked_rows']} | {d['theory']['blanked_rows']} "
            f"| {d['empirical']['blanked_rows']} |")
    add("")
    add("The theory test flags half the integrations; the empirical one flags none.")
    add("")
    add("## Why")
    add("")
    add("`excise_rfi`'s row test asks whether an integration's band median is further")
    add("from the median of medians than `sigma * 1.25/sqrt(n_f) * typical_std`. The")
    add("`1.25/sqrt(n_f)` term is the standard error of a median *if the channels were")
    add("independent*. On a real OVRO-LWA solar product the channel gains are")
    add("correlated - the Sun is one source across the band - so the row medians")
    add("scatter far more than that, and a two-sided 6 sigma threshold on the")
    add("theoretical scale flags integrations that are perfectly ordinary. It rejects")
    add("both bright *and* faint rows, keeping a subset selected on the same quantity")
    add("the search then measures, and the dedoppler sum of that subset is a")
    add("different - much quieter - statistic.")
    add("")
    add("The empirical test keeps the structure and replaces the theoretical scale")
    add("with `1.4826 * MAD(row medians)`, which is what \"one integration far")
    add("brighter than the rest\" has to be measured against. It is opt-in")
    add("(`row_noise=\"empirical\"`); the default is unchanged, because flipping it")
    add("would silently rewrite every sigma this repository has published.")
    add("")
    add("## The 2x2")
    add("")
    add("Declared design: frequency span x integration duration, with the drift window")
    add("guard printed beside each cell.")
    add("")
    add("| window | duration | integrations | drift window held | candidates |")
    add("| --- | --- | --- | --- | --- |")
    for row in two_by_two():
        add(f"| {row['window']} | {row['duration']} | {row['integrations']} | "
            f"{row['drift_window_held']} | {row['candidates']} |")
    add("")
    add("Candidate counts in this table are counts produced by the theory-scaled mask,")
    add("so they are an upper bound on what the uncorrected pipeline reports, not a")
    add("measurement of the sky. The wide-band cells truncate the drift window because")
    add("a Nyquist drift needs more channels than the band has at these durations; that")
    add("is the stabilizer's warning, not a null result.")
    add("")
    add("## Every stage, both windows")
    add("")
    add("| day | window | mask | A raw | C sigma (map) | D rows/ch blanked | E sigma | F factor | G local sigma | H recovered |")
    add("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for day, windows in by_mask.items():
        for window, traced in windows.items():
            for mask in ("none", "theory", "empirical"):
                d = traced[mask]
                for c in d["cells"][:1]:
                    add(f"| {day} | {window} | {mask} | {c['A_raw_power']:.3e} | "
                        f"{c['C_sigma_map_units']:.2f} | {d['blanked_rows']}/"
                        f"{d['blanked_channels']} | {c['E_sigma_calibrated']:.2f} | "
                        f"{c['F_calibration_factor']:.3f} | "
                        f"{c['G_sigma_local'] if c['G_sigma_local'] is not None else '-'} | "
                        f"{c['H_recovered']} |")
    add("")
    add("One line per configuration shown at its first cell; every cell is in")
    add("`results/light_path_ledger.csv`, 42 rows.")
    add("")
    add("## What the stages say")
    add("")
    add("- **A to B** (bandpass) is reversible and moves nothing much: the cells sit at")
    add("  the band's own level.")
    add("- **C** (dedoppler, nothing removed) puts every cell inside 1.5 sigma. This is")
    add("  the number that does not depend on a mask.")
    add("- **D** is where the story is. In the narrow window the theory test removes")
    add("  half the integrations and no channels; in the wide window it also removes")
    add("  a quarter to half the channels. Which set survives is decided by the band's")
    add("  brightness at that moment, which correlates with what the Sun was doing.")
    add("- **E to F** (calibration) is a factor of 0.38 to 0.87 - real, but not the")
    add("  effect. The effect is in D's selection, which C never sees.")
    add("- **G** (local scale) is undefined where a neighbourhood has been blanked to")
    add("  nothing, which is itself a diagnostic: on 20240325 the local reference has")
    add("  no live data to measure.")
    add("- **H** (peak extraction) returns the cells under the theory mask and not")
    add("  otherwise.")
    add("")
    add("## Verdict, corrected")
    add("")
    add("The four frequencies at 45.98-48.99 MHz are **not validated candidates**. The")
    add("earlier `INDETERMINATE DUE TO LIMITED SPAN` was a verdict about whether a")
    add("regularity is a comb or a scale family; it presupposed that the four")
    add("frequencies were selected by a valid detection. They were selected by a mask")
    add("that manufactures 954 sigma out of 1.44. Under an honestly scaled mask they")
    add("are at 1.44, 0.54, -0.15 and -0.87 sigma, which is not a detection at all.")
    add("")
    add("The mathematical result stands on its own terms and keeps no evidential")
    add("weight: four frequencies 0.988-1.056 MHz apart are regular at p = 0.002")
    add("against random lines, and a comb and a scale family fit them equally. Four")
    add("numbers chosen by an invalid selection can be regular for the same reason any")
    add("four chosen numbers can.")
    add("")
    add("## Polarization")
    add("")
    add("The product carries two Stokes planes and the reader kept one. Both are now")
    add("carried through the same aligned sum, and the ratio is in the ledger. It is")
    add("reported, not interpreted: without the header's CTYPE cards there is no")
    add("guarantee which two Stokes parameters these are, so a ratio is a measurement")
    add("whose meaning is not established - which is the standard this record holds")
    add("everything else to.")
    add("")
    add("## Reproduce")
    add("")
    add("```")
    add("python3 scripts/light_path.py          # ~11 s, both days, 2x2 x 3 masks")
    add("```")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(chr(10).join(lines) + chr(10))
    return path


def main() -> int:
    rows, by_mask = [], {}
    for day, target in TARGETS.items():
        plane0, plane1, meta = read_planes(os.path.join(ROOT, "data", f"{day}.fits"))
        by_mask[day] = {}
        for window, n_ch in WINDOWS.items():
            c0 = max(0, (meta["nfreq"] - n_ch) // 2)
            traced = trace(day, target["mhz"], target["drift"], plane0, plane1, n_ch, c0)
            by_mask[day][window] = traced
            for mask, doc in traced.items():
                rows.extend(doc["cells"])
        del plane0, plane1
        print(f"{day}: traced {len(WINDOWS) * len(MASKS)} configurations")

    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "light_path_stages.json"), "w",
              encoding="utf-8") as fh:
        json.dump({"channel_hz": CHANNEL_HZ, "nyquist_hz_s": NYQUIST,
                   "configurations": by_mask, "two_by_two": two_by_two()},
                  fh, indent=2, sort_keys=True)
        fh.write(chr(10))
    with open(os.path.join(ROOT, "results", "light_path_ledger.csv"), "w", newline="",
              encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    for day, windows in by_mask.items():
        for window, traced in windows.items():
            print(f"   {day} {window}")
            for mask in ("none", "theory", "empirical"):
                d = traced[mask]
                sig = " ".join(f"{c['E_sigma_calibrated']:8.2f}" for c in d["cells"])
                print(f"      {mask:9s} rows {d['blanked_rows']:4d} "
                      f"ch {d['blanked_channels']:3d}  sigma: {sig}")
    doc = write_report(by_mask, rows, os.path.join(ROOT, "LIGHT_PATH.md"))
    print(f"wrote results/light_path_stages.json and results/light_path_ledger.csv "
          f"({len(rows)} rows)")
    print(f"wrote {os.path.basename(doc)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
