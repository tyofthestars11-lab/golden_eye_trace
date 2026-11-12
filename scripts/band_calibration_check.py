"""Why the same candidates appear in a narrow window and vanish in a wide band.

The survey dedopplered 96 channels and found candidates in 45.4-54.5 MHz. The
expanded search dedopplers all 731 channels at the same threshold and integration
count and finds none of them. Nothing about the data changed, so the difference is
in the processing, and the only processing that depends on the map's extent is the
noise scale: `calibrate_map` measures one spread for the whole map and divides
every cell by it.

This script measures that, at the known candidates' own cells, in both windows, and
reports the raw and calibrated significance side by side. It does not decide which
calibration is right - the record already holds a rejected attempt at per-trial
calibration - it records the size of the effect so the paper trail states it.

    python3 scripts/band_calibration_check.py
"""

from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

from ovro_lwa_population import dynamic  # noqa: E402
from ovro_lwa_test import TSAMP_S  # noqa: E402
from seti_search.search import (  # noqa: E402
    SearchConfig,
    calibrate_map,
    control_bandpass,
    dedoppler,
    excise_rfi,
)

# frequency of every candidate the narrow-window survey reported, by day
KNOWN = {
    "20240325": [45.9843, 46.9828, 47.9380, 48.9937],
    "20240920": [45.4206, 46.5890, 47.5358],
    "20240207": [54.4530, 55.9010],
}
NARROW = 96
INTS = 1000


def map_for(block, df_mhz, cfg):
    normalised, _ = control_bandpass(block)
    clean, _ = excise_rfi(normalised, sigma=cfg.rfi_sigma, carrier_sigma=cfg.carrier_sigma)
    snr_map, drifts = dedoppler(clean, TSAMP_S, df_mhz, max_drift_hz_s=cfg.max_drift_hz_s,
                                drift_step_hz_s=cfg.drift_step_hz_s, boxcar=cfg.boxcar)
    live = ~(clean == 0).all(axis=0)
    return calibrate_map(snr_map, live), drifts, live


def cell(snr_map, drifts, freqs_mhz, channel_hz, freq_mhz):
    ch = int(round((freq_mhz - freqs_mhz[0]) * 1e6 / channel_hz))
    row = int(np.argmax(snr_map[:, ch]))
    return ch, row, float(drifts[row]), float(snr_map[row, ch])


def main() -> int:
    out = {"narrow_channels": NARROW, "integrations": INTS, "days": {}}
    for day, known_mhz in KNOWN.items():
        data, meta = dynamic(os.path.join(ROOT, "data", f"{day}.fits"))
        block_all = data.T.astype(np.float32)[:INTS]
        del data
        ch_hz = meta["channel_hz"]
        base = meta.get("low_mhz", 13.4)
        c0 = max(0, (block_all.shape[1] - NARROW) // 2)
        cfg = SearchConfig(boxcar=3, max_drift_hz_s=0.5 * ch_hz / TSAMP_S,
                           snr_threshold=5.5, max_hits=400)

        rows = {}
        for name, block, first in (("narrow", block_all[:, c0: c0 + NARROW], c0),
                                  ("wide", block_all, 0)):
            freqs = base + np.arange(first, first + block.shape[1]) * ch_hz / 1e6
            (cal, factor), drifts, live = map_for(block, ch_hz / 1e6, cfg)
            rows[name] = {
                "calibration_factor": round(factor, 4),
                "channels": int(block.shape[1]),
                "cells": {},
            }
            for f in known_mhz:
                if not (freqs[0] <= f <= freqs[-1]):
                    rows[name]["cells"][f"{f:.4f}"] = {"outside_window": True}
                    continue
                ch, row, rate, value = cell(cal, drifts, freqs, ch_hz, f)
                rows[name]["cells"][f"{f:.4f}"] = {
                    "channel": ch, "drift_hz_s": round(rate, 1),
                    "significance": round(value, 2),
                    "detected": bool(value >= 5.5),
                    "live": bool(live[ch]),
                }
        out["days"][day] = rows
        del block_all
        print(f"{day}: narrow scale {rows['narrow']['calibration_factor']:.3f}  "
              f"wide scale {rows['wide']['calibration_factor']:.3f}")
        for f in known_mhz:
            n = rows["narrow"]["cells"][f"{f:.4f}"]
            w = rows["wide"]["cells"][f"{f:.4f}"]
            if n.get("outside_window") or w.get("outside_window"):
                print(f"   {f:9.4f} MHz  outside one of the windows, skipped")
                continue
            print(f"   {f:9.4f} MHz  narrow {n['significance']:6.2f} "
                  f"({'hit' if n['detected'] else '   '})   wide {w['significance']:6.2f} "
                  f"({'hit' if w['detected'] else '   '})   drift {n['drift_hz_s']:.0f} Hz/s")

    path = os.path.join(ROOT, "prereg", "band_calibration.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write(chr(10))
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
