"""Preregister, then run, the whole-band search over every available OVRO day.

The survey searched a 96-channel window in the middle of each file. Every file is
731 channels wide - 73.5 MHz of 13.4-86.9 MHz - so 87% of the frequency coverage
had never been dedopplered. This searches all of it, over every integration in
each file, and it is preregistered before the first run so the window cannot be
chosen after seeing what is in it.

The preregistration records the configuration and the decision rules; the run
records what came out. Both are written to `prereg/`, and `--prereg` refuses to
overwrite an existing declaration, so a later attempt to re-declare after looking
fails loudly instead of quietly rewriting history.

    python3 scripts/expanded_search.py --prereg          # declare, once
    python3 scripts/expanded_search.py --days 20240207   # search, repeatable
"""

from __future__ import annotations

import argparse
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CHANNEL_HZ = 100547.19562243501
TSAMP_S = 0.256
CHANNELS = 731
BAND_MHZ = (13.4, 86.9)
THRESHOLD = 5.5
BOXCAR = 3
GROUP_RATIO = 1.05
PREREG = os.path.join(ROOT, "prereg", "expanded_prereg.json")
PREREG_BC = os.path.join(ROOT, "prereg", "expanded_prereg_bc.json")
VARIANTS = {
    # frequency coverage x integration count. A is the declared whole-band,
    # all-integration search; B and C keep one axis fixed so the two expansions
    # are not confounded with each other.
    "A": {"channels": CHANNELS, "integrations": 0, "window": "full band"},
    "B": {"channels": CHANNELS, "integrations": 1000, "window": "full band, survey duration"},
    "C": {"channels": 96, "integrations": 0, "window": "survey window, all integrations"},
}
DAYS = ("20240207", "20240322", "20240324", "20240325", "20240326", "20240920")


def declare() -> int:
    if os.path.exists(PREREG):
        print(f"refusing to overwrite {PREREG}: a preregistration is declared once")
        return 1
    payload = {
        "declared": "2026-11-12",
        "purpose": (
            "expand frequency coverage from the survey's 96-channel window to all "
            "731 channels of every available day, and from 1000 samples to every "
            "integration, before inspecting anything in the new coverage"
        ),
        "files": [f"data/{day}.fits" for day in DAYS],
        "band_mhz": list(BAND_MHZ),
        "channels": CHANNELS,
        "channel_hz": CHANNEL_HZ,
        "tsamp_s": TSAMP_S,
        "integrations": "every integration present in each file",
        "boxcar": BOXCAR,
        "threshold_sigma": THRESHOLD,
        "drift_range": "Nyquist, 0.5 channel per sample",
        "drift_step": "library default: half a boxcar of smear over the integration",
        "window_guard": "seti_search.stabilizer, 20% significance tolerance",
        "grouping_rule": (
            "sort a day's candidates by frequency; chain neighbours whose ratio is "
            "at most GROUP_RATIO; any chain of three or more is a group to fit"
        ),
        "group_ratio": GROUP_RATIO,
        "models": {
            "comb": "f_n = f_0 + n*df (two parameters)",
            "scale": "f_n = f_0 * q**n (two parameters)",
        },
        "verdict_rules": {
            "span": (
                "a group spanning under one phi rung cannot return SCALE FAMILY "
                "SUPPORTED or FREQUENCY COMB SUPPORTED; at most BOTH MODELS FIT "
                "or INDETERMINATE DUE TO LIMITED SPAN"
            ),
            "replication": (
                "SCALE FAMILY SUPPORTED additionally needs the same q on at least "
                "two separated groups or dates"
            ),
            "decisive": "delta AIC > 10 and the same model wins held-out prediction",
            "noise": "RMS residual above half a channel means neither model fits",
        },
        "no_inspection_before_declaration": True,
    }
    with open(PREREG, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"declared {PREREG}")
    return 0


def declare_bc() -> int:
    if os.path.exists(PREREG_BC):
        print(f"refusing to overwrite {PREREG_BC}")
        return 1
    payload = {
        "declared": "2026-11-12",
        "declared_after_variant_A": True,
        "why": (
            "variant A changed both frequency coverage and integration count at "
            "once, and the window guard showed the two cannot be separated from "
            "its result: A holds 21-52% of the Nyquist drift window, so candidates "
            "that vanish in A could have been lost to duration, to bandwidth, or "
            "to the drift range. B and C hold one axis at the survey's setting so "
            "each expansion is measured alone. Declaring them after seeing A is "
            "recorded here rather than hidden; neither variant's threshold, band "
            "or drift range was chosen from A's candidates."
        ),
        "variants": {
            "B": "all 731 channels, first 1000 integrations (256 s), Nyquist drift",
            "C": "survey's 96-channel centre window (317:413), every integration",
        },
        "held_constant": {
            "threshold_sigma": THRESHOLD, "boxcar": BOXCAR,
            "channel_hz": CHANNEL_HZ, "drift_range": "Nyquist",
        },
        "grouping_rule": "as declared in expanded_prereg.json",
    }
    with open(PREREG_BC, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write(chr(10))
    print(f"declared {PREREG_BC}")
    return 0


def search_day(day: str, variant: str = "A") -> dict:
    import time

    import numpy as np

    from ovro_lwa_population import dynamic
    from ovro_lwa_test import TSAMP_S as IMPLIED_TSAMP
    from seti_search.rungs import rung, rung_drift
    from seti_search.search import SearchConfig, search_band
    from seti_search.stabilizer import check_window

    path = os.path.join(ROOT, "data", f"{day}.fits")
    spec = VARIANTS[variant]
    data, meta = dynamic(path)
    block = data.T.astype(np.float32)  # (time, channel)
    n_f = block.shape[1]
    c0 = max(0, (n_f - spec["channels"]) // 2)
    block = block[:, c0: c0 + spec["channels"]]
    if spec["integrations"]:
        block = block[: spec["integrations"]]
    n_t, n_f = block.shape
    freqs_mhz = meta.get("low_mhz", BAND_MHZ[0]) + np.arange(c0, c0 + n_f) * meta["channel_hz"] / 1e6

    cfg = SearchConfig(boxcar=BOXCAR, max_drift_hz_s=0.5 * meta["channel_hz"] / IMPLIED_TSAMP,
                       snr_threshold=THRESHOLD, max_hits=400)
    started = time.time()
    hits, rfi = search_band(block, freqs_mhz, IMPLIED_TSAMP, cfg)
    elapsed = time.time() - started
    window = check_window(cfg.max_drift_hz_s, n_f, n_t, meta["channel_hz"], IMPLIED_TSAMP)

    rows = []
    for hit in hits:
        freq_hz = hit.freq_mhz * 1e6
        rows.append({
            "date": day,
            "file": f"data/{day}.fits",
            "freq_mhz": round(hit.freq_mhz, 6),
            "freq_hz": freq_hz,
            "drift_hz_s": round(hit.drift_hz_s, 3),
            "snr": round(hit.snr, 2),
            "channel_index": int(hit.channel),
            "channel_offset": round(
                (hit.freq_mhz - freqs_mhz[0]) / meta["channel_hz"] * 1e6
                - round((hit.freq_mhz - freqs_mhz[0]) / meta["channel_hz"] * 1e6), 4),
            "bandwidth_hz": n_f * meta["channel_hz"],
            "channel_hz": meta["channel_hz"],
            "duration_s": n_t * IMPLIED_TSAMP,
            "on_source_status": "on-source",
            "off_source_status": "unavailable: one disk-integrated beam",
            "phi_freq_rung": round(rung(freq_hz), 6),
            "phi_freq_nearest": int(round(rung(freq_hz))),
            "phi_freq_drift": round(rung_drift(freq_hz), 6),
            "phi_drift_rung": (round(rung(abs(hit.drift_hz_s)), 6)
                               if hit.drift_hz_s else None),
        })
    del data, block
    return {
        "date": day,
        "variant": variant,
        "window": spec["window"],
        "first_channel": c0,
        "n_channels": n_f,
        "n_integrations": n_t,
        "duration_s": n_t * IMPLIED_TSAMP,
        "seconds": round(elapsed, 1),
        "n_candidates": len(rows),
        "blanked_channels": rfi["blanked_channels"],
        "snr_calibration": round(float(rfi.get("snr_calibration", 1.0)), 4),
        "window_check": window.as_dict(),
        "candidates": rows,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prereg", action="store_true")
    ap.add_argument("--days", nargs="*", default=None)
    ap.add_argument("--variant", default="A", choices=sorted(VARIANTS))
    ap.add_argument("--prereg-bc", action="store_true")
    ap.add_argument("--max-ints", type=int, default=0,
                    help="cap integrations (0 = every integration, as declared)")
    args = ap.parse_args(argv)
    if args.prereg:
        return declare()
    if args.prereg_bc:
        return declare_bc()

    result_path = os.path.join(
        ROOT, "prereg", f"expanded_search{'' if args.variant == 'A' else '_' + args.variant}.json")
    existing = {"days": [], "config": {}}
    if os.path.exists(result_path):
        with open(result_path, encoding="utf-8") as fh:
            existing = json.load(fh)
    done = {d["date"]: d for d in existing.get("days", [])}
    wanted = args.days or [d for d in DAYS if d not in done]

    print(f"EXPANDED SEARCH variant {args.variant}  ({len(wanted)} day(s), "
          f"{VARIANTS[args.variant]['channels']} channels, "
          f"{VARIANTS[args.variant]['window']}, threshold {THRESHOLD:g} sigma)")
    for day in wanted:
        print(f"   {day} ...", flush=True)
        entry = search_day(day, args.variant)
        if args.max_ints:
            entry["note"] = f"integrations capped at {args.max_ints} (not as declared)"
        done[day] = entry
        print(f"   {day}: {entry['n_integrations']} integrations, "
              f"{entry['n_candidates']} candidates, {entry['seconds']}s, "
              f"window {entry['window_check']['fraction_held']:.1%} held "
              f"({entry['window_check']['verdict']})", flush=True)

    payload = {
        "variant": args.variant,
        "window": VARIANTS[args.variant]["window"],
        "prereg": os.path.basename(PREREG if args.variant == "A" else PREREG_BC),
        "channel_hz": CHANNEL_HZ, "threshold": THRESHOLD, "boxcar": BOXCAR,
        "days_searched": len(done),
        "days": [done[d] for d in DAYS if d in done],
    }
    with open(result_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    total = sum(len(d["candidates"]) for d in payload["days"])
    print(f"   wrote {result_path}: {len(payload['days'])} days, {total} candidates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
