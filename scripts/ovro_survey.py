"""Survey every OVRO-LWA day in the archive's size window, with a full table.

The population run searched three days. This searches all day files that fall in a
size window which is both fetchable inside the time budget and long enough to
carry a drift measurement (20-140 MB: at 12,288 bytes per 0.256 s sample that is
roughly 20 minutes to an hour of data). Each day is searched with the identical
configuration, so the resulting table is comparable row to row.

Every row records the date, frequency, drift rate, SNR, bandwidth, duration,
resolution, threshold and whether the observation has an off-source counterpart.
For `spec_fits` it never does: the product is one disk-integrated beam per file,
which is a hard limit on the confirmation gates and is recorded as such rather
than quietly ignored.

    python3 scripts/ovro_survey.py --samples 1000 --chans 96
"""

from __future__ import annotations

import argparse
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ovro_lwa_population import ensure_file, list_days, search_day  # noqa: E402
from ovro_lwa_test import TSAMP_S  # noqa: E402
from phi_echo import THETA  # noqa: E402

KNOWN = [
    {"date": "20240207", "freq_mhz": 54.452957, "drift_hz_s": 8194.7, "snr": 7.68},
    {"date": "20240207", "freq_mhz": 55.901001, "drift_hz_s": 7866.9, "snr": 7.00},
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--year", default="2024")
    ap.add_argument("--min-mb", type=float, default=20.0)
    ap.add_argument("--max-mb", type=float, default=140.0)
    ap.add_argument("--samples", type=int, default=1000)
    ap.add_argument("--chans", type=int, default=96)
    ap.add_argument("--threshold", type=float, default=5.5)
    ap.add_argument("--max-days", type=int, default=0)
    args = ap.parse_args()

    days = list_days(args.year)
    window = [(k, s) for k, s in days if args.min_mb * 1e6 <= s <= args.max_mb * 1e6]
    window.sort()
    if args.max_days:
        window = window[: args.max_days]
    print("OVRO-LWA MULTI-DAY SURVEY")
    print(f"   {len(days)} files listed for {args.year}; {len(window)} inside "
          f"{args.min_mb:g}-{args.max_mb:g} MB, all searched")
    print(f"   {args.samples} samples x {args.chans} channels, "
          f"{args.threshold:g} sigma, drift window to Nyquist")

    results, failures = [], []
    for index, (key, size) in enumerate(window, 1):
        date = os.path.basename(key).split(".")[0]
        print(f"\n   [{index}/{len(window)}] {date}  {size / 1e6:.1f} MB")
        try:
            day = search_day(ensure_file(key), args.samples, args.chans,
                             args.threshold)
        except Exception as exc:  # noqa: BLE001 - recorded, not swallowed
            print(f"      failed: {exc}")
            failures.append({"key": key, "error": str(exc)})
            continue
        day.update({"key": key, "date": date, "size_mb": size / 1e6,
                    "on_source": True, "off_source": None,
                    "bandwidth_hz": 3.0 * day["channel_hz"],
                    "duration_s": day["samples_used"] * TSAMP_S,
                    "threshold": args.threshold})
        results.append(day)
        print(f"      {day['candidates']} candidates, {day['resolved']} resolved, "
              f"{day['locked']} locked; resolution "
              f"{day['resolution_hz_s']:.1f} Hz/s, theta "
              f"{day['theta_drift_hz_s']:,.0f} Hz/s")
        for row in sorted(day["rows"], key=lambda r: -r["snr"])[:4]:
            print(f"         {row['freq_mhz']:9.4f} MHz  {row['drift_hz_s']:>10.1f} "
                  f"Hz/s  snr {row['snr']:5.2f}  rate {row['rate']:.4f}")

    rows = [dict(row, date=day["date"], resolution_hz_s=day["resolution_hz_s"],
                 channel_hz=day["channel_hz"], duration_s=day["duration_s"],
                 threshold=day["threshold"], on_source=True, off_source=None,
                 bandwidth_hz=day["bandwidth_hz"])
            for day in results for row in day["rows"]]

    print("\n   REPRODUCTION OF THE TWO KNOWN LINES")
    reproduced = []
    for known in KNOWN:
        match = [r for r in rows
                 if r["date"] == known["date"]
                 and abs(r["freq_mhz"] - known["freq_mhz"]) < 0.3
                 and abs(r["drift_hz_s"] - known["drift_hz_s"]) < 400.0]
        ok = bool(match)
        reproduced.append({"known": known, "found": match[:1], "reproduced": ok})
        print(f"      {known['date']} {known['freq_mhz']:.3f} MHz "
              f"{known['drift_hz_s']:,.1f} Hz/s: "
              f"{'reproduced' if ok else 'NOT reproduced at this configuration'}")
    print("      (exact-configuration reproduction is checked separately by "
          "scripts/ovro_confirm.py at 1600 x 192)")

    payload = {"year": args.year, "window_mb": [args.min_mb, args.max_mb],
               "config": {"samples": args.samples, "chans": args.chans,
                          "threshold": args.threshold, "tsamp_s": TSAMP_S,
                          "theta": THETA},
               "days_searched": len(results), "failures": failures,
               "days": [{k: v for k, v in day.items() if k != "rows"}
                        for day in results],
               "candidates": rows, "known": KNOWN, "reproduction": reproduced}
    path = os.path.join(ROOT, "prereg", "ovro_survey.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")

    lines = ["# OVRO-LWA candidate table", "",
             f"Survey of {len(results)} day files in {args.year} between "
             f"{args.min_mb:g} and {args.max_mb:g} MB, each searched at "
             f"{args.samples} samples x {args.chans} channels to Nyquist with a "
             f"{args.threshold:g} sigma threshold. `on-source` is always true and "
             f"`off-source` is always absent: `spec_fits` is one disk-integrated "
             f"beam per file, so no off-source data exist in this product.", "",
             "| date | freq MHz | drift Hz/s | SNR | bandwidth Hz | duration s | "
             "resolution Hz/s | threshold | on-source | off-source |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in sorted(rows, key=lambda r: (r["date"], -r["snr"])):
        lines.append(f"| {r['date']} | {r['freq_mhz']:.4f} | "
                     f"{r['drift_hz_s']:,.1f} | {r['snr']:.2f} | "
                     f"{r['bandwidth_hz'] / 1e3:.1f} | {r['duration_s']:.0f} | "
                     f"{r['resolution_hz_s']:.1f} | {r['threshold']:g} | yes | "
                     f"none available |")
    lines += ["", f"Rows: {len(rows)} candidates across {len(results)} days. "
                  f"Failures: {len(failures)}."]
    with open(os.path.join(ROOT, "CANDIDATES.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"\n   {len(rows)} candidates across {len(results)} days, "
          f"{len(failures)} failures")
    print("   wrote prereg/ovro_survey.json and CANDIDATES.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
