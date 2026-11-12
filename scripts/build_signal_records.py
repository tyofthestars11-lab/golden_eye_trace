"""Build typed SignalRecords from the survey, confirmation and local-noise runs.

The join is the point. Each source answers a different question, and the record
keeps them apart: the survey supplies the raw detection, the local-noise re-run
supplies the calibrated significance in the candidate's own neighbourhood, and
the confirmation build supplies the disposition and the gate matrix. Nothing is
averaged together, and no record is written for a candidate that no run found.

    python3 scripts/build_signal_records.py
"""

from __future__ import annotations

import csv
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from seti_search.signal_record import (  # noqa: E402
    DISPOSITIONS,
    SignalRecord,
    channel_rate,
)

TRANSFORM_VERSION = "ovro-local-noise-v1"
CHANNEL_HZ = 100547.19562243501
TSAMP_S = 0.256


def load(name: str) -> dict:
    path = os.path.join(ROOT, "prereg", name)
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def match(rows: list, date: str, freq_mhz: float) -> dict | None:
    for row in rows:
        if row["date"] == date and abs(row["freq_mhz"] - freq_mhz) < 0.01:
            return row
    return None


def main() -> int:
    survey = load("ovro_survey.json")
    confirm = load("ovro_confirm.json")
    local = load("local_noise.json")

    records: list[SignalRecord] = []
    skipped: list[dict] = []
    for entry in local["candidates"]:
        date, freq = entry["date"], entry["freq_mhz"]
        base = match(survey.get("candidates", []), date, freq)
        conf = match(confirm["candidates"], date, freq)
        if base is None or conf is None:
            skipped.append({"date": date, "freq_mhz": freq,
                            "reason": "no matching row in survey or confirmation"})
            continue
        disposition = conf.get("disposition", "UNRESOLVED")
        # the confirmation build writes free-form parentheticals; the record
        # vocabulary is closed, so map to it rather than store prose
        if disposition.startswith("ASTROPHYSICAL"):
            disposition = "ASTROPHYSICAL"
        elif disposition.startswith("INSTRUMENTAL"):
            disposition = "INSTRUMENTAL ALIAS"
        elif disposition.startswith("TERRESTRIAL"):
            disposition = "TERRESTRIAL RFI"
        elif disposition.startswith("UNRESOLVED"):
            disposition = "UNRESOLVED"
        if disposition not in DISPOSITIONS:
            disposition = "UNRESOLVED"
        records.append(SignalRecord(
            source_file=os.path.join("data", f"{date}.fits"),
            observation_date=date,
            instrument="OVRO-LWA (spec_fits, single disk-integrated beam)",
            frequency_hz=float(freq) * 1e6,
            drift_hz_per_s=float(base["drift_hz_s"]),
            snr_global=float(conf["snr"]),
            snr_at_cell=float(entry["at_candidate"]["local_snr"]),
            local_scale=float(entry["at_candidate"]["scale"]),
            resolution_hz_per_s=float(conf["resolution_hz_s"]),
            channel_rate=channel_rate(float(base["drift_hz_s"]), CHANNEL_HZ,
                                      TSAMP_S),
            duration_s=float(conf["duration_s"]),
            on_source_status="on-source",
            transform_version=TRANSFORM_VERSION,
            disposition=disposition,
            metadata={
                "drift_resolution_hz_s": float(conf["resolution_hz_s"]),
                "survey_threshold_sigma": float(conf["threshold"]),
                "off_source_available": False,
                "barycentric_reducible": False,
                "independent_telescope": False,
                "rate_cycles_per_sample": float(base.get("rate", 0.0)),
                "phi_lock_distance": conf.get("phi", {}).get("lock_distance"),
                "phi_locked": conf.get("phi", {}).get("locked"),
            },
        ))

    payload = {"transform_version": TRANSFORM_VERSION,
               "channel_width_hz": CHANNEL_HZ, "tsamp_s": TSAMP_S,
               "n_records": len(records),
               "skipped": skipped,
               "records": [r.as_dict() for r in records]}
    with open(os.path.join(ROOT, "signal_records.json"), "w",
              encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    with open(os.path.join(ROOT, "signal_records.csv"), "w", newline="",
              encoding="utf-8") as fh:
        fields = [f for f in SignalRecord.__dataclass_fields__ if f != "metadata"]
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow(record.as_dict())

    print(f"SIGNAL RECORDS  ({len(records)} written, {len(skipped)} skipped)")
    for record in sorted(records, key=lambda r: (r.observation_date,
                                                 r.frequency_hz)):
        print(f"   {record.observation_date} {record.frequency_hz / 1e6:8.4f} MHz  "
              f"global {record.snr_global:7.2f}  at-cell {record.snr_at_cell:7.2f}  "
              f"scale {record.local_scale:5.2f}  "
              f"rate {record.channel_rate:+.4f}  {record.disposition}")
    if skipped:
        print("   skipped:")
        for row in skipped:
            print(f"      {row['date']} {row['freq_mhz']:.3f}: {row['reason']}")
    print("   wrote signal_records.json, signal_records.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
