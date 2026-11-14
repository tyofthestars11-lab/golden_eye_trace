"""What a diamond NV-centre sensor could do for this survey, and what it cannot.

    python3 scripts/nv_lane.py

Writes ``NV_LANE.md`` and ``results/nv_lane.json``. The lane is a local-field
witness: a millimetre of diamond at the telescope sees the magnetic part of any
radio wave arriving there, and ``E = c B`` makes the two units interchangeable.
It exists to let a radio candidate be *detracted* by a local emitter - never to
support one.

The numbers are taken from the published literature and cited in the report:

* ensemble NV magnetometry reaches ~210 fT/sqrt(Hz) narrowband and ~460
  fT/sqrt(Hz) broadband in the best reported devices, and ~344 pT/sqrt(Hz) in a
  fully integrated portable unit;
* dynamical-decoupling sequences that reach those figures run from kHz to a few
  MHz, with the upper limit set by the Rabi frequency;
* the NV ground-state splitting is 2.87 GHz and tunes at 2.8 MHz/G, so reaching
  the LWA band (13.4-86.9 MHz) needs a bias field near 100 mT, and the practical
  route there is RF-dressed CW-ODMR with a demonstrated bandwidth of a few MHz.

None of this is a measurement made here. It is a feasibility table, and it says
plainly which of its entries are literature values and which are this survey's.
"""

from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.nv import (  # noqa: E402
    SENSING_MODES,
    band_covered,
    field_from_magnetic,
    magnetic_from_field,
    significance,
    verdict,
    zeeman_bias_t,
)
from seti_search.signal_record import SignalRecord  # noqa: E402

#: The survey's own band, from the data files.
BAND_LOW_HZ = 13.4e6
BAND_HIGH_HZ = 86.9e6

#: Literature values, each with its source. ``band`` is the frequency range the
#: technique is demonstrated over, not a specification anyone promised.
MODE_TABLE = [
    {
        "mode": "cw-odmr",
        "sensitivity_t_sqrt_hz": 1.0e-9,
        "bandwidth_hz": 1.0e3,
        "band": (0.0, 1.0e7),
        "note": "broadband continuous readout, 7.5 MHz demonstrated bandwidth "
                "with RF-dressed states (J. Appl. Phys. 135, 044401, 2024)",
    },
    {
        "mode": "ramsey",
        "sensitivity_t_sqrt_hz": 5.0e-10,
        "bandwidth_hz": 1.0e4,
        "band": (1.0e3, 1.0e5),
        "note": "broadband DC-sensitive mode; 460 fT/sqrt(Hz) reported for the "
                "best ensemble (Phys. Rev. Applied 22, 044069, 2024)",
    },
    {
        "mode": "hahn-echo",
        "sensitivity_t_sqrt_hz": 2.1e-13,
        "bandwidth_hz": 1.0e5,
        "band": (1.0e3, 1.0e6),
        "note": "narrowband AC; 210 fT/sqrt(Hz) reported for the best ensemble "
                "(Phys. Rev. Applied 22, 044069, 2024)",
    },
    {
        "mode": "cpmg",
        "sensitivity_t_sqrt_hz": 2.3e-10,
        "bandwidth_hz": 1.0e5,
        "band": (1.0e4, 1.0e7),
        "note": "229 pT/sqrt(Hz) at 1 MHz with XY4-4 (Sensors 23, 7566, 2023)",
    },
    {
        "mode": "rf-dressed",
        "sensitivity_t_sqrt_hz": 1.0e-9,
        "bandwidth_hz": 7.5e6,
        "band": (1.0e6, 1.0e8),
        "note": "frequency-tunable CW-ODMR, 7.5 MHz bandwidth (J. Appl. Phys. "
                "135, 044401, 2024)",
    },
    {
        "mode": "quantum-heterodyne",
        "sensitivity_t_sqrt_hz": 1.0e-12,
        "bandwidth_hz": 1.0e6,
        "band": (1.0e6, 1.0e9),
        "note": "representative order of magnitude for heterodyne readout at "
                "high frequency; treat as indicative",
    },
]

#: Frequencies to price, including the band edges and the candidates' cells.
PRICED_HZ = [13.4e6, 45.9843e6, 47.2416e6, 53.7833e6, 86.9e6]

#: A field strength worth comparing with: strong local RFI near a telescope is
#: of order 1 V/m at the antenna. It is a reference point, not a threshold of
#: this survey.
REFERENCE_E_V_M = 1.0


def feasibility() -> list[dict]:
    rows = []
    for spec in MODE_TABLE:
        low, high = spec["band"]
        threshold_t = (spec["sensitivity_t_sqrt_hz"]
                       * (spec["bandwidth_hz"] ** 0.5))
        rows.append({
            "mode": spec["mode"],
            "sensitivity_t_sqrt_hz": spec["sensitivity_t_sqrt_hz"],
            "sensitivity_v_m_sqrt_hz": field_from_magnetic(
                spec["sensitivity_t_sqrt_hz"]),
            "bandwidth_hz": spec["bandwidth_hz"],
            "band_low_hz": low,
            "band_high_hz": high,
            "coverage_of_lwa_band": band_covered(BAND_LOW_HZ, BAND_HIGH_HZ,
                                                 low, high),
            "threshold_b_t": threshold_t,
            "threshold_e_v_m": field_from_magnetic(threshold_t),
            "snr_for_1v_per_m": significance(
                magnetic_from_field(REFERENCE_E_V_M),
                spec["sensitivity_t_sqrt_hz"], spec["bandwidth_hz"]),
            "note": spec["note"],
        })
    return rows


def bias_table() -> list[dict]:
    return [{"frequency_hz": f, "frequency_mhz": f / 1e6,
             "bias_field_t": zeeman_bias_t(f),
             "bias_field_mt": zeeman_bias_t(f) * 1e3}
            for f in PRICED_HZ]


def record_rows(records) -> list[dict]:
    """Every radio record against the lane, with no sensor available."""
    rows = []
    for record in records:
        stamp = record.metadata.get("start_time_utc")
        radio_time = None
        if stamp:
            import datetime as dt
            try:
                radio_time = dt.datetime.fromisoformat(
                    str(stamp).replace("Z", "+00:00"))
            except ValueError:
                radio_time = None
        result = verdict(None, frequency_hz=record.frequency_hz,
                         radio_time_utc=radio_time)
        rows.append({
            "observation_date": record.observation_date,
            "freq_mhz": round(record.frequency_hz / 1e6, 4),
            "disposition": record.disposition,
            **result,
        })
    return rows


def write_report(ledger: dict) -> None:
    lines = [
        "# The diamond NV-centre lane",
        "",
        "Framework: \u03c6\u00b2 = \u03c6 + 1. Attribution: Tyree Jones (tyofthestarz), "
        "prior art from December 2025.",
        "",
        "A nitrogen-vacancy ensemble in diamond measures magnetic field, and a",
        "radio wave carries one: for a plane wave E = cB, so 1 pT is 0.3 mV/m. A",
        "millimetre of diamond at the telescope is therefore a **local field",
        "witness** - it can state that a magnetic transient happened at the site",
        "at the same frequency and second as a radio candidate.",
        "",
        "## The asymmetry, which is the whole point",
        "",
        "| outcome | what it means | allowed to |",
        "| --- | --- | --- |",
        "| a local transient in band and in window | the simplest explanation is a",
        "local emitter | **detract** |",
        "| no local transient | compatible with a real emitter polarised away from",
        "the NV axis, or coupled to the antenna and not the diamond | nothing |",
        "| no sensor, no band, no window | the question was not asked | nothing |",
        "",
        "`guard_support` raises `NVAsEvidence`, so the lane cannot be cited as",
        "support even by accident, and the OMNI gate `local_field` is required by",
        "no disposition. This is the same rule the phi layer carries, for the same",
        "reason: a check that can only ever remove belief must never add it.",
        "",
        "## Feasibility: can a diamond see 13.4-86.9 MHz at all?",
        "",
        "The NV ground-state splitting is 2.87 GHz and tunes at 2.8 MHz/G, so the",
        "band is reachable only with a substantial bias field:",
        "",
        "| frequency | bias field needed |",
        "| --- | --- |",
    ]
    for row in ledger["bias"]:
        lines.append(f"| {row['frequency_mhz']:.2f} MHz | "
                     f"{row['bias_field_mt']:.1f} mT |")
    lines += [
        "",
        "Sensitivity per technique, with the bandwidth that makes the number",
        "mean something (threshold = sensitivity x sqrt(bandwidth)):",
        "",
        "| mode | sensitivity | as E-field | bandwidth | LWA band covered | "
        "threshold | SNR for 1 V/m |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in ledger["feasibility"]:
        lines.append(
            f"| {row['mode']} | {row['sensitivity_t_sqrt_hz']:.1e} T/rtHz | "
            f"{row['sensitivity_v_m_sqrt_hz'] * 1e3:.2f} mV/m/rtHz | "
            f"{row['bandwidth_hz']:.1e} Hz | "
            f"{row['coverage_of_lwa_band']:.0%} | "
            f"{row['threshold_e_v_m'] * 1e3:.2f} mV/m | "
            f"{row['snr_for_1v_per_m']:.1f} |")
    lines += [
        "",
        "Read the two columns together and the conclusion is narrower than the",
        "sensitivity figures suggest, because band and sensitivity are mutually",
        "exclusive with present hardware. The sequences that reach pT sensitivity",
        "(echo, CPMG) stop at a few MHz and cover **none** of the survey band; the",
        "technique that does reach the band (RF-dressed CW-ODMR) is nT-sensitive",
        "over 7.5 MHz, so its threshold is ~0.8 V/m. A diamond can therefore see",
        "**strong local RFI** across 13.4-86.9 MHz, or pT-level fields in narrow",
        "slices below 10 MHz, but not both. The 20 mV/m row with 0% coverage is",
        "the trap: the impressive number belongs to a band the survey never",
        "observes.",
        "",
        "## This survey's records",
        "",
        "| date | MHz | verdict | reason |",
        "| --- | --- | --- | --- |",
    ]
    for row in ledger["records"]:
        reason = "; ".join(row.get("reasons") or [])[:110]
        lines.append(f"| {row['observation_date']} | {row['freq_mhz']:.4f} | "
                     f"{row['verdict']} | {reason} |")
    lines += [
        "",
        "**No NV sensor data exists for these observations**, so every record is",
        "`NOT TESTABLE` and none is detracted. That is the current state and not a",
        "result: the lane needs a sensor at the site, running during the",
        "observation, with a bias field and a calibration - three things this",
        "product does not have on disk.",
        "",
        "## What integration would require",
        "",
        "1. A sensor record per night: band, sensitivity, bandwidth, bias field,",
        "   axis orientation, and a UTC start with an uncertainty. The schema",
        "   refuses a record missing any of them, because each is a way for a null",
        "   to look like evidence.",
        "2. A calibration against a known local emitter, so the sensitivity is a",
        "   measurement and not a datasheet number.",
        "3. A stated coupling model: the NV axis responds to one field component,",
        "   so an unpolarised local emitter can be invisible to it. Without that",
        "   model a null cannot even be interpreted as a weak null.",
        "",
        "Until then this lane subtracts nothing and adds nothing, which is the",
        "correct behaviour for an instrument that does not exist yet.",
        "",
    ]
    with open(os.path.join(ROOT, "NV_LANE.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> int:
    payload = json.load(open(os.path.join(ROOT, "signal_records.json"),
                            encoding="utf-8"))
    records = [SignalRecord.from_dict(row) for row in payload["records"]]
    ledger = {
        "band_hz": [BAND_LOW_HZ, BAND_HIGH_HZ],
        "bias": bias_table(),
        "feasibility": feasibility(),
        "records": record_rows(records),
        "sensing_modes": list(SENSING_MODES),
        "note": ("a local-field record may only detract; absence of a transient "
                 "is compatible with a real local emitter, so it is never "
                 "evidence for an astrophysical origin"),
        "sources": [spec["note"] for spec in MODE_TABLE],
    }
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "nv_lane.json"), "w",
              encoding="utf-8") as fh:
        json.dump(ledger, fh, indent=2)
    write_report(ledger)

    print("NV LANE  (local field witness; can only detract)")
    print(f"   {'frequency':>12} {'bias field':>12}")
    for row in ledger["bias"]:
        print(f"   {row['frequency_mhz']:9.2f} MHz {row['bias_field_mt']:9.1f} mT")
    print()
    print(f"   {'mode':20} {'sensitivity':>14} {'bandwidth':>10} "
          f"{'coverage':>9} {'threshold':>12}")
    for row in ledger["feasibility"]:
        print(f"   {row['mode']:20} "
              f"{row['sensitivity_t_sqrt_hz']:9.1e} T/rHz "
              f"{row['bandwidth_hz']:9.1e} "
              f"{row['coverage_of_lwa_band']:8.0%} "
              f"{row['threshold_e_v_m'] * 1e3:9.2f} mV/m")
    print()
    print(f"   records testable: "
          f"{sum(1 for r in ledger['records'] if r['verdict'] != 'NOT TESTABLE')}"
          f" of {len(ledger['records'])}")
    print("   wrote results/nv_lane.json and NV_LANE.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
