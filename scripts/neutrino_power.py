"""What this neutrino lane could detect, as a function of its assumptions.

    python3 scripts/neutrino_power.py

Writes ``NEUTRINO_POWER.md`` and ``results/neutrino_power.json``. The point is
to replace "the lane is uninformative" with the numbers that make it so, and to
show which of the assumptions would have to change for it to stop being true.

Three curves, kept apart on purpose:

* how often a random alert lands in the box (this is the coincidence null),
* how much flux puts one neutrino in the box (this is the sensitivity), and
* how far away a flare of a given energy could be seen (this is the reach).

They move in opposite directions with the window length, which is the trap this
script exists to make visible: a longer window improves the flux limit and
destroys the coincidence.
"""

from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.sensitivity import (  # noqa: E402
    background_events,
    chance_probability,
    energy_fluence_limit_erg_cm2,
    fluence_limit_cm2,
    horizon_mpc,
    power,
    window_for_chance_probability,
)

#: Assumptions, all disputable on purpose. The alert rate is the whole-sky rate
#: of alerts an observer could act on; the effective area is a mid-range value
#: for a cubic-kilometre detector at a PeV, where the cross-section has grown;
#: the flare energy is a bright-but-not-extraordinary neutrino flare.
ASSUMPTIONS = {
    "alert_rate_per_year": 100,
    "alert_rate_range_per_year": [30, 100, 1000],
    "effective_area_m2": 1.0e4,
    "effective_area_range_m2": [1.0e2, 1.0e3, 1.0e4],
    "pivot_energy_gev": 1.0e6,
    "flare_energy_erg": 1.0e47,
    "flare_duration_s": 10.0,
    "note": ("order-of-magnitude values: the alert rate is a published-alert "
             "rate, the effective area is a PeV track value for a "
             "cubic-kilometre detector, and the flare energy is a bright "
             "transient. Each scales the answer linearly and is stated so it "
             "can be disputed."),
}

WINDOWS = [
    ("1 s", 1.0),
    ("256 s (survey window)", 256.0),
    ("1 h", 3600.0),
    ("1 day", 86400.0),
    ("30 days", 2.592e6),
]
CONES_DEG = [1.0, 3.0, 4.68]
CONE_LABELS = ["1.0 deg (alert-grade)", "3.0 deg (3 sigma of 1 deg)",
               "4.68 deg (3 sigma of the demo geometry)"]
DURATIONS = [("1 s", 1.0), ("10 s", 10.0), ("256 s", 256.0), ("1 day", 86400.0)]


def background_table() -> list[dict]:
    rows = []
    for label, window in WINDOWS:
        row = {"window": label, "window_s": window}
        for cone, cone_label in zip(CONES_DEG, CONE_LABELS):
            mu = background_events(ASSUMPTIONS["alert_rate_per_year"], cone, window)
            row[cone_label] = {"mu": mu, "p": chance_probability(mu)}
        rows.append(row)
    return rows


def timing_table() -> list[dict]:
    rows = []
    for rate in ASSUMPTIONS["alert_rate_range_per_year"]:
        for cone in CONES_DEG:
            rows.append({
                "alert_rate_per_year": rate,
                "cone_deg": cone,
                "max_window_s": window_for_chance_probability(rate, cone, 0.01),
            })
    return rows


def fluence_table() -> list[dict]:
    rows = []
    for area in ASSUMPTIONS["effective_area_range_m2"]:
        for label, duration in DURATIONS:
            number = fluence_limit_cm2(0.9, area, duration)
            energy = energy_fluence_limit_erg_cm2(0.9, area, duration,
                                                 ASSUMPTIONS["pivot_energy_gev"])
            rows.append({
                "effective_area_m2": area,
                "duration": label,
                "duration_s": duration,
                "fluence_cm2": number,
                "energy_fluence_erg_cm2": energy,
                "horizon_mpc_1e47erg": horizon_mpc(
                    ASSUMPTIONS["flare_energy_erg"], energy),
            })
    return rows


def alert_far_table() -> list[dict]:
    """The null a real alert already carries: its own false alarm rate."""
    path = os.path.join(ROOT, "results", "icecube_ingest.json")
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as fh:
        ledger = json.load(fh)
    rows = []
    for alert in ledger.get("alerts", []):
        far = alert.get("far_per_year")
        if far is None:
            continue
        rows.append({
            "event": alert["event"]["event_id"],
            "far_per_year": far,
            "p_in_1s": chance_probability(far * 1.0 / 3.1557e7),
            "p_in_256s": chance_probability(far * 256.0 / 3.1557e7),
            "p_in_1day": chance_probability(far * 86400.0 / 3.1557e7),
        })
    return rows


def write_report(ledger: dict) -> None:
    a = ASSUMPTIONS
    lines = [
        "# What the neutrino lane could detect",
        "",
        "Framework: \u03c6\u00b2 = \u03c6 + 1. Attribution: Tyree Jones (tyofthestarz), "
        "prior art from December 2025.",
        "",
        "This document exists so that `UNINFORMATIVE` is a measured statement and",
        "not an excuse. It does not change any disposition: no number here can",
        "promote a record, and the seven mask-manufactured cells stay excluded",
        "whatever the flux limit turns out to be.",
        "",
        "## Assumptions",
        "",
        f"- whole-sky alert rate: {a['alert_rate_per_year']} per year "
        f"(range {a['alert_rate_range_per_year']})",
        f"- effective area: {a['effective_area_m2']:.0e} m^2 "
        f"(range {a['effective_area_range_m2']})",
        f"- pivot energy: {a['pivot_energy_gev']:.0e} GeV",
        f"- flare: {a['flare_energy_erg']:.0e} erg over "
        f"{a['flare_duration_s']:.0f} s",
        "",
        f"{a['note']}",
        "",
        "## 1. How often a random alert lands in the box",
        "",
        "| window | " + " | ".join(CONE_LABELS) + " |",
        "| --- | " + " | ".join("---" for _ in CONE_LABELS) + " |",
    ]
    for row in ledger["background"]:
        cells = " | ".join(
            f"mu {row[label]['mu']:.2e}, p {row[label]['p']:.2e}"
            for label in CONE_LABELS)
        lines.append(f"| {row['window']} | {cells} |")
    lines += [
        "",
        "The background is tiny even for a one-degree cone and a day, because the",
        "alert rate is only ~100 per year: 100/yr over 41,253 square degrees is one",
        "event per 413 square degrees per year, and a one-degree cone is 3.1 square",
        "degrees. **A coincidence is therefore easy and means little.** The",
        "limiting quantity is not the background, it is the signal, which is why",
        "the next two tables matter more.",
        "",
        "## 2. Timing precision needed to keep the null under 1%",
        "",
        "| alert rate (/yr) | cone | longest window |",
        "| --- | --- | --- |",
    ]
    for row in ledger["timing"]:
        lines.append(f"| {row['alert_rate_per_year']} | {row['cone_deg']:.2f} deg "
                     f"| {row['max_window_s']:.3g} s "
                     f"({row['max_window_s'] / 86400.0:.1f} days) |")
    lines += [
        "",
        "The radio product has **no timestamp at all** - the files are named by",
        "date - so its window is not 256 s and not a day, it is undefined. That is",
        "the whole of the `UNINFORMATIVE` verdict: not a null result, an",
        "unbounded box.",
        "",
        "## 3. Sensitivity: flux and reach",
        "",
        "| A_eff (m^2) | duration | fluence (cm^-2) | energy fluence (erg/cm^2) | horizon for 1e47 erg |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in ledger["fluence"]:
        lines.append(
            f"| {row['effective_area_m2']:.0e} | {row['duration']} | "
            f"{row['fluence_cm2']:.2e} | {row['energy_fluence_erg_cm2']:.2e} | "
            f"{row['horizon_mpc_1e47erg']:.0f} Mpc |")
    short = next(r for r in ledger["fluence"]
                 if r["duration_s"] == 1.0
                 and r["effective_area_m2"] == a["effective_area_m2"])
    survey = next(r for r in ledger["fluence"]
                  if r["duration_s"] == 256.0
                  and r["effective_area_m2"] == a["effective_area_m2"])
    lines += [
        "",
        f"At {a['effective_area_m2']:.0e} m^2, the same detector is "
        f"**{short['energy_fluence_erg_cm2'] / survey['energy_fluence_erg_cm2']:.0f} "
        "times more sensitive** over a 256 s window than over a 1 s window, and a",
        "correspondingly smaller share of the sky is excluded per second. This is",
        "the trade the lane cannot avoid: the window that maximises the chance of",
        "catching a flare is the window that maximises the background.",
        "",
        "## 4. The null a real alert already carries",
        "",
    ]
    if ledger["alert_far"]:
        lines += [
            "| event | FAR (/yr) | p within 1 s | p within 256 s | p within 1 day |",
            "| --- | --- | --- | --- | --- |",
        ]
        for row in ledger["alert_far"]:
            lines.append(
                f"| {row['event']} | {row['far_per_year']:.3g} | "
                f"{row['p_in_1s']:.2e} | {row['p_in_256s']:.2e} | "
                f"{row['p_in_1day']:.2e} |")
        lines += [
            "",
            "These come from the ingested alerts' own false alarm rates, which is a",
            "better null than the published rate used above: it already accounts",
            "for the alert threshold. Note what it buys - even a day-long window",
            "gives a small probability - and what it does not: a small probability",
            "for a *background event*, which says nothing about whether the radio",
            "feature is real.",
        ]
    else:
        lines += ["No ingested alerts with a published false alarm rate.", ]
    lines += [
        "",
        "## What would make the lane informative",
        "",
        "1. A radio pointing: RA/Dec per record, from a beam-formed or targeted",
        "   follow-up. Without it the cone does not exist.",
        "2. A UTC start time and an integration length. That turns the window from",
        "   undefined into 256 s, at which point the background is ~1e-7 and the",
        "   background stops being the problem.",
        "3. Enough effective area or enough flare energy to put one event in the",
        "   box - the table above says how much, and the answer scales linearly in",
        "   both.",
        "",
        "Until (1) and (2) exist, the honest report is that the lane cannot be",
        "run, and this document is the measurement of that sentence.",
        "",
    ]
    with open(os.path.join(ROOT, "NEUTRINO_POWER.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> int:
    ledger = {
        "assumptions": ASSUMPTIONS,
        "background": background_table(),
        "timing": timing_table(),
        "fluence": fluence_table(),
        "alert_far": alert_far_table(),
        "power_at_90pct_events": power(2.302585),
    }
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "neutrino_power.json"), "w",
              encoding="utf-8") as fh:
        json.dump(ledger, fh, indent=2)
    write_report(ledger)

    print("NEUTRINO POWER")
    print(f"   {'window':24} " + "  ".join(f"{c:>22}" for c in CONE_LABELS))
    for row in ledger["background"]:
        cells = "  ".join(f"p {row[label]['p']:.2e}".rjust(22)
                          for label in CONE_LABELS)
        print(f"   {row['window']:24} {cells}")
    print()
    for row in ledger["fluence"]:
        if row["effective_area_m2"] != ASSUMPTIONS["effective_area_m2"]:
            continue
        print(f"   A={row['effective_area_m2']:.0e} m2 {row['duration']:>6}: "
              f"{row['energy_fluence_erg_cm2']:.2e} erg/cm2 -> "
              f"{row['horizon_mpc_1e47erg']:.0f} Mpc for 1e47 erg")
    print("   wrote results/neutrino_power.json and NEUTRINO_POWER.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
