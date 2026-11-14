"""Crossmatch the radio records against a public neutrino catalog.

Run it as ``python3 scripts/neutrino_crossmatch.py``. It reads
``data/neutrino_events.json`` when a real catalog has been placed there, and
falls back to the labelled synthetic fixture otherwise. The two are never mixed:
the ledger records which one it used, and every synthetic event carries
``alert_status: "test"`` and a ``data_release`` that says so.

Three things this script refuses to do, each because it is the way this lane
would otherwise go wrong:

1. Treat an empty query as a result. With no pointing and no UTC timestamp the
   query cannot be drawn at all, and the ledger says ``UNINFORMATIVE`` with the
   missing ingredient named rather than "no neutrino found".
2. Let a coincidence change a disposition. ``guard_support`` is called on the
   records the radio pipeline invalidated, and the refusal is written into the
   ledger as evidence that the rule is enforced rather than promised.
3. Report a chance probability without its exposure. The worked example prints
   the assumed rate beside it, because a p-value is only as good as the null
   behind it.
"""

from __future__ import annotations

import datetime as dt
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.neutrino import (  # noqa: E402
    Coincidence,
    DetectorExposure,
    NeutrinoAsEvidence,
    NeutrinoEvent,
    guard_support,
    match_catalog,
)
from seti_search.signal_record import INVALID_CANDIDATE, SignalRecord  # noqa: E402

REAL_CATALOG = os.path.join(ROOT, "data", "neutrino_events.json")
SYNTHETIC_CATALOG = os.path.join(ROOT, "data", "neutrino_events_synthetic.json")

#: A worked example of the machinery, not an observation: the radio geometry is
#: hypothetical and says so in the ledger. It exists to show what a coincidence
#: would look like and how small its chance probability would have to be, so that
#: nobody mistakes the real lane's UNINFORMATIVE rows for a null result.
DEMO = {
    "observation_date": "20240322",
    "frequency_mhz": 47.2416,
    "ra_deg": 100.0,
    "dec_deg": -20.0,
    "time_utc": dt.datetime(2024, 3, 22, 8, 15, 0, tzinfo=dt.timezone.utc),
    "radio_time_uncertainty_s": 256.0,
    "radio_localisation_deg": 1.0,
    "note": "HYPOTHETICAL geometry, chosen to be inside the synthetic cone",
}

#: An order-of-magnitude archival rate for a cubic-kilometre detector: of order
#: 1e5 events per year over 4 pi at all energies, with O(100) alert-grade ones.
#: It is stated here so the worked example's probability can be recomputed and
#: disputed; it is not a measurement of any detector's exposure.
DEMO_EXPOSURE = DetectorExposure(
    livetime_s=3.156e7, n_events=100, sky_fraction=1.0,
    energy_threshold_gev=1e14,
    note="order-of-magnitude alert rate, stated for the worked example only")


def load_catalog() -> tuple[list[NeutrinoEvent], str, str]:
    if os.path.exists(REAL_CATALOG):
        path, kind = REAL_CATALOG, "real"
    else:
        path, kind = SYNTHETIC_CATALOG, "synthetic"
    with open(path, encoding="utf-8") as fh:
        rows = json.load(fh)
    events = [NeutrinoEvent.from_dict(row) for row in rows]
    return events, kind, os.path.relpath(path, ROOT)


def load_from(path: str, kind: str):
    with open(path, encoding="utf-8") as fh:
        rows = json.load(fh)
    return ([NeutrinoEvent.from_dict(row) for row in rows], kind,
            os.path.relpath(path, ROOT))


def record_time(record: SignalRecord) -> dt.datetime | None:
    """The record's start time, if the product carries one.

    These files are named by date and carry no WCS and no UTC stamp in the record
    metadata, so this returns ``None`` and the crossmatch is uninformative. It is
    written as a function so that a product which *does* carry a stamp can be
    matched without touching the rest of this script.
    """
    raw = record.metadata.get("start_time_utc")
    if not raw:
        return None
    return dt.datetime.fromisoformat(str(raw).replace("Z", "+00:00"))


def record_position(record: SignalRecord) -> tuple[float | None, float | None]:
    ra = record.metadata.get("ra_deg")
    dec = record.metadata.get("dec_deg")
    return (ra, dec) if ra is not None and dec is not None else (None, None)


def crossmatch_records(records, events) -> list[dict]:
    rows = []
    for record in records:
        ra, dec = record_position(record)
        stamp = record_time(record)
        if record.disposition in INVALID_CANDIDATE:
            rows.append({
                "observation_date": record.observation_date,
                "freq_mhz": round(record.frequency_hz / 1e6, 4),
                "disposition": record.disposition,
                "status": "EXCLUDED_FROM_CROSSMATCH",
                "note": ("the radio event was manufactured by the excision mask; "
                         "a neutrino coincidence is not evidence about an "
                         "artifact, so the lane does not run on it"),
            })
            continue
        result = match_catalog(events,
                               ra_deg=ra, dec_deg=dec, time_utc=stamp,
                               radio_time_uncertainty_s=record.duration_s)
        rows.append({
            "observation_date": record.observation_date,
            "freq_mhz": round(record.frequency_hz / 1e6, 4),
            "disposition": record.disposition,
            "status": result.status,
            "separation_deg": result.separation_deg,
            "dt_s": result.dt_s,
            "chance_probability": result.chance_probability,
            "limits": list(result.limits),
            "note": result.note,
        })
    return rows


def demo_crossmatch(records, events) -> dict:
    record = next((r for r in records
                   if r.observation_date == DEMO["observation_date"]
                   and abs(r.frequency_hz / 1e6 - DEMO["frequency_mhz"]) < 0.05),
                  None)
    result = match_catalog(events,
                           ra_deg=DEMO["ra_deg"], dec_deg=DEMO["dec_deg"],
                           time_utc=DEMO["time_utc"],
                           radio_time_uncertainty_s=DEMO["radio_time_uncertainty_s"],
                           radio_localisation_deg=DEMO["radio_localisation_deg"],
                           exposure=DEMO_EXPOSURE)
    return {
        "synthetic_geometry": True,
        "note": DEMO["note"],
        "record": None if record is None else {
            "observation_date": record.observation_date,
            "freq_mhz": round(record.frequency_hz / 1e6, 4),
            "disposition": record.disposition,
        },
        "exposure": {
            "livetime_s": DEMO_EXPOSURE.livetime_s,
            "n_events": DEMO_EXPOSURE.n_events,
            "sky_fraction": DEMO_EXPOSURE.sky_fraction,
            "density_per_sr_s": DEMO_EXPOSURE.density_per_sr_s,
            "note": DEMO_EXPOSURE.note,
        },
        "result": result.to_dict(),
    }


def restoration_refusals(records, events) -> list[dict]:
    """Call the guard on every invalid record and record the refusal.

    This is the point of the lane: the rule is exercised, not documented. If the
    guard ever stops raising, this list empties and the run should be treated as
    a failure, which is why the script prints its length.
    """
    refusals = []
    coincidence = match_catalog(events,
                                ra_deg=DEMO["ra_deg"], dec_deg=DEMO["dec_deg"],
                                time_utc=DEMO["time_utc"],
                                radio_time_uncertainty_s=DEMO["radio_time_uncertainty_s"],
                                radio_localisation_deg=DEMO["radio_localisation_deg"],
                                exposure=DEMO_EXPOSURE)
    for record in records:
        if record.disposition not in INVALID_CANDIDATE:
            continue
        try:
            guard_support(record, coincidence,
                          claimed="CONFIRMED TECHNOSIGNATURE CANDIDATE")
        except NeutrinoAsEvidence as exc:
            refusals.append({
                "observation_date": record.observation_date,
                "freq_mhz": round(record.frequency_hz / 1e6, 4),
                "disposition": record.disposition,
                "refused": str(exc),
            })
    return refusals


def write_report(ledger: dict) -> None:
    rows = ledger["records"]
    fmt_p = lambda p: "no null available" if p is None else f"{p:.2e}"
    lines = [
        "# The neutrino lane",
        "",
        "Framework: \u03c6\u00b2 = \u03c6 + 1. Attribution: Tyree Jones (tyofthestarz), "
        "prior art from December 2025.",
        "",
    ]
    if ledger["catalog"]["kind"] == "synthetic":
        lines += [
            "> **SYNTHETIC CATALOG.** Every event below comes from",
            "> `data/neutrino_events_synthetic.json`, whose events carry",
            "> `alert_status: test` and a `data_release` saying they are invented.",
            "> Nothing in this document is a statement about the sky. Place a real",
            "> release at `data/neutrino_events.json` to run the same lane on it.",
            "",
        ]
    lines += [
        "A neutrino is a second messenger: it leaves the source without being",
        "absorbed the way radio is, and it arrives with no dispersion to speak of,",
        "so a coincidence in time and angle is a coincidence between two",
        "independent channels. That is what makes it worth asking about. It is",
        "also why the answer has to be stated carefully, because a neutrino",
        "detector is a poor instrument for ruling things out:",
        "",
        "| limitation | consequence |",
        "| --- | --- |",
        "| energy threshold (~TeV and up) | a source that is radio-bright but",
        "neutrino-quiet sits below it and is invisible |",
        "| angular resolution (0.5-15 deg) | a one-degree cone contains many",
        "square degrees of sky, so chance coincidences are normal |",
        "| sky coverage (up-going vs down-going) | half the sky is seen through",
        "the Earth, at a different effective area |",
        "| duty cycle and livetime | an event can be missed while the detector is",
        "not taking data |",
        "| flux | a single source may simply not produce a detectable event in",
        "the exposure |",
        "",
        "Every one of those is a reason that **an absent neutrino is not evidence",
        "that a radio feature is artificial**, and this module refuses to represent",
        "it otherwise: `Coincidence.absence_is_informative` is a field fixed at",
        "`False`, and an empty query must name the limit that makes it empty.",
        "",
        "## The rule that outranks the rest",
        "",
        "A neutrino coincidence cannot restore an invalid radio candidate. If the",
        "radio pipeline found that a signal was manufactured by its own mask, then",
        "the radio event does not exist, and a real particle that happens to point",
        "there is a coincidence with an artifact. The run below calls the guard on",
        f"every invalid record and records {len(ledger['refusals'])} refusals, so",
        "the rule is exercised rather than promised.",
        "",
        "## The records",
        "",
        "| date | MHz | disposition | status | chance p | limits |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        limits = "; ".join(row.get("limits") or []) or row.get("note", "")
        lines.append(
            f"| {row['observation_date']} | {row['freq_mhz']:.4f} | "
            f"{row['disposition']} | {row['status']} | "
            f"{fmt_p(row.get('chance_probability'))} | {limits} |")
    lines += [
        "",
        "## What the answer currently is",
        "",
        "**UNINFORMATIVE, and the reason is the product, not the sky.** These",
        "records have no world coordinate system and no UTC timestamp: they are a",
        "single disk-integrated beam, named by date. There is therefore no position",
        "to put a cone around and no time to draw a window on, and the lane reports",
        "that instead of reporting a null. This is the same limitation that makes",
        "the off-source and barycentric gates untestable, and it is not repaired by",
        "having a catalog.",
        "",
        "What would make it informative: a beam-formed or targeted follow-up with a",
        "pointing, a UTC start time, and a localisation of a degree or better. Then",
        "the same function decides, and the number it reports is a chance",
        "probability against the stated exposure.",
        "",
        "## Worked example, hypothetical geometry",
        "",
        f"- catalog: {ledger['demo'].get('catalog', {}).get('kind', 'unknown')} "
        f"({ledger['demo'].get('catalog', {}).get('path', '?')}), used only for "
        "this example",
        f"- radio position: RA {DEMO['ra_deg']:.2f} deg, Dec {DEMO['dec_deg']:.2f} deg "
        "(**invented to sit inside the synthetic cone**)",
        f"- radio time: {DEMO['time_utc'].isoformat()}, "
        f"+-{DEMO['radio_time_uncertainty_s']:.0f} s",
        f"- exposure: {DEMO_EXPOSURE.n_events} events over "
        f"{DEMO_EXPOSURE.livetime_s:.3e} s and the whole sky, i.e. "
        f"{DEMO_EXPOSURE.density_per_sr_s:.3e} per sr per s ({DEMO_EXPOSURE.note})",
        "",
    ]
    result = ledger["demo"]["result"]
    if result["status"] == "COINCIDENT":
        lines += [
            f"Nearest event: `{result['event']['event_id']}` at "
            f"{result['separation_deg']:.2f} deg and {result['dt_s']:.0f} s, inside "
            f"{result['allowed_separation_deg']:.2f} deg and "
            f"{result['allowed_dt_s']:.0f} s.",
            "",
            f"Chance probability: **{result['chance_probability']:.2e}**. Read it as "
            "the probability that one random event from this exposure lands in that",
            "box, which is what makes it small: the box is small in solid angle,",
            "not because the event is remarkable.",
        ]
    else:
        lines += [
            f"Status: {result['status']}. "
            + "; ".join(result.get("limits") or [result.get("note", "")]),
        ]
    lines += [
        "",
        "The example is not a candidate and is not written to the record list. It",
        "demonstrates the arithmetic, and it shows the asymmetry: a coincidence is",
        "an invitation to point a telescope, never a disposition.",
        "",
    ]
    with open(os.path.join(ROOT, "NEUTRINO_LANE.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> int:
    payload = json.load(open(os.path.join(ROOT, "signal_records.json"),
                            encoding="utf-8"))
    records = [SignalRecord.from_dict(row) for row in payload["records"]]
    events, kind, path = load_catalog()
    demo_events, demo_kind, demo_path = events, kind, path
    if kind == "real" and os.path.exists(SYNTHETIC_CATALOG):
        demo_events, demo_kind, demo_path = load_from(SYNTHETIC_CATALOG,
                                                     "synthetic")

    rows = crossmatch_records(records, events)
    demo = demo_crossmatch(records, demo_events)
    demo["catalog"] = {"kind": demo_kind, "path": demo_path}
    refusals = restoration_refusals(records, demo_events)

    ledger = {
        "catalog": {
            "kind": kind,
            "path": path,
            "n_events": len(events),
            "detectors": sorted({e.detector for e in events}),
            "data_release": sorted({e.data_release for e in events}),
            "synthetic": any(e.synthetic for e in events),
        },
        "records": rows,
        "demo": demo,
        "refusals": refusals,
        "absence_is_informative": False,
        "note": ("a neutrino coincidence may add a counterpart to a valid record; "
                 "it may never restore an invalid one, and absence proves nothing"),
    }
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "neutrino_crossmatch.json"), "w",
              encoding="utf-8") as fh:
        json.dump(ledger, fh, indent=2)
    write_report(ledger)

    print(f"NEUTRINO LANE  ({len(events)} events, {kind} catalog: {path})")
    if kind == "synthetic":
        print("   SYNTHETIC CATALOG: no claim below is about the sky")
    print(f"   {'date':8} {'MHz':>9} {'status':24} {'chance p':>12}  note")
    for row in rows:
        p = row.get("chance_probability")
        print(f"   {row['observation_date']:8} {row['freq_mhz']:9.4f} "
              f"{row['status']:24} {'-' if p is None else f'{p:.2e}':>12}  "
              f"{(row.get('limits') or [row.get('note', '')])[0][:70]}")
    result = demo["result"]
    print(f"\n   worked example (hypothetical geometry): {result['status']}"
          + (f", chance p {result['chance_probability']:.2e}"
             if result["status"] == "COINCIDENT" else ""))
    print(f"   restoration refused on {len(refusals)} invalid records")
    print("   wrote results/neutrino_crossmatch.json and NEUTRINO_LANE.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
