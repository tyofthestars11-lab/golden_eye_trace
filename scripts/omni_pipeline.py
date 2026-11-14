"""One entry point for all six layers, with the disposition audited at the end.

Every other script in this repository answers one question and writes one
artifact. This one runs the questions in order, hands each candidate's claimed
disposition to the omni layer together with the gates that produced it, and
reports the strongest disposition the gates actually permit. It writes
OMNI.md and omnipipeline.json; it does not re-search anything, so it runs in a
second and cannot drift from the runs it reads.

The layers, in the order they are allowed to speak:

  1. measurement    the stabilizer (window vs trail), the local significance,
                    and the alias signature
  2. calibration    which noise reference the significance was measured against
  3. timedomain     whether the nonlinear coordinate beats its warp control
  4. phi            the golden-ratio rung, recorded and never supporting
  5. interpretation the mechanism: parallel lanes at a coronal exciter's speed
  6. confirmation   off-source, second telescope, barycentric, repeat

A disposition is a statement about layer 6, so it is checked last and it can
only be as strong as the layers under it. Nothing here can raise a disposition:
the audit's only powers are to leave a claim alone or to lower it.

    python3 scripts/omni_pipeline.py            # exit 1 if a confirmation is claimed
    python3 scripts/omni_pipeline.py --strict   # exit 1 on any downgrade
"""

from __future__ import annotations

import argparse
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.omni import LAYERS, REQUIRED, Gate, audit  # noqa: E402
from seti_search.signal_record import SignalRecord  # noqa: E402
from seti_search.stabilizer import check_window, fastest_rate_hz_s  # noqa: E402

CHANNEL_HZ = 100547.19562243501
TSAMP_S = 0.256
THRESHOLD = 5.5
WINDOW_CHANNELS = 96
SAMPLES = 1000
ALIAS_GAP = 0.49          # |rate| this close to 0.5 cycles/sample is the fold
TOL_MHZ = 0.01
PARALLEL_TOLERANCE = 0.1  # lanes must agree within a tenth of a resolution element

MARKS = {"PASS": ".", "FAIL": "F", "UNTESTABLE": "U", "NOT RUN": "-"}


def load(name: str) -> dict:
    """Read a run artifact: prereg/ holds the registered runs, the root the joins."""
    for folder in ("prereg", ""):
        path = os.path.join(ROOT, folder, name)
        if os.path.exists(path):
            with open(path, encoding="utf-8") as fh:
                return json.load(fh)
    raise SystemExit(f"missing artifact {name}; run the run that writes it first")


def find(rows: list, date: str, freq_mhz: float, tol_mhz: float = TOL_MHZ):
    for row in rows:
        if row.get("date") == date and abs(row["freq_mhz"] - freq_mhz) <= tol_mhz:
            return row
    return None


def gate_mask_dependence(record: SignalRecord) -> Gate:
    """Does the candidate's significance survive an honestly scaled mask?

    The gate passes when the dependence it is auditing is *demonstrated*: the
    theory-scaled row test and the empirical one disagree by a factor, and the
    empirical value is below threshold. A candidate with no mask measurement is
    NOT RUN rather than PASS, because absence of the test is not absence of the
    artifact.
    """
    theory = record.metadata.get("mask_theory_sigma")
    empirical = record.metadata.get("mask_empirical_sigma")
    if theory is None or empirical is None:
        return Gate(id="mask_dependence", layer="calibration", status="NOT RUN",
                    evidence="no mask measurement was carried for this record")
    ratio = float(theory) / max(float(empirical), 1e-9)
    rows = (f"{record.metadata.get('mask_blanked_rows_theory', '?')} of "
            f"{record.metadata.get('mask_total_rows', '?')} integrations blanked "
            f"by the theory test")
    passes = ratio >= 3.0 and float(empirical) < THRESHOLD
    return Gate(
        id="mask_dependence", layer="calibration",
        status="PASS" if passes else "FAIL",
        evidence=(f"{float(theory):.2f} sigma through the theory mask against "
                  f"{float(empirical):.2f} through the empirical one, a factor "
                  f"{ratio:.1f}; {rows}"))


def gate_measurement(record: SignalRecord) -> list[Gate]:
    window = check_window(record.drift_hz_per_s, WINDOW_CHANNELS, SAMPLES,
                          CHANNEL_HZ, TSAMP_S)
    detected = record.snr_at_cell >= THRESHOLD
    gap = 0.5 - abs(record.channel_rate)
    return [
        Gate(id="window", layer="measurement",
             status="PASS" if window.fits else "FAIL",
             evidence=window.note),
        Gate(id="detected", layer="measurement",
             status="PASS" if detected else "FAIL",
             evidence=(f"{record.snr_at_cell:.2f} sigma at the reported cell "
                       f"against a {THRESHOLD:g} sigma threshold")),
        Gate(id="alias_signature", layer="measurement",
             status="PASS" if gap <= 0.5 - ALIAS_GAP else "FAIL",
             evidence=(f"{abs(record.channel_rate):.4f} cycles per sample, "
                       f"{gap:.4f} from the fold")),
        gate_mask_dependence(record),
        Gate(id="noise_reference", layer="calibration",
             status="PASS" if record.transform_version.startswith(
                 "ovro-local-noise") else "NOT RUN",
             evidence=(f"{record.transform_version}: local scale "
                       f"{record.local_scale:.2f}, "
                       f"{record.snr_global:.2f} -> {record.snr_at_cell:.2f} sigma")),
    ]


def gate_timedomain(nonlinear: dict | None) -> Gate:
    if nonlinear is None:
        return Gate(id="nonlinear_beats_control", layer="timedomain",
                    status="NOT RUN",
                    evidence="this candidate was not in the nonlinear run")
    beats = bool(nonlinear.get("beats_warp_gain"))
    return Gate(
        id="nonlinear_beats_control", layer="timedomain",
        status="PASS" if beats else "FAIL",
        evidence=(f"gain {nonlinear.get('gain_over_linear', 0.0):.2f} over the "
                  f"linear coordinate, best score "
                  f"{nonlinear.get('best_score', 0.0):.2f}, warp control "
                  f"{'beaten' if beats else 'not beaten'}"))


def gate_phi(record: SignalRecord, confirm: dict) -> Gate:
    phi = confirm.get("phi", {}) or {}
    locked = bool(phi.get("locked"))
    return Gate(
        id="phi_lock", layer="phi",
        status="PASS" if locked else "FAIL",
        evidence=(f"rung {phi.get('rung', float('nan')):.3f}, k={phi.get('k')}, "
                  f"lock distance {phi.get('lock_distance', float('nan')):.4f} "
                  f"({record.channel_rate:+.4f} cycles per sample)"))


def gate_interpretation(record: SignalRecord, confirm: dict,
                        siblings: list[float]) -> Gate:
    """Does an astrophysical mechanism actually explain this lane?"""
    if not confirm.get("resolution_hz_s"):
        return Gate(id="astrophysical_explanation", layer="interpretation",
                    status="NOT RUN",
                    evidence="no confirmation row for this candidate, so no "
                             "resolution element to compare lanes against")
    if len(siblings) < 2:
        return Gate(id="astrophysical_explanation", layer="interpretation",
                    status="FAIL",
                    evidence=("one lane on this date: no parallel structure to "
                              "test the mechanism against"))
    spread = max(siblings) - min(siblings)
    resolution = float(confirm["resolution_hz_s"])
    coronal = float(confirm["coronal_20kms_drift_hz_s"])
    speed_kms = abs(record.drift_hz_per_s) * 20.0 / coronal
    parallel = spread / resolution <= PARALLEL_TOLERANCE
    plausible = 30.0 <= speed_kms <= 300.0
    status = "PASS" if parallel and plausible else "FAIL"
    return Gate(
        id="astrophysical_explanation", layer="interpretation", status=status,
        evidence=(f"{len(siblings)} lanes on this date, spread {spread:.1f} Hz/s "
                  f"({spread / resolution:.3f} of a resolution element); "
                  f"{abs(record.drift_hz_per_s):,.0f} Hz/s implies a "
                  f"{speed_kms:.0f} km/s coronal exciter"))


def gate_neutrino(record: SignalRecord) -> Gate:
    """The neutrino lane's verdict, which can only ever be a counterpart.

    A coincidence is not a required gate for any disposition: the radio gates
    are checks on the radio, and a particle does not perform them. What this
    gate must never become is a PASS used to carry a claim, and it must never
    report an empty query as a null, so both are refused by construction here.
    """
    path = os.path.join(ROOT, "results", "neutrino_crossmatch.json")
    if not os.path.exists(path):
        return Gate(id="neutrino_counterpart", layer="confirmation",
                    status="NOT RUN",
                    evidence=("no neutrino catalog was supplied; absence of a "
                              "neutrino is not evidence about this record"))
    ledger = json.load(open(path, encoding="utf-8"))
    freq_mhz = record.frequency_hz / 1e6
    row = next((r for r in ledger.get("records", [])
                if r.get("observation_date") == record.observation_date
                and abs(r.get("freq_mhz", -1.0) - freq_mhz) <= 0.05), None)
    if row is None:
        return Gate(id="neutrino_counterpart", layer="confirmation",
                    status="NOT RUN",
                    evidence="this record is not in the neutrino ledger")
    status = row.get("status")
    if status == "EXCLUDED_FROM_CROSSMATCH":
        return Gate(id="neutrino_counterpart", layer="confirmation",
                    status="NOT RUN",
                    evidence=("the record is excluded from discovery claims, so "
                              "the lane does not run on it: a coincidence cannot "
                              "restore an artifact"))
    if status == "COINCIDENT":
        p = row.get("chance_probability")
        ok = p is not None and p < 0.01
        return Gate(id="neutrino_counterpart", layer="confirmation",
                    status="PASS" if ok else "UNTESTABLE",
                    evidence=(f"chance p {p:.2e} against the stated exposure; a "
                              "counterpart to follow up, never a disposition"
                              if p is not None else
                              "coincident but with no exposure, so no null"))
    limits = "; ".join(row.get("limits") or []) or row.get("note", "")
    return Gate(id="neutrino_counterpart", layer="confirmation",
                status="UNTESTABLE",
                evidence=(f"{status}: {limits}" if limits else status))


def gates_confirmation(record: SignalRecord, confirm: dict,
                       rfi: dict | None) -> list[Gate]:
    g = confirm.get("gates", {})

    def status(gate_id: str, flag: str, note: str) -> Gate:
        return Gate(id=gate_id, layer="confirmation",
                    status="PASS" if g.get(flag) else "UNTESTABLE",
                    evidence=note)

    off = status(
        "off_source", "off_source_available",
        "the product is one disk-integrated beam, so there is no off-source scan")
    telescope = status(
        "independent_telescope", "independent_telescope",
        "every candidate here comes from one instrument at one epoch")
    bary = status(
        "barycentric", "barycentric_reducible",
        "no world coordinate system and no sky position, so no barycentric "
        "velocity can be removed")

    repeats = bool(g.get("repeats_same_position"))
    testable = bool(g.get("repeat_testable"))
    if rfi is not None:
        repeats = False
        evidence = (f"the pair test at full integration and across windows: "
                    f"{rfi['verdict']}")
    elif not testable:
        evidence = "no second epoch at this position inside the searched band"
    else:
        evidence = (f"same position and drift on {confirm.get('repeats_other_days', 0)} "
                    "other day(s)")
    repeat = Gate(id="repeat", layer="confirmation",
                  status="PASS" if (repeats and testable) else (
                      "FAIL" if testable else "UNTESTABLE"),
                  evidence=evidence)
    return [off, telescope, bary, repeat, gate_neutrino(record)]


def run(strict: bool) -> int:
    payload = load("signal_records.json")
    confirm = load("ovro_confirm.json")
    local = load("local_noise.json")
    rfi = load("rfi_pair.json")
    nonlinear = json.load(open(os.path.join(ROOT, "nonlinear_candidates.json"),
                               encoding="utf-8"))

    records = [SignalRecord.from_dict(r) for r in payload["records"]]
    by_date: dict[str, list[SignalRecord]] = {}
    for record in records:
        by_date.setdefault(record.observation_date, []).append(record)

    rfi_freqs = [row["freq_mhz"] for row in rfi.get("pair", [])]
    rows, audits, escalations, downgrades = [], [], [], []
    for record in records:
        freq_mhz = record.frequency_hz / 1e6
        conf = find(confirm["candidates"], record.observation_date, freq_mhz)
        conf = conf or {}
        nonlin = find(nonlinear["candidates"], record.observation_date, freq_mhz,
                      0.05)
        rfi_row = rfi if any(abs(f - freq_mhz) <= 0.05 for f in rfi_freqs) else None

        siblings = [r.drift_hz_per_s for r in by_date[record.observation_date]
                    if r.disposition == "ASTROPHYSICAL"]
        if record.disposition != "ASTROPHYSICAL":
            siblings = []

        gates = gate_measurement(record)
        gates.append(gate_timedomain(nonlin))
        gates.append(gate_phi(record, conf))
        gates.append(gate_interpretation(record, conf, siblings))
        gates.extend(gates_confirmation(record, conf, rfi_row))

        result = audit(record.disposition, gates)
        window = check_window(record.drift_hz_per_s, WINDOW_CHANNELS, SAMPLES,
                              CHANNEL_HZ, TSAMP_S)
        rows.append({
            "record": record,
            "audit": result,
            "window": window,
            "missing": [gate_id for gate_id in ("detected", "window", "repeat",
                                                "off_source",
                                                "independent_telescope",
                                                "barycentric")
                        if result.status(gate_id) != "PASS"],
        })
        audits.append(result)
        if result.escalation:
            escalations.append(record)
        if result.downgraded:
            downgrades.append((record, result))

    limit = fastest_rate_hz_s(WINDOW_CHANNELS, SAMPLES, CHANNEL_HZ, TSAMP_S)
    nyquist = 0.5 * CHANNEL_HZ / TSAMP_S

    print(f"OMNI PIPELINE  ({len(records)} records, {len(LAYERS)} layers)")
    print(f"   window guard: {WINDOW_CHANNELS} channels x {SAMPLES} samples "
          f"holds a lane up to {limit / 1e3:.1f} kHz/s; the dedoppler ran to "
          f"the fold at {nyquist / 1e3:.1f} kHz/s")
    print(f"   {'date':8} {'MHz':>9} {'claimed':22} {'permitted':22} "
          f"{'trail':>6} {'held':>6}  deficits")
    for row in sorted(rows, key=lambda r: (r["record"].observation_date,
                                           r["record"].frequency_hz)):
        rec, res, win = row["record"], row["audit"], row["window"]
        mark = "" if not res.downgraded else "  <- lowered"
        print(f"   {rec.observation_date:8} {rec.frequency_hz / 1e6:9.4f} "
              f"{res.claimed[:22]:22} {res.permitted[:22]:22} "
              f"{win.trail_channels:6.1f} {win.fraction_held:6.1%}  "
              f"{len(res.deficits)}{mark}")

    print(f"\n   permitted as claimed: {len(records) - len(downgrades)} of "
          f"{len(records)}")
    for record, result in downgrades:
        print(f"   LOWERED {record.observation_date} "
              f"{record.frequency_hz / 1e6:.4f} MHz: {result.claimed} -> "
              f"{result.permitted}")
        for deficit in result.deficits:
            print(f"      {deficit}")
    claimed_top = sum(1 for a in audits
                      if a.claimed == "CONFIRMED TECHNOSIGNATURE CANDIDATE")
    allowed_top = sum(1 for a in audits
                      if a.permitted == "CONFIRMED TECHNOSIGNATURE CANDIDATE")
    print(f"   confirmations: {claimed_top} claimed, {allowed_top} permitted")
    if escalations:
        for record in escalations:
            print(f"   ESCALATION {record.observation_date} "
                  f"{record.frequency_hz / 1e6:.4f} MHz claims a confirmation "
                  "no gate supports", file=sys.stderr)

    write_outputs(rows, payload, limit, nyquist, downgrades, escalations)
    print("   wrote OMNI.md and omnipipeline.json")

    if escalations or (strict and downgrades):
        return 1
    return 0


def write_outputs(rows, payload, limit, nyquist, downgrades, escalations) -> None:
    gate_ids = ["window", "detected", "alias_signature", "mask_dependence",
                "noise_reference",
                "nonlinear_beats_control", "phi_lock",
                "astrophysical_explanation", "off_source",
                "independent_telescope", "barycentric",
                "neutrino_counterpart", "repeat"]

    out = {
        "transform_version": payload["transform_version"],
        "layers": list(LAYERS),
        "window_channels": WINDOW_CHANNELS,
        "samples": SAMPLES,
        "max_drift_hz_s_for_window": limit,
        "nyquist_hz_s": nyquist,
        "n_records": len(rows),
        "n_downgraded": len(downgrades),
        "n_escalations": len(escalations),
        "required": {d: list(g) for d, g in REQUIRED.items()},
        "records": [
            {
                **row["record"].as_dict(),
                "claimed_disposition": row["record"].disposition,
                "permitted_disposition": row["audit"].permitted,
                "window_check": row["window"].as_dict(),
                "audit": row["audit"].as_dict(),
            }
            for row in sorted(rows, key=lambda r: (r["record"].observation_date,
                                                   r["record"].frequency_hz))
        ],
    }
    path = os.path.join(ROOT, "omnipipeline.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")

    claims = sum(1 for row in rows
                 if row["audit"].claimed == "CONFIRMED TECHNOSIGNATURE CANDIDATE")
    lines = ["# The omni pipeline", "",
             f"**{len(rows) - len(downgrades)} of {len(rows)} records keep the "
             f"disposition their run gave them. {claims} claim a confirmation, "
             f"and {len(escalations)} of those would be unsupported.**", "",
             "Every other script here answers one question. This one orders the "
             "questions - measurement, calibration, nonlinear time, phi, "
             "interpretation, confirmation - and then checks the claimed "
             "disposition against the gates underneath it. It cannot raise a "
             "disposition; it can leave one alone or lower it.", "",
             "## The window guard", "",
             f"The survey's window is {WINDOW_CHANNELS} channels by {SAMPLES} "
             f"samples, and the stabilizer asks whether that window can hold the "
             f"trail a candidate's drift lays down. At this geometry the honest "
             f"drift range is {limit / 1e3:.1f} kHz/s, while the dedoppler "
             f"evaluates trials out to the fold at {nyquist / 1e3:.1f} kHz/s: a "
             f"rate between those two is searched, found, and measured through a "
             f"window that loses part of it.", "",
             "| date | MHz | cycles/sample | trail | held | amplitude lost | verdict |",
             "| --- | --- | --- | --- | --- | --- | --- |"]
    for row in sorted(rows, key=lambda r: (r["record"].observation_date,
                                           r["record"].frequency_hz)):
        rec, win = row["record"], row["window"]
        lines.append(
            f"| {rec.observation_date} | {rec.frequency_hz / 1e6:.4f} | "
            f"{win.channel_rate:+.4f} | {win.trail_channels:.0f} | "
            f"{win.fraction_held:.1%} | {win.significance_loss:.1%} | "
            f"{win.verdict} |")
    lines += ["",
              "The 20240325 group holds at 91.5% of its trail, so its sigma "
              "carries a 4.3% truncation footnote and no more. The 20240920 "
              "trio sits at 0.4965 cycles per sample, holds 19% of its trail, "
              "and is the set this geometry is least able to measure - which is "
              "another reason its strength says nothing about its origin.", "",
              "The measured cost of that loss, as drift rises past the window's "
              "limit, is in DRIFT_SENSITIVITY.md: the guard's sqrt(held) form "
              "overstates the loss, so a window it passes keeps at least the "
              "significance it claims.", "",
              "## The gate matrix", "",
              "`.` pass, `F` fail, `U` untestable with this product, `-` not run.",
              "",
              "| date | MHz | " + " | ".join(gate_ids) + " | claimed | permitted |",
              "| --- | --- | " + " | ".join("---" for _ in gate_ids) + " | --- | --- |"]
    for row in sorted(rows, key=lambda r: (r["record"].observation_date,
                                           r["record"].frequency_hz)):
        rec, res = row["record"], row["audit"]
        marks = []
        for gate_id in gate_ids:
            gate = res.gate(gate_id)
            marks.append("-" if gate is None else MARKS[gate.status])
        lines.append(f"| {rec.observation_date} | {rec.frequency_hz / 1e6:.4f} | "
                     + " | ".join(marks)
                     + f" | {res.claimed} | {res.permitted} |")

    lines += ["", "## Where the claims were lowered", ""]
    if not downgrades:
        lines.append("None: every claimed disposition was supported by its own gates.")
    for record, result in downgrades:
        lines.append(f"**{record.observation_date} "
                     f"{record.frequency_hz / 1e6:.4f} MHz**: "
                     f"{result.claimed} -> {result.permitted}.")
        for deficit in result.deficits:
            lines.append(f"- {deficit}")
        lines.append("")

    lines += ["", "## phi", "",
              "The golden-ratio gate is computed for every record and is not "
              "allowed to support anything: it is not in any disposition's "
              "requirement list, and citing it raises an error. No record's "
              "disposition depends on whether its rung is close to an integer.",
              "",
              "## What would raise a disposition", "",
              "Nothing here can be repaired by more analysis of this data. A "
              "confirmation needs an off-source scan, a second telescope and a "
              "barycentric correction, and the product has one disk-integrated "
              "beam, one instrument and no world coordinate system. Those three "
              "gates are `UNTESTABLE` rather than `FAIL`, which is the difference "
              "between a method that found nothing and a method that was never "
              "run."]

    with open(os.path.join(ROOT, "OMNI.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run every layer in order and audit the disposition each "
                    "candidate claims.")
    parser.add_argument("--strict", action="store_true",
                        help="exit non-zero on any lowered disposition, not only "
                             "on an unsupported confirmation")
    args = parser.parse_args(argv)
    return run(args.strict)


if __name__ == "__main__":
    raise SystemExit(main())
