"""Merge every lane into one object, and say what that object can and cannot be.

    python3 scripts/super_tetrahedron.py

Writes ``SUPER_TETRAHEDRON.md`` and ``results/super_tetrahedron.json``.

The merge is defined explicitly, because an undefined merge is a mood rather
than a measurement. Four lanes, each reduced to the φ-rung of a declared
quantity:

    radio        median rung of the candidate frequencies        (Hz)
    neutrino     median rung of the ingested alert energies      (GeV)
    local field  median rung of the bias field the band requires (T)
    phi geometry rung of the frozen target 1/phi^2               (dimensionless)

Four vertices, six pairwise rung separations. A **regular simplex** - the
tetrahedron - requires all six to be equal, so that is the testable content of
"tetrahedral": equal bond lengths among four vertices.

The angular part of the word is not testable here, and the report says so
plainly: a tetrahedron's 109.4712 degrees is an angle between *directions*, and
the rung ladder is one-dimensional. Worse, any four bond lengths can be drawn
along four directions 109.47 degrees apart, so the tetrahedral geometry can
always be satisfied - the lengths are the only thing carrying information. A
claim that four quantities form a tetrahedral bond therefore reduces to a claim
that four numbers are equally spaced, and that is what is measured.

The null: draw four positions uniformly from the span the four lanes occupy and
ask how often their separations are at least as equal as the observed ones.
"""

from __future__ import annotations

import csv
import json
import math
import os
import random
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.rungs import rung  # noqa: E402

PHI = (1.0 + math.sqrt(5.0)) / 2.0
DRAWS = 20_000
SEED = 20261114

TETRAHEDRAL_DEG = math.degrees(math.acos(-1.0 / 3.0))   # 109.4712206
GOLDEN_ANGLE_DEG = 137.5077641
GOLDEN_PITCH_DEG = 17.0324


def load(path: str):
    with open(os.path.join(ROOT, path), encoding="utf-8") as fh:
        if path.endswith(".csv"):
            return list(csv.DictReader(fh))
        return json.load(fh)


def lane_positions() -> list[dict]:
    records = load("signal_records.json")["records"]
    radio = statistics.median(rung(r["frequency_hz"]) for r in records)

    ingest = load("results/icecube_ingest.json")
    energies = [a["event"]["energy_estimate_gev"] for a in ingest["alerts"]]
    neutrino = statistics.median(rung(e) for e in energies)

    nv = load("results/nv_lane.json")
    bias = [row["bias_field_t"] for row in nv["bias"]]
    local = statistics.median(rung(b) for b in bias)

    target = 1.0 / PHI ** 2
    geometry = rung(target)

    return [
        {"lane": "radio", "unit": "Hz", "n": len(records),
         "rung": radio, "source": "signal_records.json (10 cells)"},
        {"lane": "neutrino", "unit": "GeV", "n": len(energies),
         "rung": neutrino, "source": "results/icecube_ingest.json (2 alerts)"},
        {"lane": "local field", "unit": "T", "n": len(bias),
         "rung": local, "source": "results/nv_lane.json (5 priced cells)"},
        {"lane": "phi geometry", "unit": "dimensionless", "n": 1,
         "rung": geometry, "source": "frozen target 1/phi^2"},
    ]


def simplex(positions: list[float]) -> dict:
    seps = [abs(b - a) for i, a in enumerate(positions)
            for b in positions[i + 1:]]
    mean = statistics.mean(seps)
    sd = statistics.pstdev(seps)
    return {
        "separations": seps,
        "mean": mean,
        "sd": sd,
        "cv": sd / mean if mean else float("inf"),
        "max_min": max(seps) / min(seps) if min(seps) else float("inf"),
        "equal_edge_deviation": max(abs(s - mean) / mean for s in seps),
    }


def line_simplex_impossible(positions: list[float]) -> dict:
    """Four distinct points on a line can never be equidistant.

    Sort them; the outer pair's separation is the sum of the interior gaps and is
    therefore strictly larger than any other separation whenever the points are
    distinct. Equality of all six forces every gap to zero. So the tetrahedral
    condition is not merely unmet in this coordinate, it is unachievable in it,
    which is a stronger statement than a failed test.
    """
    ordered = sorted(positions)
    gaps = [b - a for a, b in zip(ordered, ordered[1:])]
    outer = ordered[-1] - ordered[0]
    largest_other = max(gaps + [0.0])
    return {
        "sorted": ordered,
        "gaps": gaps,
        "outer_separation": outer,
        "sum_of_gaps": sum(gaps),
        "outer_exceeds_largest_gap": outer > largest_other,
        "impossible": True,
        "reason": ("on a line the outer separation is the sum of the interior "
                   "gaps, so it exceeds every other separation unless all gaps "
                   "are zero; a regular simplex of four vertices is therefore "
                   "unreachable in a one-dimensional rung coordinate"),
    }


def null_cv(positions: list[float], draws: int = DRAWS,
            seed: int = SEED) -> dict:
    """How equal are four random positions from the same span?"""
    lo, hi = min(positions), max(positions)
    rng = random.Random(seed)
    observed = simplex(positions)["cv"]
    cvs = []
    at_least = 0
    for _ in range(draws):
        drawn = [rng.uniform(lo, hi) for _ in range(4)]
        cv = simplex(drawn)["cv"]
        cvs.append(cv)
        if cv <= observed:  # as equal or more equal than observed
            at_least += 1
    return {"observed_cv": observed, "draws": draws,
            "null_median_cv": statistics.median(cvs),
            "p_as_equal_or_more": (at_least + 1) / (draws + 1)}


def angle_facts() -> dict:
    return {
        "tetrahedral_deg": TETRAHEDRAL_DEG,
        "golden_angle_deg": GOLDEN_ANGLE_DEG,
        "difference_deg": GOLDEN_ANGLE_DEG - TETRAHEDRAL_DEG,
        "three_tetrahedral_from_full_circle_deg": 360.0 - 3.0 * TETRAHEDRAL_DEG,
        "golden_pitch_deg": GOLDEN_PITCH_DEG,
        "note": ("these are exact geometric constants; none of them is a "
                 "measurement of any dataset in this repository"),
    }


def results_ledger() -> dict:
    records = load("signal_records.json")["records"]
    dispositions: dict[str, int] = {}
    for record in records:
        dispositions[record["disposition"]] = dispositions.get(
            record["disposition"], 0) + 1

    probe = load("results/phi_probe.json")
    phase_hits = sum(1 for p in probe["probes"] if any(p["hit"].values()))

    cross = load("results/phi_cross_domain_protocol.json")
    best_phi = min(row["median_rung_error"] for row in cross["holdout"]
                   if row["predictor"] == "phi")
    best_fitted = min(row["median_rung_error"] for row in cross["holdout"]
                      if row["predictor"] == "fitted_constant")

    neutrino = load("results/neutrino_crossmatch.json")
    nv = load("results/nv_lane.json")

    comb = load("results/comb_vs_scale.csv")
    comb_row = next(r for r in comb if r["date"] == "20240325")

    return {
        "records": {"n": len(records), "dispositions": dispositions},
        "rung_law": {"round_trip_max_rel_error":
                     max(p["round_trip_relative_error"]
                         for p in probe["probes"])},
        "phi_phase": {"hits": phase_hits, "tested": len(probe["probes"]),
                      "tolerance_deg": probe["tolerances_deg"][1],
                      "null": probe["null_phase"]},
        "phi_cross_domain": {
            "hits": sum(row["hits"] for row in cross["tests"]),
            "tests": len(cross["tests"]),
            "p_values": sorted({row["p_value"] for row in cross["tests"]}),
            "cross_domain_status": cross["cross_domain"]["status"],
            "phi_median_rung_error": best_phi,
            "fitted_median_rung_error": best_fitted,
            "power": cross["power"],
        },
        "neutrino": {
            "n_records": neutrino["summary"]["n_records"],
            "n_testable": neutrino["summary"]["n_testable"],
            "restoration_refusals": len(neutrino["refusals"]),
            "catalog": neutrino["catalog"]["kind"],
        },
        "nv": {
            "n_records": len(nv["records"]),
            "n_testable": sum(1 for r in nv["records"]
                              if r["verdict"] != "NOT TESTABLE"),
            "band_low_mhz": nv["band_hz"][0] / 1e6,
            "band_high_mhz": nv["band_hz"][1] / 1e6,
            "modes": len(nv["feasibility"]),
        },
        "comb_vs_scale": {
            "delta_aic": float(comb_row["delta_aic"]),
            "boot_lo": float(comb_row["boot_lo"]),
            "boot_hi": float(comb_row["boot_hi"]),
            "null_p": float(comb_row["null_p"]),
            "verdict": comb_row["verdict"],
            "span_rungs_phi": float(comb_row["span_rungs_phi"]),
        },
    }


def write_report(payload: dict) -> None:
    lanes = payload["lanes"]
    block = payload["simplex"]
    null = payload["null"]
    angles = payload["angle_facts"]
    ledger = payload["ledger"]

    lines = [
        "# The super-tetrahedral merge",
        "",
        "Framework: φ² = φ + 1. Attribution: Tyree Jones (tyofthestarz), prior "
        "art from December 2025.",
        "",
        "One object built from every lane in this repository, defined explicitly "
        "so that it can fail: four vertices, each the φ-rung of a declared "
        "quantity, and six bond lengths between them.",
        "",
        "## The four vertices",
        "",
        "| lane | quantity | unit | n | rung | source |",
        "| --- | --- | --- | ---: | ---: | --- |",
    ]
    for lane in lanes:
        lines.append(f"| {lane['lane']} | median of the lane's measurements | "
                     f"{lane['unit']} | {lane['n']} | {lane['rung']:+.4f} | "
                     f"{lane['source']} |")
    lines += [
        "",
        "## The six bonds",
        "",
        "| pair | rung separation |",
        "| --- | ---: |",
    ]
    names = [lane["lane"] for lane in lanes]
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            sep = abs(lanes[i]["rung"] - lanes[j]["rung"])
            lines.append(f"| {a} – {names[j]} | {sep:.4f} |")
    lines += [
        "",
        f"Mean bond {block['mean']:.4f} rungs, sd {block['sd']:.4f}, "
        f"coefficient of variation **{block['cv']:.3f}**, longest/shortest "
        f"**{block['max_min']:.1f}×**.",
        "",
        f"A regular simplex needs all six equal. The worst bond deviates from the "
        f"mean by **{block['equal_edge_deviation'] * 100:.0f}%**, and the range "
        f"spans {min(block['separations']):.2f} to "
        f"{max(block['separations']):.2f} rungs. This is not a tetrahedron; it is "
        "a line with four marks on it.",
        "",
        "It could not have been one. Sort the four rungs and the outer "
        f"separation is {payload['impossibility']['outer_separation']:.4f} rungs, "
        "which is exactly the sum of the interior gaps ("
        + " + ".join(f"{g:.4f}" for g in payload['impossibility']['gaps'])
        + "). It therefore exceeds every other separation unless all gaps are "
        "zero, and **four distinct points on a line are never equidistant**. A "
        "regular simplex of four vertices is unachievable in a one-dimensional "
        "coordinate: the failure is structural, not a property of these numbers.",
        "",
        "## What was actually testable",
        "",
        f"Null: four positions drawn uniformly from the observed span "
        f"({min(l['rung'] for l in lanes):.2f} to "
        f"{max(l['rung'] for l in lanes):.2f} rungs), {null['draws']:,} draws.",
        "",
        "| quantity | value |",
        "| --- | ---: |",
        f"| observed coefficient of variation | {null['observed_cv']:.3f} |",
        f"| median CV under the null | {null['null_median_cv']:.3f} |",
        f"| p(separation ≥ as equal as observed) | "
        f"{null['p_as_equal_or_more']:.3f} |",
        "",
        "The observed spread is typical of random placement, so the bond structure "
        "is not hidden by a small sample: there is no equal-edge structure to "
        "find.",
        "",
        "## The angular half of the word",
        "",
        "| constant | degrees |",
        "| --- | ---: |",
        f"| tetrahedral angle | {angles['tetrahedral_deg']:.4f} |",
        f"| golden angle | {angles['golden_angle_deg']:.4f} |",
        f"| difference | {angles['difference_deg']:.4f} |",
        f"| 360 − 3 × tetrahedral | "
        f"{angles['three_tetrahedral_from_full_circle_deg']:.4f} |",
        f"| golden pitch | {angles['golden_pitch_deg']:.4f} |",
        "",
        "These are exact constants and none of them is a measurement. Two further "
        "facts matter: **any** four bond lengths can be drawn along four "
        "directions 109.47° apart, so tetrahedral geometry is always satisfiable "
        "and tests nothing; and a rung ladder is one-dimensional, so it has no "
        "angles at all. The only falsifiable content of \"a tetrahedral bond "
        "between the lanes\" is equality of the six lengths, which fails above by "
        f"a factor of {block['max_min']:.0f}.",
        "",
        "## Results so far, in one place",
        "",
        "| lane | result | null beside it |",
        "| --- | --- | --- |",
        f"| rung law | round-trips to "
        f"{ledger['rung_law']['round_trip_max_rel_error']:.2e} relative error | "
        "arithmetic, not evidence |",
        f"| radio records | {ledger['records']['n']} records: "
        + ", ".join(f"{k} ×{v}" for k, v in
                    sorted(ledger['records']['dispositions'].items()))
        + " | OMNI permits what is claimed; 0 confirmations |",
        f"| comb vs scale | ΔAIC {ledger['comb_vs_scale']['delta_aic']:+.2f}, "
        f"bootstrap [{ledger['comb_vs_scale']['boot_lo']:+.1f}, "
        f"{ledger['comb_vs_scale']['boot_hi']:+.1f}] | "
        f"random-line null p = {ledger['comb_vs_scale']['null_p']:.3f}; "
        f"span {ledger['comb_vs_scale']['span_rungs_phi']:.3f} φ rungs |",
        f"| φ phase test | {ledger['phi_phase']['hits']} of "
        f"{ledger['phi_phase']['tested']} within "
        f"{ledger['phi_phase']['tolerance_deg']:g}° of 137.5° or 17.03° | "
        f"chance {ledger['phi_phase']['null']['1.0']:.3f} per record |",
        f"| cross-domain | {ledger['phi_cross_domain']['hits']} hits over "
        f"{ledger['phi_cross_domain']['tests']} tests, p = "
        f"{ledger['phi_cross_domain']['p_values']} | "
        f"{ledger['phi_cross_domain']['cross_domain_status']}: one domain |",
        f"| held-out prediction | φ {ledger['phi_cross_domain']['phi_median_rung_error']:.3f} "
        f"rungs vs fitted {ledger['phi_cross_domain']['fitted_median_rung_error']:.3f} | "
        "fitted constant wins |",
        f"| neutrino lane | {ledger['neutrino']['n_testable']} of "
        f"{ledger['neutrino']['n_records']} testable, "
        f"{ledger['neutrino']['restoration_refusals']} restorations refused | "
        f"catalog: {ledger['neutrino']['catalog']} |",
        f"| NV lane | {ledger['nv']['n_testable']} of "
        f"{ledger['nv']['n_records']} testable, {ledger['nv']['modes']} sensing "
        f"modes priced across {ledger['nv']['band_low_mhz']:.1f}–"
        f"{ledger['nv']['band_high_mhz']:.1f} MHz | detraction only |",
        "",
        "## Verdict",
        "",
        "**The merge does not produce a bond.** Six lengths spanning a factor of "
        f"{block['max_min']:.0f}, a coefficient of variation typical of random "
        "placement, and no angular content to test. Every lane keeps its own "
        "verdict, and none of them changes because they were drawn in the same "
        "picture.",
        "",
        "What the merge is: a single index of what exists, so that the state of "
        "the whole repository can be read at once. What it is not: evidence for "
        "anything. Seven cells remain MASK-MANUFACTURED SIGNAL, three records "
        "remain UNRESOLVED, the two multimessenger lanes remain untestable for "
        "lack of a clock and a pointing, and φ remains a coordinate that can "
        "only detract.",
        "",
    ]
    with open(os.path.join(ROOT, "SUPER_TETRAHEDRON.md"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> int:
    lanes = lane_positions()
    positions = [lane["rung"] for lane in lanes]
    block = simplex(positions)
    null = null_cv(positions)
    payload = {
        "lanes": lanes,
        "simplex": block,
        "null": null,
        "angle_facts": angle_facts(),
        "impossibility": line_simplex_impossible(positions),
        "ledger": results_ledger(),
        "definition": ("four lanes at the phi-rung of a declared quantity; the "
                       "tetrahedral content is equality of the six separations"),
    }
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "super_tetrahedron.json"), "w",
              encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    write_report(payload)

    print("SUPER TETRAHEDRON")
    print(f"   {'lane':14} {'unit':14} {'rung':>10}")
    for lane in lanes:
        print(f"   {lane['lane']:14} {lane['unit']:14} {lane['rung']:+10.4f}")
    print()
    for i, a in enumerate(lanes):
        for j in range(i + 1, len(lanes)):
            sep = abs(lanes[i]["rung"] - lanes[j]["rung"])
            print(f"   bond {a['lane']:14} – {lanes[j]['lane']:14} {sep:8.4f}")
    print()
    print(f"   mean {block['mean']:.4f}  sd {block['sd']:.4f}  "
          f"CV {block['cv']:.3f}  max/min {block['max_min']:.1f}x")
    print(f"   worst deviation from equal edges "
          f"{block['equal_edge_deviation'] * 100:.0f}%")
    print(f"   null CV median {null['null_median_cv']:.3f}, "
          f"p(as equal) {null['p_as_equal_or_more']:.3f}")
    print("   wrote results/super_tetrahedron.json and SUPER_TETRAHEDRON.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
