"""The rung law on the FRB sample: invariants pinned, description computed.

Only the things that must never move are asserted:

  * the shape of the sample (62 cited rows, 2 cluster hosts, every row cited);
  * the round trip, phi^A = z * DM, for every row;
  * the registered anchors, target, tolerance, predicted window and control,
    read from prereg/frb_total_address_lock.json and checked against the data;
  * the telescoping identity: the distance between two anchored addresses is
    the sum of the steps between them, whatever bursts the sample gains.

Everything that depends on which bursts are in the sample -- spans, median and
mean steps, how many steps backtrack, how many land within a tolerance of a
whole rung, the permutation null's mean and p -- is computed by build_summary()
instead of being written down. A new burst moves those numbers; it must not
break the suite.

The prereg's recorded readings are the same computation. recorded_readings()
names each one by its path in prereg/frb_total_address_lock.json, so the record
cannot drift from the sample: if the two disagree the suite fails, and
`python3 tests/test_frb_lattice.py --sync-record` rewrites the record from the
computation rather than by hand. The three readings that describe the 25-row
sample at registration are historical, cannot be recomputed from the current
CSV, and are marked as such in the record instead of being tied.

The record's sentences are written the same way: prereg_prose() fills
'why_this_pair', 'reading' and 'prior_note' from the same readings, so a quoted
number cannot outlive the sample that produced it. A burst that was checked and
did not score is logged under 'external_checks' with its address recomputed
from the published redshift and DM, and its status stays open.
"""

from __future__ import annotations

import csv
import json
import math
import os
import random
import re
import statistics
import sys
import unittest
from decimal import Decimal, getcontext

getcontext().prec = 40
PHI = (1 + Decimal(5).sqrt()) / 2
LN_PHI = PHI.ln()

SAMPLES = os.path.join("results", "frb_localized_sample.csv")
PREREG = os.path.join("prereg", "frb_total_address_lock.json")
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ROWS = 62

# The sample as it stood at registration. These three numbers cannot be
# recomputed from the current CSV, so they live here once, the record's *_before
# fields are checked against them, and the prose quotes them from here.
HISTORY = {
    "rows": 25,
    "mean_step": 1.1231,
    "null_mean": 0.39,
    "null_p": 0.33,
    "first_window": 0.007,
}
BASE_ROWS = HISTORY["rows"]
CLUSTERS = 2
REGISTERED_FROM = "20190711A"
REGISTERED_TO = "20190523A"
TOLERANCES = (0.01, 0.05, 0.10)
DRAWS = 20000
SEED = 61


def rung(x: float | str | Decimal) -> Decimal:
    """r(x) = ln(x) / ln(phi), carried at 40 digits."""
    return Decimal(str(x)).ln() / LN_PHI


def rung_fast(x: float) -> float:
    """Same map in floating point, for the resampling null."""
    return math.log(x) / math.log(1.618033988749895)


def load_record(root: str) -> dict:
    """Read the prereg."""
    with open(os.path.join(root, PREREG), encoding="utf-8") as handle:
        return json.load(handle)


def load_sample(root: str) -> list[dict]:
    """Read the committed sample and attach its rung addresses."""
    path = os.path.join(root, SAMPLES)
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            z = float(row["z"])
            dm = float(row["dm_pc_cm3"])
            rows.append({
                "name": row["name"],
                "z": z,
                "dm": dm,
                "reference": row["reference"],
                "cluster": row["cluster_host"] == "1",
                "r_z": rung(z),
                "r_dm": rung(dm),
                "A": rung(z) + rung(dm),
            })
    rows.sort(key=lambda item: item["z"])
    return rows


def whole(gap: float, eps: float, floor: int = 1) -> bool:
    """Is the gap within eps of a whole rung of at least `floor`?"""
    nearest = round(gap)
    return nearest >= floor and abs(gap - nearest) <= eps


def permutation_null(rows: list[dict], eps: float, draws: int = DRAWS,
                     seed: int = SEED) -> dict:
    """Shuffle DM across the z-order and count one-rung steps."""
    generator = random.Random(seed)
    dms = [row["dm"] for row in rows]
    z_rungs = [rung_fast(row["z"]) for row in rows]
    counts = []
    for _ in range(draws):
        shuffled = dms[:]
        generator.shuffle(shuffled)
        addresses = [z_rungs[i] + rung_fast(shuffled[i]) for i in range(len(dms))]
        gaps = [addresses[i + 1] - addresses[i] for i in range(len(addresses) - 1)]
        counts.append(sum(1 for gap in gaps if abs(gap - 1.0) <= eps))
    return {
        "mean": statistics.mean(counts),
        "sd": statistics.pstdev(counts),
        "p": sum(1 for c in counts if c >= 1) / draws,
        "draws": draws,
    }


def window_hits(rows: list[dict], record: dict) -> dict:
    """How many sample rows sit inside each registered window."""

    def inside(bounds) -> int:
        low, high = (float(value) for value in bounds)
        return sum(1 for row in rows if low <= float(row["A"]) <= high)

    return {
        "prediction": inside(record["prediction"]["window_A"]),
        "control": inside(record["control_window"]["window_A"]),
    }


def build_summary(rows: list[dict], record: dict | None = None) -> dict:
    """Every sample-dependent reading, computed from the rows as given."""
    if record is None:
        record = load_record(REPO_ROOT)
    z_gaps = [float(rows[i + 1]["r_z"] - rows[i]["r_z"]) for i in range(len(rows) - 1)]
    a_gaps = [float(rows[i + 1]["A"] - rows[i]["A"]) for i in range(len(rows) - 1)]
    names = [row["name"] for row in rows]
    start, end = names.index(REGISTERED_FROM), names.index(REGISTERED_TO)
    summary = {
        "rows": len(rows),
        "gaps": len(a_gaps),
        "clusters": sum(1 for row in rows if row["cluster"]),
        "z_span": float(rows[-1]["r_z"] - rows[0]["r_z"]),
        "z_median_gap": statistics.median(z_gaps),
        "a_span": float(rows[-1]["A"] - rows[0]["A"]),
        "a_first": float(rows[0]["A"]),
        "a_last": float(rows[-1]["A"]),
        "a_mean_gap": statistics.mean(a_gaps),
        "a_median_gap": statistics.median(a_gaps),
        "a_backsteps": sum(1 for gap in a_gaps if gap < 0),
        "a_min_gap": min(a_gaps),
        "correlation": statistics.correlation(
            [float(row["r_z"]) for row in rows],
            [float(row["r_dm"]) for row in rows],
        ),
        "registered_index_gap": end - start,
        "registered_pair_adjacent": end - start == 1,
        "registered_span_sum": math.fsum(a_gaps[start:end]),
        "registered_min_gap": min(a_gaps[start:end]) if end > start else 0.0,
        "hits": {},
        "one_rung_pairs": [],
        "windows": window_hits(rows, record),
    }
    for eps in TOLERANCES:
        summary["hits"][eps] = {
            "exactly_one": sum(1 for gap in a_gaps if abs(gap - 1.0) <= eps),
            "any_whole": sum(1 for gap in a_gaps if whole(gap, eps, 1)),
            "including_zero": sum(1 for gap in a_gaps if whole(gap, eps, 0)),
        }
    summary["one_rung_pairs"] = [
        (rows[i]["name"], rows[i + 1]["name"], a_gaps[i])
        for i in range(len(a_gaps))
        if abs(a_gaps[i] - 1.0) <= TOLERANCES[1]
    ]
    summary["null"] = permutation_null(rows, TOLERANCES[1])
    tolerance = float(record["prediction"]["tolerance_rungs"])
    summary["window_width_rungs"] = 2 * tolerance
    summary["prior_percent"] = round(100 * 2 * tolerance, 2)
    return summary


def recorded_readings(summary: dict) -> dict:
    """The prereg's readings, recomputed. Keyed by their path in the record."""
    return {
        "sample_extension.sample_size": summary["rows"],
        "sample_extension.rows_added": summary["rows"] - BASE_ROWS,
        "sample_extension.findings.registered_pair_still_adjacent":
            summary["registered_pair_adjacent"],
        "sample_extension.findings.bursts_now_between_the_pair":
            summary["registered_index_gap"] - 1,
        "sample_extension.findings.intervening_gaps":
            summary["registered_index_gap"],
        "sample_extension.findings.step_recovered_as_sum_of_the_gaps":
            round(summary["registered_span_sum"], 6),
        "sample_extension.findings.largest_backstep_inside_the_step":
            round(summary["registered_min_gap"], 3),
        "sample_extension.findings.predicted_window_hits":
            summary["windows"]["prediction"],
        "sample_extension.findings.control_window_hits":
            summary["windows"]["control"],
        "sample_extension.findings.mean_consecutive_step_after":
            round(summary["a_mean_gap"], 4),
        "sample_extension.findings.null_one_rung_mean_after":
            round(summary["null"]["mean"], 2),
        "sample_extension.findings.null_one_rung_p_after":
            round(summary["null"]["p"], 2),
        "prior_per_window": round(summary["window_width_rungs"], 4),
    }


def read_path(record: dict, path: str):
    """Follow a dotted path into the prereg."""
    node = record
    for key in path.split("."):
        node = node[key]
    return node


def address(z, dm) -> dict:
    """One burst's address, from its redshift and dispersion measure."""
    redshift, dispersion = Decimal(str(z)), Decimal(str(dm))
    return {
        "A": rung(redshift) + rung(dispersion),
        "z_times_dm": redshift * dispersion,
    }


WORDS = {
    1: "one", 2: "two", 3: "three", 4: "four", 5: "five",
    6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten",
}


def words(number: int) -> str:
    """five, not 5, where a sentence needs it."""
    return WORDS.get(number, str(number))


def prereg_prose(summary: dict, record: dict) -> dict:
    """The record's sentences, written from the checked readings.

    Every number in them is formatted from build_summary() or from HISTORY, so
    a sentence cannot quote a reading the sample no longer supports.
    """
    step = summary["registered_span_sum"]
    deviation = abs(step - 1.0)
    gaps = summary["registered_index_gap"]
    backstep = abs(summary["registered_min_gap"])
    observed = summary["hits"][TOLERANCES[1]]["exactly_one"]
    width = summary["window_width_rungs"]
    percent = summary["prior_percent"]
    anchor = record["anchor"]
    return {
        "prior_note": (
            f"Each window is {width:.4f} rungs wide; uniform placement along "
            f"the lattice gives a {percent:.2f}% chance per window per event."
        ),
        "why_this_pair": [
            f"In the {HISTORY['rows']}-row localized-FRB sample ordered by redshift, the total address",
            f"A = r(z) + r(DM) = r(z*DM) steps {step:.6f} rungs from FRB {anchor['from']['name']} to",
            f"FRB {anchor['to']['name']}. That is the closest approach to the Fibonacci step r(phi) = 1.000",
            f"anywhere in the campaign: {deviation:.6f} rungs, tighter than the {HISTORY['first_window']:.3f} window that",
            f"began it. Under the permutation null for that sample, the chance of at least one such gap was {HISTORY['null_p']:.2f},",
            "so the pair is a candidate, not a detection. This document fixes it as a",
            "prediction before further data arrive.",
        ],
        "reading": (
            f"The endpoints are fixed, so the {step:.6f}-rung distance survives; what the "
            f"extension removes is its adjacency. The sequence now reaches that "
            f"distance over {words(gaps)} gaps and backtracks {backstep:.2f} rungs on the way, "
            f"and the null mean for a one-rung step rose from {HISTORY['null_mean']:.2f} to "
            f"{summary['null']['mean']:.2f} while the observed count stayed at {observed}. "
            f"The registration stands where it was put; it is simply no longer the "
            f"best candidate in the sample."
        ),
    }


def external_check_readings(entry: dict, record: dict) -> dict:
    """The computed readings of a logged burst, keyed by their field names."""
    base = address(entry["z"], entry["dm_pc_cm3"])
    top = address(
        Decimal(str(entry["z"])) + Decimal(str(entry["z_uncertainty"])),
        entry["dm_pc_cm3"],
    )  # the far edge of the quoted redshift range, same DM
    target = Decimal(str(record["prediction"]["target_A"]))
    floor = Decimal(str(record["anchor"]["to"]["A"]))
    return {
        "total_address": round(float(base["A"]), 6),
        "redshift_times_dm": round(float(base["z_times_dm"]), 3),
        "total_address_at_top_of_redshift": round(float(top["A"]), 6),
        "offset_from_target_rungs": round(float(base["A"] - target), 6),
        "offset_at_top_of_range_rungs": round(float(top["A"] - target), 6),
        "below_sample_floor_rungs": round(float(floor - base["A"]), 6),
    }


def _same(text: str, value) -> bool:
    """Does this literal already say the computed value?"""
    try:
        parsed = json.loads(text)
    except ValueError:
        return False
    return type(parsed) is type(value) and parsed == value


PROSE_PATHS = {
    "prior_note": "prior_note",
    "why_this_pair": "why_this_pair",
    "reading": "sample_extension.reading",
}


def _replace_scalar(
    text: str, key: str, value, where: str, occurrence: int = 0, expected: int = 1
):
    """Swap one number or boolean on its own line.

    Keys inside `external_checks` repeat once per logged burst, so the caller
    says how many lines it expects and which one it means.
    """
    pattern = re.compile(
        r'(^\s*"%s":\s*)([^,\n]+?)(,?)$' % re.escape(key), re.MULTILINE
    )
    found = list(pattern.finditer(text))
    if len(found) != expected:
        raise SystemExit(
            f"{key}: expected {expected} line(s) in {where}, found {len(found)}"
        )
    match = found[occurrence]
    if _same(match.group(2), value):
        return text, None
    was = json.loads(match.group(2))
    return text[:match.start(2)] + json.dumps(value) + text[match.end(2):], was


def _render_strings(items: list[str], indent: int) -> str:
    """A list of strings in the record's own hand-set layout."""
    pad = " " * indent
    body = ",\n".join(f"{pad}  {json.dumps(item, ensure_ascii=False)}" for item in items)
    return "[\n" + body + "\n" + pad + "]"


def _replace_block(text: str, key: str, value, where: str):
    """Swap a string or a list of strings, whatever it spans."""
    head = re.search(
        r'^(?P<indent>[ \t]*)"%s":' % re.escape(key), text, re.MULTILINE
    )
    if head is None:
        raise SystemExit(f"{key}: no line in {where}")
    start = head.end()
    while start < len(text) and text[start] in " \t\r\n":
        start += 1
    current, end = json.JSONDecoder().raw_decode(text, start)
    if current == value:
        return text, None
    rendered = (
        json.dumps(value, ensure_ascii=False) if isinstance(value, str)
        else _render_strings(value, len(head.group("indent")))
    )
    return text[:start] + rendered + text[end:], current


def sync_record(path: str, summary: dict, record: dict) -> list[tuple[str, object, object]]:
    """Rewrite the readings and the sentences that quote them, nothing else."""
    with open(path, encoding="utf-8") as handle:
        text = handle.read()
    changes = []

    scalars = [(dotted, value, 0, 1) for dotted, value in recorded_readings(summary).items()]
    logged = record.get("external_checks", [])
    for position, entry in enumerate(logged):
        for key, value in external_check_readings(entry, record).items():
            scalars.append(
                (f"external_checks[{entry['name']}].{key}", value, position, len(logged))
            )
    for dotted, value, occurrence, expected in scalars:
        text, was = _replace_scalar(
            text, dotted.rsplit(".", 1)[-1], value, path, occurrence, expected
        )
        if was is not None:
            changes.append((dotted, was, value))

    for key, value in prereg_prose(summary, record).items():
        text, was = _replace_block(text, key, value, path)
        if was is not None:
            changes.append((PROSE_PATHS[key], was, value))

    if changes:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
    return changes


class FrbInvariantTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cls.root = root
        cls.rows = load_sample(root)
        cls.record = load_record(root)
        cls.summary = build_summary(cls.rows, cls.record)
        cls.a_gaps = [
            float(cls.rows[i + 1]["A"] - cls.rows[i]["A"])
            for i in range(len(cls.rows) - 1)
        ]

    # ---- invariants: the shape of the sample ------------------------------

    def test_the_sample_has_the_committed_shape(self):
        self.assertEqual(len(self.rows), ROWS)
        self.assertEqual(len(self.a_gaps), ROWS - 1)
        self.assertEqual(self.summary["clusters"], CLUSTERS)
        self.assertTrue(all(row["reference"] for row in self.rows))
        self.assertEqual(len({row["name"] for row in self.rows}), ROWS)

    def test_every_address_round_trips(self):
        for row in self.rows:
            product = Decimal(str(row["z"])) * Decimal(str(row["dm"]))
            recovered = PHI ** row["A"]
            self.assertLess(
                abs(recovered - product) / product, Decimal("1e-35"),
                msg=f"{row['name']} does not round-trip",
            )

    # ---- invariants: the registration -------------------------------------

    def test_the_registered_anchors_are_exact(self):
        pair = {row["name"]: row for row in self.rows}
        anchor = self.record["anchor"]
        self.assertIn(REGISTERED_FROM, pair)
        self.assertIn(REGISTERED_TO, pair)
        self.assertAlmostEqual(
            float(pair[REGISTERED_FROM]["A"]), anchor["from"]["A"], places=6
        )
        self.assertAlmostEqual(
            float(pair[REGISTERED_TO]["A"]), anchor["to"]["A"], places=6
        )

    def test_the_window_is_built_from_the_anchor_not_from_the_sample(self):
        pair = {row["name"]: row for row in self.rows}
        anchor_to = Decimal(str(self.record["anchor"]["to"]["A"]))
        target = Decimal(str(self.record["prediction"]["target_A"]))
        tol = Decimal(str(self.record["prediction"]["tolerance_rungs"]))
        low, high = (Decimal(str(v)) for v in self.record["prediction"]["window_A"])
        control = Decimal(str(self.record["control_window"]["target_A"]))
        self.assertEqual(target, anchor_to + 1)
        self.assertEqual(low, target - tol)
        self.assertEqual(high, target + tol)
        self.assertEqual(control, target - Decimal("0.5"))
        # the tolerance is the step's own deviation from a whole rung
        step = Decimal(str(self.record["anchor"]["observed_step_rungs"]))
        self.assertEqual(tol, abs(step - 1).quantize(Decimal("0.0001")))
        self.assertAlmostEqual(float(pair[REGISTERED_TO]["A"] + 1), float(target), places=6)

    def test_the_registered_distance_survives_any_extension(self):
        # telescoping: the distance between two anchored addresses is the sum
        # of the steps between them, however many bursts now sit in between
        names = [row["name"] for row in self.rows]
        start = names.index(REGISTERED_FROM)
        end = names.index(REGISTERED_TO)
        self.assertGreater(end, start)
        self.assertAlmostEqual(
            math.fsum(self.a_gaps[start:end]),
            float(self.rows[end]["A"] - self.rows[start]["A"]),
            places=9,
        )
        self.assertAlmostEqual(
            float(self.rows[end]["A"] - self.rows[start]["A"]),
            self.record["anchor"]["observed_step_rungs"],
            places=6,
        )

    def test_the_record_keeps_its_shape(self):
        self.assertEqual(self.record["status"], "open, no event scored")
        self.assertIn("sample_extension", self.record)
        self.assertIn("reporting_addendum", self.record)
        self.assertEqual(
            self.record["reporting_addendum"]["worked_example"]["z_times_dm"],
            [309.15942, 502.128],
        )

    # ---- invariants that hold of the computed statistics -------------------

    def test_the_steps_sum_to_the_span(self):
        self.assertAlmostEqual(
            math.fsum(self.a_gaps), self.summary["a_span"], places=9
        )

    def test_the_tolerances_order_the_hits(self):
        hits = self.summary["hits"]
        self.assertLessEqual(hits[0.01]["exactly_one"], hits[0.05]["exactly_one"])
        self.assertLessEqual(hits[0.05]["exactly_one"], hits[0.10]["exactly_one"])
        self.assertLessEqual(hits[0.05]["any_whole"], hits[0.10]["any_whole"])
        for eps in TOLERANCES:
            self.assertLessEqual(hits[eps]["any_whole"], self.summary["gaps"])

    def test_the_summary_is_internally_consistent(self):
        summary = self.summary
        self.assertEqual(summary["rows"], len(self.rows))
        self.assertLessEqual(-1.0, summary["correlation"])
        self.assertGreaterEqual(1.0, summary["correlation"])
        self.assertLessEqual(summary["a_min_gap"], summary["a_median_gap"])
        self.assertAlmostEqual(
            summary["registered_span_sum"],
            self.record["anchor"]["observed_step_rungs"],
            places=6,
        )

    # ---- the record's readings follow the computation ---------------------

    def test_the_record_matches_the_computed_readings(self):
        for dotted, computed in recorded_readings(self.summary).items():
            self.assertEqual(
                read_path(self.record, dotted), computed,
                msg=f"{dotted} has drifted from the sample; run "
                    f"python3 tests/test_frb_lattice.py --sync-record",
            )

    def test_the_prose_quotes_the_computed_readings(self):
        prose = prereg_prose(self.summary, self.record)
        self.assertEqual(self.record["prior_note"], prose["prior_note"])
        self.assertEqual(self.record["why_this_pair"], prose["why_this_pair"])
        self.assertEqual(
            self.record["sample_extension"]["reading"], prose["reading"]
        )

    def test_the_prior_follows_the_registered_window(self):
        width = 2 * float(self.record["prediction"]["tolerance_rungs"])
        self.assertAlmostEqual(self.summary["window_width_rungs"], width, places=12)
        self.assertAlmostEqual(self.record["prior_per_window"], width, places=12)
        self.assertIn(f"{width:.4f} rungs wide", self.record["prior_note"])

    def test_history_is_recorded_once(self):
        findings = self.record["sample_extension"]["findings"]
        self.assertEqual(findings["mean_consecutive_step_before"], HISTORY["mean_step"])
        self.assertEqual(findings["null_one_rung_mean_before"], HISTORY["null_mean"])
        self.assertEqual(findings["null_one_rung_p_before"], HISTORY["null_p"])
        self.assertEqual(
            self.record["sample_extension"]["sample_size"], ROWS
        )
        self.assertEqual(self.record["sample_extension"]["rows_added"], ROWS - BASE_ROWS)

    def test_every_logged_burst_is_checked_and_unscored(self):
        # The count is not pinned: logging another checked burst is allowed, so
        # long as it too is recomputed, unscored and outside the sample.
        entries = self.record["external_checks"]
        self.assertGreaterEqual(len(entries), 1)
        logged = [entry["name"] for entry in entries]
        self.assertEqual(len(logged), len(set(logged)), msg="two bursts logged as one")
        sample = {row["name"] for row in self.rows}
        window = self.record["prediction"]["window_A"]
        self.assertEqual(self.record["status"], "open, no event scored")
        for entry in entries:
            self.assertNotIn(entry["name"], sample)  # checked, not admitted
            self.assertIs(entry["scores"], False)
            self.assertEqual(entry["z_kind"], "photometric")
            self.assertLess(entry["total_address_at_top_of_redshift"], window[0])
            for key, computed in external_check_readings(entry, self.record).items():
                self.assertEqual(
                    entry[key], computed,
                    msg=f"external_checks[{entry['name']}].{key} has drifted from the computation",
                )

    def test_the_record_marks_what_it_does_not_check(self):
        findings = self.record["sample_extension"]["findings"]
        self.assertIn("prior_sample_note", findings)
        for key in findings:
            if key.endswith("_before"):
                self.assertIn("cannot be recomputed", findings["prior_sample_note"])

    def test_the_null_is_deterministic(self):
        again = permutation_null(self.rows, TOLERANCES[1], draws=500, seed=SEED)
        once = permutation_null(self.rows, TOLERANCES[1], draws=500, seed=SEED)
        self.assertEqual(again, once)
        self.assertLessEqual(0.0, self.summary["null"]["p"])
        self.assertGreaterEqual(1.0, self.summary["null"]["p"])
        self.assertGreaterEqual(self.summary["null"]["mean"], 0.0)

    # ---- the coding that must not drift -----------------------------------

    def test_the_coding_discards_the_cluster_hosts(self):
        kept = [row for row in self.rows if not row["cluster"]]
        self.assertEqual(len(kept), ROWS - CLUSTERS)
        gaps = [
            float(kept[i + 1]["r_z"] - kept[i]["r_z"])
            for i in range(len(kept) - 1)
        ]
        self.assertEqual(len(gaps), ROWS - CLUSTERS - 1)


if __name__ == "__main__":
    if "--sync-record" in sys.argv:
        readings = build_summary(load_sample(REPO_ROOT))
        changed = sync_record(
            os.path.join(REPO_ROOT, PREREG), readings, load_record(REPO_ROOT)
        )
        for dotted, was, now in changed:
            print(f"{dotted}: {was} -> {now}")
        print(
            f"{len(changed)} reading(s) rewritten" if changed
            else "the record already matches the computation"
        )
        raise SystemExit(0)
    summary = build_summary(load_sample(REPO_ROOT))
    print(json.dumps(summary, indent=2, default=float))
    unittest.main()
