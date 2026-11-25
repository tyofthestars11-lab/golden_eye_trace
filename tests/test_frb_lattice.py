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
"""

from __future__ import annotations

import csv
import json
import math
import os
import random
import statistics
import unittest
from decimal import Decimal, getcontext

getcontext().prec = 40
PHI = (1 + Decimal(5).sqrt()) / 2
LN_PHI = PHI.ln()

SAMPLES = os.path.join("results", "frb_localized_sample.csv")
PREREG = os.path.join("prereg", "frb_total_address_lock.json")

ROWS = 62
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


def build_summary(rows: list[dict]) -> dict:
    """Every sample-dependent reading, computed from the rows as given."""
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
        "registered_span_sum": math.fsum(a_gaps[start:end]),
        "hits": {},
        "one_rung_pairs": [],
        "registered_rows_in_windows": 0,
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
    return summary


class FrbInvariantTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cls.root = root
        cls.rows = load_sample(root)
        with open(os.path.join(root, PREREG), encoding="utf-8") as handle:
            cls.record = json.load(handle)
        cls.summary = build_summary(cls.rows)
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
        self.assertLessEqual(summary["a_median_gap"], summary["a_mean_gap"] * 2)
        self.assertAlmostEqual(
            summary["registered_span_sum"],
            self.record["anchor"]["observed_step_rungs"],
            places=6,
        )

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
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    summary = build_summary(load_sample(root))
    print(json.dumps(summary, indent=2, default=float))
    unittest.main()
