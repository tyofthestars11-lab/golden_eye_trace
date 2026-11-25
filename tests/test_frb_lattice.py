"""The rung law on a published FRB sample: addresses, the Fibonacci-step null, and the registered pair.

This is a fixture, not a detection.  It holds the 25-row localized-FRB sample as
data, recomputes every address, and re-runs the two tests that were reported:
the redshift stream alone and the total address A = r(z) + r(DM).  The null is
part of the test so the numbers cannot drift away from the reading they came
from.  The registered pair (20190711A -> 20190523A, one step of 1.007868 rungs)
is pinned here and in prereg/frb_total_address_lock.json; the pair is a
candidate with a 0.33 chance of arising under the permutation null, and it is
committed so that a later event can settle it.
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

# The published readings this fixture must keep reproducing.
Z_SPAN = 20.0294
A_SPAN = 26.9534
A_MEAN_GAP = 1.1231
REGISTERED_STEP = 1.007868003
REGISTERED_FROM = "20190711A"
REGISTERED_TO = "20190523A"


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


def near_integer(gap: float, eps: float, floor: int = 1) -> bool:
    """Is the gap within eps of a whole rung of at least `floor`?"""
    whole = round(gap)
    return whole >= floor and abs(gap - whole) <= eps


class FrbTotalAddressTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cls.root = root
        cls.rows = load_sample(root)
        cls.gaps = [
            float(cls.rows[i + 1]["r_z"] - cls.rows[i]["r_z"])
            for i in range(len(cls.rows) - 1)
        ]
        cls.a_gaps = [
            float(cls.rows[i + 1]["A"] - cls.rows[i]["A"])
            for i in range(len(cls.rows) - 1)
        ]

    def test_the_sample_holds_twenty_five_cited_rows(self):
        self.assertEqual(len(self.rows), 25)
        self.assertEqual(len(self.gaps), 24)
        self.assertTrue(all(row["reference"] for row in self.rows))
        self.assertEqual(sum(1 for row in self.rows if row["cluster"]), 2)

    def test_the_redshift_stream_spans_twenty_rungs(self):
        span = float(self.rows[-1]["r_z"] - self.rows[0]["r_z"])
        self.assertAlmostEqual(span, Z_SPAN, places=3)
        self.assertAlmostEqual(float(self.rows[0]["r_z"]), -18.4407, places=3)
        self.assertAlmostEqual(float(self.rows[-1]["r_z"]), 1.5888, places=3)

    def test_no_consecutive_redshift_pair_makes_a_fibonacci_step(self):
        exact = sum(1 for gap in self.gaps if abs(gap - 1.0) <= 0.05)
        anyn = sum(1 for gap in self.gaps if near_integer(gap, 0.05, 1))
        anyn_all = sum(1 for gap in self.gaps if near_integer(gap, 0.05, 0))
        self.assertEqual(exact, 0)
        self.assertEqual(anyn, 1)
        self.assertEqual(anyn_all, 4)
        # the single n>=1 hit is an n=2 step into a cluster-host burst
        hit = [
            (self.rows[i]["name"], self.rows[i + 1]["name"], gap)
            for i, gap in enumerate(self.gaps)
            if near_integer(gap, 0.05, 1)
        ]
        self.assertEqual(hit[0][0], "20180916B")
        self.assertEqual(hit[0][1], "20220509G")
        self.assertAlmostEqual(hit[0][2], 2.0274, places=3)

    def test_the_total_address_is_the_rung_of_the_product(self):
        for row in self.rows:
            product = Decimal(str(row["z"])) * Decimal(str(row["dm"]))
            recovered = PHI ** row["A"]
            self.assertLess(
                abs(recovered - product) / product, Decimal("1e-35"),
                msg=f"{row['name']} does not round-trip",
            )

    def test_the_total_address_sequence_is_not_monotone(self):
        backsteps = [gap for gap in self.a_gaps if gap < 0]
        self.assertEqual(len(backsteps), 4)
        self.assertAlmostEqual(min(backsteps), -2.8875, places=3)

    def test_the_total_address_stream_steps_about_one_rung(self):
        span = float(self.rows[-1]["A"] - self.rows[0]["A"])
        self.assertAlmostEqual(span, A_SPAN, places=3)
        self.assertAlmostEqual(statistics.mean(self.a_gaps), A_MEAN_GAP, places=3)
        # the mean spacing matches the step being looked for, which is why a
        # 1-rung hit here is not by itself informative
        self.assertAlmostEqual(statistics.mean(self.a_gaps), 1.0, delta=0.2)

    def test_the_registered_pair_is_the_only_one_rung_step(self):
        exact = sum(1 for gap in self.a_gaps if abs(gap - 1.0) <= 0.01)
        self.assertEqual(exact, 1)
        found = [
            (self.rows[i]["name"], self.rows[i + 1]["name"], gap)
            for i, gap in enumerate(self.a_gaps)
            if abs(gap - 1.0) <= 0.01
        ]
        self.assertEqual(found[0][0], REGISTERED_FROM)
        self.assertEqual(found[0][1], REGISTERED_TO)
        self.assertAlmostEqual(found[0][2], REGISTERED_STEP, places=6)

    def test_the_two_streams_are_mostly_one_stream(self):
        r_z = [float(row["r_z"]) for row in self.rows]
        r_dm = [float(row["r_dm"]) for row in self.rows]
        self.assertAlmostEqual(statistics.correlation(r_z, r_dm), 0.838, places=2)

    def test_the_permutation_null_does_not_favour_the_registered_pair(self):
        generator = random.Random(61)
        dms = [row["dm"] for row in self.rows]
        z_rungs = [rung_fast(row["z"]) for row in self.rows]
        draws = 4000
        counts = []
        for _ in range(draws):
            shuffled = dms[:]
            generator.shuffle(shuffled)
            addresses = [z_rungs[i] + rung_fast(shuffled[i]) for i in range(len(dms))]
            gaps = [addresses[i + 1] - addresses[i] for i in range(len(addresses) - 1)]
            counts.append(sum(1 for gap in gaps if abs(gap - 1.0) <= 0.05))
        observed = 1
        mean = statistics.mean(counts)
        p_value = sum(1 for count in counts if count >= observed) / draws
        self.assertAlmostEqual(mean, 0.39, delta=0.12)
        self.assertGreater(p_value, 0.20)
        self.assertLess(p_value, 0.50)

    def test_the_registered_window_matches_the_sample(self):
        path = os.path.join(self.root, PREREG)
        with open(path, encoding="utf-8") as handle:
            record = json.load(handle)
        anchor = record["anchor"]
        pair = {row["name"]: row for row in self.rows}
        self.assertAlmostEqual(
            float(pair[REGISTERED_FROM]["A"]), anchor["from"]["A"], places=6
        )
        self.assertAlmostEqual(
            float(pair[REGISTERED_TO]["A"]), anchor["to"]["A"], places=6
        )
        self.assertAlmostEqual(
            float(pair[REGISTERED_TO]["A"] - pair[REGISTERED_FROM]["A"]),
            anchor["observed_step_rungs"], places=6,
        )
        target = record["prediction"]["target_A"]
        tol = record["prediction"]["tolerance_rungs"]
        self.assertAlmostEqual(
            target, float(pair[REGISTERED_TO]["A"]) + 1.0, places=6
        )
        low, high = record["prediction"]["window_A"]
        self.assertAlmostEqual(low, target - tol, places=6)
        self.assertAlmostEqual(high, target + tol, places=6)
        self.assertEqual(record["status"], "open, no event scored")

    def test_the_coding_discards_the_cluster_hosts(self):
        kept = [row for row in self.rows if not row["cluster"]]
        gaps = [
            float(kept[i + 1]["r_z"] - kept[i]["r_z"])
            for i in range(len(kept) - 1)
        ]
        self.assertEqual(len(kept), 23)
        self.assertEqual(sum(1 for gap in gaps if near_integer(gap, 0.05, 1)), 0)


if __name__ == "__main__":
    unittest.main()
