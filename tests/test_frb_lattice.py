"""The rung law on the extended FRB sample: addresses, the Fibonacci-step null, and the registered pair.

This is a fixture, not a detection.  It holds the 62-row localized-FRB sample
(Table 1 of Perez-Casares et al. 2024 plus the DSA-110 and CHIME/KKO
localizations published after it), recomputes every address, and re-runs the
two tests that were reported: the redshift stream alone and the total address
A = r(z) + r(DM).  The permutation null is part of the test so the numbers
cannot drift away from the reading they came from.

The 2026-11-25 registration (prereg/frb_total_address_lock.json) is kept at the
values it was registered with.  The extended sample was published before that
registration, so it cannot score it; what it does is sharpen the null and show
that the registered pair is no longer adjacent in the z-ordered sequence.
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
Z_MEDIAN_GAP = 0.0718
A_SPAN = 26.9534
A_MEAN_GAP = 0.4419
A_MEDIAN_GAP = 0.3800
A_BACKSTEPS = 20
A_MIN_GAP = -2.2010
CORRELATION = 0.8110
NULL_MEAN = 1.17
REGISTERED_STEP = 1.007868003
REGISTERED_FROM = "20190711A"
REGISTERED_TO = "20190523A"
CANDIDATE_FROM = "20231123B"
CANDIDATE_TO = "20230307A"
CANDIDATE_STEP = 0.9609


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

    def test_the_extension_holds_sixty_two_cited_rows(self):
        self.assertEqual(len(self.rows), 62)
        self.assertEqual(len(self.gaps), 61)
        self.assertTrue(all(row["reference"] for row in self.rows))
        self.assertEqual(sum(1 for row in self.rows if row["cluster"]), 2)
        sources = {row["reference"] for row in self.rows}
        self.assertIn("DSA-110 2026", sources)
        self.assertIn("CHIME/KKO 2025", sources)

    def test_the_redshift_stream_spans_twenty_rungs(self):
        span = float(self.rows[-1]["r_z"] - self.rows[0]["r_z"])
        self.assertAlmostEqual(span, Z_SPAN, places=3)
        self.assertAlmostEqual(float(self.rows[0]["r_z"]), -18.4407, places=3)
        self.assertAlmostEqual(float(self.rows[-1]["r_z"]), 1.5888, places=3)
        # the extension quadruples the density but cannot move the extremes
        self.assertAlmostEqual(statistics.median(self.gaps), Z_MEDIAN_GAP, places=3)

    def test_the_redshift_stream_still_has_no_fibonacci_step(self):
        self.assertEqual(sum(1 for g in self.gaps if abs(g - 1.0) <= 0.01), 0)
        self.assertEqual(sum(1 for g in self.gaps if near_integer(g, 0.01, 1)), 0)
        # five pairs sit within 3.4% in z of their neighbour: crowding, not a step
        self.assertEqual(sum(1 for g in self.gaps if near_integer(g, 0.01, 0)), 5)
        self.assertEqual(sum(1 for g in self.gaps if near_integer(g, 0.05, 1)), 1)
        self.assertEqual(sum(1 for g in self.gaps if near_integer(g, 0.05, 0)), 26)

    def test_the_total_address_is_the_rung_of_the_product(self):
        for row in self.rows:
            product = Decimal(str(row["z"])) * Decimal(str(row["dm"]))
            recovered = PHI ** row["A"]
            self.assertLess(
                abs(recovered - product) / product, Decimal("1e-35"),
                msg=f"{row['name']} does not round-trip",
            )

    def test_the_total_address_sequence_is_mostly_backwards(self):
        backsteps = [gap for gap in self.a_gaps if gap < 0]
        self.assertEqual(len(backsteps), A_BACKSTEPS)
        self.assertAlmostEqual(min(backsteps), A_MIN_GAP, places=3)
        span = float(self.rows[-1]["A"] - self.rows[0]["A"])
        self.assertAlmostEqual(span, A_SPAN, places=3)
        self.assertAlmostEqual(float(self.rows[0]["A"]), -9.1406, places=3)
        self.assertAlmostEqual(float(self.rows[-1]["A"]), 17.8128, places=3)

    def test_the_extension_crowds_the_total_address(self):
        # 62 points in the same 26.95-rung span: the mean step halves, and a
        # 1-rung step is no longer the natural spacing of the sequence
        self.assertAlmostEqual(statistics.mean(self.a_gaps), A_MEAN_GAP, places=3)
        self.assertAlmostEqual(statistics.median(self.a_gaps), A_MEDIAN_GAP, places=3)
        self.assertLess(statistics.mean(self.a_gaps), 0.5)

    def test_the_registered_pair_is_no_longer_adjacent(self):
        names = [row["name"] for row in self.rows]
        start = names.index(REGISTERED_FROM)
        end = names.index(REGISTERED_TO)
        self.assertEqual(end - start, 5)
        # the endpoint addresses are unchanged, so the distance survives as a
        # sum over five gaps even though the single step does not
        self.assertAlmostEqual(
            math.fsum(self.a_gaps[start:end]), REGISTERED_STEP, places=6
        )
        self.assertEqual(sum(1 for g in self.a_gaps if abs(g - 1.0) <= 0.01), 0)

    def test_the_new_candidate_is_a_different_pair(self):
        hits = [
            (self.rows[i]["name"], self.rows[i + 1]["name"], gap)
            for i, gap in enumerate(self.a_gaps)
            if abs(gap - 1.0) <= 0.05
        ]
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0][0], CANDIDATE_FROM)
        self.assertEqual(hits[0][1], CANDIDATE_TO)
        self.assertAlmostEqual(hits[0][2], CANDIDATE_STEP, places=3)

    def test_the_two_streams_are_mostly_one_stream(self):
        r_z = [float(row["r_z"]) for row in self.rows]
        r_dm = [float(row["r_dm"]) for row in self.rows]
        self.assertAlmostEqual(statistics.correlation(r_z, r_dm), CORRELATION, places=3)

    def test_the_permutation_null_rose_with_the_sample(self):
        # 24 gaps gave a null mean of 0.39 and p = 0.33; 61 gaps give 1.17 and
        # p = 0.70.  More data made the one-rung step less remarkable, not more.
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
        observed = sum(1 for gap in self.a_gaps if abs(gap - 1.0) <= 0.05)
        mean = statistics.mean(counts)
        p_value = sum(1 for count in counts if count >= observed) / draws
        self.assertAlmostEqual(mean, NULL_MEAN, delta=0.20)
        self.assertGreater(p_value, 0.55)
        self.assertLess(p_value, 0.85)

    def test_the_registration_survived_the_extension_unchanged(self):
        path = os.path.join(self.root, PREREG)
        with open(path, encoding="utf-8") as handle:
            record = json.load(handle)
        pair = {row["name"]: row for row in self.rows}
        anchor = record["anchor"]
        self.assertAlmostEqual(
            float(pair[REGISTERED_FROM]["A"]), anchor["from"]["A"], places=6
        )
        self.assertAlmostEqual(
            float(pair[REGISTERED_TO]["A"]), anchor["to"]["A"], places=6
        )
        target = Decimal(str(record["prediction"]["target_A"]))
        tol = Decimal(str(record["prediction"]["tolerance_rungs"]))
        self.assertAlmostEqual(target, pair[REGISTERED_TO]["A"] + 1, places=6)
        # no burst in the extended sample lands in either registered window
        near_target = [row for row in self.rows if abs(row["A"] - target) <= tol]
        control = Decimal(str(record["control_window"]["target_A"]))
        near_control = [row for row in self.rows if abs(row["A"] - control) <= tol]
        self.assertEqual(near_target, [])
        self.assertEqual(near_control, [])
        self.assertEqual(record["status"], "open, no event scored")
        self.assertIn("sample_extension", record)

    def test_the_coding_discards_the_cluster_hosts(self):
        kept = [row for row in self.rows if not row["cluster"]]
        gaps = [
            float(kept[i + 1]["r_z"] - kept[i]["r_z"])
            for i in range(len(kept) - 1)
        ]
        self.assertEqual(len(kept), 60)
        self.assertEqual(sum(1 for gap in gaps if near_integer(gap, 0.01, 1)), 0)
        self.assertEqual(sum(1 for gap in gaps if near_integer(gap, 0.01, 0)), 5)


if __name__ == "__main__":
    unittest.main()
