"""The merge must be able to fail: equal edges, and an exact rung for 1/phi^2."""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from super_tetrahedron import (  # noqa: E402
    GOLDEN_ANGLE_DEG,
    TETRAHEDRAL_DEG,
    angle_facts,
    line_simplex_impossible,
    simplex,
)

from seti_search.rungs import rung  # noqa: E402


class SimplexTests(unittest.TestCase):
    def test_four_equally_spaced_points_are_not_a_simplex(self):
        # equal adjacent gaps do not make equal pairwise distances: the
        # separations are 1, 2, 3, 1, 2, 1
        block = simplex([0.0, 1.0, 2.0, 3.0])
        self.assertAlmostEqual(block["max_min"], 3.0)
        self.assertGreater(block["equal_edge_deviation"], 0.5)

    def test_four_points_on_a_line_can_never_be_equidistant(self):
        block = line_simplex_impossible([-4.7681, -2.0, 24.3258, 36.7230])
        self.assertTrue(block["impossible"])
        self.assertTrue(block["outer_exceeds_largest_gap"])
        self.assertAlmostEqual(block["outer_separation"],
                               block["sum_of_gaps"], places=9)

    def test_the_ladder_is_not_a_simplex(self):
        block = simplex([-4.7681, -2.0, 24.3258, 36.7230])
        self.assertGreater(block["equal_edge_deviation"], 0.8)
        self.assertGreater(block["max_min"], 10.0)

    def test_six_bonds_from_four_vertices(self):
        self.assertEqual(len(simplex([0.0, 1.0, 2.0, 3.0])["separations"]), 6)


class AngleTests(unittest.TestCase):
    def test_the_tetrahedral_angle_is_arccos_minus_one_third(self):
        self.assertAlmostEqual(TETRAHEDRAL_DEG, math.degrees(math.acos(-1 / 3)))
        self.assertAlmostEqual(TETRAHEDRAL_DEG, 109.4712, places=3)

    def test_the_constants_are_not_measurements(self):
        facts = angle_facts()
        self.assertIn("none of them is a measurement", facts["note"])
        self.assertAlmostEqual(facts["difference_deg"],
                               GOLDEN_ANGLE_DEG - TETRAHEDRAL_DEG)


class TargetRungTests(unittest.TestCase):
    def test_one_over_phi_squared_sits_exactly_two_rungs_below_unity(self):
        self.assertAlmostEqual(rung(1.0 / ((1 + math.sqrt(5)) / 2) ** 2), -2.0,
                               places=12)


if __name__ == "__main__":
    unittest.main()
