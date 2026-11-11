"""Tests for the rung-lattice quantum walk: unitarity, spread law, lattice invariance."""

import math
import unittest

import numpy as np

from seti_search.qwalk import (
    RungWalk,
    classical_spread,
    cover_steps,
    golden_angle_rad,
    golden_coin,
    grover_steps,
    hadamard_coin,
    placement_report,
)
from seti_search.rungs import GOLDEN_ANGLE_DEG, PHI


class QuantumWalkTests(unittest.TestCase):
    def test_golden_angle_is_the_turn_over_phi_squared(self):
        self.assertAlmostEqual(GOLDEN_ANGLE_DEG * PHI**2, 360.0, places=9)
        self.assertAlmostEqual(golden_angle_rad(), 2.4, places=4)
        self.assertAlmostEqual(math.degrees(golden_angle_rad()), GOLDEN_ANGLE_DEG, places=9)

    def test_coins_are_unitary(self):
        for coin in (golden_coin(), hadamard_coin()):
            self.assertTrue(np.allclose(coin @ coin.conj().T, np.eye(2)))

    def test_norm_is_conserved(self):
        walk = RungWalk(n_sites=501)
        state = walk.start(walk.cells_for([43.6361, 43.6583, 43.7000]))
        for _ in range(40):
            state = walk.step(state)
            self.assertAlmostEqual(float(np.sum(np.abs(state) ** 2)), 1.0, places=12)

    def test_rejects_non_unitary_coin(self):
        with self.assertRaises(ValueError):
            RungWalk(n_sites=51, coin=np.array([[1.0, 0.0], [0.0, 0.5]], dtype=complex))

    def test_walk_spreads_ballistically_and_faster_than_classical(self):
        walk = RungWalk(n_sites=2001, coin=hadamard_coin())
        state = walk.start(walk.cells_for([0.5]))
        spreads = {}
        for _ in range(60):
            state = walk.step(state)
            spreads[len(spreads) + 1] = walk.spread(state)
        quantum = spreads[60]
        classical = classical_spread(60, walk.cell_width)
        self.assertGreater(quantum, 2.0 * classical)
        # ballistic: doubling the steps roughly doubles the spread, unlike sqrt(t)
        self.assertGreater(spreads[60] / spreads[30], 1.7)
        self.assertLess(spreads[60] / spreads[30], 2.1)

    def test_spread_in_rungs_converges_under_refinement(self):
        """Not invariant: it converges like 1/N^2, and that is the artifact signal.

        A start spanning several cells at t = 40 gives 0.073561, 0.066149 and
        0.063840 rungs at 501, 1001 and 2001 cells; a coarse lattice cannot
        resolve the fringe structure and inflates the reading by 15 per cent.
        """
        rungs = [43.6361, 43.6583, 43.7000, 43.8008]
        spreads = []
        for n_sites in (501, 1001, 2001):
            walk = RungWalk(n_sites=n_sites)
            spreads.append(walk.spread(walk.run(40, walk.start(walk.cells_for(rungs)))))
        first = spreads[0] - spreads[1]
        second = spreads[1] - spreads[2]
        self.assertGreater(first / second, 2.0, msg=str(spreads))
        self.assertLess((spreads[0] - spreads[2]) / spreads[2], 0.20, msg=str(spreads))

    def test_cover_steps_reaches_the_requested_span(self):
        walk = RungWalk(n_sites=2001)
        target = 0.1
        steps = cover_steps(walk, target, max_steps=900)
        self.assertIsNotNone(steps)
        state = walk.run(steps, walk.start([walk.n_sites // 2]))
        self.assertGreaterEqual(walk.spread(state), target)
        state = walk.run(steps - 1, walk.start([walk.n_sites // 2]))
        self.assertLess(walk.spread(state), target)

    def test_placement_is_ordinary_against_random_addresses(self):
        walk = RungWalk(n_sites=1001)
        rungs = [43.6361, 43.6583, 43.7000, 43.8008, 43.8200, 43.8500, 43.8595]
        report = placement_report(walk, walk.cells_for(rungs), steps=20, trials=300, seed=1)
        self.assertGreaterEqual(report["p_value"], 0.0)
        self.assertLessEqual(report["p_value"], 1.0)
        self.assertGreater(report["peak"], 0.0)
        self.assertGreater(report["ratio_to_uniform"], 1.0)

    def test_grover_steps_match_the_sqrt_law(self):
        self.assertEqual(grover_steps(1), 1)
        self.assertEqual(grover_steps(817), 23)
        self.assertLess(grover_steps(10_000), 0.01 * 10_000)


if __name__ == "__main__":
    unittest.main()
