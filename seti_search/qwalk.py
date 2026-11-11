"""A discrete-time quantum walk over the rung ladder, and the steps it takes.

Framework: phi^2 = phi + 1. Attribution: Tyree Jones (tyofthestarz), phi^2 =
phi + 1 framework, prior art from December 2025.

The ladder is cyclic: two quantities are the same address when their rungs
differ by an integer, so the natural graph is a ring of ``n_sites`` cells
spanning one rung period ``[0, 1)``. A walk on that ring moves in *address*
space, not in frequency space. That distinction is the whole point of this
module: the walk never reads the data amplitudes, so it cannot make a
candidate stronger. What it can do is state a spread law, and show that the
ladder reading of a set of measurements does or does not depend on how finely
the ring is cut.

The steps, in order
-------------------
1. **Lattice.** Fix one rung period and ``n_sites`` cells of width
   ``1/n_sites``; positions wrap, because rung + 1 is the same address.
2. **Coin.** A 2x2 unitary ``C``. ``golden_coin()`` is a rotation by the golden
   angle, ``2*pi/phi**2 = 2.399963... rad = 137.507764... deg`` - the same
   constant that makes a turn exactly ``phi**2`` golden angles.
   ``hadamard_coin()`` is the standard walk coin, kept for comparison so the
   golden coin is never mistaken for something special merely by being used.
3. **Initial superposition.** ``|+>_coin`` tensor a uniform superposition over
   the measured rung cells (or over every cell, if none are given). This is
   the "superposition" step: the seven candidate addresses are entered at once,
   with equal amplitudes, instead of one at a time.
4. **Step.** Apply the coin, then the conditional shift: coin up moves +1 cell,
   coin down moves -1 cell. ``psi <- S (C tensor I) psi``.
5. **Iterate** ``t`` steps. ``t`` is the only tunable, and it is set by a
   target: the number of steps whose spread reaches the span being compared.
6. **Measure.** ``p(x) = |psi_up(x)|^2 + |psi_down(x)|^2``; report spread in
   rung units, ``sigma_rung``, so the answer is comparable across lattices.
7. **Control.** Re-run step 6 with the data cells replaced by random cells and
   compare the peak probability at the data's own addresses against that null.
   Without step 7 a walk started on seven addresses will always "find" them.

Why the walk cannot be evidence
-------------------------------
The evolution is unitary and data-independent: nothing in steps 1-6 depends on
the measured powers, only on their rung addresses. Two different skies with the
same set of addresses produce the identical walk. A peak in ``p(x)`` therefore
says what a unitary reshuffle of the address set does, never what the sky does,
and the step-7 control shows the data's addresses are ordinary among random
ones: measured p = 0.48 at t = 20 and 0.26 at t = 40 on the seven candidates,
with the data peak no higher than the 95th percentile of a random placement.

The one usable output is the refinement test, and it fails honestly. Spread
quoted in rungs should be lattice-independent; measured from a four-address
start at t = 40 it is 0.073561, 0.066149 and 0.063840 rungs at 501, 1001 and
2001 cells, i.e. ``sigma_N = sigma_inf + a/N**2``, so a coarse lattice inflates
the reading by 15 per cent. The number tracks the ruler rather than the sky,
which is what an artifact looks like by this framework's own definition. Quote
sigma with its lattice, or do not quote it.

The refinement test only means something once the spread is many cells wide. A
single-address start at 40 steps spans about 19 cells, so its sigma in rungs
gives 0.038337, 0.019188 and 0.009599 at 501, 1001 and 2001 cells - halving with
the cell width, exactly as units require and with no information in it.

One measured surprise, and the reason not to call the coin a superposition
generator: ``R(137.507764 deg)|+>`` is 99.8 per cent one coin state, because the
golden angle is only 2.499 degrees from 135 degrees, where a rotation maps
``|+>`` exactly onto a single state. The standard ``|+>`` start therefore does
not fan out; the coin is reapplied each step, the amplitude re-rotates, and what
remains is a biased ballistic walk - measured drift 0.295 cells per step with
sigma 0.48 cells per step. A symmetric spreading superposition needs ``|+i>``
or a coin far from 135 degrees. That is a property of 137.507764 degrees, not of
the sky, and it is the opposite of what a walk is usually invoked to show.

Quantum speed-up, stated honestly: searching ``n`` drift trials by amplitude
amplification takes about ``(pi/4)*sqrt(n)`` oracle calls
(:func:`grover_steps`), which is fewer *queries* than the ``n`` of a classical
trial loop. Each oracle still needs the whole band coherently prepared - about
``2.9e8`` amplitudes for a Breakthrough Listen fine channelisation - and the
readout returns one drift. The query count improves; nothing else does, and no
experiment here is close to being able to run it. The walk is reported for what
it is: a way to move through the ladder, not a way to search it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .rungs import GOLDEN_ANGLE_DEG, PHI

__all__ = [
    "RungWalk",
    "classical_spread",
    "coin_rotation",
    "cover_steps",
    "golden_angle_rad",
    "golden_coin",
    "grover_steps",
    "hadamard_coin",
    "placement_report",
]


def golden_angle_rad() -> float:
    """The golden angle in radians: 2*pi/phi^2 = 2.3999632297..."""
    return math.radians(GOLDEN_ANGLE_DEG)


def coin_rotation(theta_rad: float) -> np.ndarray:
    """A 2x2 rotation - the coin used for the golden-angle walk."""
    cos, sin = math.cos(theta_rad), math.sin(theta_rad)
    return np.array([[cos, -sin], [sin, cos]], dtype=complex)


def golden_coin() -> np.ndarray:
    """Coin rotation by the golden angle (2.399963... rad)."""
    return coin_rotation(golden_angle_rad())


def hadamard_coin() -> np.ndarray:
    """The standard Hadamard coin, for comparison."""
    return np.array([[1.0, 1.0], [1.0, -1.0]], dtype=complex) / math.sqrt(2.0)


@dataclass
class RungWalk:
    """Discrete-time quantum walk on a cyclic lattice of rung addresses."""

    n_sites: int = 1001
    period: float = 1.0
    coin: np.ndarray | None = None
    coin_name: str = field(default="golden", compare=False)

    def __post_init__(self) -> None:
        if self.n_sites < 3:
            raise ValueError(f"need at least 3 sites, got {self.n_sites}")
        if not self.period > 0.0:
            raise ValueError(f"period must be positive, got {self.period!r}")
        coin = golden_coin() if self.coin is None else np.asarray(self.coin, dtype=complex)
        if coin.shape != (2, 2):
            raise ValueError(f"coin must be 2x2, got {coin.shape}")
        if not np.allclose(coin @ coin.conj().T, np.eye(2), atol=1e-12):
            raise ValueError("coin is not unitary")
        self.coin = coin

    # -- geometry ---------------------------------------------------------
    @property
    def cell_width(self) -> float:
        """Rung width of one lattice cell."""
        return self.period / self.n_sites

    def cells_for(self, rungs) -> np.ndarray:
        """Map rung values onto cell indices (wrapping at each period)."""
        width = self.cell_width
        out = [int(round((float(r) % self.period) / width)) % self.n_sites for r in rungs]
        return np.array(out, dtype=int)

    def rung_of(self, cells) -> np.ndarray:
        return np.asarray(cells, dtype=float) * self.cell_width

    # -- state ------------------------------------------------------------
    def start(self, cells=None) -> np.ndarray:
        """|+> coin tensor a uniform superposition over ``cells`` (or all)."""
        amps = np.zeros(self.n_sites, dtype=complex)
        if cells is None:
            amps[:] = 1.0 / math.sqrt(self.n_sites)
        else:
            unique = np.unique(np.asarray(cells, dtype=int))
            if unique.size == 0:
                raise ValueError("no cells given")
            amps[unique] = 1.0 / math.sqrt(unique.size)
        plus = np.full(2, 1.0 / math.sqrt(2.0), dtype=complex)
        return np.outer(plus, amps)

    def step(self, psi: np.ndarray) -> np.ndarray:
        """Coin, then conditional shift by one cell."""
        turned = self.coin @ psi
        out = np.empty_like(turned)
        out[0] = np.roll(turned[0], 1)
        out[1] = np.roll(turned[1], -1)
        return out

    def run(self, steps: int, psi: np.ndarray | None = None) -> np.ndarray:
        state = self.start() if psi is None else np.asarray(psi, dtype=complex)
        for _ in range(int(steps)):
            state = self.step(state)
        return state

    # -- measurement ------------------------------------------------------
    def probability(self, psi: np.ndarray) -> np.ndarray:
        prob = np.abs(psi[0]) ** 2 + np.abs(psi[1]) ** 2
        total = prob.sum()
        return prob / total if total > 0.0 else prob

    def spread(self, psi: np.ndarray) -> float:
        """RMS displacement in rung units (wrapping ignored, so keep t small)."""
        prob = self.probability(psi)
        x = np.arange(self.n_sites) * self.cell_width
        mean = float((prob * x).sum())
        return float(math.sqrt(float((prob * (x - mean) ** 2).sum())))


def classical_spread(steps: int, cell_width: float) -> float:
    """RMS displacement of an unbiased classical walk, in the same units."""
    return float(cell_width) * math.sqrt(max(int(steps), 0))


def cover_steps(
    walk: RungWalk,
    target_rung: float,
    cells=None,
    max_steps: int = 2000,
) -> int | None:
    """First t whose spread reaches ``target_rung``; None if it never does.

    With no cells given the walk starts at one address (the centre cell), which
    is the question this answers: how many steps does a single address need to
    smear across a span that wide. The state is stepped in place, so the cost is
    linear in ``max_steps`` rather than quadratic.
    """
    if cells is None:
        cells = [walk.n_sites // 2]
    state = walk.start(cells)
    for steps in range(1, int(max_steps) + 1):
        state = walk.step(state)
        if walk.spread(state) >= float(target_rung):
            return steps
    return None


def placement_report(
    walk: RungWalk,
    cells,
    steps: int,
    trials: int = 2000,
    seed: int = 0,
) -> dict:
    """Is the walk concentrated on the data's addresses, or on any addresses?

    The control is the same walk with the data cells swapped for random cells:
    a walk started on seven sites will always return probability at those seven
    sites, so the peak only means something relative to that null.
    """
    cells = np.unique(np.asarray(cells, dtype=int))
    prob = walk.probability(walk.run(steps, walk.start(cells)))
    observed = float(prob[cells].max())
    uniform = 1.0 / walk.n_sites
    rng = np.random.default_rng(seed)
    null = np.array([
        float(prob[rng.choice(walk.n_sites, size=cells.size, replace=False)].max())
        for _ in range(int(trials))
    ])
    return {
        "cells": int(cells.size),
        "steps": int(steps),
        "peak": observed,
        "uniform": uniform,
        "ratio_to_uniform": observed / uniform,
        "null_median": float(np.median(null)),
        "null_p95": float(np.percentile(null, 95)),
        "p_value": float((null >= observed).mean()),
        "trials": int(trials),
    }


def grover_steps(drift_trials: int) -> int:
    """Oracle calls for amplitude amplification over ``drift_trials`` rates."""
    return int(math.ceil((math.pi / 4.0) * math.sqrt(max(int(drift_trials), 1))))
