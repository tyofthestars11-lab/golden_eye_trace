"""phi^2 = phi + 1 read as a recurrence, not an equation to be judged.

Framework: phi^2 = phi + 1. Attribution: Tyree Jones (tyofthestarz), phi^2 =
phi + 1 framework, prior art from December 2025.

The identity is usually quoted as a fact about a number. It is more useful read
as a *process*: from any positive x, add one and invert,

    T(x) = 1 + 1/x ,

and repeat. The sequence never settles on a rational because the "+1" is never
carried - it is an instruction, not a quantity. ``T(phi) = phi``, so phi is the
fixed point of the map, and everything the identity says about measurement
follows from the derivative at that fixed point:

    T'(x) = -1 / x^2        =>      T'(phi) = -1/phi^2 = -0.3819660...

That single line is why ``1/phi^2`` keeps appearing in this project. It is not a
tolerance, a band, or a threshold: it is **the ratio by which the recurrence
shrinks its own error each step**. An orbit starting at any x0 converges to phi
geometrically, losing a factor ``1/phi^2`` per iteration - 2.41 steps per
decimal place. A measurement whose structure is phi-like should therefore show
that *decay*, and a measurement whose structure is arithmetic cannot, because
arithmetic spacing changes the map's input at every scale.

What replaces a threshold
-------------------------
A test of the form ``|x - 1/phi^2| < 0.01`` rounds a rate into a window, and a
window in rungs is not scale free: 0.01 rungs is 0.14 MHz of physical width at
13.4 MHz and 0.90 MHz at 86.9 MHz, a factor of 6.49. The same test therefore
means six different things across one band. That is the "whole-ing" failure -
adding the +1 back as a *quantity* and carrying it, when it should have stayed
in the map.

Two non-binary readings are implemented here:

* :func:`orbit` - iterate the map and report the residual at every step with its
  per-step contraction ratio. No pass/fail; the trajectory is the output.
* :func:`continued_fraction` / :func:`agreement_depth` / :func:`excess` - write
  a ratio as ``a0 + 1/(a1 + 1/(a2 + ...))`` and keep *every* partial quotient.
  Phi is ``[1; 1, 1, 1, ...]``, so its excess ``sum(a_i - 1)`` is zero at every
  depth and its agreement depth is unbounded. Every other number pays: the
  channel grid ``(k+1)/k = [1; k]`` agrees for exactly one level and then hands
  over a partial quotient of 318. Depth is a number, not a verdict.

Floating point limit, declared: a double holds ~15 exact levels of a continued
fraction near phi, after which the digits are IEEE 754 noise and not evidence.
:func:`continued_fraction` stops at ``|r - phi| < 1e-13`` and says so.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction

__all__ = [
    "PHI",
    "effective_depth",
    "is_phi_shaped",
    "LN_PHI",
    "CONTRACTION",
    "STEPS_PER_DECADE",
    "Orbit",
    "agreement_depth",
    "continued_fraction",
    "convergents",
    "derivative",
    "excess",
    "fibonacci",
    "orbit",
    "step",
]

PHI = (1.0 + math.sqrt(5.0)) / 2.0
LN_PHI = math.log(PHI)

#: T'(phi) = -1/phi^2, the per-step error ratio of the recurrence.
CONTRACTION = 1.0 / PHI ** 2

#: iterations needed to lose one decimal place: 1/|log10(CONTRACTION)|
STEPS_PER_DECADE = 1.0 / abs(math.log10(CONTRACTION))

_FLOAT_DEPTH_LIMIT = 1e-13


def step(x: float) -> float:
    """One turn of the recurrence: add one, invert. Refuses x <= 0."""
    if x <= 0.0:
        raise ValueError("T(x) = 1 + 1/x needs x > 0")
    return 1.0 + 1.0 / x


def derivative(x: float) -> float:
    """T'(x) = -1/x^2. At phi this is -CONTRACTION."""
    if x <= 0.0:
        raise ValueError("T'(x) = -1/x^2 needs x > 0")
    return -1.0 / (x * x)


def fibonacci(n: int) -> int:
    """F(0) = 0, F(1) = 1. Used for the exact convergent ladder."""
    if n < 0:
        raise ValueError("n must be >= 0")
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


@dataclass(frozen=True)
class Orbit:
    """One trajectory of the recurrence, with no verdict attached."""

    start: float
    values: list = field(default_factory=list)
    residuals: list = field(default_factory=list)
    ratios: list = field(default_factory=list)

    def __post_init__(self):  # pragma: no cover - frozen dataclass guard
        if not self.values:
            raise ValueError("an orbit needs at least one value")

    @property
    def settled(self) -> bool:
        return self.residuals[-1] < 1e-12

    @property
    def contraction(self) -> float:
        """Median |residual ratio| over the last third of the orbit."""
        if not self.ratios:
            return float("nan")
        tail = self.ratios[-max(1, len(self.ratios) // 3):]
        return float(sorted(abs(r) for r in tail)[len(tail) // 2])

    def stages_to(self, decades: float) -> int:
        """How many turns to shrink the remaining error by ``decades``."""
        return int(math.ceil(decades * STEPS_PER_DECADE))


def orbit(start: float, turns: int = 40) -> Orbit:
    """Iterate T from ``start`` and record residual and contraction each turn."""
    if turns < 1:
        raise ValueError("turns must be >= 1")
    values = [float(start)]
    residuals = [abs(start - PHI)]
    ratios = []
    x = float(start)
    for _ in range(turns):
        x = step(x)
        values.append(x)
        residuals.append(abs(x - PHI))
        if residuals[-2] > 0.0:
            ratios.append((x - PHI) / (values[-2] - PHI))
    return Orbit(start=float(start), values=values,
                 residuals=residuals, ratios=ratios)


def continued_fraction(x: float, depth: int = 24) -> list:
    """Partial quotients of ``x``, every level kept, none rounded away."""
    if x <= 0.0:
        raise ValueError("continued_fraction needs x > 0")
    if depth < 1:
        raise ValueError("depth must be >= 1")
    out = []
    exact = isinstance(x, Fraction)
    r = Fraction(x) if exact else float(x)
    for _ in range(depth):
        a = math.floor(r)
        out.append(int(a))
        frac = r - a
        if frac <= 0.0:
            break
        r = (1 / frac) if exact else (1.0 / frac)
        if (not exact) and abs(r - PHI) < _FLOAT_DEPTH_LIMIT and len(out) >= 2:
            # every further level is 1 by the identity; the doubles are gone
            break
    return out


def excess(partials) -> int:
    """``sum(a_i - 1)`` over the fractional recursion, i.e. from a1 on.

    The integer part a0 is the whole that the identity says *not* to carry, so
    it is excluded: ``excess`` is zero exactly when the recursion is "add one"
    at every level, which is phi. Values in (0, 1) are handled naturally -
    their a0 is 0 and carries no information either way.

    Ambiguity, declared: a ratio that terminates immediately (an integer, e.g.
    ``2 = [2]``) has an empty fractional recursion, so it also reports depth 0
    and excess 0. Neither statistic can tell that apart from phi; callers must
    carry :func:`is_phi_shaped` beside them and read neither alone.
    """
    return int(sum(int(a) - 1 for a in list(partials)[1:]))


def agreement_depth(partials) -> int:
    """How many levels of the fractional recursion equal 1, from a1 on."""
    depth = 0
    for a in list(partials)[1:]:
        if int(a) != 1:
            break
        depth += 1
    return depth


def convergents(x: float, depth: int = 12) -> list:
    """Exact convergent ladder p_n/q_n of ``x``, with each level's error."""
    if x <= 0.0:
        raise ValueError("x must be > 0")
    if abs(x - PHI) < 1e-12:
        # phi's own ladder is the Fibonacci ladder; a float cannot hold it
        out = []
        prev_err = None
        for n in range(1, depth + 1):
            p, q = fibonacci(n + 1), fibonacci(n)
            err = abs(p / q - x)
            out.append((p, q, err if prev_err is None else err))
            prev_err = err
        return out
    partials = continued_fraction(x, depth=depth)
    p_prev, p = 1, partials[0]
    q_prev, q = 0, 1
    out = []
    for a in partials[1:]:
        p_prev, p = p, a * p + p_prev
        q_prev, q = q, a * q + q_prev
        out.append((p, q, abs(p / q - x)))
        if q > 10 ** 15:
            break
    return out


def is_phi_shaped(x: float, levels: int = 12, tol: float = 1e-9) -> bool:
    """True when the fractional recursion keeps returning to phi.

    A double cannot carry "add one forever", so an all-ones tail is detected
    structurally instead of being counted: apply the map to the fractional part
    ``levels`` times and check it never leaves phi's neighbourhood. phi, 1/phi
    and phi^2 all satisfy this, because subtracting 1 or adding 1 are exact.
    """
    if x <= 0.0:
        raise ValueError("x must be > 0")
    r = float(x) - math.floor(x)
    if r <= 0.0:
        return False
    for _ in range(levels):
        r = 1.0 / r
        r = r - math.floor(r)
        if r <= 0.0 or abs(r - (PHI - 1.0)) > tol:
            return False
    return True


def effective_depth(x: float) -> float:
    """Continuous depth: how many doublings of phi-separation would fit.

    ``|x/phi - 1| = phi ** (-2 * depth)``, so a value that agrees with phi to
    ``n`` levels reports about ``n``. Zero agreement is not a verdict, it is
    ``0.0``; a value far from phi reports a negative number and the caller must
    say so rather than clamping it to "pass" or "fail".
    """
    if x <= 0.0:
        raise ValueError("x must be > 0")
    miss = abs(float(x) / PHI - 1.0)
    if miss < 1e-16:
        return math.inf
    return -math.log(miss) / math.log(PHI ** 2)
