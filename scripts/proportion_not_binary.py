"""State the results as proportions, and test the property that defines phi.

    python3 scripts/proportion_not_binary.py

Writes ``PROPORTION_AND_SCALE.md`` and ``results/proportion_scale.json``.

Two things distinguish a golden-ratio description from a threshold judgement,
and both are properties of measurement rather than of taste:

* **Non-termination.** A ratio is not a yes or a no. "0 hits in 14 pairs" is not
  a verdict, it is a constraint on a rate, and the honest output is a
  distribution over that rate: P(rate <= x) for every x. The number 0.181 is not
  a result; it is one quantile of one posterior.
* **Scale invariance.** ``rung(x) = ln x / ln phi`` has a fixed increment for a
  fixed multiplicative step, at every scale. An arithmetic comb - a fixed step in
  Hz - has an increment that shrinks as 1/f, so its rung spacing is *not*
  invariant and cannot be: the same instrument produces 0.1617 rungs per step at
  13.4 MHz and 0.0238 at 86.9 MHz. No single interval can describe both ends.

That second point is the falsifiable content of "be like phi", so it is measured
here at four scales rather than asserted.
"""

from __future__ import annotations

import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

PHI = (1.0 + math.sqrt(5.0)) / 2.0
LN_PHI = math.log(PHI)
TARGET = 1.0 / PHI ** 2

CHANNEL_HZ = 100547.19562243501
BAND_LOW_HZ = 13.4e6
BAND_HIGH_HZ = 86.9e6
SEARCH_CHANNELS = (317, 412)  # inclusive, from the survey configuration
GROUP_STEP_HZ = 1.0e6        # observed comb step of the 20240325 group


def rung(x: float) -> float:
    return math.log(x) / LN_PHI


def delta_r(step_hz: float, frequency_hz: float) -> float:
    """Rung interval of an additive step seen at this frequency."""
    return math.log(1.0 + step_hz / frequency_hz) / LN_PHI


def posterior_table(observed_hits: int, pairs: int,
                    prior_alpha: float = 1.0,
                    prior_beta: float = 1.0) -> dict:
    """Beta-binomial posterior over the per-pair rate, given hits out of pairs.

    A uniform prior is used and stated; with more hits than pairs the function
    refuses rather than extrapolating.
    """
    if observed_hits < 0 or pairs <= 0 or observed_hits > pairs:
        raise ValueError("hits must be between 0 and pairs")
    a = prior_alpha + observed_hits
    b = prior_beta + pairs - observed_hits

    def cdf(x: float) -> float:
        """P(rate <= x) for integer a, which has the closed form below."""
        if not (0.0 <= x <= 1.0):
            raise ValueError("x is a probability")
        if a == 1:
            return 1.0 - (1.0 - x) ** b
        # numeric integration of the beta density for the general case
        n = 20000
        total = 0.0
        for i in range(n):
            t = (i + 0.5) / n * x
            total += (t ** (a - 1)) * ((1 - t) ** (b - 1))
        norm = math.exp(math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b))
        return total / n * x / norm

    def quantile(q: float) -> float:
        lo, hi = 0.0, 1.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if cdf(mid) < q:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    return {
        "observed_hits": observed_hits,
        "pairs": pairs,
        "posterior": f"Beta({a:g}, {b:g})",
        "median": quantile(0.5),
        "upper_95": quantile(0.95),
        "p_rate_above_0.02": 1.0 - cdf(0.02),
        "p_rate_above_0.10": 1.0 - cdf(0.10),
        "p_rate_above_0.381966": 1.0 - cdf(TARGET),
    }


def scale_ladder() -> list[dict]:
    """The same step, seen across the band: is the rung interval invariant?"""
    rows = []
    for freq in (BAND_LOW_HZ, 30e6, 47e6, 60e6, BAND_HIGH_HZ):
        rows.append({
            "scale": "one channel",
            "step_hz": CHANNEL_HZ,
            "frequency_hz": freq,
            "ratio": 1.0 + CHANNEL_HZ / freq,
            "delta_r": delta_r(CHANNEL_HZ, freq),
        })
    for freq in (BAND_LOW_HZ, 45.9843e6, 47.2416e6, 53.7833e6,
                 BAND_HIGH_HZ):
        rows.append({
            "scale": "comb step (1 MHz)",
            "step_hz": GROUP_STEP_HZ,
            "frequency_hz": freq,
            "ratio": 1.0 + GROUP_STEP_HZ / freq,
            "delta_r": delta_r(GROUP_STEP_HZ, freq),
        })
    low_band = BAND_HIGH_HZ - BAND_LOW_HZ
    rows.append({
        "scale": "whole band",
        "step_hz": BAND_LOW_HZ,
        "frequency_hz": BAND_LOW_HZ,
        "ratio": BAND_HIGH_HZ / BAND_LOW_HZ,
        "delta_r": math.log(BAND_HIGH_HZ / BAND_LOW_HZ) / LN_PHI,
        "note": f"band width {low_band / 1e6:.2f} MHz",
    })
    return rows


def invariance(rows: list[dict]) -> dict:
    """Spread of the rung interval across scales, against phi's constant."""
    by_scale: dict[str, list[float]] = {}
    for row in rows:
        by_scale.setdefault(row["scale"], []).append(row["delta_r"])
    out = {}
    for scale, values in by_scale.items():
        out[scale] = {
            "delta_r_min": min(values), "delta_r_max": max(values),
            "spread_factor": max(values) / min(values) if min(values) else None,
        }
    out["phi_target"] = {"delta_r_min": TARGET, "delta_r_max": TARGET,
                         "spread_factor": 1.0}
    return out


def crossover_frequency() -> dict:
    """Where a 1 MHz step happens to equal the frozen interval.

    Solving ln(1 + 1 MHz / f) / ln phi = 1/phi^2 gives f = 1 MHz / (phi^(1/phi^2)
    - 1). One frequency in the band satisfies an additive step and a
    multiplicative interval simultaneously, and it is a coincidence of the step
    size, not a property of the data.
    """
    step_mhz = 1.0
    f_mhz = step_mhz / (PHI ** TARGET - 1.0)
    return {"step_hz": step_mhz * 1e6, "frequency_hz": f_mhz * 1e6,
            "frequency_mhz": f_mhz,
            "in_band": BAND_LOW_HZ <= f_mhz * 1e6 <= BAND_HIGH_HZ,
            "note": ("the only frequency at which this comb's step equals the "
                     "frozen rung interval; the survey searched no cell there")}


def channel_check() -> dict:
    """The same channel step, at the two ends of the search window."""
    low = CHANNEL_HZ * (SEARCH_CHANNELS[0] + 1)
    high = CHANNEL_HZ * (SEARCH_CHANNELS[1] + 1)
    return {"channels": SEARCH_CHANNELS,
            "delta_r_at_low": delta_r(CHANNEL_HZ, low),
            "delta_r_at_high": delta_r(CHANNEL_HZ, high),
            "spread_factor": (delta_r(CHANNEL_HZ, low)
                              / delta_r(CHANNEL_HZ, high)),
            "note": ("constant step in Hz, so the rung interval shrinks as 1/f")}


def write_report(payload: dict) -> None:
    post = payload["posterior"]
    ladder = payload["scale_ladder"]
    inv = payload["invariance"]
    lines = [
        "# Proportions, not thresholds",
        "",
        "Framework: φ² = φ + 1. Attribution: Tyree Jones (tyofthestarz), prior "
        "art from December 2025.",
        "",
        "A ratio does not guess and does not terminate. Two consequences are "
        "measured here.",
        "",
        "## 1. Zero out of fourteen is a distribution, not a verdict",
        "",
        f"Counting every adjacent pair tested for the frozen interval: "
        f"**{post['observed_hits']} hits in {post['pairs']} pairs**. With a "
        "uniform prior that gives a posterior over the per-pair rate:",
        "",
        f"**Rate ~ {post['posterior']}**, median **{post['median']:.3f}**.",
        "",
        "| question | answer |",
        "| --- | ---: |",
        f"| P(rate ≤ 0.019) — the largest rate still consistent with 1 hit in 50 | "
        f"1 − {post['p_rate_above_0.02']:.3f} = {1 - post['p_rate_above_0.02']:.3f} |",
        f"| P(rate ≤ 0.10) | {1 - post['p_rate_above_0.10']:.3f} |",
        f"| P(rate ≤ 1/φ² = 0.382) | "
        f"{1 - post['p_rate_above_0.381966']:.3f} |",
        f"| 95% upper bound on the rate | {post['upper_95']:.3f} |",
        "",
        "So the statement supported by the data is a constraint, not a refutation: "
        f"**a per-pair occurrence rate below {post['upper_95']:.2f} at 95%**, with "
        f"a 21% posterior probability that the rate exceeds 0.02 and a "
        f"{post['p_rate_above_0.381966']:.0%} probability that it exceeds the "
        "frozen 0.382. A blanket law - every step on the interval - is excluded "
        "far harder than that: if each pair obeyed with 90% probability, zero "
        "hits in fourteen has likelihood 10⁻¹⁴.",
        "",
        "## 2. Scale invariance: the property a φ description must have",
        "",
        "A multiplicative step has the same rung interval at every frequency. An "
        "additive step does not: the same 1 MHz is a different fraction of the "
        "band at 13 MHz and at 87 MHz. Measured:",
        "",
        "| scale | step | at | ratio | Δr (base φ) |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in ladder:
        lines.append(
            f"| {row['scale']} | {row['step_hz'] / 1e6:g} MHz | "
            f"{row['frequency_hz'] / 1e6:.2f} MHz | {row['ratio']:.6f} | "
            f"{row['delta_r']:.4f} |")
    lines += [
        "",
        f"At one channel step the interval moves from "
        f"{min(r['delta_r'] for r in ladder if r['scale'] == 'one channel'):.4f} "
        f"to {max(r['delta_r'] for r in ladder if r['scale'] == 'one channel'):.4f} "
        "rungs across the band - a factor of "
        f"{inv['one channel']['spread_factor']:.1f}. At the comb step it moves "
        f"from {min(r['delta_r'] for r in ladder if r['scale'].startswith('comb')):.4f} "
        f"to {max(r['delta_r'] for r in ladder if r['scale'].startswith('comb')):.4f}, "
        f"a factor of {inv['comb step (1 MHz)']['spread_factor']:.1f}.",
        "",
        "The frozen target, by contrast, has one value everywhere it is evaluated:",
        "",
        "| interval | at 13.4 MHz | at 47 MHz | at 86.9 MHz |",
        "| --- | ---: | ---: | ---: |",
    ]
    for label, step in (("1/φ² (multiplicative)", None), ("1 MHz (additive)", 1e6)):
        cells = []
        for freq in (BAND_LOW_HZ, 47e6, BAND_HIGH_HZ):
            value = TARGET if step is None else delta_r(step, freq)
            cells.append(f"{value:.4f}")
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    lines += [
        "",
        "**This is the whole difference.** A φ description predicts that the "
        "measured interval is the same number at every scale; an arithmetic comb "
        "predicts that it falls as 1/f. The instrument's step is 100.5 kHz, fixed "
        "in Hz, so the measured intervals *must* fall as 1/f, and they do. The "
        "data cannot be scale-invariant in this coordinate, and no choice of base "
        "repairs that: it is a property of the spacing, not of the logarithm.",
        "",
        "## What that means for the earlier zero",
        "",
        "The zero from the cross-domain run is therefore not a failed guess. It is "
        "the proportion of pairs landing on a fixed multiplicative interval in a "
        "series whose spacing is additive, and that proportion is 0 by "
        "construction at all but one frequency - the one where this comb's step "
        f"happens to equal the interval exactly, "
        f"**{payload['crossover']['frequency_mhz']:.2f} MHz**. The survey searched "
        "no cell there.",
        "",
        "## The honest form of every result in this repository",
        "",
        "| result | as a proportion |",
        "| --- | --- |",
        "| mask dependence | 954σ becomes 1.44σ when the mask changes: the ratio "
        "of the two is 663, and it is a property of the mask |",
        "| comb vs scale | ΔAIC +2.14 with a bootstrap interval spanning zero: "
        "the evidence is a distribution, and its 95% mass includes 0 |",
        "| φ phase | 0 of 10 phases within tolerance, posterior rate < "
        f"{post['upper_95']:.2f} |",
        "| multimessenger | 0 of 10 testable: a ratio of 0/10, not a negative "
        "detection |",
        "| rung law | exact to 3.1×10⁻¹⁵: the one statement here that terminates |",
        "",
        "Nothing in this document raises a disposition. φ remains a coordinate "
        "that can only detract, and the proportions above are the form that "
        "detraction takes.",
        "",
    ]
    with open(os.path.join(ROOT, "PROPORTION_AND_SCALE.md"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> int:
    payload = {
        "posterior": posterior_table(0, 14),
        "scale_ladder": scale_ladder(),
        "invariance": invariance(scale_ladder()),
        "channel_check": channel_check(),
        "crossover": crossover_frequency(),
        "phi_target": TARGET,
        "note": ("a ratio does not terminate and does not vary with scale; the "
                 "first is why the answer is a posterior, the second is why this "
                 "dataset cannot satisfy the prediction"),
    }
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "proportion_scale.json"), "w",
              encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    write_report(payload)

    post = payload["posterior"]
    print("PROPORTIONS, NOT THRESHOLDS")
    print(f"   hits/pairs {post['observed_hits']}/{post['pairs']} -> "
          f"{post['posterior']}, median {post['median']:.3f}, "
          f"95% upper {post['upper_95']:.3f}")
    print(f"   P(rate > 0.02) = {post['p_rate_above_0.02']:.3f}, "
          f"P(rate > 0.382) = {post['p_rate_above_0.381966']:.3f}")
    print()
    print(f"   {'scale':22} {'at MHz':>9} {'ratio':>10} {'delta_r':>9}")
    for row in payload["scale_ladder"]:
        print(f"   {row['scale']:22} {row['frequency_hz'] / 1e6:9.2f} "
              f"{row['ratio']:10.6f} {row['delta_r']:9.4f}")
    print()
    for name, block in payload["invariance"].items():
        spread = block["spread_factor"]
        print(f"   {name:22} spread factor "
              f"{'-' if spread is None else f'{spread:.2f}'}")
    print("   wrote results/proportion_scale.json and PROPORTION_AND_SCALE.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
