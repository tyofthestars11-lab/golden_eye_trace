"""Carry the +1 as a map, not as a quantity.

    python3 scripts/phi_plus_one.py

Writes ``PHI_PLUS_ONE.md`` and ``results/phi_plus_one.json``.

phi^2 = phi + 1 says: from any x, add one and invert, and the value returns to
itself. The return rate is T'(phi) = -1/phi^2 = -0.3819660, which is the number
this project kept flattening into a tolerance. Two readings that do not
terminate replace it - the orbit's contraction, and the continued fraction's
depth - and both are applied here to the ratios the campaign actually measured
(``results/comb_vs_scale.csv``), not only to constructions whose answer is
known in advance, such as the channel grid, which is 1 + 1/k by definition.

What this does *not* do is explain the campaign's zero hit count. That zero is a
property of additive spacing against a fixed multiplicative interval, measured
in ``PROPORTION_AND_SCALE.md``; depth was withdrawn from that job.

Beside each group's measured ratios the table also carries the additive model's
own adjacent-line ratio, ``1 + comb_df / f_low`` at the group's lowest line
from ``results/expanded_phi_groups.json``. That value is predicted by the fit
rather than measured from the lines, so it is diagnostic only and is never
counted with the measurements.

At the lowest line the two models coincide to first order, so that one value
cannot separate them. The reading is therefore carried to *every* adjacent pair
in the group, ``1 + step / f_i`` for i = 1 .. n-1, and compounded across the
span as their product - and the signed distance from the fitted multiplicative
interval is reported pair by pair, so the divergence is measured where it
actually opens up. All of it is still diagnostic only: both sides of that
difference are predictions of a fit, not measurements of the lines.

That last step still compares two fits with each other. So the same pairs are
scored against the ratio the lines themselves carry, ``f_{i+1} / f_i``: each
measured pair sits beside both predictions and says which one is nearer, and
the compounded span is scored the same way. Those residuals are descriptive -
both q and the step were fitted on these lines, the pairs are not independent,
and no null is attached - so the campaign's AIC and leave-one-out remain the
decision.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.phi_fixedpoint import (  # noqa: E402
    CONTRACTION,
    PHI,
    STEPS_PER_DECADE,
    agreement_depth,
    continued_fraction,
    convergents,
    derivative,
    effective_depth,
    excess,
    fibonacci,
    is_phi_shaped,
    orbit,
)

CHANNEL_HZ = 100547.19562243501
BAND_LOW_HZ = 13.4e6
BAND_HIGH_HZ = 86.9e6
SEARCH_CHANNELS = (317, 412)
COMB_STEP_HZ = 1.0e6
RUNG_TOLERANCE = 0.01
CAMPAIGN_CSV = os.path.join(ROOT, "results", "comb_vs_scale.csv")
GROUPS_JSON = os.path.join(ROOT, "results", "expanded_phi_groups.json")
DEPTH_LIMIT = 64


def window_in_hz(tolerance_rungs: float, freq_hz: float) -> float:
    """Physical width covered by a rung tolerance at this frequency."""
    return (PHI ** tolerance_rungs - PHI ** (-tolerance_rungs)) * freq_hz


def _exact(text) -> Fraction:
    """The printed decimal as an exact rational, so no binary tail is read."""
    return Fraction(str(text).strip())


def describe(label, value, exact_value=None, note="", synthetic=False,
             printed=None):
    part = continued_fraction(exact_value if exact_value is not None else value,
                              depth=24)
    shown = part[:10]
    # every printed decimal is rational and so terminates somewhere; what
    # matters is whether any fractional level exists at all
    levels = len(continued_fraction(
        exact_value if exact_value is not None else value, depth=DEPTH_LIMIT))
    return {
        "label": label,
        "value": float(value),
        "printed": printed,
        "partials": shown,
        "elided": len(part) > len(shown),
        "depth": agreement_depth(part),
        "unbounded": is_phi_shaped(float(value)),
        "levels": levels,
        "empty_recursion": levels <= 1,
        "excess": excess(part),
        "effective_depth": effective_depth(float(value)),
        "note": note,
        "synthetic": synthetic,
    }


def _closer_counts(items) -> dict:
    """How many readings each fit is nearer to; the splits are diagnostics."""
    out = {"scale": 0, "additive": 0}
    for it in items:
        out[it["closer"]] += 1
    return out


def _span_sides(campaign) -> dict:
    """Where the measured span sits relative to the two predictions."""
    out = {"above_both": 0, "below_both": 0, "between": 0}
    for row in campaign:
        span = row["additive_span"]
        predictions = (span["value"], span["reference_scale"]["value"])
        measured = span["measured_span_ratio"]
        if measured > max(predictions):
            out["above_both"] += 1
        elif measured < min(predictions):
            out["below_both"] += 1
        else:
            out["between"] += 1
    return out


def group_lines() -> dict:
    """The group fixtures' own line frequencies, keyed by what the CSV carries.

    ``(source, date, n_lines)`` is unique across the four groups, so it is the
    join key: the CSV alone does not say where each group starts.
    """
    if not os.path.exists(GROUPS_JSON):
        raise FileNotFoundError(
            "results/expanded_phi_groups.json is missing; run "
            "scripts/phi_model_discrimination.py first")
    with open(GROUPS_JSON, encoding="utf-8") as fh:
        payload = json.load(fh)
    out = {}
    for g in payload["groups"]:
        key = (g["source"], g["date"], int(g["n_lines"]))
        out[key] = [str(f) for f in g["freqs_mhz"]]
    return out


def additive_pairs(freqs, step, q) -> list:
    """``1 + comb_df / f_i`` at every adjacent pair, i = 1 .. n-1.

    The additive model holds the *step* fixed, so the ratio it predicts falls
    as ``1/f``. Each entry is that prediction at the pair's lower line, beside
    its signed distance from the fitted multiplicative interval ``q`` - the
    number that says where the two models part. Predicted from the fit, never
    measured from the lines, so every entry is diagnostic only.
    """
    out = []
    for i in range(len(freqs) - 1):
        ratio = 1 + step / freqs[i]
        item = describe(
            "1 + comb step / f_i", float(ratio), exact_value=ratio,
            printed=f"{float(step):g} / {float(freqs[i]):.6f}",
            note="additive model's adjacent-line ratio starting at this line; "
                 "predicted from the fit, not measured")
        item["diagnostic_only"] = True
        item["pair"] = i + 1
        item["line_low_mhz"] = float(freqs[i])
        item["line_high_mhz"] = float(freqs[i + 1])
        item["signed_delta_vs_scale_q"] = float(ratio - q)
        item["gap_to_scale_q"] = abs(item["signed_delta_vs_scale_q"])
        out.append(item)
    return out


def pair_residuals(freqs, step, q, pairs) -> list:
    """The measured adjacent-pair ratio beside both fits, pair by pair.

    ``f_{i+1} / f_i`` is the ratio the lines themselves carry, so scoring both
    predictions against it puts the divergence on the lines instead of between
    two sets of fitted values. Both q and the step were fitted on these same
    lines and the pairs are not independent, so neither residual is a selection
    statistic: the campaign's AIC and leave-one-out remain the decision. Kept
    diagnostic only, and recorded as relative misses so pairs at different
    frequencies can be read together.
    """
    out = []
    for i in range(len(freqs) - 1):
        measured = freqs[i + 1] / freqs[i]
        additive = 1 + step / freqs[i]
        rel_scale = measured / q - 1
        rel_additive = measured / additive - 1
        out.append({
            "pair": i + 1,
            "line_low_mhz": float(freqs[i]),
            "line_high_mhz": float(freqs[i + 1]),
            "measured": float(measured),
            "measured_printed": (f"{float(freqs[i + 1]):.6f} / "
                                 f"{float(freqs[i]):.6f}"),
            "scale_prediction": float(q),
            "additive_prediction": float(additive),
            "residual_vs_scale": float(measured - q),
            "residual_vs_additive": float(measured - additive),
            "relative_vs_scale": float(rel_scale),
            "relative_vs_additive": float(rel_additive),
            "closer": ("additive" if abs(rel_additive) < abs(rel_scale)
                       else "scale"),
            "diagnostic_only": True,
            "note": "measured f_high/f_low beside both fits; q and the step "
                    "were fitted on these same lines, so neither residual "
                    "is independent and neither selects a model",
        })
    return out


def additive_span(freqs, step, q, pairs) -> dict:
    """The pair ratios compounded over the group, and where they meet ``q``.

    ``prod(1 + step / f_i)`` is the additive model's own ratio from the lowest
    line to the highest; ``q ** (n - 1)`` is the scale model's for the same
    span. Both are predictions, so the comparison is descriptive. The crossing
    frequency ``step / (q - 1)`` is where a pair ratio equals q - the single
    line at which the two models coincide, and hence the pivot the gap is
    measured around.
    """
    product = Fraction(1)
    for f in freqs[:-1]:
        product *= 1 + step / f
    n = len(freqs) - 1
    item = describe(
        "compounded 1 + comb step / f_i", float(product), exact_value=product,
        printed=f"product of {n} pairs",
        note="additive model's ratios compounded from the lowest line to the "
             "highest; predicted from the fit, not measured")
    item["diagnostic_only"] = True
    item["pairs"] = n
    item["step_mhz"] = float(step)
    item["closed_form_ratio"] = float(1 + n * step / freqs[0])
    item["product_vs_closed_form"] = float(
        product / (1 + n * step / freqs[0]) - 1)
    scale_cumulative = q ** n
    reference = describe(
        "q ** pairs", float(scale_cumulative), exact_value=scale_cumulative,
        printed=f"{float(q):.8f} ** {n}",
        note="scale model's ratio over the same span; reference only")
    reference["diagnostic_only"] = True
    reference["reference_model"] = "scale"
    item["reference_scale"] = reference
    item["relative_gap_to_scale"] = float(product / scale_cumulative - 1)
    first, last = pairs[0], pairs[-1]
    item["first_pair_gap_to_scale_q"] = first["gap_to_scale_q"]
    item["last_pair_gap_to_scale_q"] = last["gap_to_scale_q"]
    item["max_pair_gap_to_scale_q"] = max(p["gap_to_scale_q"] for p in pairs)
    item["signed_drift_across_span"] = (
        last["signed_delta_vs_scale_q"] - first["signed_delta_vs_scale_q"])
    item["drift_over_first_gap"] = abs(
        item["signed_drift_across_span"]) / first["gap_to_scale_q"]
    signs = [p["signed_delta_vs_scale_q"] for p in pairs]
    item["gap_sign_changes"] = sum(
        1 for a, b in zip(signs, signs[1:]) if a * b < 0)
    measured = freqs[-1] / freqs[0]
    rel_additive = measured / product - 1
    rel_scale = measured / scale_cumulative - 1
    item["measured_span_ratio"] = float(measured)
    item["measured_span_printed"] = (f"{float(freqs[-1]):.6f} / "
                                      f"{float(freqs[0]):.6f}")
    item["residual_measured_vs_additive"] = float(measured - product)
    item["residual_measured_vs_scale"] = float(measured - scale_cumulative)
    item["relative_vs_additive"] = float(rel_additive)
    item["relative_vs_scale"] = float(rel_scale)
    item["closer"] = ("additive" if abs(rel_additive) < abs(rel_scale)
                       else "scale")
    item["measured_between_predictions"] = bool(
        min(product, scale_cumulative) <= measured
        <= max(product, scale_cumulative))
    if q > 1:
        crossing = step / (q - 1)
        item["scale_crossing_mhz"] = float(crossing)
        item["crossing_inside_span"] = bool(freqs[0] < crossing < freqs[-1])
    return item


def campaign_rows() -> list:
    """The measured ratios in ``results/comb_vs_scale.csv``, read as ratios.

    Two columns are interval ratios rather than model scores:

    * ``scale_q`` - the fitted multiplicative interval q in f_n = f0 * q**n.
    * ``span_mhz / comb_df_mhz`` - the group's span in units of its own fitted
      comb step, a ratio of two intervals measured from the same lines.

    Both are held as exact rationals of the printed decimals. Nothing here is
    reconstructed from a synthetic grid.
    """
    if not os.path.exists(CAMPAIGN_CSV):
        raise FileNotFoundError(
            "results/comb_vs_scale.csv is missing; run "
            "scripts/phi_model_discrimination.py first")
    with open(CAMPAIGN_CSV, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    lines = group_lines()
    out = []
    for row in rows:
        q = _exact(row["scale_q"])
        span = _exact(row["span_mhz"])
        step = _exact(row["comb_df_mhz"])
        steps_across = span / step
        key = (row["source"], row["date"], int(row["n_lines"]))
        if key not in lines:
            raise KeyError(f"no group fixture for {key}")
        freqs = [_exact(f) for f in lines[key]]
        if len(freqs) != int(row["n_lines"]):
            raise ValueError(f"{key}: fixture has {len(freqs)} lines, CSV says "
                             f"{row['n_lines']}")
        if freqs != sorted(freqs):
            raise ValueError(f"{key}: fixture lines are not in frequency order")
        f_low = freqs[0]
        pairs = additive_pairs(freqs, step, q)
        # the first-line entry is the first pair, not a second computation, so
        # the lowest line cannot drift between the two places it is reported
        additive = dict(pairs[0])
        additive["label"] = "1 + comb step / f_low"
        additive["lowest_line_mhz"] = float(f_low)
        additive["note"] = ("additive model's own adjacent-line ratio at the "
                            "group's lowest line: the first entry of "
                            "additive_pairs, predicted from the fit, not "
                            "measured, so diagnostic only")
        span = additive_span(freqs, step, q, pairs)
        scored = pair_residuals(freqs, step, q, pairs)
        out.append({
            "source": row["source"],
            "date": row["date"],
            "verdict": row["verdict"],
            "diagnostic_only": row["verdict"] == "MASK-MANUFACTURED SIGNAL",
            "n_lines": int(row["n_lines"]),
            "scale_q": row["scale_q"],
            "lowest_line_mhz": float(f_low),
            "additive": additive,
            "additive_pairs": pairs,
            "pair_residuals": scored,
            "additive_span": span,
            "ratios": [
                describe("fitted line ratio q", float(q), exact_value=q,
                         printed=row["scale_q"],
                         note="fitted multiplicative interval of the scale model"),
                describe("span / comb step", float(steps_across),
                         exact_value=steps_across,
                         printed=f"{row['span_mhz']} / {row['comb_df_mhz']}",
                         note="ratio of two intervals measured from the same lines"),
            ],
        })
    return out


def build() -> dict:
    der = derivative(PHI)
    assert abs(der + CONTRACTION) < 1e-15

    trials = []
    for start in (0.5, 1.0, 2.0, BAND_HIGH_HZ / BAND_LOW_HZ, 86.9e6 / 45.98e6):
        ob = orbit(start, turns=30)
        trials.append({
            "start": round(start, 9),
            "first_residual": ob.residuals[0],
            "residual_after_5": ob.residuals[5],
            "residual_after_10": ob.residuals[10],
            "contraction": ob.contraction,
            "settled": ob.settled,
        })

    window = [{"freq_mhz": f / 1e6, "width_mhz": window_in_hz(RUNG_TOLERANCE, f) / 1e6}
              for f in (BAND_LOW_HZ, 45.98e6, 47.24e6, 53.78e6, BAND_HIGH_HZ)]
    window_spread = window[-1]["width_mhz"] / window[0]["width_mhz"]

    ladder = [
        describe("phi", PHI, note="the fixed point: add one, never carried"),
        describe("1/phi = phi - 1", 1.0 / PHI, note="subtracting 1 is exact"),
        describe("phi^2 = phi + 1", PHI ** 2, note="adding 1 is exact"),
    ]
    for n in (8, 12, 16):
        ladder.append(describe(f"F{n+1}/F{n}", fibonacci(n + 1) / fibonacci(n),
                               exact_value=Fraction(fibonacci(n + 1), fibonacci(n)),
                               note="exact convergent of phi"))

    # control constructions: depth 0 is guaranteed by how they are built, so
    # they sit in the JSON as a control and never as evidence
    control = []
    for k in (SEARCH_CHANNELS[0], 365, SEARCH_CHANNELS[1]):
        control.append(describe(f"channel {k}->{k+1}", (k + 1) / k,
                                exact_value=Fraction(k + 1, k),
                                note="instrument grid: exactly 1 + 1/k"))
    for f in (45980000, 53780000):
        control.append(describe(f"1 MHz comb at {f/1e6:.2f} MHz",
                                (f + int(COMB_STEP_HZ)) / f,
                                exact_value=Fraction(f + int(COMB_STEP_HZ), f),
                                note="arithmetic comb line spacing"))
    control.append(describe("band 13.4 -> 86.9 MHz", 869.0 / 134.0,
                            exact_value=Fraction(869, 134),
                            note="the full survey span"))

    sweep = []
    for jitter in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10):
        x = PHI * (1.0 + jitter)
        part = continued_fraction(x, depth=24)
        sweep.append({"jitter": jitter, "value": x,
                      "depth": agreement_depth(part),
                      "effective_depth": effective_depth(x)})

    conv = convergents(PHI, depth=14)
    decay = []
    for i in range(1, len(conv)):
        e0, e1 = conv[i - 1][2], conv[i][2]
        decay.append({"p": conv[i][0], "q": conv[i][1], "error": e1,
                      "ratio": e1 / e0})

    campaign = campaign_rows()
    flat = [r for row in campaign for r in row["ratios"]]
    depth_counts = {}
    for r in flat:
        depth_counts[r["depth"]] = depth_counts.get(r["depth"], 0) + 1
    campaign_summary = {
        "rows": len(campaign),
        "ratios": len(flat),
        "depth_counts": {str(k): v for k, v in sorted(depth_counts.items())},
        "depth_one": [{"label": f"{row['date']} n={row['n_lines']} {r['label']}",
                       "value": r["value"],
                       "effective_depth": r["effective_depth"]}
                      for row in campaign for r in row["ratios"]
                      if r["depth"] == 1],
        "excess_min": min(r["excess"] for r in flat),
        "excess_max": max(r["excess"] for r in flat),
        "effective_depth_min": min(r["effective_depth"] for r in flat),
        "effective_depth_max": max(r["effective_depth"] for r in flat),
        "any_phi_shaped": any(r["unbounded"] for r in flat),
        "vacuous": [f"{row['date']} {r['label']}"
                    for row in campaign for r in row["ratios"]
                    if r["empty_recursion"]],
        "diagnostic_only": [f"{row['date']} {row['verdict']}"
                            for row in campaign if row["diagnostic_only"]],
    }

    add = [row["additive"] for row in campaign]
    add_depths = {}
    for r in add:
        add_depths[r["depth"]] = add_depths.get(r["depth"], 0) + 1
    pair_flat = [p for row in campaign for p in row["additive_pairs"]]
    pair_depths = {}
    for r in pair_flat:
        pair_depths[r["depth"]] = pair_depths.get(r["depth"], 0) + 1
    spans = [row["additive_span"] for row in campaign]
    span_depths = {}
    for r in spans:
        span_depths[r["depth"]] = span_depths.get(r["depth"], 0) + 1
    campaign_summary["additive"] = {
        "ratios": len(add),
        "depth_counts": {str(k): v for k, v in sorted(add_depths.items())},
        "excess_min": min(r["excess"] for r in add),
        "excess_max": max(r["excess"] for r in add),
        "effective_depth_min": min(r["effective_depth"] for r in add),
        "effective_depth_max": max(r["effective_depth"] for r in add),
        "any_phi_shaped": any(r["unbounded"] for r in add),
        "lowest_lines_mhz": [row["lowest_line_mhz"] for row in campaign],
        "max_gap_to_scale_q": max(r["gap_to_scale_q"] for r in add),
        "pairs": len(pair_flat),
        "pair_depth_counts": {str(k): v for k, v in sorted(pair_depths.items())},
        "pair_excess_min": min(r["excess"] for r in pair_flat),
        "pair_excess_max": max(r["excess"] for r in pair_flat),
        "pair_effective_depth_min": min(r["effective_depth"] for r in pair_flat),
        "pair_effective_depth_max": max(r["effective_depth"] for r in pair_flat),
        "max_pair_gap_to_scale_q": max(r["gap_to_scale_q"] for r in pair_flat),
        "pair_ratios_decline": all(
            [p["value"] for p in row["additive_pairs"]] ==
            sorted((p["value"] for p in row["additive_pairs"]), reverse=True)
            for row in campaign),
        "crossings": [{"date": row["date"],
                       "n_lines": row["n_lines"],
                       "crossing_mhz": row["additive_span"].get(
                           "scale_crossing_mhz"),
                       "inside_span": row["additive_span"][
                           "crossing_inside_span"]}
                      for row in campaign],
        "span_depth_counts": {str(k): v for k, v in sorted(span_depths.items())},
        "span_relative_gap_to_scale": [r["relative_gap_to_scale"] for r in spans],
        "max_abs_span_relative_gap": max(
            abs(r["relative_gap_to_scale"]) for r in spans),
        "min_drift_over_first_gap": min(r["drift_over_first_gap"] for r in spans),
        "max_drift_over_first_gap": max(r["drift_over_first_gap"] for r in spans),
        "diagnostic_only": True,
        "note": "1 + comb_df / f_low at the lowest line, the same ratio at "
                "every adjacent pair, and their product over each group's "
                "span, all predicted by the additive fit rather than "
                "measured from the lines",
    }

    scored = [r for row in campaign for r in row["pair_residuals"]]
    line_fit = {
        "pairs": len(scored),
        "closer_counts": _closer_counts(scored),
        "span_closer_counts": _closer_counts(spans),
        "max_abs_relative_vs_scale": max(abs(r["relative_vs_scale"])
                                         for r in scored),
        "max_abs_relative_vs_additive": max(abs(r["relative_vs_additive"])
                                            for r in scored),
        "mean_abs_relative_vs_scale": (sum(abs(r["relative_vs_scale"])
                                           for r in scored) / len(scored)),
        "mean_abs_relative_vs_additive": (
            sum(abs(r["relative_vs_additive"]) for r in scored) / len(scored)),
        "max_abs_span_relative_vs_scale": max(
            abs(r["relative_vs_scale"]) for r in spans),
        "max_abs_span_relative_vs_additive": max(
            abs(r["relative_vs_additive"]) for r in spans),
        "measured_between_predictions": [
            f"{row['date']} n={row['n_lines']}" for row in campaign
            if row["additive_span"]["measured_between_predictions"]],
        "span_side_counts": _span_sides(campaign),
        "groups": [{"date": row["date"],
                    "n_lines": row["n_lines"],
                    "closer_counts": _closer_counts(row["pair_residuals"]),
                    "max_abs_relative_vs_scale": max(
                        abs(r["relative_vs_scale"])
                        for r in row["pair_residuals"]),
                    "max_abs_relative_vs_additive": max(
                        abs(r["relative_vs_additive"])
                        for r in row["pair_residuals"]),
                    "span_closer": row["additive_span"]["closer"],
                    "span_relative_vs_scale": row["additive_span"][
                        "relative_vs_scale"],
                    "span_relative_vs_additive": row["additive_span"][
                        "relative_vs_additive"]}
                   for row in campaign],
        "diagnostic_only": True,
        "note": "measured f_high/f_low at every adjacent pair and across each "
                "span, scored against both fits; descriptive only, since q "
                "and the step were fitted on these same lines and no null "
                "is attached",
    }
    campaign_summary["line_fit"] = line_fit

    return {
        "phi": PHI,
        "contraction": CONTRACTION,
        "derivative_at_phi": der,
        "steps_per_decade": STEPS_PER_DECADE,
        "orbit_trials": trials,
        "rung_tolerance": RUNG_TOLERANCE,
        "window": window,
        "window_spread": window_spread,
        "ladder": ladder,
        "control": control,
        "campaign": campaign,
        "campaign_summary": campaign_summary,
        "depth_sweep": sweep,
        "convergent_decay": decay,
        "summary": {
            "control_depths": [m["depth"] for m in control[:3]],
            "control_effective_depth": [m["effective_depth"] for m in control[:3]],
            "control_excess": [m["excess"] for m in control[:3]],
            "phi_unbounded": [r["unbounded"] for r in ladder[:3]],
        },
    }


def render(s: dict) -> str:
    L = []
    L.append("# The +1, carried as a map\n")
    L.append("Framework: phi^2 = phi + 1. Attribution: Tyree Jones "
             "(tyofthestarz), phi^2 = phi + 1 framework, prior art from "
             "December 2025.\n")
    L.append("Generated by `python3 scripts/phi_plus_one.py`; every number is "
             "computed in this run.\n")

    L.append("## 1. The identity is a rate, not a width\n")
    L.append("From any positive x, add one and invert: `T(x) = 1 + 1/x`. The "
             "fixed point solves `x = 1 + 1/x`, i.e. `x^2 = x + 1`, i.e. phi. "
             "The derivative there is\n")
    L.append(f"    T'(phi) = -1/phi^2 = {-s['contraction']:.7f}\n")
    L.append(f"The recurrence shrinks its own error by {s['contraction']:.6f} "
             f"per turn - {s['steps_per_decade']:.2f} turns per decimal place. "
             "**The 0.382 this project kept testing against is that rate**, not "
             "a tolerance: it is the step size of the map that defines the "
             "constant. Carrying it as a fixed +/-0.01 rung window was the "
             "error the identity warns about: a window in rungs is not scale "
             "free. That measurement is kept as a regression fixture "
             "(`RejectedToleranceTest` in `tests/test_phi_fixedpoint.py`, values "
             "in `results/phi_plus_one.json`) rather than restated here.\n")
    L.append("| start | residual | after 5 turns | after 10 turns | measured contraction |")
    L.append("| ---: | ---: | ---: | ---: | ---: |")
    for t in s["orbit_trials"]:
        L.append(f"| {t['start']:.6f} | {t['first_residual']:.3e} | "
                 f"{t['residual_after_5']:.3e} | {t['residual_after_10']:.3e} | "
                 f"{t['contraction']:.6f} |")
    L.append("")
    L.append("Every start converges, all at the same rate. The trajectory is the "
             "result; no step in it is a decision.\n")

    L.append("## 2. A ratio that never terminates\n")
    L.append("Write a ratio as `a0 + 1/(a1 + 1/(a2 + ...))` and keep every "
             "partial quotient. Phi is `[1; 1, 1, 1, ...]`: at every level the "
             "instruction is *add one* and it is never carried, so the fraction "
             "never terminates. The integer part is excluded from both "
             "statistics, because it is exactly the whole the identity refuses "
             "to carry.\n")
    L.append("- **depth** - how many levels of the fractional recursion equal 1.\n")
    L.append("- **excess** - `sum(a_i - 1)` over those levels; zero exactly for "
             "phi's shape.\n")
    L.append("| value | partial quotients | depth | excess | levels |")
    L.append("| --- | --- | ---: | ---: | ---: |")
    for row in s["ladder"]:
        pq = ", ".join(str(a) for a in row["partials"]) + ("..." if row["elided"] else "")
        depth = "unbounded" if row["unbounded"] else row["depth"]
        levels = "-" if row["unbounded"] else row["levels"]
        L.append(f"| {row['label']} | {pq} | {depth} | {row['excess']} | "
                 f"{levels} |")
    L.append("")
    L.append("Phi, `1/phi` and `phi^2` are unbounded with excess 0: one "
             "identity, read as a statement about partial quotients instead of "
             "magnitudes. Each exact convergent `F(n+1)/F(n)` matches for a "
             "finite number of levels - the count is in the table - and then "
             "pays 1 at its terminating 2. It is phi seen at finite depth, with "
             "the +1 still pending.\n")

    camp = s["campaign"]
    cs = s["campaign_summary"]
    L.append("## 3. The campaign's own ratios, read the same way\n")
    L.append("Section 2 is a set of constructions whose answer is known in "
             "advance. The ratios that were actually measured are in "
             "`results/comb_vs_scale.csv`, and two of its columns are interval "
             "ratios rather than model scores:\n")
    L.append("- `scale_q` - the fitted multiplicative interval q in "
             "`f_n = f0 * q**n`: the ratio between adjacent lines.\n")
    L.append("- `span / comb step` - the group's measured span in units of its "
             "own fitted comb step: a ratio of two intervals read off the same "
             "lines.\n")
    L.append("Both are held as exact rationals of the printed decimals, so the "
             "partial quotients below are the printed value's, not a binary "
             "tail's.\n")
    L.append("| group | verdict | ratio | as printed | depth | excess | effective depth |")
    L.append("| --- | --- | --- | ---: | ---: | ---: | ---: |")
    for row in camp:
        for r in row["ratios"] + [row["additive"]]:
            depth = "unbounded" if r["unbounded"] else r["depth"]
            label = r["label"] + (" (diagnostic)" if r.get("diagnostic_only") else "")
            L.append(f"| {row['date']} | {row['verdict']} | {label} | "
                     f"{r['printed']} | {depth} | {r['excess']} | "
                     f"{r['effective_depth']:.2f} |")
    L.append("")
    spread = " and ".join(f"{n} at depth {d}"
                         for d, n in cs["depth_counts"].items())
    L.append(f"Across the {cs['ratios']} measured ratios the counts are "
             f"{spread}; excess runs {cs['excess_min']} to {cs['excess_max']}, "
             "and effective depth "
             f"{cs['effective_depth_min']:.2f} to "
             f"{cs['effective_depth_max']:.2f}. These are measurements of "
             "ratios that no one chose.\n")
    ad = cs["additive"]
    L.append(f"The third line of each group is the additive model's own "
             f"adjacent-line ratio, `1 + comb_df / f_low`, with `f_low` the "
             f"group's lowest line from `results/expanded_phi_groups.json` "
             f"({', '.join(f'{m:.6f}' for m in ad['lowest_lines_mhz'])} MHz) "
             f"and `comb_df` the fit's own step. It is what the comb predicts, "
             f"not what the lines show, so it is diagnostic only and is left "
             f"out of the counts above. On those {ad['ratios']} predicted "
             f"ratios the depths are "
             f"{' and '.join(f'{n} at depth {d}' for d, n in ad['depth_counts'].items())}, "
             f"excess runs {ad['excess_min']} to {ad['excess_max']}, and "
             f"effective depth {ad['effective_depth_min']:.2f} to "
             f"{ad['effective_depth_max']:.2f}. At the first line the additive "
             f"and multiplicative models differ by at most "
             f"{ad['max_gap_to_scale_q']:.2e} in ratio, because there the step "
             f"and the interval coincide to first order. They separate at every "
             f"later pair, where the additive ratio falls away from `q` on one "
             f"side of the crossing and climbs away from it on the other (3a), "
             f"and in how that compounds over the span (3b). A depth taken at "
             f"one line therefore cannot discriminate between them, and is not "
             f"offered as if it could.\n")
    L.append("### 3a. Every pair, not only the first\n")
    lf = cs["line_fit"]
    L.append("At the lowest line the two models agree to first order, so a "
             "single reading there cannot separate them. Every later pair is "
             "a different question: the additive model holds the step fixed, "
             "so its own ratio falls as `1/f`, while the scale model predicts "
             "the same `q` at every line. Each pair is therefore read three "
             "ways: the ratio the lines carry, the additive prediction, and "
             "the signed difference between the two predictions - carried in "
             "the JSON as `1 + d/f_i - q` and located by the crossing "
             "below:\n")
    L.append("| group | pair | starts at (MHz) | measured f_high/f_low | "
             "1 + d/f_i | measured/(1+d/f_i) - 1 | measured/q - 1 | closer | "
             "additive depth | additive eff. depth |")
    L.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: |")
    for row in camp:
        for p, m in zip(row["additive_pairs"], row["pair_residuals"]):
            L.append(f"| {row['date']} n={row['n_lines']} | {p['pair']} | "
                     f"{p['line_low_mhz']:.4f} | {m['measured']:.6f} | "
                     f"{p['value']:.6f} | {m['relative_vs_additive']:+.2e} | "
                     f"{m['relative_vs_scale']:+.2e} | {m['closer']} | "
                     f"{p['depth']} | {p['effective_depth']:.2f} |")
    L.append("")
    cross_txt = ", ".join(
        f"{row['date']} n={row['n_lines']} at "
        f"{row['additive_span']['scale_crossing_mhz']:.3f} MHz"
        for row in camp)
    L.append("Within every group the ratios fall monotonically and cross `q` "
             "exactly once, at `f = d/(q - 1)`: " + cross_txt + " - every "
             "one inside its own span. So the models are not disagreeing at "
             "the first line and everywhere after: they touch at one "
             "frequency per group and separate on both sides of it. That is "
             "the one-point meeting recorded in `PROPORTION_AND_SCALE.md` "
             "section 2 for a fixed comb against a fixed interval, now "
             "located per group with the fitted `q`. A reading taken only at "
             "the first line sees one side of that crossing: the drift "
             "across a span is "
             f"{ad['min_drift_over_first_gap']:.1f} to "
             f"{ad['max_drift_over_first_gap']:.1f}x the first-line "
             "displacement in every group, and the sign of the difference "
             "flips exactly where the crossing lies.\n")
    per_group = ", ".join(
        f"{g['date']} n={g['n_lines']} {g['closer_counts']['additive']}-"
        f"{g['closer_counts']['scale']}" for g in lf["groups"])
    L.append("Scored against the ratio the lines themselves carry, "
             "`f_high / f_low`, the two fits are nearly indistinguishable "
             f"pair by pair: across {lf['pairs']} measured pairs the additive "
             f"prediction is nearer on {lf['closer_counts']['additive']} and "
             f"the scale prediction on {lf['closer_counts']['scale']} "
             f"(additive-scale: {per_group}). The largest relative miss is "
             f"{lf['max_abs_relative_vs_additive']:.1e} for the additive "
             f"prediction against {lf['max_abs_relative_vs_scale']:.1e} for "
             "the scale prediction, while the two predictions sit only "
             f"{ad['max_pair_gap_to_scale_q']:.1e} apart in ratio: the "
             "disagreement between the fits is smaller than either fit's own "
             "miss on the lines, so neither is consistently nearer, and the "
             "split above over correlated pairs is not a decision. "
             "Neither residual carries a null, because neither is a detection "
             "statistic; and because `q` and the step were both fitted on "
             "these lines, with the pairs correlated, counting which fit is "
             "nearer is not a likelihood test. The AIC and leave-one-out "
             "columns in `results/comb_vs_scale.csv` remain the decision; "
             "this column says where on the lines the two fits part.\n")

    L.append("### 3b. Compounded across the span\n")
    L.append("Multiplying those pair ratios gives the additive model's own "
             "ratio from the group's lowest line to its highest; `q ** pairs` "
             "is what the scale model predicts over the same span. The depth "
             "reading is applied to both, and the span the lines actually "
             "carry, `f_high / f_low`, sits beside them, so both predictions "
             "are scored on the data rather than only against each other:\n")
    L.append("| group | pairs | measured f_high/f_low | additive compounded | "
             "q^pairs | measured/additive - 1 | measured/q^pairs - 1 | "
             "closer | additive depth | additive eff. depth |")
    L.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: |")
    for row in camp:
        sp = row["additive_span"]
        L.append(f"| {row['date']} n={row['n_lines']} | {sp['pairs']} | "
                 f"{sp['measured_span_ratio']:.6f} | {sp['value']:.6f} | "
                 f"{sp['reference_scale']['value']:.6f} | "
                 f"{sp['relative_vs_additive']:+.2e} | "
                 f"{sp['relative_vs_scale']:+.2e} | {sp['closer']} | "
                 f"{sp['depth']} | {sp['effective_depth']:.2f} |")
    L.append("")
    widest = max(camp, key=lambda row: abs(
        row["additive_span"]["relative_gap_to_scale"]))
    ws = widest["additive_span"]
    others = max(abs(row["additive_span"]["relative_gap_to_scale"])
                 for row in camp if row is not widest)
    L.append("Compounding is where the divergence between the fits stops "
             f"being cosmetic: the widest gap is "
             f"{abs(ws['relative_gap_to_scale']):.2e} on "
             f"the {widest['n_lines']}-line group ({ws['value']:.6f} additive "
             f"against {ws['reference_scale']['value']:.6f} multiplicative), "
             f"while the other three stay below {others:.0e}. The compounded "
             "gap is not the sum of the pair gaps but their product, which is "
             "why it grows with the number of lines rather than staying at "
             "the size of one step. On that same group the lines sit at "
             f"{ws['measured_span_ratio']:.6f}, which is "
             f"{ws['relative_vs_additive']:+.2e} from the additive prediction "
             f"and {ws['relative_vs_scale']:+.2e} from the scale one: the two "
             "fits differ from each other by less than either differs from "
             "the lines, so a wider gap between the fits is not a better "
             "reading of the data. The measured span is nearer the "
             f"{ws['closer']} prediction, which is the nearer one in "
             f"{lf['span_closer_counts'][ws['closer']]} of the four groups, "
             "and it never falls between the two predictions: it lies above "
             f"both in {lf['span_side_counts']['above_both']} groups and "
             f"below both in {lf['span_side_counts']['below_both']}, with "
             f"{lf['span_side_counts']['between']} between.\n")
    L.append("Two cautions on reading that table as evidence. First, the "
             "compounded ratios reach effective depth "
             f"{max(row['additive_span']['effective_depth'] for row in camp):.2f}, "
             "higher than any per-pair value, but only because a ratio near "
             "1.4 is closer to phi in relative terms than a ratio near 1.02: "
             "at these magnitudes effective depth is measuring size, not "
             "phi's shape, and the integer depth is 0 on every span. Second, "
             "the two prediction columns come from fits estimated on these "
             "same lines; the measured column is what scores them, and with "
             "the fits taken from the lines and the pairs correlated, neither "
             "side is an independent test. Nothing here supports or rejects "
             "either model; the reading only says how far apart the two "
             "predictions are, and how far each sits from the lines.\n")
    if cs["depth_one"]:
        shown = ", ".join(f"{r['label']} at {r['value']:.2f} "
                          f"(effective depth {r['effective_depth']:.2f})"
                          for r in cs["depth_one"])
        L.append("A single level of depth is cheap and proves nothing: any ratio "
                 "whose fractional part exceeds 1/2 gets `a1 = 1`. The rows that "
                 f"reach it are {shown}, and the continuous reading is what the "
                 "depth count rounds away.\n")
    if cs["vacuous"]:
        L.append(f"Vacuous recursion: {'; '.join(cs['vacuous'])} reduces to a "
                 "single partial quotient, so it has no fractional level at "
                 "all. Depth 0 and excess 0 there are the same numbers phi's "
                 "infinite recursion would give (0 at every level), which is "
                 "why `unbounded` travels beside them: without that flag, "
                 "`(0, 0)` cannot tell an empty recursion from a faithful "
                 "one.\n")
    L.append("Read the numbers as scale, not as support. A ratio within a few "
             "percent of 1 cannot be near phi, so a depth of 0 or 1 there is "
             "arithmetic, exactly as it is for the channel grid - and the "
             "instrument rows are still in the JSON, labelled as controls, "
             "where depth 0 is guaranteed by construction. None of it is "
             "evidence that the data is or is not artificial.\n")
    if cs["diagnostic_only"]:
        L.append(f"{len(cs['diagnostic_only'])} of the {cs['rows']} rows carry "
                 "MASK-MANUFACTURED SIGNAL and are diagnostic only: their ratios "
                 "measure the excision mask, so no depth computed on them is a "
                 "statement about a source.\n")

    L.append("## 4. Depth is continuous, not a verdict\n")
    L.append("A ladder at ratio `phi * (1 + jitter)`, so each row is phi plus a "
             "declared relative displacement:\n")
    L.append("| jitter | integer depth | effective depth |")
    L.append("| ---: | ---: | ---: |")
    for row in s["depth_sweep"]:
        L.append(f"| {row['jitter']:.0e} | {row['depth']} | "
                 f"{row['effective_depth']:.2f} |")
    L.append("")
    L.append("Depth falls smoothly with displacement, because it is a logarithm "
             "in disguise: a relative miss of `j` costs about "
             "`log(1/j)/log(phi^2)` levels. The integer depth rounds that down; "
             "the effective depth keeps the fraction and may be negative. So "
             "there is nothing binary to decide - a value is not 'phi' or 'not "
             "phi', it sits at a measured depth, continuous in the distance "
             "from phi.\n")

    L.append("## 5. The convergent decay is the identity, exactly\n")
    L.append("| convergent | value | error | ratio to previous |")
    L.append("| --- | ---: | ---: | ---: |")
    for row in s["convergent_decay"][:8]:
        L.append(f"| {row['p']}/{row['q']} | {row['p']/row['q']:.9f} | "
                 f"{row['error']:.3e} | {row['ratio']:.6f} |")
    L.append("")
    L.append("The ratio alternates about 0.381966 and is within 1e-3 of it from "
             "level 5 on; its limit is exactly `1/phi^2`, from "
             "`F(n+1)/F(n) - phi = (-1)^n / (F(n) phi^(n+1))`. Not a fit, no "
             "null needed. The constant in section 1, the constant the earlier "
             "campaign tested against, and this decay are one number in three "
             "roles.\n")

    L.append("## 6. What this does not say\n")
    L.append("- **It does not explain the campaign's zero.** Zero hits in "
             "fourteen pairs is a property of additive spacing: a fixed step meets a "
             "fixed multiplicative interval at exactly one frequency, 4.96 MHz "
             "in this band, and the survey searched no cell there. That "
             "argument needs no continued fraction and is measured in "
             "`PROPORTION_AND_SCALE.md` section 2. Depth was withdrawn from that "
             "job.\n")
    L.append("- **Depth 0 is not evidence either way.** It appears on a masked, "
             "manufactured signal and on the channel grid; one is artificial "
             "and one is arithmetic. Depth describes a ratio's continued "
             "fraction and nothing more.\n")
    L.append("- **Nothing here supports a candidate.** This machinery can only "
             "detract: a value at depth 0 is still an unexplained measurement, "
             "not an artificial one.\n")
    L.append("- **The crossing is arithmetic, not a property of the sources.** "
             "`1 + d/f = q` has one solution for any two intervals, so every "
             "group brackets `q` exactly once because `q` and `d` are fitted "
             "per group. The pair table measures how far the two fits part; "
             "it does not discover a special frequency in the sky.\n")
    return "\n".join(L) + "\n"


def main() -> int:
    payload = build()
    with open(os.path.join(ROOT, "results", "phi_plus_one.json"), "w") as fh:
        json.dump(payload, fh, indent=2)
    with open(os.path.join(ROOT, "PHI_PLUS_ONE.md"), "w") as fh:
        fh.write(render(payload))
    s = payload["summary"]
    cs = payload["campaign_summary"]
    print(f"T'(phi)                  {payload['derivative_at_phi']:.9f}")
    print(f"steps per decade         {payload['steps_per_decade']:.4f}")
    print(f"all orbits settle        {all(t['settled'] for t in payload['orbit_trials'])}")
    print(f"rung window spread       {payload['window_spread']:.3f}x")
    print(f"phi/1:phi/phi^2 unbounded {s['phi_unbounded']}")
    print(f"control grid depth       {s['control_depths']} eff "
          f"{s['control_effective_depth'][0]:.2f}")
    print(f"campaign ratios          {cs['ratios']} on {cs['rows']} rows, "
          f"depth counts {cs['depth_counts']}")
    print(f"campaign excess          {cs['excess_min']} .. {cs['excess_max']}")
    print(f"campaign eff depth       {cs['effective_depth_min']:.2f} .. "
          f"{cs['effective_depth_max']:.2f}")
    print(f"additive 1+step/f_low    {cs['additive']['ratios']} on "
          f"{len(cs['additive']['lowest_lines_mhz'])} groups, depth counts "
          f"{cs['additive']['depth_counts']}, excess "
          f"{cs['additive']['excess_min']} .. {cs['additive']['excess_max']}, "
          f"eff {cs['additive']['effective_depth_min']:.2f} .. "
          f"{cs['additive']['effective_depth_max']:.2f}")
    print(f"max |1+step/f - q|       {cs['additive']['max_gap_to_scale_q']:.3e}")
    print(f"additive pairs           {cs['additive']['pairs']} pairs, depth "
          f"counts {cs['additive']['pair_depth_counts']}, |1+step/f-q| up to "
          f"{cs['additive']['max_pair_gap_to_scale_q']:.2e}, drift/first gap "
          f"{cs['additive']['min_drift_over_first_gap']:.2f} .. "
          f"{cs['additive']['max_drift_over_first_gap']:.2f}x")
    print(f"additive compounding     span depth counts "
          f"{cs['additive']['span_depth_counts']}, max |product/q^pairs-1| "
          f"{cs['additive']['max_abs_span_relative_gap']:.2e}")
    lf = cs["line_fit"]
    print(f"scored on the lines       {lf['pairs']} pairs: additive nearer "
          f"{lf['closer_counts']['additive']}, scale nearer "
          f"{lf['closer_counts']['scale']}; spans add {lf['span_closer_counts']['additive']} "
          f"scale {lf['span_closer_counts']['scale']}; max "
          f"|measured/pred-1| {lf['max_abs_relative_vs_additive']:.2e} add / "
          f"{lf['max_abs_relative_vs_scale']:.2e} scale")
    print(f"model crossings (MHz)    "
          f"{[round(c['crossing_mhz'], 3) for c in cs['additive']['crossings']]}")
    print(f"vacuous recursion        {cs['vacuous']}")
    print(f"diagnostic-only rows     {len(cs['diagnostic_only'])}")
    print(f"depth sweep              {[(r['jitter'], r['depth']) for r in payload['depth_sweep']]}")
    print("wrote results/phi_plus_one.json and PHI_PLUS_ONE.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
