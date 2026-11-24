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

A residual in ratio units still does not say whether it is large. The campaign
table records the scale each fit was judged by: ``rms_comb_channels`` and
``rms_scale_channels`` are the rms of each fit's own residuals on these same
lines, in channels. A ratio miss at a pair displaces that pair's upper line,
so multiplying gives an offset in frequency and one rms divides it - the
largest pair miss, the compounded span miss and the gap between the two fits
all read as fractions of the scatter each fit leaves. Still diagnostic only:
the rms is one number per group, estimated from the same lines, so no
chi-square and no null is formed, and none of it enters the model-selection
counts.

That divisor is still one in-sample number, so the campaign's own leave-one-out
columns are read beside it: ``loo_comb_mhz`` and ``loo_scale_mhz`` are each
fit's mean error on the line it never saw, and dividing that by the same
group's rms says whether one rms was a fair summary of the scatter. It is on
the group where the lines far outnumber the two fitted parameters, and it is
optimistic by roughly a factor two on the groups of four. The group of three
has no held-out fold - three lines leave two, and two points fit either model
exactly - so it is reported as not applicable and no claim is made about it.
The nearer-fit split is checked against the fold majority the same way.

A mean over folds can still hide a single line, so each fold's error already
recorded in ``results/expanded_phi_groups.json`` is read and set beside the
mean the campaign table prints. The worst fold is reported as its share of the
group's total, and the multiple of an even split and the mean-without-worst
ratio are left as arithmetic on that share and the fold count rather than run
as further columns. The share carries a null of its own - the largest of n
exchangeable shares, uniform on the simplex with the total fixed - so a share
can be read as unusual or not, instead of only as large. Both the held-out
error, its per-fold spread, the share's null and the rms stay diagnostic: the
columns are read as printed, never recomputed, the fold file is only read, and
the AIC and leave-one-out decisions are untouched.
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
# sqrt(2/pi): the ratio a mean absolute error has to the rms of the same
# Gaussian scatter, used to put the held-out column on the rms's footing
GAUSSIAN_MAE = math.sqrt(2.0 / math.pi)


def window_in_hz(tolerance_rungs: float, freq_hz: float) -> float:
    """Physical width covered by a rung tolerance at this frequency."""
    return (PHI ** tolerance_rungs - PHI ** (-tolerance_rungs)) * freq_hz


def _exact(text) -> Fraction:
    """The printed decimal as an exact rational, so no binary tail is read."""
    return Fraction(str(text).strip())


def _median(values) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("no values to take a median of")
    return float(ordered[len(ordered) // 2])


def fit_noise(row) -> dict:
    """The scatter each fit leaves on these lines, read from the table.

    ``rms_comb_channels`` and ``rms_scale_channels`` are the rms of each
    fit's own residuals, in channels, so one channel of miss is one of these
    units by definition. Read as printed and never re-derived from the
    frequencies. A blank value raises: a residual cannot be read against a
    scatter that was not recorded.
    """
    out = {}
    for model in ("comb", "scale"):
        key = f"rms_{model}_channels"
        text = str(row.get(key, "")).strip()
        if not text:
            raise ValueError(
                f"{row.get('date')}: {key} is empty; a miss cannot be "
                "read against a scatter that was not recorded")
        value = _exact(text)
        if value <= 0:
            raise ValueError(
                f"{row.get('date')}: {key} must be positive, got {text}")
        out[model] = value
    return out


def line_offset(relative, predicted_high_mhz, rms) -> dict:
    """A ratio residual as a displacement of the pair's upper line.

    ``relative = measured / predicted - 1`` is a miss in ratio, and the line
    it moves is the predicted upper line of the pair: their product is the
    miss in frequency, which ``rms`` - in channels - divides. ``over_comb_rms``
    and ``over_scale_rms`` carry both denominators; callers take the one
    belonging to the fit that made the prediction.
    """
    offset_mhz = predicted_high_mhz * relative
    channels = offset_mhz * 1e6 / CHANNEL_HZ
    return {
        "offset_mhz": float(offset_mhz),
        "offset_channels": float(channels),
        "over_comb_rms": float(channels / rms["comb"]),
        "over_scale_rms": float(channels / rms["scale"]),
    }


def held_out_noise(row, rms) -> dict:
    """What a line the fit never saw costs, beside the rms it was judged by.

    ``loo_comb_mhz`` and ``loo_scale_mhz`` in the campaign table are each
    fit's mean error on the line omitted from the fit, averaged over the
    folds, in MHz. They are read as printed and never re-derived. A blank
    value is the table's way of saying leave-one-out was not defined on that
    group - three lines leave two, and two points fit either model exactly -
    so it is reported as ``applicable: False`` rather than as a zero.

    The quotient with the group's own rms is the cross-check: the rms is an
    in-sample summary of the same lines, and the held-out error is those
    lines' scatter with the fit removed. One rms of held-out error means the
    single number was a fair summary; more means it was optimistic. Because
    the held-out column is a mean absolute error and the rms a root mean
    square, ``sigma_over_rms`` divides out ``sqrt(2/pi)``, the ratio those
    two summaries have for Gaussian scatter, so 1.0 there means the held-out
    line sat as far out as the in-sample rms predicts.
    """
    texts = {model: str(row.get(f"loo_{model}_mhz", "")).strip()
             for model in ("comb", "scale")}
    defined = [model for model in ("comb", "scale") if texts[model]]
    if len(defined) == 1:
        raise ValueError(
            f"{row.get('date')}: only one leave-one-out column is filled "
            f"({texts}); both fits are held out together or neither is")
    majority = str(row.get("loo_majority", "")).strip()
    if not majority:
        raise ValueError(f"{row.get('date')}: loo_majority is empty")
    out = {
        "as_printed": [texts["comb"], texts["scale"]],
        "loo_majority": majority,
        "applicable": bool(defined),
        "loo_comb_mhz": None,
        "loo_scale_mhz": None,
        "loo_comb_channels": None,
        "loo_scale_channels": None,
        "held_out_over_rms": {"comb": None, "scale": None},
        "sigma_over_rms": {"comb": None, "scale": None},
        "gaussian_mae": GAUSSIAN_MAE,
        "note": "each fit's mean error on the held-out line, read as "
                "printed and divided by the same group's in-sample rms; "
                "diagnostic only, and never recomputed from the lines",
    }
    if not defined:
        return out
    for model in ("comb", "scale"):
        value = _exact(texts[model])
        if value <= 0:
            raise ValueError(
                f"{row.get('date')}: loo_{model}_mhz must be positive, "
                f"got {texts[model]}")
        channels = float(value) * 1e6 / CHANNEL_HZ
        over = channels / float(rms[model])
        out[f"loo_{model}_mhz"] = float(value)
        out[f"loo_{model}_channels"] = channels
        out["held_out_over_rms"][model] = over
        out["sigma_over_rms"][model] = over / GAUSSIAN_MAE
    return out


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


def group_folds() -> dict:
    """Each leave-one-out fold's own error, keyed like the campaign table.

    ``results/expanded_phi_groups.json`` already records the runs that made
    the ``loo_*_mhz`` columns: every fold's held-out line and the error both
    fits made on it. Those per-fold numbers are read here and never
    recomputed; ``fold_spread`` only checks that their mean is the mean the
    campaign table prints, so the two files are one measurement seen twice.
    A group with no recorded leave-one-out block is returned as ``None`` and
    reported as not applicable rather than as a zero.
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
        out[key] = (g.get("comparison") or {}).get("loo")
    return out


def share_max_tail(share, n) -> float:
    """Probability the largest of ``n`` exchangeable shares reaches ``share``.

    The null assumes exchangeability and nothing else: the group's held-out
    error could have fallen on any of its folds with no fold favoured, and the
    total is fixed, so the share vector is uniform on the simplex,
    Dirichlet(1, ..., 1). The largest of those shares is not pinned at 1/n -
    even an even scatter of n numbers has a largest member - so the question
    the share cannot answer alone is how large the largest gets by chance.

    Its upper tail follows by inclusion-exclusion over the shares that exceed
    x, ``P(max >= x) = sum_{k=1}^{floor(1/x)} (-1)^(k+1) C(n, k) (1 - k x)^
    (n-1)``, which returns exactly 1 at x = 1/n because the maximum cannot fall
    below an even split. Evaluated in Fractions at the share as the table
    prints it, so the two columns reproduce by hand: at n = 14 the alternating
    terms nearly cancel, and a float round-off there could decide the answer.
    """
    if n < 2:
        raise ValueError("a share null needs at least two folds")
    x = Fraction(float(share)).limit_denominator(10 ** 9)
    if not 0.0 < x <= 1.0:
        raise ValueError(f"a share must lie in (0, 1], got {share}")
    tail = Fraction(0)
    k = 1
    while k * x <= 1:
        term = Fraction(math.comb(n, k)) * (1 - k * x) ** (n - 1)
        tail += term if k % 2 else -term
        k += 1
    return min(1.0, max(0.0, float(tail)))


def fold_spread(held, loo, n_lines) -> dict:
    """Whether one held-out line carries a group's mean error.

    ``held`` is the group's block from ``held_out_noise``, whose mean error
    was read as printed from the campaign table. ``loo`` is the same group's
    ``comparison.loo`` block from ``results/expanded_phi_groups.json``, which
    holds one error per fold. The two describe the same runs, so the mean of
    the folds must agree with the printed mean to the last printed decimal;
    that is checked rather than assumed, and a disagreement raises instead of
    being averaged away.

    The worst fold is reported as its share of the group's total error. A mean
    of n folds gives each line 1/n of the total, so the share is also the
    multiple of an even split divided by n: a share of 1/n is an even spread
    and a share approaching 1 is one line accounting for the whole mean. The
    share is what travels, because the multiple and the mean-without-worst
    ratio are arithmetic on it and on n - ``n share`` and
    ``n (1 - share) / (n - 1)`` - and tabulating three forms of one number
    invites reading them as three findings. Each share carries a null, the
    largest of n exchangeable shares, so a value can be read as unusual or
    not; no threshold, no outlier test and no verdict is attached, and no fold
    is removed from anything.
    """
    out = {"n_folds": 0, "models": {}, "same_worst_line": None,
           "share_null": {
               "name": "largest of n exchangeable shares",
               "distribution": "Dirichlet(1, ..., 1)",
               "tail": "P(max >= x) = sum_{k=1}^{floor(1/x)} (-1)^(k+1) "
                       "C(n, k) (1 - k x)^(n-1)",
               "diagnostic_only": True,
               "note": "the fold shares are uniform on the simplex, so the "
                       "total is held fixed and no fold is favoured; the "
                       "null is on the share alone, says nothing about "
                       "which fit is right, and enters no selection count"},
           "diagnostic_only": True,
           "note": "per-fold errors read from "
                   "results/expanded_phi_groups.json and set beside the mean "
                   "printed in results/comb_vs_scale.csv; the worst fold is "
                   "reported as its share of the group total, with that "
                   "share's own null, and nothing enters the model-selection "
                   "counts"}
    if not held["applicable"]:
        return out
    folds = (loo or {}).get("folds") or []
    out["n_folds"] = len(folds)
    if len(folds) != int(n_lines):
        raise ValueError(
            f"{held.get('date')}: {len(folds)} folds recorded for "
            f"{n_lines} lines; leave-one-out has one fold per line")
    worst_index = {}
    for model in ("comb", "scale"):
        errors = [float(f["error_hz"][model]) for f in folds]
        if any(not math.isfinite(e) or e < 0 for e in errors):
            raise ValueError(f"{held.get('date')}: a {model} fold error is "
                             "negative or not finite")
        mean_hz = float((loo.get(model) or {})["mean_error_hz"])
        if mean_hz <= 0:
            raise ValueError(f"{held.get('date')}: the {model} fold mean "
                             "must be positive")
        printed = held[f"loo_{model}_mhz"]
        if abs(mean_hz / 1e6 - printed) > 1e-6:
            raise ValueError(
                f"{held.get('date')}: the {model} fold mean "
                f"{mean_hz / 1e6} MHz disagrees with the printed "
                f"loo_{model}_mhz {printed}; the two files must describe "
                "the same folds")
        n = len(errors)
        worst = max(range(n), key=lambda i: errors[i])
        total = sum(errors)
        worst_index[model] = folds[worst]["index"]
        out["models"][model] = {
            "errors_mhz": [e / 1e6 for e in errors],
            "mean_mhz": printed,
            "worst_fold": int(folds[worst]["index"]),
            "worst_line_mhz": float(folds[worst]["actual_hz"]) / 1e6,
            "worst_mhz": errors[worst] / 1e6,
            "worst_over_mean": errors[worst] / mean_hz,
            "worst_share_of_total": errors[worst] / total,
            "worst_share_printed": float(f"{errors[worst] / total:.2f}"),
            "share_max_p": share_max_tail(
                float(f"{errors[worst] / total:.2f}"), n),
            "even_share": 1.0 / n,
            "mean_without_worst_over_mean":
                (total - errors[worst]) / ((n - 1) * mean_hz),
        }
    out["same_worst_line"] = (worst_index["comb"]
                              == worst_index["scale"])
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


def pair_residuals(freqs, step, q, pairs, rms) -> list:
    """The measured adjacent-pair ratio beside both fits, pair by pair.

    ``f_{i+1} / f_i`` is the ratio the lines themselves carry, so scoring both
    predictions against it puts the divergence on the lines instead of between
    two sets of fitted values. Both q and the step were fitted on these same
    lines and the pairs are not independent, so neither residual is a selection
    statistic: the campaign's AIC and leave-one-out remain the decision. Kept
    diagnostic only, and recorded as relative misses so pairs at different
    frequencies can be read together.

    Each relative miss is also carried as the displacement of the pair's upper
    line, in channels, and divided by that fit's own scatter - the comb's rms
    for the additive prediction, the scale's rms for the scale one - so the
    size of the miss is stated in the units the campaign judged the fits by.
    The separation between the two predictions is carried the same way.
    """
    out = []
    for i in range(len(freqs) - 1):
        measured = freqs[i + 1] / freqs[i]
        additive = 1 + step / freqs[i]
        rel_scale = measured / q - 1
        rel_additive = measured / additive - 1
        high_add = additive * freqs[i]      # = freqs[i] + step, exactly
        high_scale = q * freqs[i]
        miss_add = line_offset(rel_additive, high_add, rms)
        miss_scale = line_offset(rel_scale, high_scale, rms)
        gap_channels = (high_add - high_scale) * 1e6 / CHANNEL_HZ
        out.append({
            "pair": i + 1,
            "line_low_mhz": float(freqs[i]),
            "line_high_mhz": float(freqs[i + 1]),
            "measured": float(measured),
            "measured_printed": (f"{float(freqs[i + 1]):.6f} / "
                                 f"{float(freqs[i]):.6f}"),
            "scale_prediction": float(q),
            "additive_prediction": float(additive),
            "predicted_high_mhz": {"additive": float(high_add),
                                   "scale": float(high_scale)},
            "residual_vs_scale": float(measured - q),
            "residual_vs_additive": float(measured - additive),
            "relative_vs_scale": float(rel_scale),
            "relative_vs_additive": float(rel_additive),
            "offset_mhz_vs_additive": miss_add["offset_mhz"],
            "offset_channels_vs_additive": miss_add["offset_channels"],
            "offset_mhz_vs_scale": miss_scale["offset_mhz"],
            "offset_channels_vs_scale": miss_scale["offset_channels"],
            "over_comb_rms": miss_add["over_comb_rms"],
            "over_scale_rms": miss_scale["over_scale_rms"],
            "fit_gap_channels": float(gap_channels),
            "fit_gap_over_comb_rms": float(gap_channels / rms["comb"]),
            "fit_gap_over_scale_rms": float(gap_channels / rms["scale"]),
            "closer": ("additive" if abs(rel_additive) < abs(rel_scale)
                       else "scale"),
            "diagnostic_only": True,
            "note": "measured f_high/f_low beside both fits; q and the step "
                    "were fitted on these same lines, so neither residual "
                    "is independent and neither selects a model; the "
                    "channel offsets use each fit's own rms as the unit",
        })
    return out


def additive_span(freqs, step, q, pairs, rms) -> dict:
    """The pair ratios compounded over the group, and where they meet ``q``.

    ``prod(1 + step / f_i)`` is the additive model's own ratio from the lowest
    line to the highest; ``q ** (n - 1)`` is the scale model's for the same
    span. Both are predictions, so the comparison is descriptive. The crossing
    frequency ``step / (q - 1)`` is where a pair ratio equals q - the single
    line at which the two models coincide, and hence the pivot the gap is
    measured around.

    The measured span is scored against both predictions, and each residual is
    also carried as a displacement of the group's highest line, in channels,
    over that fit's own rms - the same conversion the pair table uses, applied
    to the compounded ratio.
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
    high_add = product * freqs[0]
    high_scale = scale_cumulative * freqs[0]
    miss_add = line_offset(rel_additive, high_add, rms)
    miss_scale = line_offset(rel_scale, high_scale, rms)
    gap_channels = (high_add - high_scale) * 1e6 / CHANNEL_HZ
    item["predicted_high_mhz"] = {"additive": float(high_add),
                                  "scale": float(high_scale)}
    item["offset_mhz_vs_additive"] = miss_add["offset_mhz"]
    item["offset_channels_vs_additive"] = miss_add["offset_channels"]
    item["offset_mhz_vs_scale"] = miss_scale["offset_mhz"]
    item["offset_channels_vs_scale"] = miss_scale["offset_channels"]
    item["over_comb_rms"] = miss_add["over_comb_rms"]
    item["over_scale_rms"] = miss_scale["over_scale_rms"]
    item["fit_gap_channels"] = float(gap_channels)
    item["fit_gap_over_comb_rms"] = float(gap_channels / rms["comb"])
    item["fit_gap_over_scale_rms"] = float(gap_channels / rms["scale"])
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
    fold_runs = group_folds()
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
        rms = fit_noise(row)
        held = held_out_noise(row, rms)
        held["folds"] = fold_spread(held, fold_runs.get(key), key[2])
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
        span = additive_span(freqs, step, q, pairs, rms)
        scored = pair_residuals(freqs, step, q, pairs, rms)
        out.append({
            "source": row["source"],
            "date": row["date"],
            "verdict": row["verdict"],
            "diagnostic_only": row["verdict"] == "MASK-MANUFACTURED SIGNAL",
            "n_lines": int(row["n_lines"]),
            "scale_q": row["scale_q"],
            "lowest_line_mhz": float(f_low),
            "noise": {
                "rms_comb_channels": float(rms["comb"]),
                "rms_scale_channels": float(rms["scale"]),
                "as_printed": [row["rms_comb_channels"],
                               row["rms_scale_channels"]],
                "channel_khz": CHANNEL_HZ / 1e3,
                "note": "rms of each fit's own residuals on these same "
                        "lines, in channels: the unit the pair and span "
                        "misses are divided by, never a weight",
            },
            "held_out": held,
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
                    "max_abs_over_comb_rms": max(
                        abs(r["over_comb_rms"])
                        for r in row["pair_residuals"]),
                    "max_abs_over_scale_rms": max(
                        abs(r["over_scale_rms"])
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

    # --- the misses in the units each fit was judged by ---------------------
    pairs_grouped = [(row, r) for row in campaign for r in row["pair_residuals"]]

    def _largest(items, key) -> dict:
        row, r = max(items, key=lambda it: abs(it[1][key]))
        return {"date": row["date"], "n_lines": row["n_lines"],
                "pair": r["pair"], "line_high_mhz": r["line_high_mhz"],
                "value": abs(r[key]),
                "relative_vs_additive": r["relative_vs_additive"],
                "relative_vs_scale": r["relative_vs_scale"],
                "offset_channels": abs(
                    r["offset_channels_vs_additive" if key == "over_comb_rms"
                      else "offset_channels_vs_scale"])}

    pair_comb = [abs(r["over_comb_rms"]) for _, r in pairs_grouped]
    pair_scale = [abs(r["over_scale_rms"]) for _, r in pairs_grouped]
    span_comb = [abs(r["over_comb_rms"]) for r in spans]
    span_scale = [abs(r["over_scale_rms"]) for r in spans]
    # --- the same misses against the line each fit never saw -----------------
    def _largest_loo(items, model, key=None):
        """The widest miss in held-out widths, skipping groups with no fold."""
        best = None
        channel_key = key or ("offset_channels_vs_"
                              + ("additive" if model == "comb" else "scale"))
        for row, item in items:
            loo = row["held_out"][f"loo_{model}_channels"]
            if loo is None or channel_key not in item:
                continue
            value = abs(item[channel_key]) / loo
            if best is None or value > best["value"]:
                best = {"date": row["date"], "n_lines": row["n_lines"],
                        "pair": item.get("pair"), "value": value}
        return best

    def _largest_fold(rows, model):
        best = max(rows, key=lambda r: r[model]["worst_over_mean"])
        return {"date": best["date"], "n_lines": best["n_lines"],
                "value": best[model]["worst_over_mean"],
                "worst_line_mhz": best[model]["worst_line_mhz"]}

    def _smallest_share_p(rows, model):
        """The least likely worst-fold share under the exchangeable null."""
        best = min(rows, key=lambda r: r[model]["share_max_p"])
        return {"date": best["date"], "n_lines": best["n_lines"],
                "value": best[model]["share_max_p"],
                "share": best[model]["worst_share_of_total"]}

    held_groups = []
    for row in campaign:
        held = row["held_out"]
        closer = _closer_counts(row["pair_residuals"])
        nearer = ("additive" if closer["additive"] > closer["scale"]
                  else "scale" if closer["scale"] > closer["additive"]
                  else "tie")
        vote = {"comb": "additive", "scale": "scale"}.get(
            held["loo_majority"])
        comparable = held["applicable"] and vote is not None
        held_groups.append({
            "date": row["date"],
            "n_lines": row["n_lines"],
            "rms_comb_channels": row["noise"]["rms_comb_channels"],
            "rms_scale_channels": row["noise"]["rms_scale_channels"],
            "loo_comb_channels": held["loo_comb_channels"],
            "loo_scale_channels": held["loo_scale_channels"],
            "held_out_over_rms": dict(held["held_out_over_rms"]),
            "sigma_over_rms": dict(held["sigma_over_rms"]),
            "pairs_closer_counts": closer,
            "pairs_nearer": nearer,
            "loo_majority": held["loo_majority"],
            "agrees": (nearer == vote) if comparable else None,
        })
    fold_rows = []
    for row in campaign:
        spread = row["held_out"]["folds"]
        fold_rows.append({
            "date": row["date"],
            "n_lines": row["n_lines"],
            "n_folds": spread["n_folds"],
            "same_worst_line": spread["same_worst_line"],
            "comb": spread["models"].get("comb"),
            "scale": spread["models"].get("scale"),
        })
    with_folds = [f for f in fold_rows if f["comb"] is not None]
    fold_null = campaign[0]["held_out"]["folds"]["share_null"]
    defined_groups = [g for g in held_groups
                      if g["sigma_over_rms"]["comb"] is not None]
    scored_groups = [g for g in held_groups if g["agrees"] is not None]
    combo = [g["sigma_over_rms"]["comb"] for g in defined_groups]
    scalo = [g["sigma_over_rms"]["scale"] for g in defined_groups]
    spans_grouped = [(row, row["additive_span"]) for row in campaign]

    line_fit["noise_scale"] = {
        "rms_comb_channels": [row["noise"]["rms_comb_channels"]
                              for row in campaign],
        "rms_scale_channels": [row["noise"]["rms_scale_channels"]
                               for row in campaign],
        "pairs": len(pairs_grouped),
        "median_abs_over_comb_rms": _median(pair_comb),
        "median_abs_over_scale_rms": _median(pair_scale),
        "max_abs_over_comb_rms": max(pair_comb),
        "max_abs_over_scale_rms": max(pair_scale),
        "pairs_within_comb_rms": sum(1 for v in pair_comb if v <= 1.0),
        "pairs_within_scale_rms": sum(1 for v in pair_scale if v <= 1.0),
        "largest_over_comb_rms": _largest(pairs_grouped, "over_comb_rms"),
        "largest_over_scale_rms": _largest(pairs_grouped, "over_scale_rms"),
        "max_abs_fit_gap_over_comb_rms": max(
            abs(r["fit_gap_over_comb_rms"]) for _, r in pairs_grouped),
        "max_abs_fit_gap_over_scale_rms": max(
            abs(r["fit_gap_over_scale_rms"]) for _, r in pairs_grouped),
        "max_abs_span_over_comb_rms": max(span_comb),
        "max_abs_span_over_scale_rms": max(span_scale),
        "max_abs_span_fit_gap_over_comb_rms": max(
            abs(r["fit_gap_over_comb_rms"]) for r in spans),
        "max_abs_span_fit_gap_over_scale_rms": max(
            abs(r["fit_gap_over_scale_rms"]) for r in spans),
        "held_out": {
            "groups": held_groups,
            "groups_defined": len(defined_groups),
            "groups_without_a_fold": [
                f"{g['date']} n={g['n_lines']}" for g in held_groups
                if g["sigma_over_rms"]["comb"] is None],
            "min_sigma_over_rms_comb": min(combo) if combo else None,
            "max_sigma_over_rms_comb": max(combo) if combo else None,
            "min_sigma_over_rms_scale": min(scalo) if scalo else None,
            "max_sigma_over_rms_scale": max(scalo) if scalo else None,
            "max_abs_pair_over_comb_loo": _largest_loo(pairs_grouped, "comb"),
            "max_abs_pair_over_scale_loo": _largest_loo(pairs_grouped,
                                                         "scale"),
            "max_abs_span_over_comb_loo": _largest_loo(spans_grouped, "comb"),
            "max_abs_span_over_scale_loo": _largest_loo(spans_grouped,
                                                        "scale"),
            "max_abs_fit_gap_over_comb_loo": _largest_loo(
                pairs_grouped, "comb", "fit_gap_channels"),
            "max_abs_fit_gap_over_scale_loo": _largest_loo(
                pairs_grouped, "scale", "fit_gap_channels"),
            "agreement": {
                "comparable": len(scored_groups),
                "agree": sum(1 for g in scored_groups if g["agrees"]),
                "not_scored": [
                    f"{g['date']} n={g['n_lines']} ({g['loo_majority']})"
                    for g in held_groups if g["agrees"] is None],
            },
            "fold_spread": {
                "groups": fold_rows,
                "groups_with_folds": len(with_folds),
                "largest_worst_over_mean": {
                    model: _largest_fold(with_folds, model)
                    for model in ("comb", "scale")},
                "max_worst_share_of_total": {
                    model: max(f[model]["worst_share_of_total"]
                               for f in with_folds)
                    for model in ("comb", "scale")},
                "share_null": dict(fold_null),
                "smallest_share_p": {
                    model: _smallest_share_p(with_folds, model)
                    for model in ("comb", "scale")},
                "same_worst_line_groups": [
                    f"{f['date']} n={f['n_lines']}" for f in with_folds
                    if f["same_worst_line"]],
                "different_worst_line_groups": [
                    f"{f['date']} n={f['n_lines']}" for f in with_folds
                    if not f["same_worst_line"]],
                "diagnostic_only": True,
                "note": "each fold's error read from "
                        "results/expanded_phi_groups.json and set beside "
                        "the mean printed in results/comb_vs_scale.csv, "
                        "which is only checked against it and never "
                        "replaced; the worst fold is quoted as its share of "
                        "the group total and carries the exchangeable-share "
                        "null, and nothing here enters the AIC, BIC or "
                        "verdict counts or the campaign's own nulls",
            },
            "gaussian_mae": GAUSSIAN_MAE,
            "diagnostic_only": True,
            "note": "each group's held-out mean error set beside the rms "
                    "that group was judged by; the held-out columns are read "
                    "as printed and never recomputed, and nothing here "
                    "replaces them or enters the AIC, BIC or verdict counts "
                    "or the campaign's own nulls",
        },
        "diagnostic_only": True,
        "note": "the pair and span misses divided by each fit's own rms "
                "residual in channels, one rms per group estimated from the "
                "same lines, so no chi-square or null is formed and none of "
                "it enters the model-selection counts",
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
             "1 + d/f_i | measured/(1+d/f_i) - 1 | measured/q - 1 | "
             "miss / comb rms | miss / scale rms | closer | "
             "additive depth | additive eff. depth |")
    L.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | "
             "--- | ---: | ---: |")
    for row in camp:
        for p, m in zip(row["additive_pairs"], row["pair_residuals"]):
            L.append(f"| {row['date']} n={row['n_lines']} | {p['pair']} | "
                     f"{p['line_low_mhz']:.4f} | {m['measured']:.6f} | "
                     f"{p['value']:.6f} | {m['relative_vs_additive']:+.2e} | "
                     f"{m['relative_vs_scale']:+.2e} | "
                     f"{abs(m['over_comb_rms']):.2f} | "
                     f"{abs(m['over_scale_rms']):.2f} | {m['closer']} | "
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
    ns = lf["noise_scale"]
    wa = ns["largest_over_comb_rms"]
    wsc = ns["largest_over_scale_rms"]
    L.append("A residual in ratio units still does not say whether it is "
             "large, and the campaign table already records the scale each "
             "fit was judged by: `rms_comb_channels` and `rms_scale_channels` "
             "are the rms of each fit's own residuals on these same lines, in "
             f"channels ({min(ns['rms_comb_channels']):.3f}-"
             f"{max(ns['rms_comb_channels']):.2f} for the comb, "
             f"{min(ns['rms_scale_channels']):.3f}-"
             f"{max(ns['rms_scale_channels']):.2f} for the scale). A miss in "
             "ratio at a pair displaces that pair's upper line, so "
             "multiplying by it gives an offset in frequency, and one rms "
             "turns that offset into a fraction of the scatter each fit left "
             "on those lines.\n")
    L.append(f"The largest pair miss in its own fit's scatter is "
             f"{wa['value']:.2f} rms, at {wa['date']} n={wa['n_lines']} "
             f"pair {wa['pair']}: {wa['offset_channels']:.1f} channels, from "
             f"a {abs(wa['relative_vs_additive']):.1e} ratio miss on a line "
             f"at {wa['line_high_mhz']:.2f} MHz. For the scale fit the "
             f"largest is {wsc['value']:.2f} rms "
             f"({wsc['offset_channels']:.1f} channels, "
             f"{abs(wsc['relative_vs_scale']):.1e}) at {wsc['date']} "
             f"n={wsc['n_lines']} pair {wsc['pair']}. The medians are "
             f"{ns['median_abs_over_comb_rms']:.2f} and "
             f"{ns['median_abs_over_scale_rms']:.2f} rms, and "
             f"{ns['pairs_within_comb_rms']} of {ns['pairs']} additive misses "
             f"and {ns['pairs_within_scale_rms']} of {ns['pairs']} scale "
             "misses sit inside one rms.\n")
    L.append(f"The {ad['max_pair_gap_to_scale_q']:.1e} by which the two "
             "predictions differ at their widest, read the same way, is never "
             f"more than {ns['max_abs_fit_gap_over_comb_rms']:.2f} rms of "
             f"the comb fit or {ns['max_abs_fit_gap_over_scale_rms']:.2f} rms "
             "of the scale fit, while the worst individual misses on the "
             f"lines reach {ns['max_abs_over_comb_rms']:.1f} and "
             f"{ns['max_abs_over_scale_rms']:.1f} rms. The two numbers are "
             "not the same size once each is divided by the scatter of the "
             "fit that produced it: the disagreement between the fits stays "
             "near one scatter width and the misses on the lines run to "
             "several, so it is the miss that sets the scale here and the "
             "fit-to-fit gap that is the smaller quantity.\n")
    L.append("That ordering is as far as the rescaling goes. Each rms is a "
             "single number per group, estimated from the same four lines (or "
             "fourteen) the residuals are measured on, so numerator and "
             "denominator share one small sample; a pair miss is one line's "
             "displacement rather than an independent draw; and no "
             "chi-square, p-value or null is formed here. Dividing by the "
             "scatter says how big these misses are next to the noise each "
             "fit left behind, not which fit the lines prefer - that remains "
             "with the AIC and leave-one-out columns, and none of these "
             "ratios is counted with them.\n")

    L.append("### 3b. Compounded across the span\n")
    L.append("Multiplying those pair ratios gives the additive model's own "
             "ratio from the group's lowest line to its highest; `q ** pairs` "
             "is what the scale model predicts over the same span. The depth "
             "reading is applied to both, and the span the lines actually "
             "carry, `f_high / f_low`, sits beside them, so both predictions "
             "are scored on the data rather than only against each other:\n")
    L.append("| group | pairs | measured f_high/f_low | additive compounded | "
             "q^pairs | measured/additive - 1 | measured/q^pairs - 1 | "
             "miss / comb rms | miss / scale rms | closer | "
             "additive depth | additive eff. depth |")
    L.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | "
             "--- | ---: | ---: |")
    for row in camp:
        sp = row["additive_span"]
        L.append(f"| {row['date']} n={row['n_lines']} | {sp['pairs']} | "
                 f"{sp['measured_span_ratio']:.6f} | {sp['value']:.6f} | "
                 f"{sp['reference_scale']['value']:.6f} | "
                 f"{sp['relative_vs_additive']:+.2e} | "
                 f"{sp['relative_vs_scale']:+.2e} | "
                 f"{abs(sp['over_comb_rms']):.2f} | "
                 f"{abs(sp['over_scale_rms']):.2f} | {sp['closer']} | "
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
    span_ns = widest["additive_span"]
    worst_span = max(camp, key=lambda row: abs(
        row["additive_span"]["over_comb_rms"]))["additive_span"]
    L.append("The span misses read in the same units. On the "
             f"{widest['n_lines']}-line group the "
             f"{abs(span_ns['relative_vs_additive']):.1e} additive miss over "
             f"a span ending at {span_ns['predicted_high_mhz']['additive']:.2f} "
             f"MHz is {abs(span_ns['offset_channels_vs_additive']):.1f} "
             f"channels, or {abs(span_ns['over_comb_rms']):.2f} rms of the "
             "comb fit, and its scale miss is "
             f"{abs(span_ns['over_scale_rms']):.2f} rms of the scale fit. "
             "Across the four spans the largest misses are "
             f"{lf['noise_scale']['max_abs_span_over_comb_rms']:.2f} comb "
             f"rms and {lf['noise_scale']['max_abs_span_over_scale_rms']:.2f} "
             "scale rms, while the gap between the two compounded predictions "
             "is at most "
             f"{lf['noise_scale']['max_abs_span_fit_gap_over_comb_rms']:.2f} "
             "and "
             f"{lf['noise_scale']['max_abs_span_fit_gap_over_scale_rms']:.2f} "
             "rms. Compounding does not lift either fit's disagreement with "
             "the lines past "
             f"{max(lf['noise_scale']['max_abs_span_over_comb_rms'], lf['noise_scale']['max_abs_span_over_scale_rms']):.1f} "
             "scatter widths, and the fits still disagree with each other by "
             "less than they disagree with the data. The same caution "
             "applies as in 3a: one rms per group, estimated from these "
             "lines, and no null attached.\n")
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

    L.append("### 3c. One rms per group, checked against the held-out line\n")
    ho = ns["held_out"]
    L.append("Sections 3a and 3b divide every miss by one rms per group, and "
             "that rms is an in-sample number: it is computed from the same "
             "lines the fit was estimated on, so by construction it understates "
             "how far a line sits that the fit never saw. The campaign table "
             "already records that second measurement - `loo_comb_mhz` and "
             "`loo_scale_mhz`, each fit's mean error on the line held out of "
             "the fit, averaged over the folds. Set beside the rms, it says "
             "whether one rms is a fair summary of the line scatter:\n")
    L.append("| group | lines | rms comb (ch) | rms scale (ch) | "
             "held-out comb (ch) | held-out scale (ch) | "
             "held-out sigma / rms comb | held-out sigma / rms scale | "
             "nearer fit (pairs) | held-out majority | same way |")
    L.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | "
             "--- | --- | --- |")
    for g in ho["groups"]:
        cells = []
        for key in ("loo_comb_channels", "loo_scale_channels"):
            v = g[key]
            cells.append("-" if v is None else f"{v:.2f}")
        for model in ("comb", "scale"):
            v = g["sigma_over_rms"][model]
            cells.append("-" if v is None else f"{v:.2f}")
        same = ("-" if g["agrees"] is None
                else ("yes" if g["agrees"] else "no"))
        L.append(f"| {g['date']} n={g['n_lines']} | {g['n_lines']} | "
                 f"{g['rms_comb_channels']:.3f} | "
                 f"{g['rms_scale_channels']:.3f} | "
                 + " | ".join(cells)
                 + f" | {g['pairs_nearer']} "
                 f"({g['pairs_closer_counts']['additive']}-"
                 f"{g['pairs_closer_counts']['scale']}) | "
                 f"{g['loo_majority']} | {same} |")
    L.append("")
    defined = [g for g in ho["groups"]
               if g["sigma_over_rms"]["comb"] is not None]
    big = max(defined, key=lambda g: g["n_lines"])
    small = [g for g in defined if g["n_lines"] < big["n_lines"]]
    small_sigma = [g["sigma_over_rms"][m] for g in small
                   for m in ("comb", "scale")]
    small_text = " and ".join(sorted({str(g["n_lines"]) for g in small}))
    no_fold = [g for g in ho["groups"]
               if g["sigma_over_rms"]["comb"] is None]
    no_fold_text = " and ".join(sorted({str(g["n_lines"]) for g in no_fold}))
    L.append("The columns are not the same kind of average, so the sigma "
             "columns carry the comparison. The held-out number is a mean "
             "absolute error and the rms a root mean square; for Gaussian "
             "scatter the two stand in the ratio sqrt(2/pi) = "
             f"{ho['gaussian_mae']:.3f}, so their raw quotient would read low "
             "by that factor whatever the data did; the table carries only "
             "the divided form, the held-out line in in-sample rms widths. "
             "1.0 means the fit's own "
             "scatter already covered it, and anything above 1.0 is the part "
             "an in-sample rms cannot see.\n")
    L.append(f"The answer differs by group size. On the {big['n_lines']}-line "
             "group the held-out line sits "
             f"{big['sigma_over_rms']['comb']:.2f} comb rms and "
             f"{big['sigma_over_rms']['scale']:.2f} scale rms from its fit, "
             f"{big['sigma_over_rms']['comb'] - 1:.0%} and "
             f"{big['sigma_over_rms']['scale'] - 1:.0%} above the in-sample "
             "reading, so one rms is within a quarter of the held-out scatter "
             f"there. On the {small_text}-line groups it sits "
             f"{min(small_sigma):.2f}-{max(small_sigma):.2f} rms out, close to "
             f"double: {small_text} lines supporting two fitted parameters "
             "leave a residual that is small partly because the fit was bent "
             "to them, and the held-out error is what that bending costs. "
             f"The {no_fold_text}-line group has no fold to read - "
             f"{no_fold_text} lines leave two, and two points fit either "
             "model exactly - so nothing is claimed about it. A single rms is "
             "therefore a fair summary only where the lines outnumber the "
             "parameters; where they do not, it is a lower bound on the "
             "scatter, and the 3a and 3b multiples are understatements by "
             "about that factor.\n")
    fs = ho["fold_spread"]
    folded = [f for f in fs["groups"] if f["comb"] is not None]
    small_f = [f for f in folded if f["n_lines"] < big["n_lines"]]

    def _fold_range(rows, model, key):
        values = [r[model][key] for r in rows]
        return min(values), max(values)

    c_share = _fold_range(folded, "comb", "worst_share_of_total")
    s_share = _fold_range(folded, "scale", "worst_share_of_total")
    c_drop = _fold_range(folded, "comb", "mean_without_worst_over_mean")
    s_drop = _fold_range(folded, "scale", "mean_without_worst_over_mean")
    c_p = _fold_range(folded, "comb", "share_max_p")
    s_p = _fold_range(folded, "scale", "share_max_p")
    L.append("The mean over folds can still hide a single line, and "
             "`results/expanded_phi_groups.json` already carries each fold's "
             "own error, so the question is put to that file rather than to a "
             "new fit. Every fold's error is read as recorded, and the mean of "
             "the folds is checked against the mean printed in "
             "`results/comb_vs_scale.csv` to the last printed decimal; a "
             "disagreement raises rather than being averaged away. The table "
             "gives the fold count and the share of the total held-out error "
             "that the worst fold carries; the multiple of the mean and the "
             "mean without it are arithmetic on those two - n share and "
             "n (1 - share) / (n - 1) - so neither is tabulated. A share of "
             "1/n is an even spread, a share approaching 1 is one line "
             "carrying the whole mean, and the last two share columns give "
             "what such a share is worth under the null below:\n")
    L.append("| group | folds | worst fold comb | worst fold scale | "
             "worst line comb (MHz) | worst line scale (MHz) | "
             "share of total comb | share of total scale | "
             "share null p comb | share null p scale | same line worst |")
    L.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | "
             "---: | ---: | --- |")
    for f in fs["groups"]:
        c, sc = f["comb"], f["scale"]
        head = f"| {f['date']} n={f['n_lines']} | {f['n_folds']} |"
        if c is None or sc is None:
            L.append(head + " - | - | - | - | - | - | - | - | - |")
            continue
        same = "yes" if f["same_worst_line"] else "no"
        L.append(head
                 + f" {c['worst_fold']} | {sc['worst_fold']} |"
                 + f" {c['worst_line_mhz']:.4f} | {sc['worst_line_mhz']:.4f} |"
                 + f" {c['worst_share_printed']:.2f} | "
                 + f"{sc['worst_share_printed']:.2f} |"
                 + f" {c['share_max_p']:.2f} | {sc['share_max_p']:.2f} | {same} |")
    L.append("")
    L.append(f"Across the {len(folded)} groups with folds the worst line "
             f"carries {c_share[0]:.0%}-{c_share[1]:.0%} of the comb fit's "
             f"total held-out error and {s_share[0]:.0%}-{s_share[1]:.0%} of "
             "the scale fit's, against the even share at each group's own "
             "fold count - "
             + ", ".join(f"{n} folds {1 / n:.0%}" for n in
                         sorted({f["n_folds"] for f in folded}))
             + ". No share reaches half, so no group's mean is one line's "
             "error in disguise; dropping that worst line lowers the mean by "
             f"{1 - c_drop[1]:.0%}-{1 - c_drop[0]:.0%} for the comb fit and "
             f"{1 - s_drop[1]:.0%}-{1 - s_drop[0]:.0%} for the scale fit, "
             "which is n share - 1 over n - 1, the same two columns again. "
             "Size alone does not say whether a share is large, so it gets a "
             "null.\n")
    null = fs["share_null"]
    L.append("If a group's total held-out error could have fallen on any of "
             "its folds with no fold favoured, the total being fixed, the "
             "share vector is "
             f"{null['distribution']}, uniform on the simplex, and the "
             "largest of n such shares is not pinned at 1/n: even an even "
             "scatter of n numbers has a largest member, and the question is "
             "how large the largest gets by chance. Its upper tail is "
             f"{null['tail']}, evaluated in exact fractions for the printed "
             "shares. The tail probabilities run "
             f"{min(c_p[0], s_p[0]):.2f}-{max(c_p[1], s_p[1]):.2f} over the "
             f"{len(folded)} groups and both fits, and the smallest tail in "
             "each fit is "
             f"p = {fs['smallest_share_p']['comb']['value']:.2f} on "
             f"{fs['smallest_share_p']['comb']['date']} "
             f"n={fs['smallest_share_p']['comb']['n_lines']} for the comb "
             f"fit and p = {fs['smallest_share_p']['scale']['value']:.2f} "
             f"on {fs['smallest_share_p']['scale']['date']} "
             f"n={fs['smallest_share_p']['scale']['n_lines']} for the scale "
             "fit: the largest of n exchangeable folds reaches the observed "
             "share in most null draws, so no group's worst fold is unusual "
             "under this reading. The null assumes nothing about either fit, "
             "knows only that the error could have fallen on any fold, and "
             "says nothing about which model is right; it is attached to the "
             "share alone, no threshold is drawn from it, and it stays out "
             "of the model-selection counts.\n")
    small_fold_mult = _fold_range(small_f, "comb", "worst_over_mean")
    small_fold_mult_hi = max(small_fold_mult[1],
                             _fold_range(small_f, "scale",
                                         "worst_over_mean")[1])
    small_fold_mult_lo = min(small_fold_mult[0],
                             _fold_range(small_f, "scale",
                                         "worst_over_mean")[0])
    big_f = [f for f in folded if f["n_lines"] == big["n_lines"]]
    big_fold_mult = [_fold_range(big_f, m, "worst_over_mean")[1]
                     for m in ("comb", "scale")]
    same_text = ", ".join(f"{f['date']} n={f['n_lines']}"
                          for f in folded if f["same_worst_line"]) or "none"
    other_text = ", ".join(f"{f['date']} n={f['n_lines']}"
                           for f in folded
                           if not f["same_worst_line"]) or "none"
    L.append("What survives is a reading about fold count, not outliers: the "
             "largest of n folds takes a growing multiple of an even share as "
             "n rises even when every fold comes from the same scatter, which "
             f"is why the {big['n_lines']}-line group shows "
             f"{min(big_fold_mult):.1f}-{max(big_fold_mult):.1f} times an even "
             f"share against {small_fold_mult_lo:.1f}-{small_fold_mult_hi:.1f} "
             f"on the {small_text}-line groups without its folds being any "
             "less even. Read as a proportion, what 3c reports on the "
             f"{small_text}-line groups is the fold count and the two fitted "
             "parameters, not one droppable line. The two fits are worst at "
             f"the same line on {same_text} and at different lines on "
             f"{other_text}: where they agree the bad fold belongs to that "
             "line, and where they part it belongs to that fit's "
             "extrapolation. Fold indices in the table are the ones recorded "
             "in the file, counted from zero. No fold is removed from "
             "anything - the fold file supplies the errors and the printed "
             "means stay the campaign's own.\n")
    wpair = ho["max_abs_pair_over_comb_loo"]
    L.append("Carried through, the largest pair miss is "
             f"{wpair['value']:.2f} held-out comb widths "
             f"({wpair['date']} n={wpair['n_lines']} pair {wpair['pair']}) "
             f"where the rms reading was {ns['max_abs_over_comb_rms']:.2f}, "
             "and the largest span miss is "
             f"{ho['max_abs_span_over_comb_loo']['value']:.2f} where it was "
             f"{ns['max_abs_span_over_comb_rms']:.2f}. On the scale side the "
             "same two move from "
             f"{ns['max_abs_over_scale_rms']:.2f} and "
             f"{ns['max_abs_span_over_scale_rms']:.2f} rms to "
             f"{ho['max_abs_pair_over_scale_loo']['value']:.2f} and "
             f"{ho['max_abs_span_over_scale_loo']['value']:.2f} held-out "
             "widths. The gap between the two fits, read the same way, reaches "
             f"{ho['max_abs_fit_gap_over_comb_loo']['value']:.2f} comb and "
             f"{ho['max_abs_fit_gap_over_scale_loo']['value']:.2f} scale "
             "held-out widths. The rescaling changes the size of the multiple "
             "and not the ordering: the misses on the lines still run wider "
             "than the gap between the two fits.\n")
    ag = ho["agreement"]
    scored_groups = [g for g in ho["groups"] if g["agrees"] is not None]
    agree_text = ", ".join(f"{g['date']} n={g['n_lines']}"
                           for g in scored_groups if g["agrees"]) or "none"
    differ_text = ", ".join(
        f"{g['date']} n={g['n_lines']} ({g['pairs_nearer']} pairs against "
        f"{g['loo_majority']} folds)"
        for g in scored_groups if not g["agrees"]) or "none"
    if scored_groups and ag["agree"] == len(scored_groups):
        compare = f"it matches on {agree_text} and differs nowhere"
    elif scored_groups and ag["agree"] == 0:
        compare = f"it differs on {differ_text} and matches nowhere"
    else:
        compare = f"it matches on {agree_text} and differs on {differ_text}"
    L.append(f"The nearer-fit split from 3a is the other thing the held-out "
             "columns can check. The pairs put one model nearer in "
             f"{sum(g['pairs_closer_counts']['additive'] for g in ho['groups'])} "
             "additive and "
             f"{sum(g['pairs_closer_counts']['scale'] for g in ho['groups'])} "
             "scale readings; the fold vote names a model in "
             f"{ag['comparable']} of {len(ho['groups'])} groups, and "
             f"{compare}. The {big['n_lines']}-line group votes "
             f"{big['loo_majority']} and the {no_fold_text}-line group cannot "
             "be held out, so neither is scored either way and both are shown "
             "above. Where the two disagree, the pair split is still a set of "
             "correlated in-sample nearnesses and the fold vote is an "
             "out-of-sample count on the same lines; reporting both is the "
             "check, and neither can be used to shore up the other.\n")
    L.append("Nothing here is a selection statistic, and nothing here makes "
             "the rms wrong - it dates it. The held-out columns are read from "
             "`results/comb_vs_scale.csv` as printed and never recomputed, and "
             "the per-fold errors that sit beside them are read from "
             "`results/expanded_phi_groups.json`; the "
             "rms on the other side is the same in-sample number 3a and 3b "
             "used; no AIC, BIC, null or verdict is touched, and the campaign's "
             "counts are unchanged. What the comparison adds is the size of "
             "the correction: a miss quoted as three rms is about "
             f"{3.0 / max(small_sigma):.1f}-{3.0 / min(small_sigma):.1f} "
             f"held-out widths on the {small_text}-line "
             f"groups and about {3.0 / big['sigma_over_rms']['comb']:.1f} on "
             f"the {big['n_lines']}-line one, so the multiple should be quoted "
             "with "
             "the group size beside it. The held-out error is "
             "itself a mean over n folds of the same small group, with the "
             "pairs correlated and no null attached, so neither denominator "
             "is a likelihood. And none of it moves a ratio towards phi: a "
             "smaller scatter makes a measurement sharper, not more golden, "
             "and phi stays a coordinate.\n")
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
    ns = lf["noise_scale"]
    print(f"noise-scaled misses      max pair "
          f"{ns['max_abs_over_comb_rms']:.2f} comb / "
          f"{ns['max_abs_over_scale_rms']:.2f} scale rms; fit gap "
          f"{ns['max_abs_fit_gap_over_comb_rms']:.2f} / "
          f"{ns['max_abs_fit_gap_over_scale_rms']:.2f} rms; spans "
          f"{ns['max_abs_span_over_comb_rms']:.2f} / "
          f"{ns['max_abs_span_over_scale_rms']:.2f} rms")
    ho = ns["held_out"]
    if ho["groups_defined"]:
        print(f"held-out vs rms          sigma "
              f"{ho['min_sigma_over_rms_comb']:.2f}-"
              f"{ho['max_sigma_over_rms_comb']:.2f} comb, "
              f"{ho['min_sigma_over_rms_scale']:.2f}-"
              f"{ho['max_sigma_over_rms_scale']:.2f} scale; pair split vs "
              f"fold vote {ho['agreement']['agree']}/"
              f"{ho['agreement']['comparable']}")
    if ho["fold_spread"]["groups_with_folds"]:
        fs = ho["fold_spread"]
        print(f"worst held-out fold      comb "
              f"{fs['largest_worst_over_mean']['comb']['value']:.2f}x mean, "
              f"scale {fs['largest_worst_over_mean']['scale']['value']:.2f}x; "
              f"max share of a group's total "
              f"{fs['max_worst_share_of_total']['comb']:.2f} comb / "
              f"{fs['max_worst_share_of_total']['scale']:.2f} scale, "
              "smallest share tail p "
              f"{min(fs['smallest_share_p']['comb']['value'], fs['smallest_share_p']['scale']['value']):.2f}")
    print(f"model crossings (MHz)    "
          f"{[round(c['crossing_mhz'], 3) for c in cs['additive']['crossings']]}")
    print(f"vacuous recursion        {cs['vacuous']}")
    print(f"diagnostic-only rows     {len(cs['diagnostic_only'])}")
    print(f"depth sweep              {[(r['jitter'], r['depth']) for r in payload['depth_sweep']]}")
    print("wrote results/phi_plus_one.json and PHI_PLUS_ONE.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
