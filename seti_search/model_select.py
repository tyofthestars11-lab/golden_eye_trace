"""Comb against scale family: fit both, then say when the data cannot decide.

Two ways to make a group of lines that looks regular.

    comb      f_n = f_0 + n * df        a fixed frequency step: a clock, a
                                        channel grid, an intermodulation product
    scale     f_n = f_0 * q**n          one thing measured at several scales

Both have two free parameters, so on the same data the same number of lines
cannot prefer either by parameter count alone; the comparison has to be about
residuals, and about what the residuals would look like if the other model were
true. Everything here is built to make that comparison honest:

  * both models are fitted to the *same* quantity - frequency in Hz - so their
    residuals are on one scale and AIC can be compared at all
  * the scale fit is profiled: for each candidate q the best f_0 is closed form,
    so there is one parameter to search and no optimiser to trust
  * held-out prediction is leave-one-out, refitting from scratch each fold, and
    reports the error in Hz and in phi rungs separately
  * a sign-flip bootstrap resamples the better model's residuals to say whether
    the AIC difference survives noise
  * `synthetic_calibration` builds groups from a *known* generator and reports how
    often this pipeline recovers the truth at that line count and span. That
    number, not the fit, is what licenses a verdict.

The verdict is deliberately hard to earn. A four-line chain inside one phi rung
returns INDETERMINATE DUE TO LIMITED SPAN however clean its residuals look,
because a narrow span is exactly where the two models coincide.
"""

from __future__ import annotations

import math

import numpy as np

from .rungs import PHI, rung

__all__ = [
    "BAND_SPAN_REQUIRED",
    "DECISIVE_AIC",
    "NOISE_FRACTION",
    "compare",
    "fit_comb",
    "fit_scale",
    "LABELS",
    "structure_controls",
    "synthetic_calibration",
    "verdict",
]

DECISIVE_AIC = 10.0        # a model has to win by this much to be named
BAND_SPAN_REQUIRED = 1.0   # phi rungs a group must span before a name is allowed
NOISE_FRACTION = 0.5       # residual RMS above this many channels = neither fits
BASES = {"phi": PHI, "2": 2.0, "e": math.e, "10": 10.0, "sqrt2": math.sqrt(2.0)}


def fit_comb(freqs) -> dict:
    """Least squares f_n = f0 + n*df, coefficients exact."""
    f = np.asarray(freqs, dtype=float)
    n = np.arange(f.size, dtype=float)
    if f.size < 3:
        raise ValueError("a comb needs at least three lines")
    df, f0 = np.polyfit(n, f, 1)
    model = f0 + n * df
    resid = f - model
    return {"model": "comb", "f0": float(f0), "df": float(df),
            "residuals_hz": resid, "rss": float(resid @ resid)}


def _profile(f: np.ndarray, n: np.ndarray, a_grid: np.ndarray):
    """For each ln q on the grid, the best f0 in closed form and its RSS.

    Vectorised over the grid: with f0 profiled out there is one parameter left, so
    a dense grid plus one refinement pass is an exhaustive search of the only
    nonlinear parameter, and no optimiser has to be trusted.
    """
    w = np.exp(np.outer(n, a_grid))            # (n_lines, grid)
    denom = (w * w).sum(axis=0)
    f0 = (w.T @ f) / denom
    resid = f[:, None] - f0[None, :] * w
    return (resid * resid).sum(axis=0), f0


def fit_scale(freqs) -> dict:
    """Least squares f_n = f0*q**n, profiled over log q."""
    f = np.asarray(freqs, dtype=float)
    n = np.arange(f.size, dtype=float)
    if f.size < 3:
        raise ValueError("a scale family needs at least three lines")
    if np.any(f <= 0):
        raise ValueError("frequencies must be positive")

    per_step = float(np.log(f[-1] / f[0]) / max(f.size - 1, 1))
    grid = np.linspace(per_step * 0.2, per_step * 5.0, 4001)
    rss, f0s = _profile(f, n, grid)
    k = int(np.argmin(rss))
    a, f0, rss_best = float(grid[k]), float(f0s[k]), float(rss[k])
    step = grid[1] - grid[0]
    # three refinements: each pass shrinks the bracket by 200x, which takes the
    # residual from grid-limited to floating-point limited
    for _ in range(3):
        fine = np.linspace(a - step, a + step, 401)
        rss_f, f0_f = _profile(f, n, fine)
        j = int(np.argmin(rss_f))
        a, f0, rss_best = float(fine[j]), float(f0_f[j]), float(rss_f[j])
        step /= 200.0
    model = f0 * np.exp(a * n)
    return {"model": "scale", "f0": f0, "q": float(math.exp(a)),
            "residuals_hz": f - model, "rss": rss_best}


def _information(fit: dict, n_lines: int) -> dict:
    """AIC and BIC for a two-parameter model with a Gaussian residual scale."""
    n = n_lines
    rss = max(fit["rss"], 1e-30)
    sigma2 = rss / n
    ln_like = -0.5 * n * (math.log(2.0 * math.pi * sigma2) + 1.0)
    k = 2
    return {"aic": 2 * k - 2 * ln_like, "bic": k * math.log(n) - 2 * ln_like,
            "rms_hz": math.sqrt(sigma2), "ln_like": ln_like}


def _predict(fit: dict, freqs, index: int) -> float:
    n = float(index)
    if fit["model"] == "comb":
        return fit["f0"] + n * fit["df"]
    return fit["f0"] * fit["q"] ** n


def leave_one_out(freqs) -> dict:
    """Refit on N-1 lines and predict the omitted one, for every omission."""
    f = np.asarray(freqs, dtype=float)
    if f.size < 4:
        # holding out one of three lines leaves two, and two points fit either
        # model exactly: leave-one-out is not defined, so it is not reported
        return {"applicable": False,
                "why": "three lines cannot be held out one at a time; two points "
                       "fit a comb and a scale family exactly",
                "majority": "undefined", "comb_wins": 0, "scale_wins": 0,
                "comb": {"mean_error_hz": None, "max_error_hz": None,
                         "mean_error_rung": None},
                "scale": {"mean_error_hz": None, "max_error_hz": None,
                          "mean_error_rung": None},
                "folds": []}
    folds = []
    for i in range(f.size):
        keep = np.delete(f, i)
        n_keep = np.delete(np.arange(f.size, dtype=float), i)
        preds = {}
        for name, fn in (("comb", fit_comb), ("scale", fit_scale)):
            if keep.size < 3:
                preds[name] = float("nan")
                continue
            # index the fit on the surviving positions, then extrapolate to i
            shifted = _fit_at(keep, n_keep)
            preds[name] = _predict(shifted[name], None, float(i))
        err = {name: abs(preds[name] - f[i]) for name in preds}
        err_rung = {name: abs(math.log(max(abs(preds[name]), 1e-9) / f[i]) / math.log(PHI))
                    for name in preds}
        folds.append({"index": i, "actual_hz": float(f[i]),
                      "predicted": {k: float(v) for k, v in preds.items()},
                      "error_hz": {k: float(v) for k, v in err.items()},
                      "error_rung": {k: float(v) for k, v in err_rung.items()},
                      "winner": min(err, key=err.get)})
    out = {"folds": folds}
    for name in ("comb", "scale"):
        e = [fo["error_hz"][name] for fo in folds]
        r = [fo["error_rung"][name] for fo in folds]
        out[name] = {"mean_error_hz": float(np.mean(e)),
                     "max_error_hz": float(np.max(e)),
                     "mean_error_rung": float(np.mean(r))}
    wins = [fo["winner"] for fo in folds]
    out["comb_wins"] = int(sum(1 for w in wins if w == "comb"))
    out["scale_wins"] = int(sum(1 for w in wins if w == "scale"))
    out["majority"] = ("comb" if out["comb_wins"] > out["scale_wins"]
                       else "scale" if out["scale_wins"] > out["comb_wins"] else "tie")
    return out


def _fit_at(values, positions) -> dict:
    """Fit both models on lines at the given (non-consecutive) indices."""
    f = np.asarray(values, dtype=float)
    n = np.asarray(positions, dtype=float)
    df, f0 = np.polyfit(n, f, 1)
    comb = {"model": "comb", "f0": float(f0), "df": float(df)}

    per_step = float(np.log(f[-1] / f[0]) / max(n[-1] - n[0], 1.0))
    grid = np.linspace(per_step * 0.2, per_step * 5.0, 2001)
    rss, f0s = _profile(f, n, grid)
    k = int(np.argmin(rss))
    a, f0 = float(grid[k]), float(f0s[k])
    step = grid[1] - grid[0]
    for _ in range(3):
        fine = np.linspace(a - step, a + step, 201)
        rss_f, f0_f = _profile(f, n, fine)
        j = int(np.argmin(rss_f))
        a, f0 = float(fine[j]), float(f0_f[j])
        step /= 100.0
    return {"comb": comb, "scale": {"model": "scale", "f0": f0, "q": math.exp(a)}}


def spans(freqs) -> dict:
    """How much of a rung ladder the group covers, in each base."""
    f = np.asarray(freqs, dtype=float)
    low, high = float(f.min()), float(f.max())
    out = {}
    for name, base in BASES.items():
        out[name] = math.log(high / low) / math.log(base)
    out["phi_rungs"] = out["phi"]
    return out


def compare(freqs, channel_hz: float | None = None, bootstrap: int = 2000,
            seed: int = 20261112) -> dict:
    """Fit both models to one group and report every number the verdict needs."""
    f = np.sort(np.asarray(freqs, dtype=float))
    comb = fit_comb(f)
    scale = fit_scale(f)
    out = {"n_lines": int(f.size), "freqs_hz": [float(x) for x in f]}
    for name, fit in (("comb", comb), ("scale", scale)):
        info = _information(fit, f.size)
        entry = {k: v for k, v in fit.items() if k != "residuals_hz"}
        entry["residuals_hz"] = [float(x) for x in fit["residuals_hz"]]
        entry.update(info)
        out[name] = entry
    out["delta_aic"] = out["comb"]["aic"] - out["scale"]["aic"]  # >0 favours scale
    out["delta_bic"] = out["comb"]["bic"] - out["scale"]["bic"]
    out["loo"] = leave_one_out(f)
    out["span"] = spans(f)
    out["span_rungs"] = out["span"]["phi"]
    if channel_hz:
        out["rms_in_channels"] = {
            "comb": out["comb"]["rms_hz"] / channel_hz,
            "scale": out["scale"]["rms_hz"] / channel_hz,
        }

    # sign-flip bootstrap of the better model's residuals: does the AIC gap
    # survive the noise the fit itself leaves behind?
    better = "scale" if out["delta_aic"] > 0 else "comb"
    base_fit = scale if better == "scale" else comb
    rng = np.random.default_rng(seed)
    n = np.arange(f.size, dtype=float)
    base_model = (base_fit["f0"] * base_fit["q"] ** n) if better == "scale" \
        else (base_fit["f0"] + n * base_fit["df"])
    resid = base_fit["residuals_hz"]
    deltas = np.empty(bootstrap)
    for b in range(bootstrap):
        draw = base_model + rng.choice([-1.0, 1.0], size=f.size) * resid
        c = fit_comb(draw)
        s = fit_scale(draw)
        deltas[b] = _information(c, f.size)["aic"] - _information(s, f.size)["aic"]
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    out["bootstrap"] = {"draws": int(bootstrap), "resampled_model": better,
                        "delta_aic_ci95": [float(lo), float(hi)],
                        "excludes_zero": bool(lo > 0 or hi < 0)}

    # sensitivity: drop one line at a time, refit, see if the preference flips
    drop = []
    for i in range(f.size if f.size >= 4 else 0):
        keep = np.delete(f, i)
        c = _information(fit_comb(keep), keep.size)
        s = _information(fit_scale(keep), keep.size)
        drop.append({"dropped_hz": float(f[i]), "delta_aic": c["aic"] - s["aic"],
                     "prefers": "scale" if c["aic"] > s["aic"] else "comb"})
    out["drop_one"] = drop
    out["drop_one_stable"] = (len({d["prefers"] for d in drop}) == 1) if drop else None
    return out


def synthetic_calibration(n_lines: int = 4, span_hz: float = 3.0e6,
                          trials: int = 400, seed: int = 7,
                          jitter_hz: float = 2000.0) -> dict:
    """How often can this pipeline name the right model at this size and span?

    Groups are generated from a known model with small jitter, fitted both ways,
    and scored. The result is a confusion matrix: if true combs are called scale
    families (or the reverse) about as often as they are called correctly, then a
    verdict of either name on real data of this size is a coin toss, whatever the
    residuals say.
    """
    rng = np.random.default_rng(seed)
    out = {"n_lines": n_lines, "span_hz": span_hz, "trials": trials,
           "jitter_hz": jitter_hz, "confusion": {}}
    for truth in ("comb", "scale"):
        calls = {"comb": 0, "scale": 0, "tie": 0}
        for _ in range(trials):
            centre = rng.uniform(45e6, 55e6)
            if truth == "comb":
                df = span_hz / (n_lines - 1)
                f = (centre + df * np.arange(n_lines)
                     + rng.uniform(0, jitter_hz) * rng.standard_normal(n_lines))
            else:
                q = (1.0 + span_hz / centre) ** (1.0 / (n_lines - 1))
                f = centre * q ** np.arange(n_lines)
                f = f + rng.uniform(0, jitter_hz) * rng.standard_normal(n_lines)
            c = _information(fit_comb(f), n_lines)
            s = _information(fit_scale(f), n_lines)
            d = c["aic"] - s["aic"]
            if abs(d) <= DECISIVE_AIC:
                calls["tie"] += 1
            else:
                calls["scale" if d > 0 else "comb"] += 1
        out["confusion"][truth] = calls
    out["correct_rate"] = {
        "comb": out["confusion"]["comb"]["comb"] / trials,
        "scale": out["confusion"]["scale"]["scale"] / trials,
    }
    return out


# --- the seven labels -------------------------------------------------------

LABELS = (
    # not one of the seven: a group whose members were selected by an artefact of
    # the processing is disqualified before any model is fitted to it
    "MASK-MANUFACTURED SIGNAL",
    "SCALE FAMILY SUPPORTED",
    "FREQUENCY COMB SUPPORTED",
    "BOTH MODELS FIT",
    "NEITHER MODEL FIT",
    "INDETERMINATE DUE TO LIMITED SPAN",
    "ASTROPHYSICAL STRUCTURE",
    "INSTRUMENTAL STRUCTURE",
)

GRID_TOLERANCE = 0.02   # channels: a spacing this close to an integer is the grid
COHERENT_SPAN = 10.0    # span in units of the median spacing
COHERENT_DRIFT = 0.05   # fractional spread of drift rate allowed for one event


def structure_controls(freqs, drifts=None, channel_hz=None) -> dict:
    """Two physical readings that outrank a model fit when they fire.

    A comb on the channel grid is instrumental whatever its residuals; a set of
    lines that share one drift across a span many times their own spacing is one
    broadband event, not a scale family. Both are measured, both are reported even
    when they do not fire, and `GRID_TOLERANCE` is fixed in code rather than per
    group.
    """
    f = np.sort(np.asarray(freqs, dtype=float))
    steps = np.diff(f)
    out = {"steps_hz": [float(x) for x in steps]}
    if channel_hz:
        in_channels = steps / channel_hz
        out["steps_in_channels"] = [float(x) for x in in_channels]
        off_integer = np.abs(in_channels - np.round(in_channels))
        out["steps_off_integer_channels"] = [float(x) for x in off_integer]
        out["integer_grid_steps"] = int((off_integer <= GRID_TOLERANCE).sum())
        out["half_integer_steps"] = int(
            (np.abs(in_channels - np.round(in_channels + 0.5) - 0.5) <= GRID_TOLERANCE).sum())
        out["grid_comb"] = bool(out["integer_grid_steps"] >= max(2, len(steps) - 1))
        out["half_integer_suspected"] = bool(out["half_integer_steps"] >= 2)
    out["median_step_hz"] = float(np.median(steps))
    out["span_over_spacing"] = float((f[-1] - f[0]) / np.median(steps))
    if drifts is not None:
        d = np.asarray(drifts, dtype=float)
        scale = max(abs(float(np.median(d))), 1.0)
        out["drift_hz_s"] = [float(x) for x in d]
        out["drift_spread"] = float((d.max() - d.min()) / scale)
        out["one_event"] = bool(out["drift_spread"] <= COHERENT_DRIFT
                                and out["span_over_spacing"] >= COHERENT_SPAN)
    else:
        out["one_event"] = False
    return out


def verdict(compared: dict, controls: dict, replicated: bool = False,
            channel_hz: float | None = None) -> dict:
    """One of the seven labels, in the preregistered order of precedence."""
    reasons = []
    rms = compared.get("rms_in_channels") or {}
    if controls.get("grid_comb"):
        reasons.append("spacings land on the channel grid, so the regularity is "
                       "the instrument's")
        return {"label": "INSTRUMENTAL STRUCTURE", "reasons": reasons}
    if controls.get("one_event"):
        reasons.append(
            f"all {compared['n_lines']} lines share one drift to within "
            f"{100 * controls['drift_spread']:.1f}% across "
            f"{controls['span_over_spacing']:.0f} times their own spacing: one "
            "event crossing the band, not several lines at a scale")
        return {"label": "ASTROPHYSICAL STRUCTURE", "reasons": reasons}

    if channel_hz and rms:
        # the better of the two, not the worse: a group that one model describes
        # perfectly and the other badly is not a group neither model fits
        best = min(rms.values())
        if best > NOISE_FRACTION:
            reasons.append(
                f"even the better model leaves {best:.2f} channels of residual, above "
                f"the {NOISE_FRACTION:g}-channel noise floor, and neither the grid nor a "
                "single drifting event explains the group")
            return {"label": "NEITHER MODEL FIT", "reasons": reasons}

    span = compared["span_rungs"]
    delta = compared["delta_aic"]
    loo = compared["loo"].get("majority", "undefined")
    winner = "scale" if delta > 0 else "comb"
    decisive = (abs(delta) > DECISIVE_AIC and loo == winner
                and compared["bootstrap"]["excludes_zero"])
    if span < BAND_SPAN_REQUIRED:
        reasons.append(
            f"the group spans {span:.3f} phi rungs; under one rung the two models "
            "differ by less than the line spacing itself, so no fit can separate "
            "them (delta AIC {:.2f}, bootstrap CI [{:.1f}, {:.1f}])".format(
                delta, *compared["bootstrap"]["delta_aic_ci95"]))
        if decisive:
            reasons.append("the fit would favour " + winner +
                           " if the span were wider, which is not evidence")
        return {"label": "INDETERMINATE DUE TO LIMITED SPAN", "reasons": reasons}

    if decisive:
        if winner == "scale":
            if not replicated:
                reasons.append("the scale family fits better, but a single group "
                               "is one scale interval and the rule needs two")
                return {"label": "BOTH MODELS FIT", "reasons": reasons}
            reasons.append(f"scale family wins by delta AIC {delta:.1f} and "
                           "replicates on a second group")
            return {"label": "SCALE FAMILY SUPPORTED", "reasons": reasons}
        reasons.append(f"comb wins by delta AIC {-delta:.1f} with agreement from "
                       "held-out prediction")
        return {"label": "FREQUENCY COMB SUPPORTED", "reasons": reasons}

    reasons.append(
        f"delta AIC {delta:.2f} is inside the noise (bootstrap CI "
        f"[{compared['bootstrap']['delta_aic_ci95'][0]:.1f}, "
        f"{compared['bootstrap']['delta_aic_ci95'][1]:.1f}]), so both fit")
    return {"label": "BOTH MODELS FIT", "reasons": reasons}
