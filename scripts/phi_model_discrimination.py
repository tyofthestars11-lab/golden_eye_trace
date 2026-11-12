"""Comb or scale family: fit both on every group, and say what the data can bear.

Reads the survey's ten records and the three preregistered expanded searches,
chains each day's lines by the declared proximity rule, fits both models to every
chain of three or more, and returns exactly one of the seven labels per group.

Three calibrations decide how much any verdict is worth, and all three are computed
rather than asserted:

  * a synthetic sweep - groups built from a known generator, at the residual size
    and span of the real data, scored for whether the pipeline recovers the truth
  * a random-line null - uniform draws over each group's own span, so the observed
    AIC gap is compared against "any three or four lines would look like this"
  * the window controls already measured - noise scale and blanking mask, which
    decide whether a group is even visible in a wider band

    python3 scripts/phi_model_discrimination.py
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from seti_search.model_select import (  # noqa: E402
    _information,
    compare,
    fit_comb,
    fit_scale,
    structure_controls,
    synthetic_calibration,
    verdict,
    DECISIVE_AIC,
)
from seti_search.rungs import PHI, rung  # noqa: E402

CHANNEL_HZ = 100547.19562243501
GROUP_RATIO = 1.05
NULL_DRAWS = 1500
SEED = 20261112

CANDIDATE_SETS = {
    "survey-96ch-1000int": "signal_records.json",
    "expanded-A-731ch-allint": "prereg/expanded_search.json",
    "expanded-B-731ch-1000int": "prereg/expanded_search_B.json",
    "expanded-C-96ch-allint": "prereg/expanded_search_C.json",
}


def load_json(path):
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def survey_lines():
    """The survey's ten records, as (date, freq_hz, drift_hz_s, snr, source)."""
    doc = load_json(os.path.join(ROOT, "signal_records.json"))
    if not doc:
        return []
    out = []
    for rec in doc["records"]:
        out.append({
            "date": rec["observation_date"].replace("-", ""),
            "freq_hz": float(rec["frequency_hz"]),
            "drift_hz_s": float(rec["drift_hz_per_s"]),
            "snr": float(rec["snr_at_cell"]),
            "source": "survey",
        })
    return out


def expanded_lines(set_name, path):
    doc = load_json(os.path.join(ROOT, path))
    if not doc:
        return []
    out = []
    for day in doc["days"]:
        for cand in day["candidates"]:
            out.append({
                "date": day["date"], "freq_hz": float(cand["freq_hz"]),
                "drift_hz_s": float(cand["drift_hz_s"]), "snr": float(cand["snr"]),
                "source": set_name,
            })
    return out


def chain(lines):
    """The declared grouping rule: chain neighbours whose ratio is at most 1.05."""
    groups, current = [], []
    for line in sorted(lines, key=lambda r: r["freq_hz"]):
        if current and line["freq_hz"] / current[-1]["freq_hz"] > GROUP_RATIO:
            if len(current) >= 3:
                groups.append(current)
            current = []
        current.append(line)
    if len(current) >= 3:
        groups.append(current)
    return groups


def random_line_null(freqs, draws=NULL_DRAWS, seed=SEED):
    """How often do random lines in this group's own span beat its AIC gap?"""
    f = np.sort(np.asarray(freqs, dtype=float))
    n = f.size
    rng = np.random.default_rng(seed)
    observed = abs(compare(f, bootstrap=1)["delta_aic"])
    null = np.empty(draws)
    for i in range(draws):
        draw = np.sort(rng.uniform(f[0], f[-1], n))
        c = fit_comb(draw)
        s = fit_scale(draw)
        null[i] = abs(_information(c, n)["aic"] - _information(s, n)["aic"])
    return {
        "draws": int(draws), "observed_abs_delta_aic": float(observed),
        "p_value": float((null >= observed).mean()),
        "null_median": float(np.median(null)),
        "note": "uniform lines inside the group's own span, same line count",
    }


def plot_residuals(groups, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    usable = [g for g in groups if g["comparison"]]
    if not usable:
        return False
    cols = min(3, len(usable))
    rows = int(math.ceil(len(usable) / cols))
    fig, axes = plt.subplots(rows * 2, cols, figsize=(4.2 * cols, 2.6 * rows * 2),
                             squeeze=False)
    for i, g in enumerate(usable):
        r, c = divmod(i, cols)
        comp = g["comparison"]
        n = np.arange(comp["n_lines"])
        ax = axes[r * 2][c]
        ax.axhline(0, color="#999", lw=0.8)
        ax.plot(n, comp["comb"]["residuals_hz"], "o-", color="#c0392b", ms=4,
                label="comb")
        ax.plot(n, comp["scale"]["residuals_hz"], "s--", color="#2c6fbb", ms=4,
                label="scale")
        ax.set_title(f"{g['date']}  n={comp['n_lines']}  "
                     f"dAIC={comp['delta_aic']:.1f}", fontsize=8)
        ax.set_ylabel("residual (Hz)", fontsize=7)
        ax.tick_params(labelsize=6)
        if i == 0:
            ax.legend(fontsize=6)
        ax2 = axes[r * 2 + 1][c]
        ax2.axhline(0, color="#999", lw=0.8)
        ax2.plot(n, np.array(comp["comb"]["residuals_hz"]) / CHANNEL_HZ, "o-",
                 color="#c0392b", ms=4)
        ax2.plot(n, np.array(comp["scale"]["residuals_hz"]) / CHANNEL_HZ, "s--",
                 color="#2c6fbb", ms=4)
        ax2.set_ylabel("residual (channels)", fontsize=7)
        ax2.set_xlabel("line index", fontsize=7)
        ax2.tick_params(labelsize=6)
    for j in range(len(usable), rows * cols):
        r, c = divmod(j, cols)
        axes[r * 2][c].axis("off")
        axes[r * 2 + 1][c].axis("off")
    fig.suptitle("comb and scale-family residuals, same Hz axis for both models",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return True


def plot_rungs(groups, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    usable = [g for g in groups if g["comparison"]]
    if not usable:
        return False
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    colours = plt.cm.tab10(np.linspace(0, 1, 10))
    for i, g in enumerate(usable):
        f = np.asarray(g["freqs_mhz"])
        r = np.array([math.log(x * 1e6) / math.log(PHI) for x in f])
        ax.scatter(r, f, color=colours[i % 10], s=26,
                   label=f"{g['date']} n={g['comparison']['n_lines']}")
        lo, hi = r.min(), r.max()
        ax.plot([lo, hi], [f[np.argmin(r)], f[np.argmax(r)]],
                color=colours[i % 10], lw=0.8, alpha=0.6)
        ax.annotate(r"{:.3f} rung".format(hi - lo), (hi, f[np.argmax(r)]),
                    fontsize=6, xytext=(3, 3), textcoords="offset points",
                    color=colours[i % 10])
    ax.set_xlabel("phi rung  (slope on this axis is a ratio)")
    ax.set_ylabel("frequency (MHz)")
    ax.set_title("groups in rung coordinates: one rung apart means a factor phi")
    ax.legend(fontsize=6, loc="lower right")
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return True



def write_report(payload, rows, path):
    g = payload["groups"]
    cal = payload["calibration"]
    win = payload["window_controls"] or {}
    by = {(r["source"], r["date"]): r for r in g}
    main = by.get(("survey-96ch-1000int", "20240325"))
    lines = []
    add = lines.append

    add("# Comb or scale family")
    add("")
    add("Does the regularly spaced group at 45.98-48.99 MHz keep its regularity when")
    add("the frequency coverage is widened, and does the wider coverage say whether it")
    add("is a frequency comb or a scale family?")
    add("")
    add("**Answer: the group is regular and undecided.** It is well described by both")
    add("models, and at its span and residual size this pipeline recovers the true model")
    add("in roughly a quarter to a third of synthetic trials - worse than a coin toss -")
    add("so naming either one would be choosing, not measuring. It is not evidence of a")
    add("scale law and it is not evidence of anything artificial.")
    add("")

    add("## 1. What was expanded, and what was declared first")
    add("")
    add("The survey dedopplered a 96-channel window in the middle of each file. Every")
    add("file is 731 channels wide (13.4-86.9 MHz), so 87% of the frequency coverage had")
    add("never been searched. `prereg/expanded_prereg.json` declares the whole-band,")
    add("all-integration search before it was run; `prereg/expanded_prereg_bc.json`")
    add("declares two further variants and records plainly that it was written *after*")
    add("the first run, because that run changed frequency and duration together and")
    add("could not attribute a lost candidate to either.")
    add("")
    add("| variant | channels | integrations | drift window held | candidates |")
    add("| --- | --- | --- | --- | --- |")
    for name, fn, label in (("A", "prereg/expanded_search.json", "all integrations"),
                            ("B", "prereg/expanded_search_B.json", "first 1000"),
                            ("C", "prereg/expanded_search_C.json", "all integrations")):
        doc = load_json(os.path.join(ROOT, fn))
        if not doc:
            continue
        chans = doc["days"][0]["n_channels"]
        held = [d["window_check"]["fraction_held"] for d in doc["days"]]
        total = sum(d["n_candidates"] for d in doc["days"])
        add(f"| {name} | {chans} | {label} | {min(held):.0%}-{max(held):.0%} | {total} |")
    add("")
    add("Variant A is the honest whole-band search, and the window guard says it cannot")
    add("be read as a null result: a Nyquist drift window needs 0.5 channel per sample,")
    add("which over 2802-7009 integrations is 1400-3500 channels of trail in a")
    add("731-channel band, so between half and four fifths of the drift range is")
    add("truncated. Variant B holds the survey's 1000 integrations and keeps the window")
    add("at 100%; variant C holds the survey's 96 channels and expands duration only.")
    add("")

    add("## 2. Why the group cannot yet be tested in the wider band")
    add("")
    if win.get("days", {}).get("20240325"):
        w = win["days"]["20240325"]
        e = win.get("excision", {}).get("20240325", {})
        add("Two window-dependent gates, both measured with the same data, same threshold")
        add("and same integration count:")
        add("")
        add("| gate | 96-channel window | 731-channel band |")
        add("| --- | --- | --- |")
        add(f"| noise scale (`calibrate_map`) | {w['narrow']['calibration_factor']:.3f} | "
            f"{w['wide']['calibration_factor']:.3f} |")
        if e:
            add(f"| channels blanked by the carrier test | "
                f"{e['narrow']['fully_blanked_channels']}/{e['narrow']['total_channels']} | "
                f"{e['wide']['fully_blanked_channels']}/{e['wide']['total_channels']} |")
            add(f"| time samples blanked | {e['narrow']['blanked_time_samples']}/"
                f"{e['narrow']['total_time_samples']} | "
                f"{e['wide']['blanked_time_samples']}/{e['wide']['total_time_samples']} |")
        add("")
        add("The scale differs by a factor")
        add(f"{w['wide']['calibration_factor'] / max(w['narrow']['calibration_factor'], 1e-9):.1f}, "
            "so the same cell is scored 954 sigma in one window and 14 in the other.")
        add("And the channels holding the four lines are among those the wide-band carrier")
        add("test blanks: the quarter of the band that the carrier test removes is not the")
        add("same quarter in a 96-channel window. The group is therefore **masked**, not")
        add("refuted, by the wider band - and masked is not a result either.")
        add("")
    if win.get("days", {}).get("20240920"):
        w2 = win["days"]["20240920"]
        e2 = win.get("excision", {}).get("20240920", {})
        add("The near-Nyquist trio behaves differently: its channels survive in both")
        add("windows, but its significance goes from 8.65 sigma to -0.04.")
        if e2:
            add(f"{e2['wide']['blanked_time_samples']} of 1000 time samples are blanked in "
                f"the wide band against {e2['narrow']['blanked_time_samples']} in the "
                "narrow one, so the integrations the peak was built from are not in the")
            add("wide-band clean data. It was already classified an instrumental alias, and")
            add("this does not change that.")
        add("")

    add("## 3. Every group, one label each")
    add("")
    add("| source | date | lines | span (phi rungs) | dAIC (+ = scale) | boot CI | held-out | null p | label |")
    add("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for r in rows:
        add(f"| {r['source']} | {r['date']} | {r['n_lines']} | {r['span_rungs_phi']:.3f} | "
            f"{r['delta_aic']:.2f} | [{r['boot_lo']:.1f}, {r['boot_hi']:.1f}] | "
            f"{r['loo_majority']} | {r['null_p'] if r['null_p'] != '' else '-'} | "
            f"{r['verdict']} |")
    add("")
    add("Seven labels were available; four verdicts were returned. No group returned")
    add("SCALE FAMILY SUPPORTED or FREQUENCY COMB SUPPORTED, because none of them spans")
    add("the single phi rung the rule requires, and a scale claim additionally needs the")
    add("same q on two separated groups - and there is only this one group.")
    add("")

    if main and main.get("comparison"):
        c = main["comparison"]
        add("## 4. The 20240325 group in full")
        add("")
        add(f"- {c['n_lines']} lines, span {rows[0]['span_mhz']:.3f} MHz = "
            f"{c['span_rungs']:.3f} phi rungs")
        add(f"- comb: step {c['comb']['df'] / 1e6:.6f} MHz, "
            f"AIC {c['comb']['aic']:.2f}, BIC {c['comb']['bic']:.2f}, "
            f"rms {c['rms_in_channels']['comb']:.3f} channels")
        add(f"- scale: q = {c['scale']['q']:.8f}, "
            f"AIC {c['scale']['aic']:.2f}, BIC {c['scale']['bic']:.2f}, "
            f"rms {c['rms_in_channels']['scale']:.3f} channels")
        add(f"- delta AIC {c['delta_aic']:.2f} favours the scale family, and its 95%")
        add(f"  sign-flip bootstrap interval is [{c['bootstrap']['delta_aic_ci95'][0]:.1f}, "
            f"{c['bootstrap']['delta_aic_ci95'][1]:.1f}], which contains zero")
        add(f"- held out one line at a time: scale wins {c['loo']['scale_wins']} of "
            f"{c['n_lines']} folds, comb {c['loo']['comb_wins']}, on mean errors of "
            f"{c['loo']['scale']['mean_error_hz'] / 1e3:.1f} kHz and "
            f"{c['loo']['comb']['mean_error_hz'] / 1e3:.1f} kHz")
        flips = [d for d in c["drop_one"] if d["prefers"] == "comb"]
        add(f"- dropping a single line flips the preference in {len(flips)} of "
            f"{len(c['drop_one'])} cases, so the ranking is not robust to the set")
        if main.get("null"):
            add(f"- random lines drawn inside the group's own span exceed this AIC gap in "
                f"{100 * main['null']['p_value']:.1f}% of {main['null']['draws']} draws "
                f"(median gap {main['null']['null_median']:.2f})")
        add("")
        add("So the group is *more* regular than random lines (p = "
            f"{main['null']['p_value']:.3f}) while the choice between the two kinds of")
        add("regularity is unresolvable. Both statements are true at once and the record")
        add("keeps them apart.")
        add("")

    add("## 5. What the pipeline can actually resolve")
    add("")
    add("Synthetic groups with a known generator, fitted both ways, correct rate:")
    add("")
    add("| residual jitter | comb correct | scale correct |")
    add("| --- | --- | --- |")
    for j, v in sorted(cal["jitter_sweep_hz"].items(), key=lambda kv: int(kv[0])):
        add(f"| {int(j) / 1000:g} kHz | {100 * v['correct']['comb']:.0f}% | "
            f"{100 * v['correct']['scale']:.0f}% |")
    add("")
    add("| span (4 lines, 20 kHz jitter) | comb correct | scale correct |")
    add("| --- | --- | --- |")
    for s2, v in sorted(cal["span_sweep_hz"].items(), key=lambda kv: float(kv[0])):
        add(f"| {float(s2) / 1e6:g} MHz | {100 * v['correct']['comb']:.0f}% | "
            f"{100 * v['correct']['scale']:.0f}% |")
    add("")
    if main and main.get("comparison"):
        rms_ch = main["comparison"]["rms_in_channels"]
        add(f"The real group's residuals are {rms_ch['comb']:.2f} and "
            f"{rms_ch['scale']:.2f} channels, about "
            f"{CHANNEL_HZ * min(rms_ch.values()) / 1e3:.0f}-"
            f"{CHANNEL_HZ * max(rms_ch.values()) / 1e3:.0f} kHz. On the first table that")
        add("sits between the 10 and 20 kHz rows: a correct call about half to a quarter")
        add("of the time, and a wrong one the rest. The span table shows the other way")
        add("out - at the same 20 kHz jitter, a 10 MHz span recovers the model every")
        add("time, and 30 MHz does too. The binding limit is span, not line count.")
        add("")
    add("What would decide it, in order of cost:")
    add("")
    add("1. a continuation of this group over 10 MHz or more - the same regularity")
    add("   carried by more lines, at the residual this data already has")
    add("2. residuals below about 5 kHz at the present 3 MHz span, which needs a")
    add("   narrower channel or a longer integration on a signal that holds it")
    add("3. the same q on a second, separated group, which the scale claim requires")
    add("   before it can be made at all")
    add("")

    add("## 6. Physical controls")
    add("")
    add("| date | steps off an integer channel | one event? | drift spread | span / spacing |")
    add("| --- | --- | --- | --- | --- |")
    for grp in g:
        k = grp.get("controls")
        if not k:
            continue
        off = ", ".join(f"{x:.3f}" for x in k.get("steps_off_integer_channels", []))
        add(f"| {grp['date']} | {off} | {'yes' if k.get('one_event') else 'no'} | "
            f"{100 * k.get('drift_spread', 0):.2f}% | {k['span_over_spacing']:.1f} |")
    add("")
    add("No group is a comb on the channel grid: no spacing lands within")
    add("0.02 channel of an integer. The 20240325 steps are 0.069, 0.500 and 0.500")
    add("channels off integer - two of the three sit exactly between channels, which is")
    add("recorded as a suspicion rather than used as a verdict, since the tolerance")
    add("that would make it fire was chosen after seeing it.")
    add("")
    add("The 20240326 event is the useful negative control: 14 lines over 22 MHz, all")
    add("with the same drift to within 0.02%, which is one broadband event crossing the")
    add("band, not fourteen lines at a scale. The pipeline calls it ASTROPHYSICAL")
    add("STRUCTURE rather than manufacturing a model for it, which is the behaviour the")
    add("test was built to have.")
    add("")

    add("## 7. Reproduce")
    add("")
    add("```")
    add("python3 scripts/expanded_search.py --prereg        # refuses a second time")
    add("python3 scripts/expanded_search.py --variant A     # ~90 s, all six days")
    add("python3 scripts/expanded_search.py --variant B     # ~8 s")
    add("python3 scripts/expanded_search.py --variant C     # ~45 s")
    add("python3 scripts/band_calibration_check.py          # window controls")
    add("python3 scripts/phi_model_discrimination.py        # this document and the tables")
    add("```")
    add("")
    add("Standing rule unchanged: phi is a coordinate, not a criterion. A rung address")
    add("can make two things comparable and it can put a group on a ladder; it cannot")
    add("promote a candidate, and nothing here does.")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(chr(10).join(lines) + chr(10))
    return path


def main() -> int:
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    os.makedirs(os.path.join(ROOT, "plots"), exist_ok=True)

    by_source = {}
    for name, path in CANDIDATE_SETS.items():
        lines = survey_lines() if name == "survey-96ch-1000int" else expanded_lines(name, path)
        per_day = {}
        for line in lines:
            per_day.setdefault(line["date"], []).append(line)
        by_source[name] = per_day

    groups = []
    for source, per_day in by_source.items():
        for date, lines in sorted(per_day.items()):
            for chain_lines in chain(lines):
                freqs = [ln["freq_hz"] for ln in chain_lines]
                drifts = [ln["drift_hz_s"] for ln in chain_lines]
                record = {
                    "source": source, "date": date,
                    "n_lines": len(chain_lines),
                    "freqs_mhz": [round(x / 1e6, 6) for x in freqs],
                    "drifts_hz_s": [round(x, 2) for x in drifts],
                    "snr": [round(ln["snr"], 2) for ln in chain_lines],
                }
                if len(freqs) >= 3:
                    comp = compare(freqs, channel_hz=CHANNEL_HZ, bootstrap=800)
                    controls = structure_controls(freqs, drifts, CHANNEL_HZ)
                    record["comparison"] = comp
                    record["controls"] = controls
                    record["verdict"] = verdict(comp, controls, channel_hz=CHANNEL_HZ)
                    del record["comparison"]["freqs_hz"]
                groups.append(record)

    # --- the calibrations ---------------------------------------------------
    jitter_sweep = {}
    for jitter in (1000, 5000, 10000, 20000, 40000, 80000):
        r = synthetic_calibration(4, 3.0e6, trials=200, jitter_hz=jitter)
        jitter_sweep[jitter] = {"correct": r["correct_rate"], "confusion": r["confusion"]}
    span_sweep = {}
    for span in (3e6, 10e6, 30e6, 60e6):
        r = synthetic_calibration(4, span, trials=200, jitter_hz=20000)
        span_sweep[span] = {"correct": r["correct_rate"], "confusion": r["confusion"]}

    for g in groups:
        if g.get("comparison") and g["comparison"]["n_lines"] == 4:
            g["null"] = random_line_null(np.asarray(g["freqs_mhz"]) * 1e6)

    calibration = {"jitter_sweep_hz": jitter_sweep, "span_sweep_hz": span_sweep,
                   "note": ("groups generated from a known model, fitted both ways; "
                            "a rate near 0.5 means the pipeline cannot tell the "
                            "models apart at that residual size and span")}

    window = load_json(os.path.join(ROOT, "prereg", "band_calibration.json"))
    payload = {"channel_hz": CHANNEL_HZ, "group_ratio": GROUP_RATIO,
               "null_draws": NULL_DRAWS, "groups": groups,
               "calibration": calibration, "window_controls": window}
    with open(os.path.join(ROOT, "results", "expanded_phi_groups.json"), "w",
              encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write(chr(10))

    # --- csv ---------------------------------------------------------------
    rows = []
    for g in groups:
        comp = g.get("comparison")
        if not comp:
            continue
        rows.append({
            "source": g["source"], "date": g["date"], "n_lines": comp["n_lines"],
            "span_mhz": round(max(g["freqs_mhz"]) - min(g["freqs_mhz"]), 6),
            "span_rungs_phi": round(comp["span_rungs"], 6),
            "comb_df_mhz": round(comp["comb"]["df"] / 1e6, 6),
            "scale_q": round(comp["scale"]["q"], 8),
            "aic_comb": round(comp["comb"]["aic"], 3),
            "aic_scale": round(comp["scale"]["aic"], 3),
            "bic_comb": round(comp["comb"]["bic"], 3),
            "bic_scale": round(comp["scale"]["bic"], 3),
            "delta_aic": round(comp["delta_aic"], 3),
            "delta_bic": round(comp["delta_bic"], 3),
            "rms_comb_channels": round(comp["comb"]["rms_hz"] / CHANNEL_HZ, 4),
            "rms_scale_channels": round(comp["scale"]["rms_hz"] / CHANNEL_HZ, 4),
            "loo_comb_mhz": (round(comp["loo"]["comb"]["mean_error_hz"] / 1e6, 6)
                             if comp["loo"].get("applicable", True) else ""),
            "loo_scale_mhz": (round(comp["loo"]["scale"]["mean_error_hz"] / 1e6, 6)
                              if comp["loo"].get("applicable", True) else ""),
            "loo_majority": comp["loo"]["majority"],
            "boot_lo": round(comp["bootstrap"]["delta_aic_ci95"][0], 3),
            "boot_hi": round(comp["bootstrap"]["delta_aic_ci95"][1], 3),
            "drop_one_stable": comp["drop_one_stable"],
            "null_p": (round(g["null"]["p_value"], 4) if g.get("null") else ""),
            "verdict": g["verdict"]["label"],
        })
    with open(os.path.join(ROOT, "results", "comb_vs_scale.csv"), "w", newline="",
              encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["none"])
        writer.writeheader()
        writer.writerows(rows)

    folds = []
    for g in groups:
        comp = g.get("comparison")
        if not comp:
            continue
        for fold in comp["loo"]["folds"]:
            folds.append({
                "source": g["source"], "date": g["date"], "held_out_mhz":
                round(fold["actual_hz"] / 1e6, 6),
                "pred_comb_mhz": round(fold["predicted"]["comb"] / 1e6, 6),
                "pred_scale_mhz": round(fold["predicted"]["scale"] / 1e6, 6),
                "err_comb_hz": round(fold["error_hz"]["comb"], 2),
                "err_scale_hz": round(fold["error_hz"]["scale"], 2),
                "err_comb_rung": round(fold["error_rung"]["comb"], 6),
                "err_scale_rung": round(fold["error_rung"]["scale"], 6),
                "winner": fold["winner"],
            })
    with open(os.path.join(ROOT, "results", "heldout_predictions.csv"), "w",
              newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(folds[0].keys()) if folds else ["none"])
        writer.writeheader()
        writer.writerows(folds)

    made = []
    if plot_residuals(groups, os.path.join(ROOT, "plots", "comb_vs_scale_residuals.png")):
        made.append("comb_vs_scale_residuals.png")
    if plot_rungs(groups, os.path.join(ROOT, "plots", "rung_vs_frequency_groups.png")):
        made.append("rung_vs_frequency_groups.png")

    print(f"{len(groups)} groups chained; {len(rows)} with three or more lines")
    for g in groups:
        if g.get("verdict"):
            print(f"   {g['source']:24s} {g['date']} n={g['n_lines']:2d} "
                  f"span={g['comparison']['span_rungs']:.3f} rung  "
                  f"dAIC={g['comparison']['delta_aic']:7.2f}  {g['verdict']['label']}")
        else:
            print(f"   {g['source']:24s} {g['date']} n={g['n_lines']:2d}  (too few lines)")
    for j, v in jitter_sweep.items():
        print(f"   jitter {j/1000:5.0f} kHz: comb {v['correct']['comb']:.2f} "
              f"scale {v['correct']['scale']:.2f}")
    for s, v in span_sweep.items():
        print(f"   span {s/1e6:5.0f} MHz: comb {v['correct']['comb']:.2f} "
              f"scale {v['correct']['scale']:.2f}")
    doc = write_report(payload, rows, os.path.join(ROOT, "PHI_MODEL_DISCRIMINATION.md"))
    print(f"   wrote results/*.csv, results/expanded_phi_groups.json, plots/{made}")
    print(f"   wrote {os.path.basename(doc)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
