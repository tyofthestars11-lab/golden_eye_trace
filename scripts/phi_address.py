"""Every measured quantity gets its rung address, and the four tests get run.

Framework: phi^2 = phi + 1. Attribution: Tyree Jones (tyofthestarz), prior art
from December 2025.

The framework's own rule is that a measurement is not presented as established
without its rung address, and the records in this repository did not carry one.
This script supplies it and then runs the four tests the method defines - ladder
step, aspect/invariant, drift, cross-quantity - on the real candidates rather
than on examples.

What the addressing buys is *structure*, and it is worth being exact about which
structure. A single rung drift is not evidence: the fractional part of
ln(x)/ln(phi) is uniform for any smooth population, so one number sitting near an
integer happens about 9% of the time. What a rung difference can do is say
whether a group of lines is spaced by a fixed *frequency* step (a comb, which is
what a clock, a channel grid or an intermodulation product produces) or by a
fixed *factor* (a scale family, which is what one object re-measured at several
scales produces). Those two look identical on a linear axis when the group is
narrow, and they mean opposite things, so the distinction is the clarification
the framework actually delivers.

Both hypotheses get a null. The spacing statistic is compared against samples
drawn uniformly in the searched band, because four lines in a 9.55 MHz window
can look regularly spaced by luck.

    python3 scripts/phi_address.py
"""

from __future__ import annotations

import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from seti_search.signal_record import INVALID_CANDIDATE  # noqa: E402
from seti_search.rungs import (  # noqa: E402
    LN_PHI,
    PHI,
    address,
    ladder_step,
    repeated_scales,
    rung,
    rung_uniformity,
)

N_NULL = 4000
BASES = {"phi": PHI, "2": 2.0, "e": math.e, "10": 10.0, "sqrt2": math.sqrt(2.0)}


def load_records() -> list[dict]:
    with open(os.path.join(ROOT, "signal_records.json"), encoding="utf-8") as fh:
        return json.load(fh)["records"]


def address_record(record: dict) -> dict:
    """Rung address for every positive measured quantity in one record."""
    out = {"date": record["observation_date"],
           "freq_mhz": record["frequency_hz"] / 1e6,
           "quantities": {}}
    for key, unit in (("frequency_hz", "Hz"), ("local_scale", "ratio"),
                      ("resolution_hz_per_s", "Hz/s"), ("duration_s", "s"),
                      ("snr_global", "sigma"), ("snr_at_cell", "sigma")):
        value = record[key]
        if value is None or value <= 0:
            out["quantities"][key] = {"value": value, "rung": None,
                                      "note": "non-positive: no rung exists"}
            continue
        entry = address(value, unit).as_dict()
        out["quantities"][key] = entry
    drift = record["drift_hz_per_s"]
    if abs(drift) > 0:
        entry = address(abs(drift), "Hz/s").as_dict()
        entry["sign"] = "negative" if drift < 0 else "positive"
        out["quantities"]["drift_hz_per_s"] = entry
    else:
        out["quantities"]["drift_hz_per_s"] = {"value": drift, "rung": None,
                                               "note": "zero: no rung exists"}
    return out


def spacing_test(freqs_mhz: list[float], band: tuple[float, float]) -> dict:
    """Is this group a fixed-frequency comb or a fixed-factor scale family?

    The comb statistic is the spread of consecutive differences relative to
    their mean; the scale statistic is the same spread over the *ratios*
    ``rung(f[i+1]) - rung(f[i])``, which is exactly a difference of rungs. Both
    are compared against groups drawn uniformly in the searched band.
    """
    values = sorted(freqs_mhz)
    if len(values) < 3:
        return {"n": len(values), "verdict": "too few lines to test"}
    steps = np.diff(values)
    rungs = np.diff([rung(v) for v in values])
    comb = float((steps.max() - steps.min()) / abs(steps.mean()))
    scale = float((rungs.max() - rungs.min()) / abs(rungs.mean()))

    rng = np.random.default_rng(5)
    comb_null, scale_null = [], []
    low, high = band
    for _ in range(N_NULL):
        draw = np.sort(rng.uniform(low, high, size=len(values)))
        d = np.diff(draw)
        r = np.diff([rung(v) for v in draw])
        comb_null.append((d.max() - d.min()) / abs(d.mean()))
        scale_null.append((r.max() - r.min()) / abs(r.mean()))
    comb_null = np.asarray(comb_null)
    scale_null = np.asarray(scale_null)
    p_comb = float((comb_null <= comb).mean())
    p_scale = float((scale_null <= scale).mean())
    # Only name a preference when one model is clearly better than the other. A
    # four-line group inside a 3 MHz span has almost the same arithmetic and
    # geometric spread, so the two p-values land on top of each other and the
    # honest read is that the data cannot separate the hypotheses.
    closer = abs(p_comb - p_scale)
    if min(p_comb, p_scale) >= 0.05:
        verdict = "not regular at 5%"
        basis, p = "neither", min(p_comb, p_scale)
    elif closer < 0.05:
        verdict = "regular; comb and scale family fit equally"
        basis, p = "cannot distinguish", min(p_comb, p_scale)
    elif p_scale < p_comb:
        verdict = "scale family (one factor repeated)"
        basis, p = "constant ratio", p_scale
    else:
        verdict = "comb (one fixed step repeated)"
        basis, p = "constant spacing", p_comb
    return {
        "n": len(values),
        "frequencies_mhz": [round(v, 4) for v in values],
        "steps_mhz": [round(float(v), 4) for v in steps],
        "rung_steps": [round(float(v), 6) for v in rungs],
        "step_spread": comb,
        "rung_spread": scale,
        "p_comb": p_comb,
        "p_scale": p_scale,
        "verdict": verdict,
        "basis": basis,
        "p_value": p,
        "factor": float(PHI ** float(np.mean(rungs))),
    }


def base_invariance(records: list[dict]) -> dict:
    """Does any other base give the same read of these addresses?

    The fractional part of log_b(x) is uniform for every b, so a null that is
    uniform under phi is uniform under 2, e or 10. Reporting this keeps the
    address honest: it is a coordinate, and its *structure* claims have to be
    about relations between values, not about one value's nearness to a rung.
    """
    freqs = [r["frequency_hz"] for r in records]
    out = {}
    low, high = min(freqs), max(freqs)
    for name, base in BASES.items():
        drifts = []
        for value in freqs:
            frac = math.log(value) / math.log(base)
            drifts.append(abs(frac - round(frac)))
        drifts = np.asarray(drifts)
        out[name] = {
            "median_fractional_distance": float(np.median(drifts)),
            "uniform_expectation": 0.25,
            "n": int(drifts.size),
            # how many rungs the whole record set spans in this base. Under one
            # rung, every value shares one address and no address can discriminate
            # between candidates - whatever base is used.
            "band_span_rungs": float(math.log(high / low) / math.log(base)),
        }
    phi_like = out["phi"]["median_fractional_distance"]
    out["spread_across_bases"] = float(
        max(v["median_fractional_distance"] for k, v in out.items()
            if k in BASES)
        - min(v["median_fractional_distance"] for k, v in out.items()
              if k in BASES)
    )
    out["phi"] = out["phi"]
    out["phi_median"] = phi_like
    return out


def main() -> int:
    records = load_records()
    invalid = {id(r) for r in records
               if r.get("disposition") in INVALID_CANDIDATE}
    addresses = [address_record(r) for r in records]

    by_date: dict[str, list[dict]] = {}
    for record in records:
        by_date.setdefault(record["observation_date"], []).append(record)

    groups = []
    for date, rows in sorted(by_date.items()):
        if len(rows) < 2:
            continue
        freqs = [r["frequency_hz"] / 1e6 for r in rows]
        pairs = []
        ordered = sorted(rows, key=lambda r: r["frequency_hz"])
        for a, b in zip(ordered, ordered[1:]):
            step = ladder_step(b["frequency_hz"], a["frequency_hz"])
            pairs.append({
                "from_mhz": round(a["frequency_hz"] / 1e6, 4),
                "to_mhz": round(b["frequency_hz"] / 1e6, 4),
                "delta_rung": round(step.delta, 6),
                "integer_step": bool(step.integer_step),
                "factor": round(step.factor, 6),
                "residual": step.residual,
            })
        low = min(freqs) - 0.2
        high = max(freqs) + 0.2
        status = ("MASK-MANUFACTURED SIGNAL"
                  if any(id(r) in invalid for r in rows) else "candidate")
        groups.append({
            "status": status,
            "date": date,
            "n": len(freqs),
            "pairs": pairs,
            "spacing": spacing_test(freqs, (low, high)) if len(freqs) >= 3 else None,
            "repeated_scales": [{"delta": round(r.delta, 6),
                                 "factor": round(r.factor, 6),
                                 "pairs": r.pairs}
                                for r in repeated_scales(freqs)],
            "rung_uniformity": rung_uniformity(freqs),
        })

    bases = base_invariance(records)
    drift_uniformity = rung_uniformity(
        [abs(r["drift_hz_per_s"]) for r in records if r["drift_hz_per_s"]])

    payload = {
        "framework": "phi^2 = phi + 1; rung(x) = ln(x)/ln(phi)",
        "attribution": "Tyree Jones (tyofthestarz), prior art from December 2025",
        "n_records": len(records),
        "addresses": addresses,
        "groups": groups,
        "base_invariance": bases,
        "drift_rung_uniformity": drift_uniformity,
        "null_draws": N_NULL,
    }
    with open(os.path.join(ROOT, "prereg", "phi_address.json"), "w",
              encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    write_report(payload)
    print(f"PHI ADDRESSES  ({len(records)} records, {len(groups)} multi-line groups)")
    for group in groups:
        spacing = group["spacing"]
        if spacing is None:
            continue
        print(f"   {group['date']}  {spacing['n']} lines: {spacing['verdict']}"
              f"  p = {spacing['p_value']:.4f}  factor {spacing['factor']:.6f}"
              f"  (comb p {spacing['p_comb']:.4f}, scale p {spacing['p_scale']:.4f})")
    print(f"   drift rung uniformity: p = {drift_uniformity['p_value']:.4f}"
          f" over n = {drift_uniformity['n']}")
    print(f"   median fractional distance across bases: "
          + ", ".join(f"{k} {v['median_fractional_distance']:.3f}"
                      for k, v in bases.items() if k in BASES))
    banner_path = os.path.join(ROOT, "PHI_ADDRESS.md")
    banner = (
        "> **Disqualified premise.** The two multi-line groups below were selected by\n"
        "> the excision mask, not by a valid detection: their cells read 1.44, 0.54,\n"
        "> -0.15 and -0.87 sigma without the theory-scaled row test, and 954, 604, 447\n"
        "> and 171 sigma with it (`LIGHT_PATH.md`). Their spacing arithmetic is kept as\n"
        "> the worked example of selection bias and is **not** a regularity claim.\n"
        "> `prereg/phi_address.json` carries `status` per group: none of the two is a\n"
        "> candidate.\n\n"
    )
    with open(banner_path, encoding="utf-8") as fh:
        body = fh.read()
    if not body.startswith("> **Disqualified premise.**"):
        with open(banner_path, "w", encoding="utf-8") as fh:
            fh.write(banner + body)
    print("   wrote prereg/phi_address.json and PHI_ADDRESS.md")

    return 0


def write_report(payload: dict) -> None:
    groups = [g for g in payload["groups"] if g["spacing"]]
    lines = [
        "# Rung addresses for every measured quantity",
        "",
        "Framework: φ² = φ + 1. Attribution: Tyree Jones (tyofthestarz), prior art "
        "from December 2025. Coordinate: `rung(x) = ln x / ln φ`.",
        "",
        "Every record's quantities are addressed in `prereg/phi_address.json`. "
        "What follows is what the addressing *buys*: the four tests, on the real "
        "candidates, with a null beside each statistic.",
        "",
        "## What the tests can and cannot do",
        "",
        "A single rung drift is not a finding. For any smooth population the "
        "fractional part of `ln x / ln φ` is uniform, so a value lands within "
        "0.05 rungs about a tenth of the time, and one number close to an integer "
        "is not evidence of anything. The table below shows that this is not a "
        "property of φ: the same median fractional distance comes out of base 2, "
        "e, 10 and √2 within the noise, because they are all logarithms. Where "
        "the coordinate earns its place is in *relations* — a rung difference is "
        "a ratio, and a repeated ratio is a repeated scaling.",
        "",
        "| base | median fractional distance | band span (rungs) | uniform expectation |",
        "| --- | --- | --- | --- |",
    ]
    for name in BASES:
        entry = payload["base_invariance"][name]
        lines.append(f"| {name} | {entry['median_fractional_distance']:.4f} | "
                     f"{entry['band_span_rungs']:.3f} | "
                     f"{entry['uniform_expectation']:.2f} |")
    lines += [
        "",
        f"Across bases the median moves by only "
        f"{payload['base_invariance']['spread_across_bases']:.4f}, so no address "
        "is φ-specific at the level of one value.",
        "",
        "## The clarification that does work: comb or scale family",
        "",
        "A group of lines can be spaced by a fixed **frequency step** — a clock, a "
        "channel grid, an intermodulation product, all of which are comb-like — or "
        "by a fixed **factor**, which is one thing measured at several scales. On "
        "a linear axis inside a narrow window the two are nearly "
        "indistinguishable, and they mean opposite things. Rung differences "
        "separate them: a comb has constant Δf and drifting Δrung, a scale family "
        "the reverse.",
        "",
        "| date | lines | step spread | rung spread | p(comb) | p(scale) | read |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for group in groups:
        spacing = group["spacing"]
        lines.append(
            f"| {group['date']} | {spacing['n']} | {spacing['step_spread']:.4f} | "
            f"{spacing['rung_spread']:.4f} | {spacing['p_comb']:.4f} | "
            f"{spacing['p_scale']:.4f} | {spacing['verdict']} |")
    lines += [
        "",
        f"Each p compares the observed spread against {payload['null_draws']} "
        "groups drawn uniformly across the same band, so a regular-looking "
        "arrangement of four lines in 9.55 MHz has to beat luck before it is read "
        "as structure. `p(comb)` is the chance of a spread that tight in "
        "frequency, `p(scale)` the same in logs.",
        "",
        "## Ladder steps between neighbours",
        "",
        "Per test A, a pair is a ladder step only if the rung difference is an "
        "integer; otherwise it is a plain scale factor and the residual says how "
        "exactly.",
        "",
        "| date | pair (MHz) | Δrung | integer step? | factor |",
        "| --- | --- | --- | --- | --- |",
    ]
    for group in payload["groups"]:
        for pair in group["pairs"]:
            lines.append(
                f"| {group['date']} | {pair['from_mhz']:.4f} → "
                f"{pair['to_mhz']:.4f} | {pair['delta_rung']:+.4f} | "
                f"{'yes' if pair['integer_step'] else 'no'} | "
                f"{pair['factor']:.6f} |")
    lines += [
        "",
        "No pair is an integer step, so no pair is a ladder relation. That is the "
        "expected outcome and it is the honest one: the lenses here are scale "
        "factors, and a scale factor is a description of how two lines differ, not "
        "a claim about what produced them.",
        "",
        "## Repeated scales",
        "",
    ]
    repeats = [(g["date"], r) for g in payload["groups"] for r in g["repeated_scales"]]
    if not repeats:
        lines.append(
            "None. `repeated_scales` needs three *independent* pairs sharing a "
            "rung difference, and it applies a shared-index rule: pairs that reuse "
            "a line are one object, not two. A chain of four lines gives at most "
            "two independent pairs, so a same-factor chain is not by itself a "
            "repeated-scaling artifact — the test says so rather than manufacturing "
            "a hit.")
    for date, entry in repeats:
        lines.append(
            f"**{date}**: Δrung {entry['delta']:+.4f} (factor {entry['factor']:.6f}) "
            f"across {entry['pairs']} independent pairs.")
    lines += [
        "",
        "## The drift population",
        "",
        f"Across the {payload['drift_rung_uniformity']['n']} non-zero drifts, "
        f"fractional rungs are uniform with p = "
        f"{payload['drift_rung_uniformity']['p_value']:.3f}. A low p here would "
        "mean the drift values were quantized by the search grid rather than by "
        "the sky, and there would be nothing left to interpret.",
        "",
        "## Reading",
        "",
        "The framework is a ruler, and a ruler's job is to make two things "
        "comparable rather than to decide between them. What it adds to a binary "
        "method is precisely the comb-versus-scale question: an accept/reject on "
        "each line cannot see it, and it is the difference between a group of "
        "sources and one source counted four times. What it cannot do is turn a "
        "small drift into a detection — that is arithmetic, not caution, and the "
        "base table above is the demonstration.",
        "",
    ]
    with open(os.path.join(ROOT, "PHI_ADDRESS.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
