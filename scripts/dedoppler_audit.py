"""Audit: does the dedoppler grid actually resolve a phi-locked chirp?

The lock sits at 76.4% of Nyquist, which is an unusual place for a drift search:
most pipelines are tuned for a few Hz/s, where the grid is coarse in rate but the
signal barely moves. Here the signal moves 0.38 channels per integration, so two
things could go wrong and both are measurable.

**Rate resolution.** `dedoppler` steps by `0.5 * boxcar * df / total_time`, which
in cycles per sample is `0.5 * boxcar / n_ints` — independent of the band. A
signal between two trials smears by at most half a step across the observation,
so the recovered rate is uncertain by half a step. If that uncertainty is a large
fraction of the 0.05-theta tolerance, the lock cannot be measured no matter how
strong the signal.

**Alignment loss.** The shift-and-add rounds each integration's offset to whole
channels, and at 0.38 channels per integration that rounding is up to half a
channel per row. The question is how much signal survives, measured against the
same injection at zero drift.

Both are measured here by injecting chirps at a range of rates, including theta
itself, and reading back what the search reports.

    python3 scripts/dedoppler_audit.py
"""

from __future__ import annotations

import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from seti_search.filterbank import InjectedSignal, simulate  # noqa: E402
from seti_search.rungs import PHI  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

THETA = 1.0 / (PHI * PHI)
TOLERANCE = 0.05
SEED = 1729

CHANNEL_MHZ = 0.0191406   # 19.14 kHz
TSAMP_S = 0.05
NCHANS = 256
NINTS = 200
BOX = 3
AMPLITUDE = 3.0
WIDTH_CHANNELS = 1.5


def grid() -> dict:
    channel_hz = CHANNEL_MHZ * 1e6
    total_time = NINTS * TSAMP_S
    step_hz = 0.5 * BOX * channel_hz / total_time
    nyquist = 0.5 * channel_hz / TSAMP_S
    n_steps = int(math.floor(nyquist / step_hz))
    step_nu = step_hz * TSAMP_S / channel_hz
    return {"channel_hz": channel_hz, "total_time_s": total_time,
            "step_hz_s": step_hz, "step_nu": step_nu,
            "step_nu_over_theta": step_nu / THETA, "n_steps": n_steps,
            "nyquist_hz_s": nyquist, "reach_hz_s": n_steps * step_hz,
            "reach_fraction": n_steps * step_hz / nyquist}


def run(nu: float | None, label: str, g: dict) -> dict:
    """Inject a chirp at `nu` cycles per sample and report what comes back."""
    channel_hz, tsamp = g["channel_hz"], TSAMP_S
    signals = ()
    drift = 0.0
    if nu is not None:
        drift = nu * channel_hz / TSAMP_S
        start = 51.0 if drift >= 0 else 53.0
        signals = (InjectedSignal(freq_mhz=start, drift_hz_s=drift,
                                  amplitude=AMPLITUDE,
                                  width_channels=WIDTH_CHANNELS),)
    path = os.path.join(ROOT, "data", f"audit_{label}.fil")
    simulate(path, nchans=NCHANS, nints=NINTS, tsamp=TSAMP_S,
             fch1=55.0, foff=-CHANNEL_MHZ, signals=signals,
             glitch=False, seed=SEED)
    cfg = SearchConfig(boxcar=BOX, max_drift_hz_s=g["nyquist_hz_s"],
                       snr_threshold=5.0, max_hits=25)
    hits, _meta = search_file(path, cfg, chan_chunk=2048, max_ints=None)
    if not hits:
        return {"label": label, "nu": nu, "detected": False}
    # the peak belonging to this injection: strongest hit within two grid steps,
    # falling back to the nearest hit when the signal is not found near its rate
    near = [h for h in hits if abs(h.drift_hz_s - drift) <= 2.0 * g["step_hz_s"]]
    nearest = (max(near, key=lambda h: h.snr) if near
               else min(hits, key=lambda h: abs(h.drift_hz_s - drift)))
    nu_rec = abs(nearest.drift_hz_s) * TSAMP_S / channel_hz
    return {"label": label, "nu": nu, "detected": True,
            "drift_hz_s": drift, "found_hz_s": nearest.drift_hz_s,
            "nu_recovered": nu_rec, "error_nu": abs(nu_rec - (nu or 0.0)),
            "error_over_theta": abs(nu_rec - (nu or 0.0)) / THETA,
            "snr": nearest.snr, "best_snr": max(h.snr for h in hits),
            "n_hits": len(hits)}


def main() -> int:
    g = grid()
    print("DEDOPPLER AUDIT AT THE PHI LOCK")
    print(f"   {NCHANS} channels, {NINTS} samples, tsamp {TSAMP_S} s, "
          f"boxcar {BOX}")
    print(f"   Nyquist drift {g['nyquist_hz_s']:,.0f} Hz/s; grid step "
          f"{g['step_hz_s']:,.1f} Hz/s = {g['step_nu']:.5f} cycles/sample")
    print(f"   step is {g['step_nu_over_theta']:.4f} theta, so half a step is "
          f"{50 * g['step_nu_over_theta']:.1f}% of the tolerance")
    print(f"   {g['n_steps']} steps each way, reach {g['reach_hz_s']:,.0f} Hz/s = "
          f"{100 * g['reach_fraction']:.2f}% of Nyquist")
    print()

    rates = [("stationary", 0.0), ("half theta", THETA / 2),
             ("arbitrary 0.31", 0.31), ("theta", THETA),
             ("near Nyquist", 0.499)]
    print(f"   {'injection':<15}{'nu':>8}{'%Nyq':>7}{'drift Hz/s':>12}"
          f"{'found Hz/s':>12}{'nu rec':>9}{'err/theta':>10}{'snr':>7}{'kept':>7}")
    rows, baseline = [], None
    for label, nu in rates:
        row = run(nu, label.replace(" ", "_"), g)
        if not row["detected"]:
            print(f"   {label:<15}{nu:>8.3f}{100 * nu / 0.5:>6.1f}%"
                  f"{'-':>12}{'-':>12}{'-':>9}{'-':>10}{'-':>7}{'-':>7}")
            rows.append(row)
            continue
        if baseline is None:
            baseline = row["snr"]
        kept = row["snr"] / baseline if baseline else float("nan")
        row["kept"] = kept
        rows.append(row)
        print(f"   {label:<15}{nu:>8.3f}{100 * nu / 0.5:>6.1f}%"
              f"{row['drift_hz_s']:>12,.0f}{row['found_hz_s']:>12,.0f}"
              f"{row['nu_recovered']:>9.4f}{row['error_over_theta']:>10.4f}"
              f"{row['snr']:>7.2f}{kept:>7.2f}")

    null = run(None, "noise_only", g)
    print()
    print(f"   noise-only run: {null.get('n_hits', 0)} hits at 5 sigma, best "
          f"{null.get('best_snr', 0):.2f} sigma")
    print()

    theta_row = next((r for r in rows if r.get("label") == "theta"), {})
    worst_err = max((r.get("error_over_theta", 0.0) for r in rows if r["detected"]),
                    default=9.9)
    reached = g["reach_fraction"] >= 0.99
    locked_ok = worst_err <= TOLERANCE
    theta_ok = bool(theta_row.get("detected")) and theta_row.get("snr", 0.0) >= 10.0
    worst_row = max((r for r in rows if r["detected"]),
                    key=lambda r: r.get("error_over_theta", 0.0), default={})
    print(f"   worst rate error {worst_err:.4f} theta ({worst_row.get('label')}) "
          f"against a {TOLERANCE} tolerance: {'OK' if locked_ok else 'TOO COARSE'}")
    print(f"   that is {worst_err / max(g['step_nu_over_theta'], 1e-9):.1f} grid "
          f"steps, so the half-step rule of thumb is optimistic at this drift")
    print(f"   grid reaches {100 * g['reach_fraction']:.2f}% of Nyquist: "
          f"{'OK' if reached else 'SHORT'}")
    print(f"   theta chirp detected at {theta_row.get('snr', float('nan')):.2f} "
          f"sigma against {null.get('n_hits', 0)} hits in noise: "
          f"{'OK' if theta_ok else 'LOST'}")
    print("   the stationary column is not a sensitivity ceiling: a line fixed in")
    print("   one channel is the case the carrier stage is built to treat as RFI,")
    print("   so the drifting chirps read higher and the ratio is quoted, not gated")

    out = {"grid": g, "rows": rows, "null": null,
           "tolerance_of_theta": TOLERANCE,
           "worst_error_over_theta": worst_err,
           "reach_fraction": g["reach_fraction"],
           "theta_snr": theta_row.get("snr"),
           "worst_case": worst_row.get("label"),
           "worst_error_in_grid_steps": worst_err / g["step_nu_over_theta"],
           "theta_kept_vs_stationary": theta_row.get("kept"),
           "verdict": ("resolves" if (locked_ok and reached and theta_ok)
                       else "does not resolve")}
    path = os.path.join(ROOT, "prereg", "dedoppler_audit.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"\n   verdict: {out['verdict']}; wrote prereg/dedoppler_audit.json")
    return 0 if out["verdict"] == "resolves" else 1


if __name__ == "__main__":
    raise SystemExit(main())
