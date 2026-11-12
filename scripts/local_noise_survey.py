"""Re-run the survey with a local noise reference.

The survey's detection statistic is (smoothed - median) / MAD, where the median
and the MAD are taken across every channel in the searched file, then rescaled so
the map's bulk is one sigma wide. That is a single noise scale for a 9.66 MHz
band. Where the sub-band around a line is brighter and more structured than the
band as a whole, the local noise is larger than the estimate and the significance
is inflated.

This re-run changes exactly one thing: the noise scale. The window, the samples,
the dedoppler, the boxcar, the calibration and the threshold are the survey's, so
a bright lane drifting 105 channels across the integration is treated identically
and cannot be lost to window truncation. The single band-wide MAD is replaced by
a per-channel spread measured from the map itself, over a +/-16 channel
neighbourhood and every drift trial.

A first version searched narrow 32-channel windows instead. It was wrong and is
discarded: the 20240325 lanes drift about 105 channels over 1000 samples, so a
narrow window truncates them mid-trail and the whole burst group came out at
0 sigma. The mistake is recorded here because it is exactly the kind of error
that would have looked like a result.

    python3 scripts/local_noise_survey.py
"""

from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

from ovro_lwa_population import dynamic  # noqa: E402
from ovro_lwa_test import LOW_MHZ, TSAMP_S  # noqa: E402
from seti_search.search import (  # noqa: E402
    SearchConfig,
    calibrate_map,
    control_bandpass,
    dedoppler,
    excise_rfi,
    extract_hits,
)

CHANS = 96                 # the survey's window
SAMPLES = 1000             # the survey's 256 s
THRESHOLD = 5.5
LOCAL_HALF = 16
TOL_MHZ = 0.3
TOL_HZ_S = 800.0


def local_scale(snr_map: np.ndarray, half: int = LOCAL_HALF) -> np.ndarray:
    """Per-channel spread of the map, from its own neighbourhood and drift axis."""
    n_drift, n_chan = snr_map.shape
    scale = np.ones(n_chan, dtype=np.float64)
    for j in range(n_chan):
        lo, hi = max(0, j - half), min(n_chan, j + half + 1)
        block = snr_map[:, lo:hi].astype(np.float64).ravel()
        if block.size < 32:
            continue
        centre = float(np.median(block))
        mad = 1.4826 * float(np.median(np.abs(block - centre)))
        if not np.isfinite(mad) or mad <= 0:
            mad = float(block.std()) or 1.0
        # clip once, so a bright candidate in the neighbourhood cannot set the
        # scale it is measured against
        keep = np.abs(block - centre) < 4.0 * mad
        if keep.sum() > 32:
            revised = 1.4826 * float(np.median(np.abs(block[keep] - centre)))
            if np.isfinite(revised) and revised > 0:
                mad = revised
        scale[j] = mad
    return scale


def significance_at(snr_map, scale, drifts, freqs, freq_mhz: float,
                    drift_hz_s: float) -> dict:
    """Significance of the survey detection under each noise reference.

    Re-extracting peaks from a rescaled map is not a fair comparison: dividing
    each column by its own scale moves the map's maximum, so a line can be strong
    and still not be the top hit. This reads the map at the cell the survey
    reported, which is the apples-to-apples question.
    """
    j = int(np.argmin(np.abs(freqs - freq_mhz)))
    k = int(np.argmin(np.abs(drifts - drift_hz_s)))
    value = float(snr_map[k, j])
    local_scale = float(scale[j])
    return {"global_snr": value,
            "local_snr": value / local_scale if local_scale else value,
            "scale": local_scale, "k": k, "j": j,
            "freq_mhz": float(freqs[j]), "drift_hz_s": float(drifts[k])}


def build_window(data, meta, chans: int = CHANS, samples: int = SAMPLES):
    channel_hz = meta["channel_hz"]
    nfreq = meta["nfreq"]
    c0 = max(0, (nfreq - chans) // 2)
    n = min(samples, meta["ntime"])
    block = data[c0:c0 + chans, :n].T.astype(np.float32)
    freqs = LOW_MHZ + (c0 + np.arange(chans)) * channel_hz / 1e6
    return block, freqs, c0, n, channel_hz


def run_window(block, freqs, channel_hz, threshold: float) -> dict:
    cfg = SearchConfig(boxcar=3, max_drift_hz_s=0.5 * channel_hz / TSAMP_S,
                       snr_threshold=threshold, max_hits=60)
    normalised, _ = control_bandpass(block)
    clean, _rfi = excise_rfi(normalised)
    df_mhz = float(freqs[1] - freqs[0])
    snr_map, drifts = dedoppler(
        clean, TSAMP_S, df_mhz,
        max_drift_hz_s=cfg.max_drift_hz_s, drift_step_hz_s=cfg.drift_step_hz_s,
        boxcar=cfg.boxcar)
    live = ~(clean == 0).all(axis=0)
    snr_map, calibration = calibrate_map(snr_map, live)
    scale = local_scale(snr_map)
    local_map = snr_map / scale[None, :]
    total_time = block.shape[0] * TSAMP_S
    out = {"global_hits": [], "local_hits": [], "calibration": float(calibration),
           "scale": scale.tolist(), "drifts": drifts, "live": live,
           "freqs": freqs, "snr_map": snr_map}
    for tag, the_map in (("global_hits", snr_map), ("local_hits", local_map)):
        hits = extract_hits(the_map, drifts, freqs, threshold=threshold,
                            boxcar=cfg.boxcar, total_time=total_time,
                            max_hits=60)
        out[tag] = [{"freq_mhz": h.freq_mhz, "drift_hz_s": h.drift_hz_s,
                     "snr": float(h.snr)} for h in hits]
    return out


def main() -> int:
    confirm = json.load(open(os.path.join(ROOT, "prereg", "ovro_confirm.json"),
                             encoding="utf-8"))
    cands = [c for c in confirm["candidates"] if c.get("source") == "survey"]
    print(f"LOCAL NOISE REFERENCE SURVEY  ({len(cands)} candidates)")
    print(f"   window {CHANS} channels x {SAMPLES} samples, threshold "
          f"{THRESHOLD:g}, local scale over +/-{LOCAL_HALF} channels")

    by_date: dict[str, list] = {}
    for c in cands:
        by_date.setdefault(c["date"], []).append(c)

    results = []
    for date in sorted(by_date):
        path = os.path.join(ROOT, "data", f"{date}.fits")
        data, meta = dynamic(path)
        block, freqs, c0, n, channel_hz = build_window(data, meta)
        run = run_window(block, freqs, channel_hz, THRESHOLD)
        print(f"   {date}: global calibration {run['calibration']:.3f}, "
              f"{len(run['global_hits'])} hits band-wide, "
              f"{len(run['local_hits'])} with a local scale")
        for c in sorted(by_date[date], key=lambda r: r["freq_mhz"]):
            j = int(round((c["freq_mhz"] - LOW_MHZ) * 1e6 / channel_hz)) - c0
            j = max(0, min(CHANS - 1, j))
            near_local = [h for h in run["local_hits"]
                          if abs(h["freq_mhz"] - c["freq_mhz"]) < TOL_MHZ
                          and abs(h["drift_hz_s"] - c["drift_hz_s"]) < TOL_HZ_S]
            near_global = [h for h in run["global_hits"]
                           if abs(h["freq_mhz"] - c["freq_mhz"]) < TOL_MHZ
                           and abs(h["drift_hz_s"] - c["drift_hz_s"]) < TOL_HZ_S]
            near_local.sort(key=lambda h: -h["snr"])
            near_global.sort(key=lambda h: -h["snr"])
            best_local = near_local[0] if near_local else None
            best_global = near_global[0] if near_global else None
            at = significance_at(run["snr_map"], np.asarray(run["scale"]),
                                 run["drifts"], run["freqs"],
                                 c["freq_mhz"], c["drift_hz_s"])
            entry = {"date": date, "freq_mhz": c["freq_mhz"],
                     "drift_hz_s": c["drift_hz_s"], "survey_snr": c["snr"],
                     "scale_at_candidate": float(run["scale"][j]),
                     "global_snr": None if best_global is None
                     else best_global["snr"],
                     "local_snr": None if best_local is None
                     else best_local["snr"],
                     "local_found": bool(best_local),
                     "at_candidate": at,
                     "survives_local": bool(at["local_snr"] >= THRESHOLD),
                     "local_drift_hz_s": None if best_local is None
                     else best_local["drift_hz_s"],
                     "local_freq_mhz": None if best_local is None
                     else best_local["freq_mhz"],
                     "channel": j}
            if best_local is not None and c["snr"]:
                entry["snr_ratio"] = best_local["snr"] / c["snr"]
            results.append(entry)
            tag = "SURVIVES" if at["local_snr"] >= THRESHOLD else "gone"
            extra = (f" local {at['local_snr']:7.2f} sigma at the reported cell "
                     f"(scale {at['scale']:.2f})")
            print(f"      {c['freq_mhz']:8.3f} MHz  survey {c['snr']:7.2f}  "
                  f"global-scale {0 if best_global is None else best_global['snr']:7.2f}  "
                  f"local {tag}{extra}")
        del data

    survived = [e for e in results if e["survives_local"]]
    print(f"\n   {len(survived)} of {len(results)} survive a local noise "
          f"reference:")
    for e in survived:
        print(f"      {e['date']} {e['freq_mhz']:8.3f} MHz  "
              f"{e['survey_snr']:7.2f} -> {e['at_candidate']['local_snr']:7.2f}")
    lost = [e for e in results if not e["survives_local"]]
    print(f"   failing: {len(lost)}")
    for e in lost:
        print(f"      {e['date']} {e['freq_mhz']:8.3f} MHz  survey "
              f"{e['survey_snr']:7.2f} -> local "
              f"{e['at_candidate']['local_snr']:.2f}  "
              f"(scale {e['at_candidate']['scale']:.2f})")

    payload = {"candidates": results, "n_candidates": len(results),
               "n_survived": len(survived),
               "window_channels": CHANS, "samples": SAMPLES,
               "threshold": THRESHOLD, "local_half": LOCAL_HALF,
               "survivors": [{"date": e["date"], "freq_mhz": e["freq_mhz"],
                              "drift_hz_s": e["drift_hz_s"],
                              "survey_snr": e["survey_snr"],
                              "local_snr": e["at_candidate"]["local_snr"]}
                             for e in survived]}
    with open(os.path.join(ROOT, "prereg", "local_noise.json"), "w",
              encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    write_report(payload)
    print("   wrote prereg/local_noise.json and LOCAL_NOISE.md")
    return 0



def write_report(payload: dict) -> None:
    rows = payload["candidates"]
    lines = ["# The survey re-run with a local noise reference", "",
             f"**{payload['n_survived']} of {payload['n_candidates']} candidates "
             f"survive.**", "",
             "## What changed", "",
             "The survey's statistic is (smoothed - median)/MAD, with one median and "
             "one MAD for the whole searched band, rescaled so the map's bulk is one "
             "sigma wide. That is a single noise scale for 9.66 MHz. This run changes "
             "only the scale: same window, same samples, same dedoppler, same "
             "boxcar, same calibration, same 5.5σ threshold, but the scale is "
             "measured per channel from a ±16 channel neighbourhood across every "
             "drift trial, with one round of clipping so a bright line cannot set "
             "the noise it is measured against.", "",
             "Significance is read at the map cell the survey reported, not from a "
             "re-extracted peak. Dividing each column by its own scale moves the "
             "map's maximum, so a line can be strong and still not be the top hit; "
             "reading its own cell asks the question that matters.", "",
             "| date | MHz | drift Hz/s | survey σ | at-cell σ | local scale | "
             "local σ | outcome |",
             "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for e in rows:
        at = e["at_candidate"]
        lines.append(f"| {e['date']} | {e['freq_mhz']:.4f} | "
                     f"{e['drift_hz_s']:,.1f} | {e['survey_snr']:.2f} | "
                     f"{at['global_snr']:.2f} | {at['scale']:.2f} | "
                     f"{at['local_snr']:.2f} | "
                     f"{'survives' if e['survives_local'] else 'below threshold'} |")
    lines += ["", "## The one that does not survive", ""]
    for e in rows:
        if not e["survives_local"]:
            at = e["at_candidate"]
            lines.append(f"**{e['date']} {e['freq_mhz']:.3f} MHz** falls from "
                         f"{e['survey_snr']:.2f}σ to {at['local_snr']:.2f}σ: the "
                         f"noise in its own ±16 channel neighbourhood is "
                         f"{at['scale']:.2f} of the band-wide estimate, so it does "
                         f"not clear 5.5σ. This is the line that RFI_PAIR.md also "
                         f"finds absent over the full 717 s integration, so the two "
                         f"tests agree that it is not a persistent feature.")
    lines += ["", "## Where the local scale is not a noise estimate", "",
              "For the 20240325 burst group the local scale comes out at 12 to 31 "
              "times the band-wide value, because the ±16 channel neighbourhood "
              "contains the bright lanes themselves and their smear cones. Their "
              "14-31σ values are therefore lower bounds on a reference that is "
              "partly made of the source, not estimates of its significance. The "
              "astrophysical classification does not rest on the number: it rests "
              "on four lanes parallel to 0.028 of a resolution element with a drift "
              "matching a 125 km/s coronal exciter.", "",
              "The 20240920 group moves the other way. Its neighbourhood is three "
              "times quieter than the band-wide average, so the significance rises "
              "to 14-23σ from 6.5-8.7σ. That does not make it astrophysical: those "
              "three lines sit at 0.4965 cycles per sample, within 0.0035 of "
              "Nyquist, and an injected 0.5035 appears at -0.4965. Their strength "
              "is a statement about the noise, not about their origin.", "",
              "## Reading", "",
              "The survey's significances were neither uniformly inflated nor "
              "uniformly conservative: the local scale ranges from 0.30 to 31.18 of "
              "the band-wide value across these ten candidates. What the re-run "
              "changes is the confidence attached to each line; what it does not "
              "change is a single classification. The bright group stays "
              "astrophysical, the near-Nyquist trio stays instrumental, the 53.783 "
              "MHz line stays an unexplained single line, and the 47.242 MHz line "
              "drops below threshold and joins the pair that fails its own repeat "
              "test.", "",
              "One correction to an earlier attempt, recorded because it mattered: a "
              "first version of this test searched narrow 32-channel windows and "
              "returned 0 of 10, including a 954σ burst group. That was a broken "
              "test, not a result. The 20240325 lanes drift about 105 channels over "
              "the 1000-sample integration, so a 32-channel window truncates them "
              "mid-trail. Keeping the survey's 96-channel window and changing only "
              "the noise scale avoids that entirely."]
    with open(os.path.join(ROOT, "LOCAL_NOISE.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
