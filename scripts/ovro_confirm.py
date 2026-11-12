"""Confirmation pipeline for the OVRO-LWA candidates.

Reproduces the two known lines at the exact configuration that found them, places
every candidate from `prereg/ovro_survey.json` through the SETI confirmation
gates, applies the corrections that can be applied without a sky position, runs
injection-recovery for completeness, scores the golden-ratio features against a
randomized control, and writes three reports.

What this pipeline cannot do is as important as what it can, and it says so in the
output rather than in a footnote:

* `spec_fits` is one disk-integrated beam per file, so there is no off-source
  observation to reject against;
* nothing here has a second telescope, so the independent-detection gate cannot
  be attempted;
* the product carries no WCS and no source coordinates, so a barycentric
  correction cannot be applied to a candidate. What can be done is a bound: the
  largest drift Earth's rotation and orbit can induce at this frequency, which is
  a few Hz/s. A candidate at thousands of Hz/s is outside that bound by three
  orders of magnitude, and no choice of sky position brings it back inside.

    python3 scripts/ovro_confirm.py
"""

from __future__ import annotations

import datetime as dt
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

from ovro_lwa_population import dynamic, ensure_file, search_day  # noqa: E402
from ovro_lwa_test import HIGH_MHZ, LOW_MHZ, TSAMP_S  # noqa: E402
from phi_echo import LOCK_TOLERANCE, THETA  # noqa: E402
from seti_search.cadence import SITES, SkyPosition, topocentric_drift  # noqa: E402
from seti_search.filterbank import (  # noqa: E402
    InjectedSignal,
    simulate,
    write_filterbank,
)
from seti_search.rungs import PHI  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

C_LIGHT = 299792458.0
GOLDEN_ANGLE_DEG = 137.50776405
KNOWN = [
    {"date": "20240207", "freq_mhz": 54.452957, "drift_hz_s": 8194.7, "snr": 7.68},
    {"date": "20240207", "freq_mhz": 55.901001, "drift_hz_s": 7866.9, "snr": 7.00},
]
WIDE = (1600, 192)
SURVEY = (1000, 96)


def mjd_from_date(date: str) -> float:
    d = dt.date(int(date[:4]), int(date[4:6]), int(date[6:8]))
    return (d - dt.date(1858, 11, 17)).days + 0.5


def band_range(nfreq: int, chans: int) -> tuple[float, float]:
    channel = (HIGH_MHZ - LOW_MHZ) / nfreq
    c0 = max(0, (nfreq - chans) // 2)
    return (LOW_MHZ + c0 * channel, LOW_MHZ + (c0 + chans) * channel)


def earth_drift_bounds(freq_mhz: float, mjd: float, draws: int = 300) -> dict:
    """Range of drift Earth's own motion induces, over the whole sky."""
    observer = SITES["ovro"]
    rng = np.random.default_rng(1729)
    values = []
    for _ in range(draws):
        ra = float(rng.uniform(0.0, 360.0))
        dec = float(np.degrees(np.arcsin(rng.uniform(-1.0, 1.0))))
        values.append(topocentric_drift(0.0, freq_mhz * 1e6, observer, mjd,
                                        SkyPosition(ra, dec)))
    values = np.array(values)
    return {"min_hz_s": float(values.min()), "max_hz_s": float(values.max()),
            "max_abs_hz_s": float(np.abs(values).max())}


def physical_bounds(freq_mhz: float) -> list[dict]:
    """Drift a moving emitter could produce, for orbits and for coronal plasma."""
    rows = []
    for name, radius_m, speed in (
        ("hot Jupiter, a = 0.05 AU", 7.5e9, 1.4e5),
        ("planet at 1 AU", 1.5e11, 3.0e4),
        ("pulsar planet, r = 1e9 m", 1.0e9, 1.0e5),
        ("ionosphere, v = 100 m/s", 6.4e6, 1.0e2),
    ):
        acc = speed ** 2 / radius_m
        rows.append({"scenario": name, "speed_m_s": speed,
                     "accel_m_s2": acc,
                     "drift_hz_s": acc * freq_mhz * 1e6 / C_LIGHT})
    for name, height_m, speed in (("corona, H = 0.1 R_sun, v = 20 km/s", 6.96e7, 2e4),
                                  ("corona, H = 0.1 R_sun, v = 125 km/s", 6.96e7, 1.25e5)):
        rows.append({"scenario": name, "speed_m_s": speed, "accel_m_s2": None,
                     "drift_hz_s": freq_mhz * 1e6 * speed / (2 * height_m)})
    return rows


def band_note(freq_mhz: float) -> str:
    if 54.0 <= freq_mhz <= 60.0:
        return "inside the 54-60 MHz TV channel 2 allocation (55.25 video, 59.75 audio)"
    if 50.0 <= freq_mhz < 54.0:
        return "inside the 50-54 MHz amateur 6 m allocation"
    if 45.0 <= freq_mhz < 50.0:
        return "below the 6 m band; no allocated broadcast service here"
    return "outside the allocations tabulated for this band"


def phi_features(freq_mhz: float, drift_hz_s: float, channel_hz: float) -> dict:
    nu = abs(drift_hz_s) * TSAMP_S / channel_hz
    # k = 0 is excluded: a stationary line is not a phase advance, and including
    # it lets a near-zero rate score an arbitrarily small "lock distance", which
    # is how a 172 Hz/s line reached 0.0012 and a spurious control p of 0.0265
    k = max(1, round(nu / THETA))
    distance = abs(nu / THETA - k)
    rung = math.log(freq_mhz) / math.log(PHI)
    phase = (rung * 360.0) % 360.0
    separation = min(abs(phase - GOLDEN_ANGLE_DEG), 360 - abs(phase - GOLDEN_ANGLE_DEG))
    return {"rate": nu, "k": int(k), "lock_distance": distance,
            "locked": bool(k != 0 and distance <= LOCK_TOLERANCE),
            "rung": rung, "golden_phase_deg": phase,
            "golden_separation_deg": separation}


def inject_and_search(path: str, nu: float, amp: float, samples: int,
                      chans: int, threshold: float, seed: int) -> dict:
    data, meta = dynamic(path)
    channel_hz, nfreq, ntime = meta["channel_hz"], meta["nfreq"], meta["ntime"]
    ntime_use, nch_use = min(samples, ntime), min(chans, nfreq)
    c0 = max(0, (nfreq - nch_use) // 2)
    block = data[c0:c0 + nch_use, :ntime_use].T.astype(np.float64).copy()
    # the solar spectra are orders of magnitude brighter than unity, so an
    # amplitude of 1.0 is invisible; scale it to the excerpt's own robust noise
    median = float(np.median(block))
    mad = 1.4826 * float(np.median(np.abs(block - median)))
    amplitude = amp * (mad if mad > 0 else 1.0)
    drift = nu * channel_hz / TSAMP_S
    # a rate past Nyquist is not lost, it folds: the search can only report the
    # aliased value, so that is what recovery has to be tested against
    folded = abs(nu) if abs(nu) <= 0.5 else 1.0 - abs(nu)
    expected = folded * channel_hz / TSAMP_S
    t = np.arange(block.shape[0])
    shift = drift * TSAMP_S * t / channel_hz
    start = block.shape[1] / 2.0 + shift[0]
    idx = np.arange(block.shape[1])[None, :] - (start + shift)[:, None]
    block += amplitude * np.exp(-0.5 * (idx / 1.5) ** 2)
    fil = os.path.join(ROOT, "data", "ovro_inject.fil")
    write_filterbank(fil, block.astype(np.float32),
                     fch1=LOW_MHZ + c0 * channel_hz / 1e6,
                     foff=channel_hz / 1e6, tsamp=TSAMP_S,
                     source="SUN + injection", nbits=32, telescope_id=6,
                     tstart=0.0)
    nyquist = 0.5 * channel_hz / TSAMP_S
    cfg = SearchConfig(boxcar=3, max_drift_hz_s=nyquist,
                       snr_threshold=threshold, max_hits=40)
    hits, _meta = search_file(fil, cfg, chan_chunk=2048, max_ints=None)
    step = 0.5 * 3 * channel_hz / (ntime_use * TSAMP_S)
    near = [h for h in hits if abs(h.drift_hz_s - expected) <= 3 * step]
    best = max(near, key=lambda h: h.snr) if near else None
    return {"nu": nu, "injected_hz_s": drift, "expected_hz_s": expected,
            "aliased": bool(abs(nu) > 0.5), "seed": seed,
            "recovered": best is not None,
            "found_hz_s": None if best is None else best.drift_hz_s,
            "snr": None if best is None else best.snr,
            "error_hz_s": None if best is None else abs(best.drift_hz_s - expected)}


def inject_synthetic(nu, amp, samples, chans, threshold, seed) -> dict:
    '''Inject into a synthetic waterfall with unit-variance exponential noise.

    The real solar spectra carry structure that the bandpass and excision stages
    act on, which confounds a completeness measurement; a synthetic file isolates
    the detection pipeline, which is what these numbers are for.
    '''
    channel_hz = 100547.19562243501
    fch1 = 40.0
    freq = fch1 + (chans // 2) * channel_hz / 1e6
    drift = nu * channel_hz / TSAMP_S
    # a rate past Nyquist folds, and it folds in sign: 0.6 cycles per sample
    # appears as 0.6 - 1 = -0.4, not +0.4. Getting that wrong made every aliased
    # injection look unrecovered when it had in fact been found.
    if nu > 0.5:
        expected_nu = nu - 1.0
    elif nu < -0.5:
        expected_nu = nu + 1.0
    else:
        expected_nu = nu
    expected = expected_nu * channel_hz / TSAMP_S
    fil = os.path.join(ROOT, 'data', 'ovro_inject.fil')
    simulate(fil, nchans=chans, nints=samples, tsamp=TSAMP_S, fch1=fch1,
             foff=channel_hz / 1e6, glitch=False, seed=seed,
             signals=(InjectedSignal(freq_mhz=freq, drift_hz_s=drift,
                                     amplitude=amp, width_channels=1.5),))
    cfg = SearchConfig(boxcar=3, max_drift_hz_s=0.5 * channel_hz / TSAMP_S,
                       snr_threshold=threshold, max_hits=40)
    hits, _meta = search_file(fil, cfg, chan_chunk=2048, max_ints=None)
    step = 0.5 * 3 * channel_hz / (samples * TSAMP_S)
    near = [h for h in hits if abs(h.drift_hz_s - expected) <= 3 * step]
    best = max(near, key=lambda h: h.snr) if near else None
    return {'nu': nu, 'injected_hz_s': drift, 'expected_hz_s': expected,
            'expected_nu': expected_nu, 'aliased': bool(abs(nu) > 0.5),
            'seed': seed,
            'recovered': best is not None,
            'found_hz_s': None if best is None else best.drift_hz_s,
            'snr': None if best is None else best.snr,
            'error_hz_s': None if best is None else abs(best.drift_hz_s - expected)}


def noise_only(samples, chans, threshold, seed) -> dict:
    channel_hz = 100547.19562243501
    fil = os.path.join(ROOT, 'data', 'ovro_noise.fil')
    simulate(fil, nchans=chans, nints=samples, tsamp=TSAMP_S, fch1=40.0,
             foff=channel_hz / 1e6, glitch=False, seed=seed, signals=())
    cfg = SearchConfig(boxcar=3, max_drift_hz_s=0.5 * channel_hz / TSAMP_S,
                       snr_threshold=threshold, max_hits=40)
    hits, _meta = search_file(fil, cfg, chan_chunk=2048, max_ints=None)
    return {'seed': seed, 'n_hits': len(hits),
            'max_snr': max((h.snr for h in hits), default=0.0)}


def main() -> int:
    survey = json.load(open(os.path.join(ROOT, "prereg", "ovro_survey.json")))
    rows = survey["candidates"]
    print("CONFIRMATION PIPELINE FOR THE OVRO-LWA CANDIDATES")
    print(f"   {len(rows)} candidates from {survey['days_searched']} days")

    # --- 1. reproduce the two known lines at the configuration that found them
    date = "20240207"
    path = ensure_file(f"spec_fits/2024/{date}.fits")
    repro = search_day(path, WIDE[0], WIDE[1], 5.5)
    matches = []
    for known in KNOWN:
        best = min(repro["rows"], key=lambda r: abs(r["freq_mhz"] - known["freq_mhz"]))
        ok = (abs(best["freq_mhz"] - known["freq_mhz"]) < 0.3
              and abs(best["drift_hz_s"] - known["drift_hz_s"]) < 400.0)
        matches.append({"known": known, "found": best, "reproduced": bool(ok)})
    print(f"   reproduction at {WIDE[0]} x {WIDE[1]}: "
          f"{sum(m['reproduced'] for m in matches)}/2 lines")
    for m in matches:
        f = m["found"]
        print(f"      {m['known']['freq_mhz']:.3f} MHz: found {f['freq_mhz']:.3f} "
              f"MHz, {f['drift_hz_s']:,.1f} Hz/s, snr {f['snr']:.2f} -> "
              f"{'reproduced' if m['reproduced'] else 'not reproduced'}")

    # --- 2. coverage of each configuration, to explain a missing line
    nfreq = repro["nfreq"]
    wide_band = band_range(nfreq, WIDE[1])
    survey_band = band_range(nfreq, SURVEY[1])
    coverage = []
    for known in KNOWN:
        coverage.append({
            "freq_mhz": known["freq_mhz"],
            "in_wide_band": wide_band[0] <= known["freq_mhz"] <= wide_band[1],
            "in_survey_band": survey_band[0] <= known["freq_mhz"] <= survey_band[1],
        })
    print(f"   wide band {wide_band[0]:.2f}-{wide_band[1]:.2f} MHz, "
          f"survey band {survey_band[0]:.2f}-{survey_band[1]:.2f} MHz")

    # --- 3. gates, corrections and phi features per candidate
    # the published pair is carried explicitly: it falls outside the survey band,
    # so a table built only from the survey would silently omit the two lines the
    # build is supposed to confirm or reject
    published = []
    for m in matches:
        published.append({
            "date": date, "freq_mhz": m["known"]["freq_mhz"],
            "drift_hz_s": m["found"]["drift_hz_s"],
            "snr": m["found"]["snr"],
            "resolution_hz_s": repro["resolution_hz_s"],
            "channel_hz": repro["channel_hz"],
            "duration_s": repro["samples_used"] * TSAMP_S,
            "bandwidth_hz": 3.0 * repro["channel_hz"],
            "threshold": 5.5, "source": "published pair"})
    for row in rows:
        row.setdefault("source", "survey")
    entries = rows + published

    enriched = []
    for row in sorted(entries, key=lambda r: (r["date"], -r["snr"])):
        mjd = mjd_from_date(row["date"])
        bounds = earth_drift_bounds(row["freq_mhz"], mjd)
        phys = physical_bounds(row["freq_mhz"])
        dv_dt = abs(row["drift_hz_s"]) * C_LIGHT / (row["freq_mhz"] * 1e6)
        others = [r for r in entries if r is not row and r["date"] != row["date"]
                  and abs(r["freq_mhz"] - row["freq_mhz"]) < 0.3
                  and abs(r["drift_hz_s"] - row["drift_hz_s"]) < 800.0]
        coronal = next(p["drift_hz_s"] for p in phys if p["scenario"].startswith("corona, H = 0.1 R_sun, v = 20"))
        enriched.append({
            "date": row["date"], "freq_mhz": row["freq_mhz"],
            "drift_hz_s": row["drift_hz_s"], "snr": row["snr"],
            "resolution_hz_s": row["resolution_hz_s"],
            "channel_hz": row["channel_hz"], "duration_s": row["duration_s"],
            "bandwidth_hz": row["bandwidth_hz"], "threshold": row["threshold"],
            "earth_drift_max_hz_s": bounds["max_abs_hz_s"],
            "required_accel_m_s2": dv_dt,
            "coronal_20kms_drift_hz_s": coronal,
            "band_note": band_note(row["freq_mhz"]),
            "source": row.get("source", "survey"),
            "repeats_other_days": len(others),
            "repeat_testable": bool(survey_band[0] <= row["freq_mhz"]
                                     <= survey_band[1]),
            "phi": phi_features(row["freq_mhz"], row["drift_hz_s"], row["channel_hz"]),
            "gates": {
                "repeat_testable": bool(survey_band[0] <= row["freq_mhz"]
                                        <= survey_band[1]),
                "repeats_same_position": bool(others),
                "off_source_available": False,
                "independent_telescope": False,
                "barycentric_reducible": bool(bounds["max_abs_hz_s"]
                                              >= abs(row["drift_hz_s"]) * 0.5),
                "physically_possible_orbital": bool(dv_dt < 1.0),
                "explained_by_rfi": bool(54.0 <= row["freq_mhz"] <= 60.0),
            },
            "survey_config": SURVEY,
        })
    print(f"   {sum(1 for e in enriched if e['phi']['locked'])} candidates pass "
          f"the registered lock (k != 0)")

    # --- 4. injections, including an aliased rate beyond Nyquist
    rates = [0.0, 0.02, THETA, 0.6]
    injections = []
    for nu in rates:
        for seed in (11, 12):
            # 384 channels, not 96: a 0.38 cycles/sample chirp crosses 306
            # channels in 800 samples, so a narrow band truncates it and the
            # injection test would measure band width rather than the pipeline
            injections.append(inject_synthetic(nu, 3.0, 800, 384, 5.0, seed))
    fp_runs = [noise_only(800, 384, 5.0, seed) for seed in (21, 22, 23)]
    fp_hits = sum(1 for r in fp_runs if r["n_hits"] > 0)
    recovered = sum(1 for i in injections if i["recovered"])
    print(f"   injections: {recovered}/{len(injections)} recovered "
          f"(rates 0, 0.02, theta, 0.6 folded)")
    by_rate: dict = {}
    for inj in injections:
        by_rate.setdefault(inj["nu"], []).append(inj["recovered"])
    completeness = {f"{nu:.4f}": sum(v) / len(v) for nu, v in by_rate.items()}

    # --- 5. phi ranking against a randomized-phase control
    observed = np.array([e["phi"]["lock_distance"] for e in enriched])
    rng = np.random.default_rng(1729)
    control = []
    for _ in range(2000):
        sample = rng.uniform(0.0, 0.5, size=observed.size)
        k = np.maximum(1.0, np.round(sample / THETA))
        dist = np.abs(sample / THETA - k)
        control.append(float(dist.min()))
    best_observed = float(observed.min())
    p_control = float((1 + np.sum(np.array(control) <= best_observed))
                      / (1 + len(control)))
    print(f"   phi: best lock distance {best_observed:.4f} against a control "
          f"p = {p_control:.4f}")

    # --- 6. verdicts
    for e in enriched:
        g = e["gates"]
        if g["explained_by_rfi"] and e["snr"] > 50:
            e["disposition"] = "ASTROPHYSICAL OR INSTRUMENTAL (bright, in a TV band)"
        elif e["snr"] > 50:
            e["disposition"] = "ASTROPHYSICAL (bright burst group)"
        elif e["phi"]["rate"] > 0.45:
            e["disposition"] = "INSTRUMENTAL ALIASING (rate at Nyquist)"
        elif g["explained_by_rfi"] and g["repeats_same_position"]:
            e["disposition"] = ("TERRESTRIAL RFI (repeats at the same frequency "
                                "and drift)")
        elif g["explained_by_rfi"] and g["repeat_testable"]:
            e["disposition"] = ("UNRESOLVED (in a TV band, but no repeat across "
                                "days)")
        elif g["explained_by_rfi"]:
            e["disposition"] = ("UNRESOLVED (in a TV band; repeat untestable, "
                                "outside survey coverage)")
        else:
            e["disposition"] = "UNRESOLVED"
    headline = "NULL RESULT"
    print(f"   headline verdict: {headline}")

    payload = {"reproduction": matches, "coverage": coverage,
               "wide_band_mhz": wide_band, "survey_band_mhz": survey_band,
               "candidates": enriched,               "injections": injections,
               "completeness": completeness,
               "false_positive_runs": len(fp_runs),
               "false_positive_hits": fp_hits,
               "false_positive_detail": fp_runs,
               "phi_control_p": p_control, "phi_best_distance": best_observed,
               "headline": headline}
    with open(os.path.join(ROOT, "prereg", "ovro_confirm.json"), "w",
              encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")

    write_reports(payload)
    print("   wrote prereg/ovro_confirm.json, REPRODUCIBILITY.md, "
          "CANDIDATES.md, CONFIRMATION.md")
    return 0


def write_reports(payload: dict) -> None:
    repro, cands = payload["reproduction"], payload["candidates"]
    lines = ["# Reproducibility report — OVRO-LWA candidates", "",
             "## Exact-configuration reproduction", "",
             f"Day 20240207 searched at {WIDE[0]} samples x {WIDE[1]} channels, "
             f"5.5 sigma, drift window to Nyquist — the same configuration that "
             f"produced the published pair.", "",
             "| known MHz | known Hz/s | found MHz | found Hz/s | SNR | verdict |",
             "| --- | --- | --- | --- | --- | --- |"]
    for m in repro:
        k, f = m["known"], m["found"]
        lines.append(f"| {k['freq_mhz']:.3f} | {k['drift_hz_s']:,.1f} | "
                     f"{f['freq_mhz']:.4f} | {f['drift_hz_s']:,.1f} | "
                     f"{f['snr']:.2f} | "
                     f"{'reproduced' if m['reproduced'] else 'NOT reproduced'} |")
    lines += ["", "## Why the pair is absent from the survey table", "",
              f"The survey uses {SURVEY[0]} samples x {SURVEY[1]} channels, which "
              f"covers {payload['survey_band_mhz'][0]:.2f}-"
              f"{payload['survey_band_mhz'][1]:.2f} MHz, while the wide run covers "
              f"{payload['wide_band_mhz'][0]:.2f}-{payload['wide_band_mhz'][1]:.2f} "
              f"MHz. Coverage per known line:", "",
              "| line MHz | inside wide band | inside survey band |", "| --- | --- | --- |"]
    for c in payload["coverage"]:
        lines.append(f"| {c['freq_mhz']:.3f} | {c['in_wide_band']} | "
                     f"{c['in_survey_band']} |")
    lines += ["", "A line outside a configuration's band cannot appear in it; that "
                  "is a coverage effect and not a failure of reproducibility. The "
                  "45-49 MHz survey band is a different part of the spectrum, so "
                  "the two tables are not directly comparable row to row."]
    with open(os.path.join(ROOT, "REPRODUCIBILITY.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    lines = ["# Candidate table — raw and corrected", "",
             "Every candidate from `prereg/ovro_survey.json`, with the corrections "
             "that can be applied without a sky position. `earth max` is the "
             "largest drift Earth's rotation and orbit can induce at that "
             "frequency, over the whole sky, for that date; `dv/dt` is the "
             "acceleration a transmitter would need to produce the measured drift; "
             "`off` is whether an off-source observation exists in this product "
             "(it does not).", "",
             "| date | MHz | drift Hz/s | SNR | res Hz/s | earth max Hz/s | dv/dt m/s² | "
             "band | repeats | off | k | lock dist | φ phase° | disposition |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for e in cands:
        lines.append(
            f"| {e['date']} | {e['freq_mhz']:.4f} | {e['drift_hz_s']:,.1f} | "
            f"{e['snr']:.2f} | {e['resolution_hz_s']:.1f} | "
            f"{e['earth_drift_max_hz_s']:.2f} | {e['required_accel_m_s2']:,.0f} | "
            f"{e['band_note'].split('(')[0].strip()} | {e['repeats_other_days']} | "
            f"none | {e['phi']['k']} | {e['phi']['lock_distance']:.3f} | "
            f"{e['phi']['golden_phase_deg']:.1f} | {e['disposition']} |")
    lines += ["", f"Phi control: best lock distance {payload['phi_best_distance']:.4f} "
                  f"against 2,000 randomized-phase draws, p = "
                  f"{payload['phi_control_p']:.4f}. The phi columns are a "
                  f"preregistered feature, not evidence."]
    with open(os.path.join(ROOT, "CANDIDATES.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    gates = ["repeat_testable", "repeats_same_position", "off_source_available",
             "independent_telescope", "barycentric_reducible",
             "physically_possible_orbital", "explained_by_rfi"]
    lines = ["# Confirmation report", "",
             f"**Headline: {payload['headline']}** — no candidate passes the "
             "confirmation gates, so no signal is claimed as evidence of "
             "intelligent life.", "",
             "## Gate matrix", "",
             "| date | MHz | " + " | ".join(gates) + " |",
             "| --- | --- | " + " | ".join("---" for _ in gates) + " |"]
    for e in cands:
        cells = " | ".join("pass" if e["gates"][g] else "fail" for g in gates)
        lines.append(f"| {e['date']} | {e['freq_mhz']:.4f} | {cells} |")
    lines += ["", "Three gates cannot be attempted with this product at all, which "
                  "is the central result: `off_source_available` and "
                  "`independent_telescope` are structurally false for one "
                  "disk-integrated beam from one array, and `barycentric_reducible` "
                  "is false for every candidate whose drift exceeds the few Hz/s "
                  "that Earth's motion can supply.", "",
                  "## Injection-recovery", "",
                  "| rate cycles/sample | recovered |", "| --- | --- |"]
    for rate, frac in sorted(payload["completeness"].items()):
        lines.append(f"| {rate} | {frac:.2f} |")
    lines += ["", "Rates 0 and 0.02 test stationary and slow lines, theta tests a "
                  "lock-rate line, and 0.6 folds past Nyquist to test aliasing "
                  "handling. Amplitude is fixed, so these are completeness at one "
                  "signal strength, not a curve.", "",
                  f"Noise-only runs: {payload['false_positive_runs']} trials at the "
                  f"same threshold, {payload['false_positive_hits']} with any hit.", "",
                  "## Dispositions", ""]
    for e in cands:
        lines.append(f"- **{e['date']} {e['freq_mhz']:.3f} MHz, "
                     f"{e['drift_hz_s']:,.0f} Hz/s, SNR {e['snr']:.1f}** — "
                     f"{e['disposition']}; {e['band_note']}; nearest phi integer "
                     f"k = {e['phi']['k']} at distance "
                     f"{e['phi']['lock_distance']:.3f}.")
    lines += ["", "## What would change this", "",
              "An off-source observation, a second telescope, and a candidate whose "
              "drift survives barycentric correction. None of the three exists in "
              "this data set, so the honest statement is that the archive contains "
              "drifting lines and none of them is a technosignature candidate."]
    with open(os.path.join(ROOT, "CONFIRMATION.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
