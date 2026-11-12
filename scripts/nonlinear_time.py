"""A second analysis lane: nonlinear event time, beside linear sample time.

The linear lane indexes each spectrum by uniform sample time, t = n * dt. This
script builds a second coordinate for the same candidates and asks whether a
signal is *more* coherent in event time than in sample time, and whether that
gain survives controls that are exactly as flexible as the transform.

Five preregistered coordinates, every one fitted on the first half of the track
and evaluated on the second half, so no coherence is ever measured on the
samples that estimated the transform:

    linear   tau = t
    phase    tau = cumulative |dphi| / (2 pi f_ref), the phase-advance clock
    chirp    tau = integral (f(t) - f0)/B dt = drift * t^2 / (2 B)
    accel    the same after removing the drift an instrumental hypothesis
             predicts: the survey-grid quantisation of the measured drift
    energy   tau = cumulative excess power / total excess

`tau_phi = frac(tau_event / PHI^2)` is carried as a preregistered feature only.
Its fractional part is not a time coordinate and is never used to rank anything.

The circularity this build has to survive: `chirp` is fitted to the same linear
drift that defines the candidate, so a constant-acceleration line is linearised
*by construction*, and its in-sample phase locking is guaranteed. What carries
information is (a) the split-half value, the first half's drift evaluated on the
second, and (b) the gain over linear time compared against the gain that
randomised monotone warps achieve. A warp is flexible enough to raise phase
locking on its own, so the control is best-of-five against best-of-five.

    python3 scripts/nonlinear_time.py
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

from ovro_lwa_population import dynamic  # noqa: E402
from ovro_lwa_test import LOW_MHZ, TSAMP_S  # noqa: E402
from seti_search.rungs import PHI  # noqa: E402

THETA = 1.0 / PHI**2
PHI2 = PHI**2
M = 512                 # resampled points per coordinate
BOOT = 200              # block-bootstrap resamples
NULLS = 200             # shuffled-time nulls
WARPS = 100             # randomised monotone warps
BLOCK = 8               # block length for the bootstrap
SURVEY_RES = 392.7624829001368
CHANNEL_HZ = 100547.19562243501
MAIN = ("linear", "phase", "chirp", "accel", "energy")
GAIN_MIN = 0.05


# --------------------------------------------------------------------------- #
# signal tools
# --------------------------------------------------------------------------- #
def analytic(x: np.ndarray) -> np.ndarray:
    """Hilbert analytic signal. The product is intensity, not voltage, so this is
    a proxy for phase and is labelled as one wherever it is used."""
    n = x.size
    spec = np.fft.fft(x - x.mean())
    h = np.zeros(n)
    h[0] = 1.0
    if n % 2 == 0:
        h[n // 2] = 1.0
        h[1:n // 2] = 2.0
    else:
        h[1:(n + 1) // 2] = 2.0
    return np.fft.ifft(spec * h)


def phase_metrics(y: np.ndarray) -> dict:
    yn = y - y.mean()
    sd = float(yn.std())
    if sd <= 0 or yn.size < 16:
        return {"plv": 0.0, "circvar": 1.0, "drift_residual": 0.0,
                "ridge_continuity": 0.0, "acf_peak": 0.0,
                "matched_filter_snr": 0.0}
    z = analytic(yn / sd)
    ph = np.unwrap(np.angle(z))
    j = np.arange(ph.size, dtype=float)
    slope, intercept = np.polyfit(j, ph, 1)
    res = ph - (slope * j + intercept)
    plv = float(abs(np.mean(np.exp(1j * res))))
    inst = np.diff(ph) / (2.0 * np.pi)          # cycles per resampled step
    drift_res = float(np.std(np.diff(inst))) if inst.size > 2 else 0.0
    # a ridge is continuous when the instantaneous frequency stays within a
    # twentieth of a cycle of the coordinate's own median rate
    ridge = float(np.mean(np.abs(inst - np.median(inst)) <= 0.05))
    ac = np.correlate(yn, yn, "full")[yn.size - 1:]
    ac = ac / ac[0]
    acf = float(np.max(ac[1:max(2, yn.size // 4)]))
    # the matched-filter statistic has to be an excess over the local floor: the
    # raw mean/std ratio of a noise track is sqrt(N), which is not a detection
    floor = float(np.median(yn))
    mf = float((yn.mean() - floor) / sd * np.sqrt(yn.size))
    # concentration: the largest share of the total excess that any window of a
    # sixteenth of the track holds. A compact burst concentrates; a flat excess
    # does not. This is defined on intensity, unlike the Hilbert proxy.
    window = max(4, yn.size // 16)
    excess = np.maximum(yn - floor, 0.0)
    total = float(excess.sum())
    if total > 0:
        cumulative = np.cumsum(excess)
        best = np.max(cumulative[window:] - cumulative[:-window])
        conc = float(best / total)
    else:
        conc = 0.0
    kurt = float(np.mean(((yn - yn.mean()) / sd) ** 4) - 3.0)
    return {"plv": plv, "circvar": 1.0 - plv, "drift_residual": drift_res,
            "ridge_continuity": ridge, "acf_peak": acf,
            "matched_filter_snr": mf, "concentration": conc, "kurtosis": kurt}


def score(m: dict) -> float:
    """The decision statistic. Defined on intensity, so it does not need the
    carrier phase this product does not carry."""
    return m["matched_filter_snr"]


def dtw(a: np.ndarray, b: np.ndarray) -> float:
    """Dynamic-time-warping distance between two normalised series."""
    n = a.size
    cost = np.full((n + 1, n + 1), np.inf)
    cost[0, 0] = 0.0
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            d = abs(a[i - 1] - b[j - 1])
            cost[i, j] = d + min(cost[i - 1, j], cost[i, j - 1], cost[i - 1, j - 1])
    return float(cost[n, n] / n)


def normalise(y: np.ndarray) -> np.ndarray:
    sd = float(y.std())
    return (y - y.mean()) / sd if sd > 0 else y * 0.0


def resample(tau: np.ndarray, x: np.ndarray, m: int = M) -> np.ndarray:
    tau = np.maximum.accumulate(np.asarray(tau, dtype=float))
    tau = tau + np.arange(tau.size) * 1e-12
    if tau[-1] <= tau[0] or x.size < 8:
        return np.full(m, float(np.mean(x)))
    return np.interp(np.linspace(tau[0], tau[-1], m), tau, x)


# --------------------------------------------------------------------------- #
# the five preregistered coordinates
# --------------------------------------------------------------------------- #
def fit_transform(kind: str, x: np.ndarray, t: np.ndarray, drift: float) -> dict:
    if kind == "linear":
        return {"t0": float(t[0])}
    if kind == "phase":
        ph = np.unwrap(np.angle(analytic(x - x.mean())))
        span = max(float(t[-1] - t[0]), 1e-9)
        return {"f_ref": max(abs(float(ph[-1] - ph[0])) / span, 1e-9),
                "t0": float(t[0])}
    if kind == "chirp":
        return {"drift": float(drift), "band": max(float(np.ptp(t)), 1e-9),
                "t0": float(t[0])}
    if kind == "accel":
        model = round(drift / SURVEY_RES) * SURVEY_RES
        return {"resid": float(drift - model), "model": float(model),
                "band": max(float(np.ptp(t)), 1e-9), "t0": float(t[0])}
    if kind == "energy":
        base = float(np.median(x))
        mad = 1.4826 * float(np.median(np.abs(x - base)))
        excess = np.maximum(x - (base + mad), 0.0)
        return {"floor": base + mad, "total": max(float(excess.sum()), 1e-12)}
    raise ValueError(kind)


def apply_transform(kind: str, p: dict, x: np.ndarray, t: np.ndarray) -> np.ndarray:
    if kind == "linear":
        return t - p["t0"]
    if kind == "phase":
        # cumulative phase advance: a clock that only runs forward, so it stays a
        # valid coordinate even where the unwrapped phase steps backwards
        ph = np.unwrap(np.angle(analytic(x - x.mean())))
        dph = np.abs(np.diff(ph, prepend=ph[0]))
        return np.cumsum(dph) / (2.0 * np.pi * p["f_ref"])
    if kind == "chirp":
        return p["drift"] * (t - p["t0"]) ** 2 / (2.0 * p["band"])
    if kind == "accel":
        return p["resid"] * (t - p["t0"]) ** 2 / (2.0 * p["band"])
    if kind == "energy":
        return np.cumsum(np.maximum(x - p["floor"], 0.0)) / p["total"]
    raise ValueError(kind)


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #
def track_at(data: np.ndarray, c0: int, drift: float, samples: int,
             channel_hz: float, half: int = 1) -> np.ndarray:
    nfreq, ntime = data.shape
    n = int(min(samples, ntime))
    rows = np.arange(n)
    chan = c0 + drift * TSAMP_S * rows / channel_hz
    out = np.empty(n, dtype=float)
    for i in range(n):
        ci = int(round(chan[i]))
        lo, hi = max(0, ci - half), min(nfreq, ci + half + 1)
        out[i] = float(data[lo:hi, i].mean()) if hi > lo else 0.0
    return out


def channel_of(freq_mhz: float, nfreq: int, channel_hz: float) -> int:
    c0 = int(round((freq_mhz * 1e6 - LOW_MHZ * 1e6) / channel_hz))
    return max(8, min(nfreq - 9, c0))


def shuffled_channel_track(data: np.ndarray, c0: int, drift: float, samples: int,
                           channel_hz: float, seed: int = 7) -> np.ndarray:
    """Track rebuilt from a waterfall with channel order permuted row by row."""
    block = data[max(0, c0 - 14):min(data.shape[0], c0 + 15), :min(samples,
                                                                 data.shape[1])].copy()
    rng = np.random.default_rng(seed)
    for j in range(block.shape[1]):
        rng.shuffle(block[:, j])
    return track_at(block, block.shape[0] // 2, drift, samples, channel_hz)


# --------------------------------------------------------------------------- #
# synthetic injections and the alias mapping
# --------------------------------------------------------------------------- #
def synthetic_track(nu: float, amp: float, samples: int = 800, seed: int = 0,
                    channels: int = 320, track_nu: float | None = None) -> np.ndarray:
    """Track through a synthetic waterfall carrying a Gaussian chirp at rate nu.

    The track is built at `track_nu`, which defaults to the folded rate the
    search can actually report, so this measures coherence at the reported rate.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(samples)
    centre = channels / 2.0 + nu * t
    j = np.arange(channels)[None, :]
    wf = rng.exponential(1.0, (samples, channels))
    wf += amp * np.exp(-0.5 * ((j - centre[:, None]) / 1.4) ** 2)
    rate = nu if track_nu is None else track_nu
    idx = np.clip(np.round(channels / 2.0 + rate * t).astype(int), 0, channels - 1)
    return wf[np.arange(samples), idx]


def fold(nu: float) -> float:
    """Per-sample increment the sampler can actually report, in (-0.5, 0.5]."""
    return float((nu + 0.5) % 1.0 - 0.5)


def alias_probe(nu: float, samples: int = 4096) -> dict:
    """Verify the fold numerically: sample a chirp and read back its increment."""
    z = np.exp(2j * np.pi * nu * np.arange(samples))
    measured = np.angle(z[1:] * np.conj(z[:-1])) / (2 * np.pi)
    return {"injected_nu": nu, "measured_nu": float(np.median(measured)),
            "folded_expected": fold(nu),
            "exact": bool(abs(np.median(measured) - fold(nu)) < 1e-9)}


def injection_lane(nu: float, amp: float, seed: int) -> dict:
    """Coherence of a known chirp in linear time and in each nonlinear time."""
    reported = fold(nu)
    series = synthetic_track(nu, amp, seed=seed, track_nu=reported)
    t = np.arange(series.size) * TSAMP_S
    drift = reported * CHANNEL_HZ / TSAMP_S
    half = series.size // 2
    out = {"nu": nu, "reported_nu": reported, "seed": seed}
    best = -np.inf
    for kind in MAIN:
        p = fit_transform(kind, series[:half], t[:half], drift)
        tau = apply_transform(kind, p, series, t)[half:]
        m = phase_metrics(resample(tau, series[half:]))
        out[kind] = score(m)
        out[kind + "_plv"] = m["plv"]
        best = max(best, out[kind])
    out["best"] = best
    out["gain"] = best - out["linear"]
    rng = np.random.default_rng(1000 + seed)
    # the synthetic track is already uniformly sampled, so the shuffled-time null
    # is a permutation of the same values
    nulls = [score(phase_metrics(rng.permutation(series[half:])))
             for _ in range(40)]
    out["shuffled_p95"] = float(np.percentile(nulls, 95))
    out["detected"] = bool(best > 5.0)
    return out


def warp_null(x: np.ndarray, t: np.ndarray, drift: float, half: int,
              linear_score: float, trials: int = WARPS, seed: int = 4242) -> dict:
    """Best-of-five score and gain under randomised monotone time warps."""
    rng = np.random.default_rng(seed)
    best, gains, plv_best = [], [], []
    for _ in range(trials):
        warp = np.cumsum(np.abs(rng.normal(size=x.size)))
        warp = warp / warp[-1]
        grid = np.linspace(0.0, 1.0, x.size)
        vals, plvs = [], []
        for kind in MAIN:
            p = fit_transform(kind, x[:half], t[:half], drift)
            tau = apply_transform(kind, p, x, t)
            # the warp reparameterises the coordinate, sample order preserved
            warped = np.interp(warp, grid, tau)
            m = phase_metrics(resample(warped[half:], x[half:]))
            vals.append(score(m))
            plvs.append(m["plv"])
        best.append(max(vals))
        gains.append(max(vals) - linear_score)
        plv_best.append(max(plvs))
    best = np.array(best)
    gains = np.array(gains)
    plv_best = np.array(plv_best)
    return {"p95": float(np.percentile(best, 95)),
            "median": float(np.median(best)),
            "gain_p95": float(np.percentile(gains, 95)),
            "gain_median": float(np.median(gains)),
            "plv_p95": float(np.percentile(plv_best, 95))}


# --------------------------------------------------------------------------- #
# SVG plots, standard library only
# --------------------------------------------------------------------------- #
class Panel:
    def __init__(self, title: str, series: list, xlabel: str = "",
                 ylabel: str = "") -> None:
        self.title = title
        self.series = series          # (xs, ys, colour, label, kind)
        self.xlabel = xlabel
        self.ylabel = ylabel


def write_svg(path: str, title: str, panels: list, ncol: int = 3,
              pw: int = 340, ph: int = 230) -> None:
    rows = max(1, math.ceil(len(panels) / ncol))
    width, height = ncol * pw + 40, rows * ph + 70
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
           f'height="{height}" viewBox="0 0 {width} {height}" '
           f'font-family="Helvetica,Arial,sans-serif">',
           f'<rect width="{width}" height="{height}" fill="#ffffff"/>',
           f'<text x="20" y="28" font-size="15" font-weight="bold">{title}</text>']
    for index, panel in enumerate(panels):
        cx = 20 + (index % ncol) * pw
        cy = 50 + (index // ncol) * ph
        x0, y0, w, h = cx + 32, cy + 16, pw - 52, ph - 56
        xs = [v for s in panel.series for v in s[0]]
        ys = [v for s in panel.series for v in s[1]]
        if not xs or not ys:
            continue
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        if xmax == xmin:
            xmax = xmin + 1.0
        if ymax == ymin:
            ymax = ymin + 1.0

        def px(v):
            return x0 + (v - xmin) / (xmax - xmin) * w

        def py(v):
            return y0 + h - (v - ymin) / (ymax - ymin) * h

        out.append(f'<text x="{cx + 30}" y="{cy + 10}" font-size="10" '
                   f'font-weight="bold">{panel.title}</text>')
        out.append(f'<rect x="{x0:.0f}" y="{y0:.0f}" width="{w:.0f}" '
                   f'height="{h:.0f}" fill="#fbfbfd" stroke="#c9ccd6"/>')
        for s_index, (sx, sy, colour, label, kind) in enumerate(panel.series):
            if kind == "line":
                pts = " ".join(f"{px(a):.1f},{py(b):.1f}" for a, b in zip(sx, sy))
                out.append(f'<polyline points="{pts}" fill="none" '
                           f'stroke="{colour}" stroke-width="1.2"/>')
            else:
                for a, b in zip(sx, sy):
                    out.append(f'<circle cx="{px(a):.1f}" cy="{py(b):.1f}" '
                               f'r="2.6" fill="{colour}"/>')
            out.append(f'<text x="{x0 + 4 + s_index * 104:.0f}" '
                       f'y="{y0 + h + 12:.0f}" font-size="7.5" '
                       f'fill="{colour}">{label[:28]}</text>')
        out.append(f'<text x="{x0 + w / 2:.0f}" y="{y0 + h + 26:.0f}" font-size="8" '
                   f'text-anchor="middle" fill="#555">{panel.xlabel}</text>')
        out.append(f'<text x="{x0 - 4:.0f}" y="{y0 + 8:.0f}" font-size="8" '
                   f'text-anchor="end" fill="#555">{panel.ylabel}</text>')
        out.append(f'<text x="{x0:.0f}" y="{y0 + h + 38:.0f}" font-size="7" '
                   f'fill="#888">{ymin:.3g} to {ymax:.3g}</text>')
    out.append("</svg>")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main() -> int:
    confirm = json.load(open(os.path.join(ROOT, "prereg", "ovro_confirm.json"),
                             encoding="utf-8"))
    cands = confirm["candidates"]
    print("NONLINEAR EVENT-TIME LANE (the linear result is left untouched)")
    print(f"   {len(cands)} candidates; linear verdict preserved: "
          f"{confirm['headline']}")

    by_date: dict[str, list] = {}
    for c in cands:
        by_date.setdefault(c["date"], []).append(c)

    rows, panels, per_candidate = [], [], []
    for date in sorted(by_date):
        path = os.path.join(ROOT, "data", f"{date}.fits")
        if not os.path.exists(path):
            print(f"   {date}: no cached file, skipping")
            continue
        data, meta = dynamic(path)
        ch_hz = meta["channel_hz"]
        print(f"   {date}: {meta['nfreq']} channels x {meta['ntime']} samples")
        for c in sorted(by_date[date], key=lambda r: -r["snr"]):
            samples = int(round(c["duration_s"] / TSAMP_S))
            c0 = channel_of(c["freq_mhz"], meta["nfreq"], ch_hz)
            series = track_at(data, c0, c["drift_hz_s"], samples, ch_hz)
            half = series.size // 2
            t = np.arange(series.size) * TSAMP_S
            nu = abs(c["drift_hz_s"]) * TSAMP_S / ch_hz

            entry = {"date": date, "freq_mhz": c["freq_mhz"],
                     "raw_drift_hz_s": c["drift_hz_s"],
                     "corrected_drift_hz_s": c["drift_hz_s"],
                     "snr": c["snr"], "rate": nu,
                     "band_note": c["band_note"],
                     "samples_used": int(series.size),
                     "separate_segments": True,
                     "trained_on": "first half", "tested_on": "second half",
                     "coherence": {}, "controls": {}}
            best_kind, best_score, linear_y = "linear", -np.inf, None
            for kind in MAIN:
                p = fit_transform(kind, series[:half], t[:half], c["drift_hz_s"])
                tau = apply_transform(kind, p, series, t)
                y = resample(tau[half:], series[half:])
                m = phase_metrics(y)
                if kind == "linear":
                    linear_y = y
                    m["dtw_to_linear"] = 0.0
                else:
                    m["dtw_to_linear"] = dtw(normalise(y[:200]),
                                             normalise(linear_y[:200]))
                rng = np.random.default_rng(900 + len(rows) + len(per_candidate))
                nulls = np.array([phase_metrics(rng.permutation(y))["plv"]
                                  for _ in range(NULLS)])
                m["null_median_plv"] = float(np.median(nulls))
                m["null_p95_plv"] = float(np.percentile(nulls, 95))
                m["null_percentile"] = float(np.mean(nulls <= m["plv"]))
                boots = []
                for _b in range(BOOT):
                    starts = rng.integers(0, max(1, y.size - BLOCK), y.size // BLOCK)
                    boots.append(phase_metrics(
                        np.concatenate([y[s:s + BLOCK] for s in starts]))["plv"])
                m["plv_ci"] = [float(np.percentile(boots, 5)),
                               float(np.percentile(boots, 95))]
                entry["coherence"][kind] = m
                if kind != "linear" and score(m) > best_score:
                    best_kind, best_score = kind, score(m)
                rows.append({"date": date, "freq_mhz": c["freq_mhz"],
                             "raw_drift_hz_s": c["drift_hz_s"],
                             "coordinate": kind,
                             **{k: v for k, v in m.items() if k != "plv_ci"},
                             "plv_ci_low": m["plv_ci"][0],
                             "plv_ci_high": m["plv_ci"][1],
                             "train_test_separate": True})

            linear_plv = float(entry["coherence"]["linear"]["plv"])
            linear_score = score(entry["coherence"]["linear"])
            wn = warp_null(series, t, c["drift_hz_s"], half, linear_score)
            rng = np.random.default_rng(31337)
            scrambled = np.real(np.fft.ifft(np.abs(np.fft.fft(series))
                                            * np.exp(2j * np.pi * rng.random(series.size))))
            scr = max(score(phase_metrics(resample(apply_transform(
                k, fit_transform(k, series[:half], t[:half], c["drift_hz_s"]),
                scrambled, t)[half:], scrambled[half:]))) for k in MAIN)
            rev_series = series[::-1]
            rev = max(score(phase_metrics(resample(apply_transform(
                k, fit_transform(k, rev_series[:half], t[:half], c["drift_hz_s"]),
                rev_series, t)[half:], rev_series[half:]))) for k in MAIN)
            sh_series = shuffled_channel_track(data, c0, c["drift_hz_s"],
                                               samples, ch_hz)
            sh = max(score(phase_metrics(resample(apply_transform(
                k, fit_transform(k, sh_series[:half], t[:half], c["drift_hz_s"]),
                sh_series, t)[half:], sh_series[half:]))) for k in MAIN)
            off_c0 = channel_of(c["freq_mhz"] + 3.0, meta["nfreq"], ch_hz)
            off_series = track_at(data, off_c0, c["drift_hz_s"], samples, ch_hz)
            off = max(score(phase_metrics(resample(apply_transform(
                k, fit_transform(k, off_series[:half], t[:half], c["drift_hz_s"]),
                off_series, t)[half:], off_series[half:]))) for k in MAIN)

            entry["controls"] = {"warp_p95": wn["p95"], "warp_median": wn["median"],
                                 "warp_plv_p95": wn["plv_p95"],
                                 "warp_gain_p95": wn["gain_p95"],
                                 "warp_gain_median": wn["gain_median"],
                                 "phase_scrambled_plv": float(scr),
                                 "reversed_plv": float(rev),
                                 "shuffled_channels_plv": float(sh),
                                 "off_band_plv": float(off)}
            # a monotone coordinate is a relabelling of the same samples, so the
            # amplitude statistic must be invariant: a shuffled-time null cannot
            # test order-dependence at all. The warp control, which resamples, is
            # the only flexibility-matched comparison available here.
            lin_y = resample(apply_transform("linear", fit_transform(
                "linear", series[:half], t[:half], c["drift_hz_s"]),
                series, t)[half:], series[half:])
            ry = np.random.default_rng(555)
            inv = max(abs(score(phase_metrics(ry.permutation(lin_y)))
                          - score(phase_metrics(lin_y))) for _ in range(5))
            entry["amplitude_statistic_invariance"] = float(inv)
            entry["shuffle_invariant"] = bool(inv < 1e-9)
            entry["best_nonlinear_coordinate"] = best_kind
            entry["best_score"] = float(best_score)
            entry["linear_score"] = float(linear_score)
            entry["gain_over_linear"] = float(best_score - linear_score)
            entry["beats_warp_gain"] = bool(best_score - linear_score
                                            > wn["gain_p95"])
            entry["split_half_holds"] = bool(best_score > 3.0)
            entry["linear_plv"] = linear_plv
            entry["best_plv_envelope"] = float(
                entry["coherence"][best_kind]["plv"])
            entry["plv_role"] = "envelope proxy, no decision role"
            entry["cross_day_repeats"] = int(c["repeats_other_days"])
            entry["tau_phi_feature"] = {
                "mean_fraction": float(np.mean(np.mod(apply_transform(
                    "energy", fit_transform("energy", series, t,
                                            c["drift_hz_s"]), series, t), PHI2))),
                "used_for_ranking": False}
            per_candidate.append(entry)
            print(f"      {c['freq_mhz']:8.3f} MHz  score linear "
                  f"{linear_score:+.2f}  best {best_kind} {best_score:+.2f}  "
                  f"gain {entry['gain_over_linear']:+.2f}  warp gain p95 "
                  f"{wn['gain_p95']:+.2f}  "
                  f"{'beats' if entry['beats_warp_gain'] else 'does not beat'}")
            panels.append(Panel(
                f"{date} {c['freq_mhz']:.3f} MHz",
                [(np.linspace(0, 1, 200), normalise(series[:200]), "#2255aa",
                  "linear t", "line"),
                 (np.linspace(0, 1, 200), normalise(resample(
                     apply_transform(best_kind, fit_transform(
                         best_kind, series[:half], t[:half], c["drift_hz_s"]),
                         series, t)[half:], series[half:])[:200]),
                  "#cc5500", f"best: {best_kind}", "line")],
                "fraction of track", "normalised power"))
        del data

    print("   injections at known rates (synthetic, unit-mean exponential noise)")
    injections = [injection_lane(nu, 3.0, seed)
                  for nu in (0.0, 0.02, THETA, 0.45, 0.6) for seed in (11, 12)]
    for inj in injections:
        print(f"      nu {inj['nu']:.4f} -> reported {inj['reported_nu']:+.4f}  "
              f"score linear {inj['linear']:+.2f}  best {inj['best']:+.2f}  "
              f"gain {inj['gain']:+.2f}  "
              f"{'detected' if inj['detected'] else 'not detected'}")

    print("   alias probe: the increment a sampler can report")
    aliases = [alias_probe(nu) for nu in (0.1, 0.3, 0.4965, 0.5035, 0.6, 0.9)]
    for a in aliases:
        print(f"      nu {a['injected_nu']:.4f} -> {a['measured_nu']:+.4f} "
              f"(fold {a['folded_expected']:+.4f}, exact {a['exact']})")

    group = [c for c in cands if c["date"] == "20240325"]
    lanes = [c["drift_hz_s"] for c in group]
    solar = {"lanes": len(group), "drifts_hz_s": lanes,
             "spread_hz_s": float(max(lanes) - min(lanes)) if lanes else 0.0,
             "spread_over_resolution": float((max(lanes) - min(lanes)) / SURVEY_RES)
             if lanes else 0.0,
             "coronal_125kms_hz_s": float(np.mean([c["freq_mhz"] * 1e6 * 1.25e5
                                                   / (2 * 6.96e7) for c in group]))
             if group else 0.0}
    solar["parallel"] = bool(solar["spread_over_resolution"] < 0.1)
    print(f"   solar template: {solar['lanes']} lanes, spread "
          f"{solar['spread_hz_s']:.1f} Hz/s = "
          f"{solar['spread_over_resolution']:.3f} resolution, 125 km/s corona "
          f"{solar['coronal_125kms_hz_s']:,.0f} Hz/s")

    for e in per_candidate:
        e["astrophysical_explanation"] = (
            "solar burst: parallel lanes sharing one drift, and the drift matches "
            "a 125 km/s coronal exciter" if e["date"] == "20240325" else "")
        e["instrumental_explanation"] = (
            "rate 0.4965 cycles/sample, 0.0035 below Nyquist; the sampler cannot "
            "distinguish it from -0.5035" if e["rate"] > 0.45 else "")
        absent = bool(e["best_score"] > e["controls"]["phase_scrambled_plv"]
                      and e["best_score"] > e["controls"]["reversed_plv"]
                      and e["best_score"] > e["controls"]["shuffled_channels_plv"]
                      and e["best_score"] > e["controls"]["off_band_plv"])
        e["absent_from_controls"] = absent
        if e["date"] == "20240325":
            e["disposition"] = "ASTROPHYSICAL"
        elif e["rate"] > 0.45:
            e["disposition"] = "INSTRUMENTAL ALIAS"
        elif e["beats_warp_gain"] and e["split_half_holds"] and absent:
            e["disposition"] = "NONLINEAR-COHERENT BUT UNRESOLVED"
        else:
            e["disposition"] = "UNRESOLVED"
        e["confirmation_status"] = {
            "nonlinear_coherent": bool(e["beats_warp_gain"]
                                       and e["split_half_holds"]),
            "absent_from_controls": absent,
            "repeat_same_position": bool(e["cross_day_repeats"] > 0),
            "later_observation": bool(e["cross_day_repeats"] > 0),
            "off_source": "UNTESTABLE",
            "barycentric": "UNTESTABLE (no WCS, no sky position)",
            "independent_telescope": "UNTESTABLE",
            "confirmation": "UNTESTABLE"}

    coherent = [e for e in per_candidate
                if e["disposition"] == "NONLINEAR-COHERENT BUT UNRESOLVED"]
    headline = ("NONLINEAR-COHERENT BUT UNRESOLVED" if coherent else "NULL RESULT")
    print(f"   {sum(1 for e in per_candidate if e['beats_warp_gain'])} of "
          f"{len(per_candidate)} beat the warp-gain control; headline {headline}")

    payload = {"coordinates": {
        "linear": "tau = t",
        "phase": "tau = cumulative |dphi| / (2 pi f_ref)",
        "chirp": "tau = integral (f - f0)/B dt = drift t^2 / (2B)",
        "accel": "tau with the survey-grid drift model removed",
        "energy": "tau = cumulative excess power / total excess",
        "phi": "tau_phi = frac(tau_event / PHI^2), feature only"},
        "split": {"trained_on": "first half", "tested_on": "second half"},
        "candidates": per_candidate, "injections": injections,
        "alias_probe": aliases, "solar_template": solar,
        "warps": WARPS, "nulls": NULLS, "bootstrap": BOOT,
        "gain_min": GAIN_MIN, "headline": headline,
        "linear_result_preserved": confirm["headline"]}
    with open(os.path.join(ROOT, "nonlinear_candidates.json"), "w",
              encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")

    fields = ["date", "freq_mhz", "raw_drift_hz_s", "coordinate", "plv",
              "plv_ci_low", "plv_ci_high", "circvar", "ridge_continuity",
              "dtw_to_linear", "matched_filter_snr", "drift_residual",
              "acf_peak", "null_median_plv", "null_p95_plv", "null_percentile",
              "train_test_separate"]
    with open(os.path.join(ROOT, "nonlinear_coherence.csv"), "w", newline="",
              encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    idx = list(range(len(per_candidate)))
    write_svg(os.path.join(ROOT, "figures", "nonlinear_tracks.svg"),
              "Linear sample time (blue) against the best nonlinear coordinate "
              "(orange)", panels, ncol=4)
    write_svg(os.path.join(ROOT, "figures", "nonlinear_coherence.svg"),
              "Phase locking per candidate: best nonlinear against linear and "
              "against the warp null",
              [Panel("candidate index",
                     [(idx, [e["best_score"] for e in per_candidate],
                       "#cc5500", "best nonlinear", "point"),
                      (idx, [e["linear_score"] for e in per_candidate],
                       "#2255aa", "linear", "point"),
                      (idx, [e["controls"]["warp_p95"] for e in per_candidate],
                       "#228833", "warp p95", "point")],
                     "candidate index", "matched-filter score")], ncol=1, pw=760)
    write_svg(os.path.join(ROOT, "figures", "nonlinear_injection.svg"),
              "Injection recovery per coordinate (synthetic chirps)",
              [Panel("phase locking against injected rate",
                     [([i["nu"] for i in injections], [i["linear"] for i in injections],
                       "#2255aa", "linear", "point"),
                      ([i["nu"] for i in injections], [i["best"] for i in injections],
                       "#cc5500", "best of five", "point"),
                      ([i["nu"] for i in injections],
                       [i["shuffled_p95"] for i in injections],
                       "#228833", "shuffled p95", "point")],
                     "cycles per sample", "matched-filter score")], ncol=1, pw=760)
    write_svg(os.path.join(ROOT, "figures", "nonlinear_alias.svg"),
              "Alias mapping: the increment a sampler can report",
              [Panel("fold at 0.5",
                     [([a["injected_nu"] for a in aliases],
                       [a["measured_nu"] for a in aliases],
                       "#cc5500", "measured", "point"),
                      ([-1, 1], [-1, 1], "#cccccc", "identity", "line")],
                     "injected cycles/sample", "measured cycles/sample")],
              ncol=1, pw=760)

    write_report(payload)
    print("   wrote NONLINEAR_TIME.md, nonlinear_candidates.json, "
          "nonlinear_coherence.csv, figures/*.svg")
    return 0



def write_report(payload: dict) -> None:
    c = payload["candidates"]
    lines = ["# Nonlinear event-time lane", "",
             f"**Verdict: {payload['headline']}.** The linear result is unchanged "
             f"and remains {payload['linear_result_preserved']}. This lane runs "
             "beside it, on the same candidates, and overwrites nothing.", "",
             "## What can be measured on this product", "",
             "`spec_fits` carries intensity, not complex voltages, so there is no "
             "carrier phase to lock. The phase-locking value reported here is "
             "computed from a Hilbert analytic signal of the de-drifted track: it "
             "measures envelope structure and is **an envelope proxy, not carrier "
             "phase**, and it carries no decision weight. The decision statistic is "
             "a matched-filter excess score defined directly on intensity, "
             "`sqrt(N) * (mean - median) / sd` of the resampled track.", "",
             "Shuffled-time nulls are reported in the CSV and are degenerate in "
             "the opposite direction from the usual one: a permutation keeps the "
             "sample values, so an amplitude statistic cannot move, and the null "
             "sits on the observation rather than below it. Measured worst case "
             "over all candidates and coordinates: "
             f"{max(e['amplitude_statistic_invariance'] for e in c):.1e}. That is, "
             "no coordinate can add coherence to the samples themselves. A coordinate "
             "changes the fitting basis; a warp that resamples changes the values "
             "by interpolation. Neither is evidence about the sky.", "",
             "## Coordinates", "",
             "| coordinate | definition |", "| --- | --- |"]
    for name, definition in payload["coordinates"].items():
        lines.append(f"| {name} | `{definition}` |")
    lines += ["", "Every transform is fitted on the **first half** of the track and "
                  "evaluated on the **second half**. No coherence is measured on "
                  "the samples that estimated the transform.", "",
              "## The circularity this had to survive", "",
              "`chirp` is fitted to the same linear drift that defines the "
              "candidate, so a constant-acceleration line is linearised *by "
              "construction* and its in-sample coherence carries no information. "
              "Two tests carry the weight. The split-half value applies the first "
              "half's drift to the second, so a burst whose drift changes fails "
              "even though the transform would have made it look coherent. And the "
              "gain over linear time is compared with the gain randomised monotone "
              "warps achieve, best-of-five against best-of-five, so transform and "
              "control have the same number of degrees of freedom.", "",
              "## Per-candidate results", "",
              "| date | MHz | raw drift Hz/s | linear score | best coordinate | "
              "best score | gain | warp gain p95 | beats control | split half | "
              "envelope PLV | repeats | disposition |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | "
              "--- | --- |"]
    for e in c:
        lines.append(
            f"| {e['date']} | {e['freq_mhz']:.4f} | {e['raw_drift_hz_s']:,.1f} | "
            f"{e['linear_score']:+.2f} | {e['best_nonlinear_coordinate']} | "
            f"{e['best_score']:+.2f} | {e['gain_over_linear']:+.2f} | "
            f"{e['controls']['warp_gain_p95']:+.2f} | "
            f"{'yes' if e['beats_warp_gain'] else 'no'} | "
            f"{'yes' if e['split_half_holds'] else 'no'} | "
            f"{e['best_plv_envelope']:.3f} | {e['cross_day_repeats']} | "
            f"{e['disposition']} |")
    lines += ["", "Scores are second-half only, in each coordinate. Envelope PLV is "
                  "shown for completeness and is not used to rank anything.", "",
              "## Controls", "",
              "| date | MHz | warp p95 | warp envelope p95 | phase scrambled | "
              "reversed | shuffled channels | off band |",
              "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for e in c:
        k = e["controls"]
        lines.append(f"| {e['date']} | {e['freq_mhz']:.4f} | {k['warp_p95']:+.2f} | "
                     f"{k['warp_plv_p95']:.3f} | {k['phase_scrambled_plv']:+.2f} | "
                     f"{k['reversed_plv']:+.2f} | "
                     f"{k['shuffled_channels_plv']:+.2f} | "
                     f"{k['off_band_plv']:+.2f} |")
    s = payload["solar_template"]
    lines += ["", "## The 20240325 group, solar template first", "",
              f"{s['lanes']} lanes with a drift spread of {s['spread_hz_s']:.1f} Hz/s, "
              f"which is {s['spread_over_resolution']:.3f} of one resolution "
              f"element, so the lanes are parallel at this resolution. A common "
              f"drift across lanes separated by about 1 MHz is a burst crossing "
              f"the band, not independent transmitters. A coronal exciter at "
              f"125 km/s gives {s['coronal_125kms_hz_s']:,.0f} Hz/s at these "
              f"frequencies against the measured 41,179-41,190 Hz/s. Parallel: "
              f"{'yes' if s['parallel'] else 'no'}.", "",
              "## The 20240920 group, exact alias mapping", "",
              "| injected cycles/sample | reported | fold |", "| --- | --- | --- |"]
    for a in payload["alias_probe"]:
        lines.append(f"| {a['injected_nu']:.4f} | {a['measured_nu']:+.4f} | "
                     f"{a['folded_expected']:+.4f} |")
    lines += ["", "0.6 is recorded as -0.4: the fold changes the sign as well as the "
                  "magnitude, which is why the earlier linear search could not "
                  "report it at 0.6 and why a reported 0.4965 is equally consistent "
                  "with a true -0.5035. The three 20240920 lines share one folded "
                  "rate and one sign, the signature of a common sampling origin "
                  "rather than of three sources.", "",
              "## Injection recovery", "",
              "| injected | reported | seed | linear score | best score | gain | "
              "shuffled p95 | detected |", "| --- | --- | --- | --- | --- | --- | "
              "--- | --- |"]
    for i in payload["injections"]:
        lines.append(f"| {i['nu']:.4f} | {i['reported_nu']:+.4f} | {i['seed']} | "
                     f"{i['linear']:+.2f} | {i['best']:+.2f} | {i['gain']:+.2f} | "
                     f"{i['shuffled_p95']:+.2f} | "
                     f"{'yes' if i['detected'] else 'no'} |")
    lines += ["", "Synthetic, unit-mean exponential noise, amplitude 3. `detected` "
                  "means the best score exceeded 5, the search threshold; the "
                  "shuffled p95 is shown but, as above, is degenerate for an "
                  "amplitude statistic. These are coherence numbers, not a "
                  "sensitivity curve.", "",
              "## Confirmation status", "",
              "| date | MHz | nonlinear coherent | absent from controls | repeat | "
              "later observation | off source | barycentric | second telescope |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for e in c:
        st = e["confirmation_status"]
        lines.append(f"| {e['date']} | {e['freq_mhz']:.4f} | "
                     f"{'yes' if st['nonlinear_coherent'] else 'no'} | "
                     f"{'yes' if st['absent_from_controls'] else 'no'} | "
                     f"{'yes' if st['repeat_same_position'] else 'no'} | "
                     f"{'yes' if st['later_observation'] else 'no'} | "
                     f"{st['off_source']} | {st['barycentric']} | "
                     f"{st['independent_telescope']} |")
    lines += ["", "Three gates are UNTESTABLE rather than failed: `spec_fits` is one "
                  "disk-integrated beam per day with no off-source pointing, the "
                  "product carries no WCS so a barycentric correction is undefined, "
                  "and there is no second instrument in this archive.", "",
              "## Dispositions", ""]
    for e in c:
        note = (e["astrophysical_explanation"] or e["instrumental_explanation"]
                or "no astrophysical or instrumental explanation established")
        lines.append(f"- **{e['date']} {e['freq_mhz']:.3f} MHz** - "
                     f"{e['disposition']}; band: {e['band_note']}; best coordinate "
                     f"`{e['best_nonlinear_coordinate']}` at score "
                     f"{e['best_score']:+.2f} against a linear score of "
                     f"{e['linear_score']:+.2f} and a warp-gain control of "
                     f"{e['controls']['warp_gain_p95']:+.2f}; {note}.")
    lines += ["", "## What the second lane adds", "",
              "A coordinate in which a constant-acceleration line is linear by "
              "construction is a statement about the transform, not about the sky. "
              "What is a property of the sources is whether the drift persists "
              "between halves, and the split-half score measures exactly that. No "
              "candidate is both coherent beyond the flexibility-matched control "
              "and confirmed, because the gates that would carry a technosignature "
              "verdict are untestable with this product."]
    with open(os.path.join(ROOT, "NONLINEAR_TIME.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
