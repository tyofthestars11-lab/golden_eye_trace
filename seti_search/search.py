"""Narrowband drift-rate search ("dedoppler") for technosignatures.

One band of channels goes through four stages:

  1. bandpass control  divide out a smooth estimate of the bandpass so every
                       channel has the same noise level
  2. RFI excision      blank carrier-dominated and broadband-glitch samples
  3. dedoppler         for each trial drift, shift each integration by the
                       accumulated Doppler offset and sum: a signal drifting at
                       that rate adds coherently, noise does not
  4. detection         boxcar-smooth the summed spectrum, normalise by a robust
                       noise estimate (median absolute deviation), rescale the
                       map by its own measured bulk width, then take peaks with
                       non-maximum suppression

The detection statistic is therefore "sigma above the band median" in units of
the robust noise level after smoothing, which is the convention used by
turboSETI and the Breakthrough Listen search. The final rescaling is what makes
the printed sigma a measurement of this file rather than of an idealised
Gaussian; `map_tail` reports the part of it the rescaling cannot fix.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .filterbank import Filterbank
from .stabilizer import WindowTooNarrow, check_window
from .trials import tail_probability

__all__ = [
    "Candidate",
    "SearchConfig",
    "control_bandpass",
    "dedoppler",
    "excise_rfi",
    "extract_hits",
    "map_tail",
    "search_band",
    "search_file",
    "smooth_bandpass",
]


@dataclass
class Candidate:
    """One narrowband candidate.

    freq_mhz is the frequency at the first integration (that is where the
    dedoppler sum peaks); drift_hz_s is positive for increasing frequency with
    time. A single-file candidate is not a detection: a confirmed signal has to
    survive an ON/OFF cadence, where it appears in the ON scan and vanishes in
    the OFF scan (plus its drift must be consistent with a sky position).
    """

    freq_mhz: float
    drift_hz_s: float
    snr: float
    channel: int = -1
    boxcar: int = 3

    def as_dict(self) -> dict:
        return {
            "freq_mhz": self.freq_mhz,
            "drift_hz_s": self.drift_hz_s,
            "snr": self.snr,
            "channel": self.channel,
        }


@dataclass
class SearchConfig:
    boxcar: int = 3
    max_drift_hz_s: float = 4.0
    drift_step_hz_s: float | None = None
    snr_threshold: float = 10.0
    max_hits: int = 50
    rfi_sigma: float = 6.0
    carrier_sigma: float = 20.0
    #: Refuse a drift range whose trail the search window cannot hold. Off by
    #: default so a fixture can deliberately search a truncated window and show
    #: what that does; the geometry is measured and recorded either way.
    require_flat_sensitivity: bool = False


def smooth_bandpass(spectrum: np.ndarray, width: int = 101) -> tuple[np.ndarray, int]:
    """Block-median bandpass estimate, interpolated back to every channel.

    A block median ignores narrowband features narrower than half a block, which
    is what we want: the bandpass should follow the receiver's smooth response,
    not the carriers sitting on top of it.
    """
    spectrum = np.asarray(spectrum, dtype=np.float64)
    n = spectrum.size
    width = int(max(9, min(width, max(9, n // 4))))
    if width % 2 == 0:
        width += 1
    nblocks = max(1, n // width)
    edges = np.linspace(0, n, nblocks + 1).astype(int)

    centres, medians = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        if b <= a:
            continue
        centres.append(0.5 * (a + b - 1))
        medians.append(float(np.median(spectrum[a:b])))

    centres_arr = np.asarray(centres, dtype=np.float64)
    medians_arr = np.asarray(medians, dtype=np.float64)
    if medians_arr.size >= 3:
        padded = np.concatenate(([medians_arr[0]], medians_arr, [medians_arr[-1]]))
        medians_arr = np.convolve(padded, np.ones(3) / 3.0, mode="valid")

    return np.interp(np.arange(n, dtype=np.float64), centres_arr, medians_arr), width


def control_bandpass(waterfall: np.ndarray, width: int = 101):
    """Divide out the bandpass; return fractional fluctuations and the estimate."""
    data = np.asarray(waterfall, dtype=np.float32)
    if data.ndim != 2:
        raise ValueError(f"expected a 2-D waterfall, got shape {data.shape}")
    if data.shape[1] < 8:
        raise ValueError(f"need at least 8 channels, got {data.shape[1]}")

    bandpass, _ = smooth_bandpass(data.mean(axis=0), width=width)
    bandpass = np.where(bandpass > 0, bandpass, 1.0)
    normalised = data / bandpass[None, :].astype(np.float32) - 1.0
    return normalised.astype(np.float32, copy=False), bandpass.astype(np.float32)


def excise_rfi(
    normalised: np.ndarray,
    sigma: float = 6.0,
    carrier_sigma: float = 20.0,
    pad: int = 1,
):
    """Blank samples that cannot hold a narrowband signal.

    Three scale-free tests, measured against the band's own spread:

    * dead or saturated channels - the noise variance is gone, so whatever is
      there is not noise plus signal. This is the one unambiguous signature;
    * overwhelming carriers - median power tens of sigma above the noise floor.
      A tone that loud is terrestrial in practice;
    * broadband glitches - one integration far brighter than the rest, applied
      along time rather than frequency.

    Deliberately *not* excised: stationary tones at ordinary amplitude. They
    look exactly like a beacon that happens not to drift, and a single file
    cannot distinguish the two - the ON/OFF follow-up is what settles it. They
    come out as zero-drift candidates instead of being silently deleted, which
    is also why a slow-drift signal survives this stage: it is bright for a
    while in a channel, not bright in every integration.

    Zeroing rather than masking keeps the arithmetic simple; blanked channels
    contribute no flux and no noise.
    """
    data = np.asarray(normalised, dtype=np.float32)
    n_t, n_f = data.shape
    if n_t < 2:
        raise ValueError("need at least two integrations")

    std = data.std(axis=0)
    typical_std = float(np.median(std))
    if not np.isfinite(typical_std) or typical_std <= 0:
        typical_std = 1.0

    dead = std < 0.5 * typical_std
    carrier = np.median(data, axis=0) > carrier_sigma * typical_std

    bad_channels = dead | carrier
    if pad > 0 and bad_channels.any():
        kernel = np.ones(2 * pad + 1)
        bad_channels = np.convolve(bad_channels.astype(np.int8), kernel, mode="same") > 0

    bad_rows = np.zeros(n_t, dtype=bool)
    if n_f > 1:
        # Relative to the band's own baseline: the normalised data is not
        # centred on zero, so comparing against zero would flag every row. The
        # baseline is the median of the row medians, not the median of the whole
        # array: for skewed total-power noise those differ by a few per cent,
        # which exceeds the row threshold whenever n_f is large enough to make
        # row_noise small, and then every row of a wide band is thrown away.
        row_median = np.median(data, axis=1)
        row_excess = row_median - float(np.median(row_median))
        row_noise = 1.25 / np.sqrt(n_f) * typical_std
        bad_rows = np.abs(row_excess) > sigma * row_noise

    clean = data.copy()
    if bad_channels.any():
        clean[:, bad_channels] = 0.0
    if bad_rows.any():
        clean[bad_rows, :] = 0.0

    return clean, {
        "blanked_channels": int(bad_channels.sum()),
        "blanked_rows": int(bad_rows.sum()),
        "total_channels": int(n_f),
        "total_rows": int(n_t),
    }


def _smoothed_snr(accumulated: np.ndarray, boxcar: int) -> np.ndarray:
    """Boxcar sum (variance-normalised) turned into robust signal-to-noise."""
    width = max(1, int(boxcar))
    half = width // 2
    n = accumulated.size

    cumulative = np.concatenate(([0.0], np.cumsum(accumulated, dtype=np.float64)))
    index = np.arange(n)
    lo = np.clip(index - half, 0, n)
    hi = np.clip(index - half + width, 0, n)
    count = np.maximum(hi - lo, 1)
    # divide by sqrt(count) so the noise level is flat across the band edges
    smoothed = (cumulative[hi] - cumulative[lo]) / np.sqrt(count)

    median = float(np.median(smoothed))
    mad = 1.4826 * float(np.median(np.abs(smoothed - median)))
    if not np.isfinite(mad) or mad <= 0:
        mad = float(smoothed.std()) or 1.0
    return ((smoothed - median) / mad).astype(np.float32)


def dedoppler(
    data: np.ndarray,
    tsamp_s: float,
    df_mhz: float,
    *,
    max_drift_hz_s: float = 4.0,
    drift_step_hz_s: float | None = None,
    boxcar: int = 3,
):
    """Incoherent dedoppler search over one band.

    Returns (snr_map, drifts) where snr_map[k] is the SNR spectrum at drifts[k].
    Frequency is assumed to increase with channel index (callers reverse bands
    that descend), and a trial drift d>0 means the signal's frequency rises with
    time.
    """
    data = np.asarray(data, dtype=np.float32)
    n_t, n_f = data.shape
    if n_t < 2 or n_f < 8:
        raise ValueError(f"waterfall too small: {data.shape}")

    df_mhz = abs(float(df_mhz))
    if df_mhz <= 0:
        raise ValueError("channel width must be positive")
    total_time = n_t * float(tsamp_s)

    # above this a signal jumps more than half a channel per integration, where
    # the shift-and-add can no longer align it
    nyquist_drift = 0.5 * df_mhz / float(tsamp_s) * 1e6
    max_drift = min(float(max_drift_hz_s), nyquist_drift)
    if drift_step_hz_s is None:
        # keep the smear across the whole observation below half a boxcar
        step = 0.5 * boxcar * df_mhz * 1e6 / total_time
    else:
        step = float(drift_step_hz_s)
    step = max(step, 1e-5)
    n_steps = int(np.floor(max_drift / step))
    drifts = np.arange(-n_steps, n_steps + 1, dtype=np.float64) * step

    times = np.arange(n_t, dtype=np.float64) * float(tsamp_s)
    snr = np.empty((drifts.size, n_f), dtype=np.float32)

    for k, drift in enumerate(drifts):
        # channel offset accumulated by integration i at this trial drift
        shifts = np.round(-drift * 1e-6 * times / df_mhz).astype(np.int64)
        acc = np.zeros(n_f, dtype=np.float32)
        for i in range(n_t):
            shift = int(shifts[i])
            row = data[i]
            if shift == 0:
                acc += row
            elif 0 < shift < n_f:
                acc[shift:] += row[: n_f - shift]
            elif -n_f < shift < 0:
                acc[: n_f + shift] += row[-shift:]
            # shifts beyond the band drop out entirely: no circular wrap-around,
            # which would otherwise fold one edge of the band onto the other
        snr[k] = _smoothed_snr(acc, boxcar)

    return snr, drifts


def calibrate_map(
    snr_map: np.ndarray,
    live=None,
) -> tuple[np.ndarray, float]:
    """Rescale a dedoppler map by its own measured spread.

    The detection statistic is power divided by 1.4826 x MAD, which equals the
    standard deviation only for Gaussian noise. Total-power data is wider and
    heavier-tailed than that: on the real break-through-listen fast-transient
    file the map's spread measured 1.33-1.44 instead of 1, so a peak printed as
    5 sigma was about 3.5 sigma and every flux limit and false-alarm budget
    inherited the optimism. Dividing by the measured spread makes the printed
    sigma the sigma of this map.

    ``live`` masks blanked channels (whose map entries are exactly zero and
    would otherwise bias the spread low). Returns (map, spread); the map is
    returned unchanged with spread 1.0 when the input is too small to measure.

    The spread is measured robustly, by sigma-clipping around the bulk: a plain
    standard deviation is driven by the very peaks being judged. On a 4,096-cell
    band holding one 219-sigma injection, the standard deviation read 3.73
    instead of 1.1, and the same source then appeared 5.4 times fainter on one
    receiver than on another. Clipping at four pilot sigmas keeps a real signal
    from setting the noise level it is measured against.
    """
    values = np.asarray(snr_map, dtype=np.float64)
    if live is not None:
        mask = np.asarray(live, dtype=bool)
        if mask.size == values.shape[-1] and mask.any():
            values = values[:, mask]
    spread = _clipped_spread(values.ravel())
    if spread is None:
        return np.asarray(snr_map, dtype=np.float32), 1.0
    return np.asarray(snr_map, dtype=np.float32) / np.float32(spread), spread


def map_tail(
    snr_map: np.ndarray,
    live=None,
    sigmas: tuple[float, ...] = (4.0, 5.0),
) -> dict[float, dict[str, float]]:
    """How far this map's measured upper tail sits from the Gaussian model.

    Calibration makes the bulk one sigma wide, but a false-alarm budget is a
    claim about the far tail, and total-power data keeps heavier tails than a
    Gaussian even after the bulk is corrected. So count what actually sits above
    each sigma and compare it with what a normal would give over the same cells.

    These are map cells, not independent trials: neighbouring drift steps share
    integrations and the boxcar smears channels, so the ratio is an upper bound
    on the excess a per-trial rate should carry, not a correction factor to
    divide by. It says which way the budget is wrong, which is the part that was
    previously invisible.

    A real signal is counted too: a bright candidate's smear cone runs across
    many drift trials, so its own map fills cells above every threshold. A large
    ratio can therefore mean "this map contains a beacon", not "this noise is
    heavy-tailed", and the two are only distinguishable by reading the map.
    """
    values = np.asarray(snr_map, dtype=np.float64)
    if live is not None:
        mask = np.asarray(live, dtype=bool)
        if mask.size == values.shape[-1] and mask.any():
            values = values[:, mask]
    finite = values[np.isfinite(values)]
    report: dict[float, dict[str, float]] = {}
    for sigma in sigmas:
        expected = tail_probability(float(sigma)) * float(finite.size)
        observed = int(np.count_nonzero(finite > float(sigma)))
        if expected > 0.0:
            ratio = float(observed) / expected
        else:
            ratio = float("inf") if observed else 0.0
        report[float(sigma)] = {
            "cells": int(finite.size),
            "observed": observed,
            "expected": expected,
            "ratio": ratio,
        }
    return report


def _clipped_spread(values: np.ndarray) -> float | None:
    """Width of the bulk, from iterated clipping around the median."""
    finite = values[np.isfinite(values)]
    if finite.size < 8:
        return None
    centre = float(np.median(finite))
    width = float(np.percentile(np.abs(finite - centre), 68.27))
    if not np.isfinite(width) or width <= 0.0:
        width = float(np.std(finite))
    for _ in range(5):
        keep = np.abs(finite - centre) <= 4.0 * width
        if int(keep.sum()) < 8:
            break
        revised = float(np.std(finite[keep]))
        if not np.isfinite(revised) or revised <= 0.0:
            break
        settled = abs(revised - width) <= 1e-3 * width
        width = revised
        if settled:
            break
    return width if np.isfinite(width) and width > 0.0 else None


def extract_hits(
    snr_map: np.ndarray,
    drifts: np.ndarray,
    freqs_mhz: np.ndarray,
    *,
    threshold: float,
    boxcar: int,
    total_time: float,
    max_hits: int,
) -> list[Candidate]:
    """Peak extraction with suppression in both frequency and drift.

    One real signal peaks at its true drift but also scores across neighbouring
    trials, so suppressing only in frequency would report the same signal many
    times. The drift suppression half-width is the range over which the drift
    mismatch smears the signal by less than a boxcar.
    """
    n_d, n_f = snr_map.shape
    freqs_mhz = np.asarray(freqs_mhz, dtype=np.float64)
    df = abs(float(freqs_mhz[1] - freqs_mhz[0])) if n_f > 1 else 0.0
    step = abs(float(drifts[1] - drifts[0])) if drifts.size > 1 else 0.0

    if step > 0 and total_time > 0:
        smear = boxcar * df * 1e6 / total_time  # Hz/s smeared by one boxcar
        drift_window = int(np.ceil(smear / step)) + 1
    else:
        drift_window = 1
    freq_window = max(2, int(boxcar) + 1)

    work = snr_map.copy()
    hits: list[Candidate] = []
    while len(hits) < max_hits:
        flat = int(np.argmax(work))
        k, j = np.unravel_index(flat, work.shape)
        peak = float(work[k, j])
        if not np.isfinite(peak) or peak < threshold:
            break

        # parabolic interpolation on the unsuppressed map for sub-step accuracy
        frac_f = 0.0
        if 0 < j < n_f - 1:
            y0, y1, y2 = float(snr_map[k, j - 1]), peak, float(snr_map[k, j + 1])
            denom = y0 - 2.0 * y1 + y2
            if denom < 0:
                frac_f = float(np.clip(0.5 * (y0 - y2) / denom, -0.5, 0.5))
        frac_d = 0.0
        if 0 < k < n_d - 1:
            y0, y1, y2 = float(snr_map[k - 1, j]), peak, float(snr_map[k + 1, j])
            denom = y0 - 2.0 * y1 + y2
            if denom < 0:
                frac_d = float(np.clip(0.5 * (y0 - y2) / denom, -0.5, 0.5))

        hits.append(
            Candidate(
                freq_mhz=float(freqs_mhz[j] + frac_f * df),
                drift_hz_s=float(drifts[k] + frac_d * step),
                snr=peak,
                channel=int(j),
                boxcar=int(boxcar),
            )
        )

        k0, k1 = max(0, k - drift_window), min(n_d, k + drift_window + 1)
        j0, j1 = max(0, j - freq_window), min(n_f, j + freq_window + 1)
        work[k0:k1, j0:j1] = -np.inf

    # A single signal also scores several drift steps away, where the mismatch
    # smears it instead of aligning it. Such an echo is displaced along the
    # smear cone: a drift error dd can move the peak by up to dd * T in
    # frequency, plus a boxcar or two. Collapse hits that fall inside each
    # other's cone, keeping the strongest of the group. Two genuine signals that
    # close together are merged - the drift axis is coarse enough at that
    # separation that they are not separable anyway.
    def same_signal(a: Candidate, b: Candidate) -> bool:
        separation = abs(a.freq_mhz - b.freq_mhz)
        spread = abs(a.drift_hz_s - b.drift_hz_s) * 1e-6 * max(total_time, 1.0)
        return separation <= spread + (2 * int(boxcar) + 1) * df

    merged: list[Candidate] = []
    for hit in sorted(hits, key=lambda h: -h.snr):
        if not any(same_signal(hit, kept) for kept in merged):
            merged.append(hit)
    return merged


def search_band(
    data: np.ndarray,
    freqs_mhz: np.ndarray,
    tsamp_s: float,
    cfg: SearchConfig,
):
    """Full pipeline for one band whose channels ascend in frequency."""
    data = np.asarray(data, dtype=np.float32)
    freqs_mhz = np.asarray(freqs_mhz, dtype=np.float64)
    if data.shape[1] != freqs_mhz.size:
        raise ValueError(
            f"band has {data.shape[1]} channels but {freqs_mhz.size} frequencies"
        )
    if freqs_mhz.size < 2:
        raise ValueError("need at least two channels")

    df_mhz = float(freqs_mhz[1] - freqs_mhz[0])
    normalised, _ = control_bandpass(data)
    clean, rfi = excise_rfi(
        normalised, sigma=cfg.rfi_sigma, carrier_sigma=cfg.carrier_sigma
    )

    snr_map, drifts = dedoppler(
        clean,
        tsamp_s,
        df_mhz,
        max_drift_hz_s=cfg.max_drift_hz_s,
        drift_step_hz_s=cfg.drift_step_hz_s,
        boxcar=cfg.boxcar,
    )
    live = ~(clean == 0).all(axis=0)
    snr_map, calibration = calibrate_map(snr_map, live)
    rfi["snr_calibration"] = float(calibration)
    rfi["snr_tail"] = map_tail(snr_map, live)
    hits = extract_hits(
        snr_map,
        drifts,
        freqs_mhz,
        threshold=cfg.snr_threshold,
        boxcar=cfg.boxcar,
        total_time=data.shape[0] * float(tsamp_s),
        max_hits=cfg.max_hits,
    )
    return hits, rfi


def search_file(
    path,
    cfg: SearchConfig | None = None,
    *,
    chan_chunk: int = 1 << 18,
    max_ints: int | None = 4096,
    progress=None,
    tsamp: float | None = None,
):
    """Search a filterbank file, band by band, and return (hits, metadata).

    tsamp overrides the header value, for files whose header omits it.
    """
    cfg = cfg or SearchConfig()
    fb = Filterbank.open(path)

    tsamp_s = float(fb.tsamp if tsamp is None else tsamp)
    if not np.isfinite(tsamp_s) or tsamp_s <= 0:
        raise ValueError(
            "this file's header has no usable tsamp; pass tsamp=<seconds>"
            " (CLI: --tsamp)"
        )

    n_ints = fb.nints if not max_ints else min(fb.nints, int(max_ints))
    if n_ints < 2:
        raise ValueError("need at least two integrations")

    file_freqs = fb.frequencies()
    flip = fb.foff < 0
    asc_freqs = file_freqs[::-1].copy() if flip else file_freqs
    if fb.nchans < 2:
        raise ValueError("need at least two channels")
    df_mhz = abs(float(asc_freqs[1] - asc_freqs[0]))
    if df_mhz <= 0:
        raise ValueError(f"channel width is zero (foff={fb.foff})")

    chunk = max(1024, min(int(chan_chunk), fb.nchans))
    # keep each band block near 150M float32 elements
    budget = 1.5e8
    if n_ints * chunk > budget:
        chunk = max(1024, int(budget / n_ints))

    total_time = n_ints * tsamp_s
    # a signal can drift up to this many channels across the observation, so the
    # overlap between chunks has to cover the whole excursion, not just the
    # boxcar, or signals near a boundary are truncated
    drift_channels = abs(cfg.max_drift_hz_s) * 1e-6 * total_time / df_mhz
    pad = max(int(cfg.boxcar) * 2 + 8, min(int(np.ceil(drift_channels)) + 2, 8192))

    # The chunk is padded by the whole drift excursion, so the window the
    # dedoppler integrates over is chunk + 2*pad, not chunk. The pad is capped at
    # 8192 channels, so a long integration at a fast drift is the case where the
    # window stops being able to hold the trail - and that is measured here
    # rather than assumed, because a sigma from a truncated window is a statement
    # about the window.
    window = check_window(
        cfg.max_drift_hz_s,
        min(fb.nchans, chunk + 2 * pad),
        n_ints,
        df_mhz * 1e6,
        tsamp_s,
    )
    if not window.fits and cfg.require_flat_sensitivity:
        raise WindowTooNarrow(window)

    hits: list[Candidate] = []
    blanked_channels = 0
    blanked_rows = 0
    calibrations: list[tuple[float, int]] = []
    tails: dict[float, dict[str, float]] = {}

    for c0 in range(0, fb.nchans, chunk):
        c1 = min(c0 + chunk, fb.nchans)
        lo, hi = max(0, c0 - pad), min(fb.nchans, c1 + pad)

        block = fb.read(0, n_ints, lo, hi)
        if flip:
            block = np.ascontiguousarray(block[:, ::-1])
            sub_freqs = asc_freqs[fb.nchans - hi:fb.nchans - lo]
        else:
            sub_freqs = asc_freqs[lo:hi]

        band_hits, rfi = search_band(block, sub_freqs, tsamp_s, cfg)
        blanked_channels += rfi["blanked_channels"]
        blanked_rows += rfi["blanked_rows"]
        live_cells = max((hi - lo) - rfi["blanked_channels"], 1)
        calibrations.append((float(rfi.get("snr_calibration", 1.0)), live_cells))
        for sigma, entry in (rfi.get("snr_tail") or {}).items():
            slot = tails.setdefault(
                float(sigma), {"cells": 0, "observed": 0, "expected": 0.0}
            )
            slot["cells"] += int(entry["cells"])
            slot["observed"] += int(entry["observed"])
            slot["expected"] += float(entry["expected"])

        for hit in band_hits:
            # map the band index back to the file's own channel order
            file_chan = (hi - 1 - hit.channel) if flip else (lo + hit.channel)
            if c0 <= file_chan < c1:
                hit.channel = int(file_chan)
                hits.append(hit)

        if progress is not None:
            progress(c0, c1, fb.nchans, len(hits))

    hits.sort(key=lambda h: -h.snr)
    hits = hits[: cfg.max_hits]

    meta = {
        "path": str(path),
        "source": fb.source,
        "tstart": fb.tstart,
        "tsamp": tsamp_s,
        "n_ints": n_ints,
        "n_ints_in_file": fb.nints,
        "nchans": fb.nchans,
        "df_mhz": df_mhz,
        "fch1": fb.fch1,
        "foff": fb.foff,
        "total_time": total_time,
        "blanked_channels": blanked_channels,
        "blanked_rows": blanked_rows,
        "truncated": bool(n_ints < fb.nints),
        "snr_threshold": cfg.snr_threshold,
        "boxcar": cfg.boxcar,
        "window_check": window.as_dict(),
        "snr_calibration": (
            float(
                np.average(
                    [value for value, _ in calibrations],
                    weights=[weight for _, weight in calibrations],
                )
            )
            if calibrations
            else 1.0
        ),
        # summed over bands, so it describes the whole searched map; the ratio is
        # recomputed here because a band's own ratio does not average
        "snr_tail": {
            sigma: {
                **entry,
                "ratio": (
                    entry["observed"] / entry["expected"]
                    if entry["expected"] > 0.0
                    else (float("inf") if entry["observed"] else 0.0)
                ),
            }
            for sigma, entry in sorted(tails.items())
        },
    }
    return hits, meta
