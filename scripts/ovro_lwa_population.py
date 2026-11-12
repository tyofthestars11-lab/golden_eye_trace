"""Measure the OVRO-LWA drift population, and run the lock test on a wider band.

Two questions the single-day run left open. Does a wider band over a longer stretch
of a day contain any drifting line at all? And what does the archive's drifting-line
population actually look like, rather than the one quiet 253-second window?

This lists `spec_fits/` from the public bucket, picks days inside a size window,
downloads each (cached in `data/`), reads the FITS cube with numpy and stdlib only,
and searches to Nyquist. The widest run uses more samples and channels than the
others, because the coverage limit of the first test was its 253 s and 13 MHz.

Solar Type III bursts drift at MHz/s, i.e. tens of channels per 256 ms sample, so
most real drifting lines here alias past Nyquist and land wherever `mod 1` puts
them. That is worth measuring rather than assuming: the reported rate distribution
is a measurement of what this archive can and cannot resolve, and the lock distance
is scored with the same rule as `prereg/phi_lock.json`.

    python3 scripts/ovro_lwa_population.py --max-days 3 --max-mb 140
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.request
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

from ovro_lwa_test import (  # noqa: E402
    BUCKET,
    HIGH_MHZ,
    LOW_MHZ,
    TSAMP_S,
    ensure_file,
    read_primary,
)
from phi_echo import LOCK_TOLERANCE, THETA, drift_lock  # noqa: E402
from seti_search.filterbank import write_filterbank  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

NS = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}


def list_days(year: str, limit: int = 400) -> list[tuple[str, int]]:
    out: list[tuple[str, int]] = []
    token = None
    while len(out) < limit:
        url = f"{BUCKET}/?list-type=2&prefix=spec_fits/{year}/&max-keys=200"
        if token:
            url += f"&continuation-token={urllib.parse.quote(token, safe='')}"
        with urllib.request.urlopen(url, timeout=60) as response:
            tree = ET.fromstring(response.read())
        for item in tree.findall("s3:Contents", NS):
            key = item.find("s3:Key", NS).text
            size = int(item.find("s3:Size", NS).text)
            if key.endswith(".fits"):
                out.append((key, size))
        truncated = tree.findtext("s3:IsTruncated", "false", NS) == "true"
        token = tree.findtext("s3:NextContinuationToken", None, NS)
        if not truncated or not token:
            break
    return out[:limit]


def dynamic(path: str) -> tuple[np.ndarray, dict]:
    """(nfreq, ntime) Stokes-I dynamic spectrum, axes mapped from the header.

    The product carries no CTYPE cards, so the layout is inferred from axis sizes
    and then checked against the documented span and cadence before use.
    """
    header, cube = read_primary(path)
    naxis = cube.ndim

    def axis_with(token: str, default: int) -> int:
        for i in range(1, naxis + 1):
            if token in str(header.get(f"CTYPE{i}", "")).upper():
                return i
        return default

    time_i, freq_i = axis_with("TIME", 1), axis_with("FREQ", 2)
    stokes_i = axis_with("STOKES", naxis)
    ntime, nfreq = int(header[f"NAXIS{time_i}"]), int(header[f"NAXIS{freq_i}"])
    named = [stokes_i, freq_i, time_i]
    rest = [i for i in range(1, naxis + 1) if i not in named]
    order = [i - 1 for i in named] + [i - 1 for i in rest]
    moved = np.transpose(cube, order).reshape(
        int(header[f"NAXIS{stokes_i}"]), nfreq, ntime)
    channel_hz = (HIGH_MHZ - LOW_MHZ) * 1e6 / nfreq
    ok = (abs(nfreq * channel_hz / 1e6 - (HIGH_MHZ - LOW_MHZ)) < 2.0
          and ntime * TSAMP_S / 3600.0 <= 24.5)
    if not ok:
        raise ValueError(f"layout check failed for {path}: "
                         f"nfreq={nfreq}, ntime={ntime}")
    data = np.nan_to_num(moved[0].astype(np.float64),
                         nan=0.0, posinf=0.0, neginf=0.0)
    return data, {"nfreq": nfreq, "ntime": ntime, "channel_hz": channel_hz,
                  "hours": ntime * TSAMP_S / 3600.0}


def search_day(path: str, samples: int, chans: int, threshold: float) -> dict:
    data, meta = dynamic(path)
    channel_hz, nfreq, ntime = meta["channel_hz"], meta["nfreq"], meta["ntime"]
    ntime_use, nch_use = min(samples, ntime), min(chans, nfreq)
    c0 = max(0, (nfreq - nch_use) // 2)
    block = data[c0:c0 + nch_use, :ntime_use].T
    fil = os.path.join(ROOT, "data", "ovro_pop.fil")
    write_filterbank(fil, block.astype(np.float32),
                     fch1=LOW_MHZ + c0 * channel_hz / 1e6,
                     foff=channel_hz / 1e6, tsamp=TSAMP_S,
                     source="SUN (OVRO-LWA spec_fits)", nbits=32,
                     telescope_id=6, tstart=0.0)
    nyquist = 0.5 * channel_hz / TSAMP_S
    cfg = SearchConfig(boxcar=3, max_drift_hz_s=nyquist,
                       snr_threshold=threshold, max_hits=60)
    hits, _meta = search_file(fil, cfg, chan_chunk=2048, max_ints=None)
    resolution = channel_hz / (ntime_use * TSAMP_S)
    rows = []
    for hit in hits:
        resolved = abs(hit.drift_hz_s) > resolution
        nu, distance = drift_lock(hit.drift_hz_s, channel_hz, TSAMP_S)
        # registered criterion: k = 0 is a stationary line and never counts
        k = round(nu / THETA)
        rows.append({"freq_mhz": hit.freq_mhz, "drift_hz_s": hit.drift_hz_s,
                     "snr": hit.snr, "resolved": bool(resolved), "rate": nu,
                     "k": int(k), "distance": distance,
                     "locked": bool(resolved and k != 0
                                    and distance <= LOCK_TOLERANCE)})
    return {"ntime": ntime, "nfreq": nfreq, "channel_hz": channel_hz,
            "samples_used": int(ntime_use), "channels_used": int(nch_use),
            "hours": meta["hours"], "nyquist_hz_s": nyquist,
            "resolution_hz_s": resolution,
            "theta_drift_hz_s": THETA * channel_hz / TSAMP_S,
            "candidates": len(hits), "rows": rows,
            "resolved": int(sum(r["resolved"] for r in rows)),
            "locked": int(sum(r["locked"] for r in rows))}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--year", default="2024")
    ap.add_argument("--max-days", type=int, default=3)
    ap.add_argument("--max-mb", type=float, default=140.0)
    ap.add_argument("--wide-samples", type=int, default=1600)
    ap.add_argument("--wide-chans", type=int, default=192)
    ap.add_argument("--threshold", type=float, default=5.5)
    args = ap.parse_args()

    print("OVRO-LWA DRIFT POPULATION AND WIDE-BAND LOCK TEST")
    days = list_days(args.year)
    print(f"   {len(days)} day files listed for {args.year}")
    pool = [(k, s) for k, s in days if 20e6 < s <= args.max_mb * 1e6]
    pool.sort(key=lambda kv: kv[1], reverse=True)
    chosen = pool[: args.max_days]
    if not chosen:
        print("   no day file inside the size window; nothing to search")
        return 2
    print(f"   {len(pool)} files inside the window; searching "
          f"{len(chosen)} of them")
    for key, size in chosen:
        print(f"      {key}  {size / 1e6:.1f} MB")

    results, failures = [], []
    for index, (key, size) in enumerate(chosen):
        wide = index == 0
        samples = args.wide_samples if wide else max(600, args.wide_samples // 2)
        chans = args.wide_chans if wide else max(64, args.wide_chans // 2)
        print(f"\n   [{index + 1}/{len(chosen)}] {os.path.basename(key)}, "
              f"{samples} samples x {chans} channels"
              f"{' (widest run)' if wide else ''}")
        try:
            path = ensure_file(key)
            day = search_day(path, samples, chans, args.threshold)
        except Exception as exc:  # noqa: BLE001 - recorded, not swallowed
            print(f"      failed: {exc}")
            failures.append({"key": key, "error": str(exc)})
            continue
        day["key"] = key
        day["size_mb"] = size / 1e6
        results.append(day)
        print(f"      {day['samples_used']} x {day['channels_used']}, "
              f"{day['samples_used'] * TSAMP_S:.0f} s, resolution "
              f"{day['resolution_hz_s']:.1f} Hz/s, theta "
              f"{day['theta_drift_hz_s']:,.0f} Hz/s")
        print(f"      {day['candidates']} candidates at {args.threshold:g} sigma, "
              f"{day['resolved']} with resolved drift, {day['locked']} locked")
        for row in sorted(day["rows"], key=lambda r: -abs(r["drift_hz_s"]))[:5]:
            print(f"         {row['freq_mhz']:9.3f} MHz  drift "
                  f"{row['drift_hz_s']:>12.1f} Hz/s  snr {row['snr']:5.2f}  "
                  f"rate {row['rate']:.4f}  dist {row['distance']:.3f}  "
                  f"{'resolved' if row['resolved'] else 'below resolution'}")

    all_rows = [r for day in results for r in day["rows"]]
    resolved = [r for r in all_rows if r["resolved"]]
    locked = [r for r in resolved if r["locked"]]
    chance = 2.0 * LOCK_TOLERANCE * THETA / 0.5
    print("\n   POPULATION")
    print(f"   {len(results)} days searched, {len(all_rows)} candidates, "
          f"{len(resolved)} resolved drifts, {len(locked)} locked")
    if resolved:
        drifts = np.array([r["drift_hz_s"] for r in resolved])
        rates = np.array([r["rate"] for r in resolved])
        print(f"   |drift| median {np.median(np.abs(drifts)):,.0f} Hz/s, "
              f"max {np.abs(drifts).max():,.0f} Hz/s")
        print(f"   rate (cycles/sample) median {np.median(rates):.4f}, "
              f"min {rates.min():.4f}, max {rates.max():.4f}")
        print(f"   expected locked by chance: {len(resolved) * chance:.2f}")
    verdict = ("NOT MEASURED: no drifting line resolved in any day searched"
               if not resolved else
               f"{len(resolved)} resolved, {len(locked)} within "
               f"{LOCK_TOLERANCE} theta, {len(resolved) * chance:.2f} expected")
    print(f"   verdict: {verdict}")

    out = {"year": args.year, "threshold": args.threshold,
           "max_mb": args.max_mb, "days": results, "failures": failures,
           "candidates": len(all_rows), "resolved": len(resolved),
           "locked": len(locked), "chance_per_candidate": chance,
           "expected_locked_by_chance": len(resolved) * chance,
           "drift_hz_s": [r["drift_hz_s"] for r in resolved],
           "rates": [r["rate"] for r in resolved], "verdict": verdict}
    path = os.path.join(ROOT, "prereg", "ovro_lwa_population.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"\n   wrote prereg/ovro_lwa_population.json")
    return 0 if resolved else 3


if __name__ == "__main__":
    raise SystemExit(main())
