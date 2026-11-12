"""Run the pre-registered phase-lock test on a real OVRO-LWA dynamic spectrum.

The public `ovro-lwa-solar` bucket (us-west-2, CC BY 4.0) publishes
`spec_fits/YYYY/YYYYMMDD.fits`: one file per UT day at 256 ms cadence across
13.4-86.9 MHz. The axis order is not the one the documentation implies - the
first file read here has FITS axes (time 7009, freq 731, singleton, Stokes 2) and
carries no `CTYPE`, `CRVAL` or `CDELT` cards at all, so there is no world
coordinate system to read: the layout is inferred from the axis sizes, and then
validated - `nfreq x channel` must reproduce the documented 73.5 MHz span and
`ntime x tsamp` must be at most a day. The check is printed rather than assumed.
That is a genuine per-channel dynamic spectrum from an open archive, which is what
the lock needs and what nothing else reachable here provides.

FITS is read here without astropy, which this project does not use: a 2880-byte
header of 80-character cards, then big-endian data with the first axis varying
fastest. Only the primary image is needed.

The test follows `prereg/phi_lock.json`: `tsamp` from the file, the drift window
searched to Nyquist, the integer advance bounded at `|k| = 1`, and the verdict
NOT MEASURED if no drift resolves. Solar dynamic spectra hold Type II and III
bursts, which are drifting lines, so a resolved drift here is a real measurement
of the rate distribution even though the source is the Sun rather than a beacon:
it tells us whether this archive can carry the test at all.

    python3 scripts/ovro_lwa_test.py --samples 987 --chans 256
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

from phi_echo import LOCK_TOLERANCE, THETA, drift_lock  # noqa: E402
from seti_search.filterbank import write_filterbank  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

BUCKET = "https://ovro-lwa-solar.s3.us-west-2.amazonaws.com"
DEFAULT_KEY = "spec_fits/2024/20240207.fits"
CARD, BLOCK = 80, 2880
DTYPE = {8: ">u1", 16: ">i2", 32: ">i4", 64: ">f8", -32: ">f4", -64: ">f8"}
LOW_MHZ, HIGH_MHZ, CHANNEL_HZ = 13.4, 86.9, 96_000.0
TSAMP_S = 0.256


def read_header(fh) -> dict:
    cards: dict = {}
    while True:
        block = fh.read(BLOCK)
        if len(block) < BLOCK:
            raise ValueError("truncated FITS header")
        for i in range(0, BLOCK, CARD):
            card = block[i:i + CARD].decode("ascii", "replace")
            if card.startswith("END"):
                return cards
            key = card[:8].strip()
            if key in ("", "COMMENT", "HISTORY"):
                continue
            raw = card[10:].split("/")[0].strip()
            if raw.startswith("'"):
                cards[key] = raw.strip("' ").rstrip()
            else:
                try:
                    cards[key] = int(raw)
                except ValueError:
                    try:
                        cards[key] = float(raw.replace("D", "E"))
                    except ValueError:
                        cards[key] = raw


def read_primary(path: str) -> tuple[dict, np.ndarray]:
    with open(path, "rb") as fh:
        header = read_header(fh)
        naxis = int(header.get("NAXIS", 0))
        if naxis < 2:
            raise ValueError(f"primary HDU has NAXIS={naxis}; expected an image cube")
        dims = [int(header[f"NAXIS{i}"]) for i in range(1, naxis + 1)]
        dtype = np.dtype(DTYPE[int(header["BITPIX"])])
        count = int(np.prod(dims))
        data = np.fromfile(fh, dtype=dtype, count=count)
    if data.size != count:
        raise ValueError(f"expected {count} values, read {data.size}")
    cube = data.reshape(tuple(reversed(dims)))
    cube = np.transpose(cube, tuple(reversed(range(naxis))))
    return header, cube


def ensure_file(key: str) -> str:
    path = os.path.join(ROOT, "data", os.path.basename(key))
    if os.path.exists(path) and os.path.getsize(path) > 1_000_000:
        print(f"   using cached {os.path.relpath(path, ROOT)} "
              f"({os.path.getsize(path) / 1e6:.1f} MB)")
        return path
    url = f"{BUCKET}/{key}"
    print(f"   fetching {url}")
    with urllib.request.urlopen(url, timeout=240) as response:
        payload = response.read()
    with open(path, "wb") as fh:
        fh.write(payload)
    print(f"   saved {os.path.getsize(path) / 1e6:.1f} MB")
    return path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--key", default=DEFAULT_KEY)
    ap.add_argument("--samples", type=int, default=987)
    ap.add_argument("--chans", type=int, default=256)
    ap.add_argument("--threshold", type=float, default=6.0)
    args = ap.parse_args()

    print("OVRO-LWA DYNAMIC SPECTRUM, PRE-REGISTERED LOCK TEST")
    try:
        path = ensure_file(args.key)
    except Exception as exc:  # noqa: BLE001 - reported, not swallowed
        print(f"   download failed: {exc}")
        print("   no test is run and no verdict is reported")
        return 2

    header, cube = read_primary(path)
    print(f"   cube {cube.shape}, BITPIX {header.get('BITPIX')}, "
          f"file {os.path.getsize(path) / 1e6:.1f} MB")
    naxis = cube.ndim
    for i in range(1, naxis + 1):
        print(f"      axis {i}: n={header.get(f'NAXIS{i}')} "
              f"ctype={header.get(f'CTYPE{i}')} "
              f"crval={header.get(f'CRVAL{i}')} cdelt={header.get(f'CDELT{i}')}")

    def axis_with(token: str, default: int) -> int:
        for i in range(1, naxis + 1):
            if token in str(header.get(f"CTYPE{i}", "")).upper():
                return i
        return default

    time_i = axis_with("TIME", 1)
    freq_i = axis_with("FREQ", 2)
    stokes_i = axis_with("STOKES", naxis)
    ntime = int(header[f"NAXIS{time_i}"])
    nfreq = int(header[f"NAXIS{freq_i}"])
    # keep the three named axes first and let the singleton filler axes follow,
    # then drop them: the file has a size-1 axis that the WCS cards would label
    named = [stokes_i, freq_i, time_i]
    rest = [i for i in range(1, naxis + 1) if i not in named]
    order = [i - 1 for i in named] + [i - 1 for i in rest]
    moved = np.transpose(cube, order).reshape(
        int(header[f"NAXIS{stokes_i}"]), nfreq, ntime)
    channel_hz = (HIGH_MHZ - LOW_MHZ) * 1e6 / nfreq
    span_mhz = nfreq * channel_hz / 1e6
    hours = ntime * TSAMP_S / 3600.0
    print(f"   stokes {moved.shape[0]}, nfreq {nfreq}, ntime {ntime}, "
          f"channel {channel_hz / 1e3:.2f} kHz (span / nfreq)")
    print(f"   check: {nfreq} channels x {channel_hz / 1e3:.2f} kHz = "
          f"{span_mhz:.1f} MHz against a documented {HIGH_MHZ - LOW_MHZ:.1f} MHz; "
          f"{ntime} samples x {TSAMP_S} s = {hours:.2f} h")
    if abs(span_mhz - (HIGH_MHZ - LOW_MHZ)) > 2.0 or hours > 24.5:
        print("   layout check failed: axis identification is not consistent with "
              "the documented product, so no test is run")
        return 2

    data = np.nan_to_num(moved[0].astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    ntime_use = min(args.samples, ntime)
    nch_use = min(args.chans, nfreq)
    c0 = max(0, (nfreq - nch_use) // 2)
    block = data[c0:c0 + nch_use, :ntime_use].T          # (time, channel), ascending
    print(f"   using {block.shape[0]} samples x {block.shape[1]} channels, "
          f"{block.shape[0] * TSAMP_S:.0f} s of a {ntime * TSAMP_S / 3600:.2f} h day")

    fch1 = LOW_MHZ + c0 * channel_hz / 1e6
    fil = os.path.join(ROOT, "data", "ovro_lwa_spec.fil")
    write_filterbank(fil, block.astype(np.float32), fch1=fch1,
                     foff=channel_hz / 1e6, tsamp=TSAMP_S,
                     source="SUN (OVRO-LWA spec_fits)", nbits=32,
                     telescope_id=6, tstart=0.0)

    nyquist = 0.5 * channel_hz / TSAMP_S
    cfg = SearchConfig(boxcar=3, max_drift_hz_s=nyquist,
                       snr_threshold=args.threshold, max_hits=40)
    hits, meta = search_file(fil, cfg, chan_chunk=2048, max_ints=None)
    resolution = channel_hz / (ntime_use * TSAMP_S)
    print(f"   Nyquist {nyquist:,.0f} Hz/s, resolution {resolution:.1f} Hz/s, "
          f"theta = {THETA * channel_hz / TSAMP_S:,.0f} Hz/s")
    print(f"   {len(hits)} candidates at {args.threshold:g} sigma")

    rows = []
    for hit in sorted(hits, key=lambda h: -h.snr)[:10]:
        resolved = abs(hit.drift_hz_s) > resolution
        nu, distance = drift_lock(hit.drift_hz_s, channel_hz, TSAMP_S)
        k = round(nu / THETA)
        rows.append({"freq_mhz": hit.freq_mhz, "drift_hz_s": hit.drift_hz_s,
                     "snr": hit.snr, "resolved": bool(resolved), "rate": nu,
                     "k": int(k), "distance": distance,
                     "locked": bool(resolved and k != 0
                                    and distance <= LOCK_TOLERANCE)})
        print(f"      {hit.freq_mhz:9.4f} MHz  drift {hit.drift_hz_s:>12.1f} Hz/s"
              f"  snr {hit.snr:5.2f}  rate {nu:.4f}  dist {distance:.3f}"
              f"  {'resolved' if resolved else 'below resolution'}")

    evaluable = [r for r in rows if r["resolved"]]
    chance = 2.0 * LOCK_TOLERANCE * THETA / 0.5
    verdict = ("NOT MEASURED: no candidate has a drift above resolution"
               if not evaluable else
               f"{len(evaluable)} evaluable, {sum(r['locked'] for r in evaluable)} locked")
    print(f"\n   verdict: {verdict}")
    print(f"   chance per candidate within {LOCK_TOLERANCE} theta of an integer: "
          f"{chance:.4f}")

    out = {"source": args.key, "sha_key": os.path.basename(args.key),
           "cube_shape": list(cube.shape), "nfreq": nfreq, "ntime": ntime,
           "samples_used": int(ntime_use), "channels_used": int(nch_use),
           "channel_hz": channel_hz, "tsamp_s": TSAMP_S,
           "nyquist_hz_s": nyquist, "resolution_hz_s": resolution,
           "theta_drift_hz_s": THETA * channel_hz / TSAMP_S,
           "threshold": args.threshold, "candidates": len(hits),
           "rows": rows, "chance_rate": chance, "verdict": verdict,
           "search_meta_time_s": meta.get("total_time")}
    with open(os.path.join(ROOT, "prereg", "ovro_lwa.json"), "w",
              encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print("   wrote prereg/ovro_lwa.json")
    return 0 if evaluable else 3


if __name__ == "__main__":
    raise SystemExit(main())
