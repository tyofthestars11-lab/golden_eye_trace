"""Minimal PSRFITS reader: primary header plus a SUBINT binary table.

Enough of PSRFITS to get a dynamic spectrum out without astropy or h5py: the
primary header, the ``SUBINT`` BINTABLE, the column layout from ``TFORMn``/
``TDIMn``, per-channel scaling from ``DAT_SCL``/``DAT_OFFS``, and the frequency
axis from ``DAT_FREQ``.

The one convention worth stating, because it is the usual source of silently
transposed data: ``TDIMn`` lists dimensions fastest-varying first, so the numpy
shape is that list reversed. The channel axis is then found by matching the
length of ``DAT_FREQ`` rather than assumed to be second.

    python3 scripts/psrfits.py --info data/f.fits
    python3 scripts/psrfits.py --to-fil data/f.fits data/f.fil
"""

from __future__ import annotations

import argparse
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from seti_search.filterbank import Filterbank, write_filterbank  # noqa: E402

CODE_SIZE = {"B": 1, "I": 2, "J": 4, "K": 8, "E": 4, "D": 8, "A": 1, "L": 1}
CODE_DTYPE = {"B": "u1", "I": "i2", "J": "i4", "K": "i8", "E": "f4", "D": "f8"}


def _value(text: str) -> str:
    text = text.strip()
    if text.startswith("'"):
        end = text.find("'", 1)
        return text[1:end] if end > 0 else text[1:]
    return text.split("/")[0].strip()


def scan_header(path: str, offset: int = 0) -> tuple[dict, int]:
    """Read one 2880-byte-blocked header; return (cards, next offset)."""
    cards: dict[str, str] = {}
    blocks = 0
    with open(path, "rb") as fh:
        fh.seek(offset)
        while True:
            block = fh.read(2880)
            if len(block) < 2880:
                raise ValueError(f"{path}: header at {offset} is truncated")
            blocks += 1
            for i in range(0, 2880, 80):
                card = block[i:i + 80].decode("ascii", "replace")
                if card.startswith("END"):
                    return cards, offset + blocks * 2880
                key = card[:8].strip()
                if not key or card[8:10] != "= ":
                    continue
                cards.setdefault(key, _value(card[10:]))
    raise ValueError("unreachable")


def as_int(cards: dict, key: str, default: int | None = None) -> int:
    if key not in cards:
        if default is None:
            raise ValueError(f"header lacks {key}")
        return default
    return int(float(cards[key]))


def column_specs(cards: dict) -> tuple[list[dict], int]:
    """Column names, dtypes, dims and byte offsets inside one table row."""
    specs: list[dict] = []
    offset = 0
    for i in range(1, as_int(cards, "TFIELDS") + 1):
        tform = cards.get(f"TFORM{i}", "")
        match = re.match(r"\s*(\d*)\s*([A-Z])", tform)
        if not match:
            raise ValueError(f"cannot parse TFORM{i}={tform!r}")
        repeat, code = int(match.group(1) or 1), match.group(2)
        size = repeat * CODE_SIZE[code]
        dims = None
        if f"TDIM{i}" in cards:
            dims = [int(v) for v in cards[f"TDIM{i}"].strip("()").split(",")]
        specs.append({
            "name": cards.get(f"TTYPE{i}", f"COL{i}").strip().upper(),
            "repeat": repeat,
            "code": code,
            "size": size,
            "dims": dims,
            "offset": offset,
        })
        offset += size
    return specs, offset


def read_column(raw: np.ndarray, spec: dict, nrows: int) -> np.ndarray:
    chunk = np.ascontiguousarray(raw[:, spec["offset"]:spec["offset"] + spec["size"]])
    # FITS tables are big-endian, whatever the host's byte order
    values = np.frombuffer(chunk.tobytes(), dtype=">" + CODE_DTYPE[spec["code"]],
                           count=nrows * spec["repeat"])
    shape = tuple(reversed(spec["dims"])) if spec["dims"] else (spec["repeat"],)
    return values.reshape((nrows,) + shape)


def load(path: str, reduce: str = "mean") -> dict:
    """Read a PSRFITS file as a (time, channel) dynamic spectrum."""
    primary, off1 = scan_header(path, 0)
    subint, off2 = scan_header(path, off1)
    if subint.get("XTENSION", "").strip() != "BINTABLE":
        raise ValueError(f"{path}: extension 1 is not a BINTABLE")
    if subint.get("EXTNAME", "").strip().upper() != "SUBINT":
        raise ValueError(f"{path}: extension is {subint.get('EXTNAME')!r}, not SUBINT")
    if as_int(subint, "PCOUNT", 0) != 0:
        raise ValueError("heap (PCOUNT != 0) is not supported")

    nrows = as_int(subint, "NAXIS2")
    na1 = as_int(subint, "NAXIS1")
    specs, row_size = column_specs(subint)
    if row_size != na1:
        raise ValueError(f"column sizes sum to {row_size}, NAXIS1 is {na1}")

    by_name = {spec["name"]: spec for spec in specs}
    for needed in ("DATA", "DAT_FREQ"):
        if needed not in by_name:
            raise ValueError(f"{path}: no {needed} column")

    raw = np.memmap(path, dtype="u1", mode="r", offset=off2,
                    shape=(nrows, na1), order="C")

    freq_spec = by_name["DAT_FREQ"]
    freqs_all = read_column(raw, freq_spec, nrows).astype(np.float64)
    if freq_spec["dims"] and len(freq_spec["dims"]) > 1:
        freqs_all = freqs_all.reshape(nrows, -1)
    nchan = freqs_all.shape[1]

    data_spec = by_name["DATA"]
    data = read_column(raw, data_spec, nrows)
    inner = data.shape[1:]
    channels = [i for i, size in enumerate(inner) if size == nchan]
    if not channels:
        raise ValueError(f"DATA {inner} has no axis matching {nchan} channels")
    chan_axis = channels[0]
    moved = np.moveaxis(data, 1 + chan_axis, -1)
    nsub = int(np.prod(moved.shape[1:-1]))
    matrix = moved.reshape(nrows, nsub, nchan).astype(np.float64)

    if "DAT_SCL" in by_name and "DAT_OFFS" in by_name:
        scl = read_column(raw, by_name["DAT_SCL"], nrows).astype(np.float64)
        offs = read_column(raw, by_name["DAT_OFFS"], nrows).astype(np.float64)
        matrix = matrix * scl[:, None, :] + offs[:, None, :]

    stacked = matrix.mean(axis=1) if reduce == "mean" else matrix[:, 0, :]

    tsamp = None
    row_tsamp = None
    if "TSUBINT" in by_name:
        tsub = read_column(raw, by_name["TSUBINT"], nrows).astype(np.float64).ravel()
        # TSUBINT is the duration of a whole row; the stored values inside the
        # row are phase bins, so a row-mean spectrum advances by TSUBINT, not by
        # a single bin. Dividing by nsub here would lose a factor of nsub.
        row_tsamp = float(np.median(tsub))
        tsamp = row_tsamp / max(nsub, 1)

    freqs = freqs_all[0]
    return {
        "path": path,
        "rows": int(nrows),
        "channels": int(nchan),
        "nsub": nsub,
        "freqs_mhz": freqs,
        "low_mhz": float(freqs.min()),
        "high_mhz": float(freqs.max()),
        "channel_khz": float(abs(freqs[1] - freqs[0]) * 1e3) if nchan > 1 else 0.0,
        "tsamp_s": tsamp,
        "row_tsamp_s": row_tsamp,
        "total_time_s": (row_tsamp * nrows) if row_tsamp else None,
        "data": stacked.astype(np.float32),
        "primary": primary,
        "subint": subint,
    }


def resolution(info: dict) -> float | None:
    if not info["total_time_s"] or info["channel_khz"] <= 0:
        return None
    return info["channel_khz"] * 1e3 / info["total_time_s"]


def to_filterbank(path: str, out_path: str, reduce: str = "mean") -> dict:
    info = load(path, reduce)
    data = info["data"]
    if info["channels"] < 2:
        raise ValueError("need at least two channels to be a band")
    tsamp = info.get("row_tsamp_s") or info["tsamp_s"]
    if not tsamp or tsamp <= 0:
        raise ValueError("no usable TSUBINT: cannot determine the sample interval")
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    write_filterbank(
        out_path,
        data,
        fch1=float(info["freqs_mhz"][0]),
        foff=float(info["freqs_mhz"][1] - info["freqs_mhz"][0]),
        tsamp=float(tsamp),
        source=str(info["primary"].get("SRC_NAME", "PSRFITS"))[:32],
        nbits=32,
    )
    fb = Filterbank.open(out_path)
    back = fb.read()
    return {
        **info,
        "out": out_path,
        "bytes": os.path.getsize(out_path),
        "round_trip_close": bool(np.allclose(back, data, rtol=1e-6, atol=1e-4)),
        "filterbank_time_s": fb.total_time,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fits", nargs="?")
    parser.add_argument("--info", action="store_true")
    parser.add_argument("--to-fil")
    parser.add_argument("--reduce", choices=("mean", "first"), default="mean")
    args = parser.parse_args()
    if not args.fits:
        parser.print_help()
        return 2

    info = load(args.fits, args.reduce)
    print(f"   {os.path.basename(args.fits)}")
    print(f"   {info['rows']} rows x {info['channels']} channels "
          f"({info['nsub']} values per channel per row)")
    print(f"   {info['low_mhz']:.4f} - {info['high_mhz']:.4f} MHz, "
          f"channel {info['channel_khz']:.4f} kHz")
    if info["tsamp_s"]:
        print(f"   tsamp {info['tsamp_s']:.6g} s, total {info['total_time_s']:.2f} s, "
              f"drift resolution {resolution(info):.3g} Hz/s")
    data = info["data"]
    print(f"   data {data.shape}, range {np.nanmin(data):.5g} .. {np.nanmax(data):.5g}")

    if args.to_fil:
        result = to_filterbank(args.fits, args.to_fil, args.reduce)
        print(f"   wrote {result['out']} ({result['bytes']:,} bytes), "
              f"{result['filterbank_time_s']:.2f} s, "
              f"round trip close: {result['round_trip_close']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
