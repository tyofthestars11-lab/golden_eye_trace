"""Read a PRESTO ``.dat``/``.inf`` pair, decide what shape it is, and convert it
to SIGPROC.

PRESTO's ``.dat`` payload is plain little-endian single-precision floats; the
text ``.inf`` beside it holds the band, the channel count and the sample
interval. The one question the format does not answer on its own is whether
those floats are a single time series or a waterfall of ``nchans`` columns per
sample, and the answer changes what the file can be used for:

* a time series has no frequency axis, so no drift can be measured in it;
* a waterfall does, and that is what a dedoppler search consumes.

This module decides from the numbers rather than from the file name: the payload
size must equal either the ``.inf`` bin count (time series) or an integer
multiple of the channel count (waterfall), and the surviving hypothesis has to be
consistent with the stated observation length. It prints both hypotheses with
their evidence so the choice is auditable, and it refuses to convert when
neither closes.

    python3 scripts/presto_dat.py --dat f.dat --inf f.inf --inspect
    python3 scripts/presto_dat.py --dat f.dat --inf f.inf --out data/lwa.fil

The conversion is checked by reading the written file back with this project's
own SIGPROC parser and comparing every sample.
"""

from __future__ import annotations

import argparse
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from seti_search.filterbank import Filterbank, write_filterbank  # noqa: E402

FLOAT = np.dtype("<f4")


def read_inf(path: str) -> dict:
    """Parse a PRESTO ``.inf`` sidecar into typed values."""
    fields: dict[str, str] = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        fields[key.strip().lower()] = value.strip()

    def number(*prefixes, default=None):
        for prefix in prefixes:
            for key, value in fields.items():
                if key.startswith(prefix):
                    try:
                        return float(value)
                    except ValueError:
                        return default
        return default

    low = number("central freq of low channel")
    span = number("total bandwidth")
    chan_bw = number("channel bandwidth")
    bins = number("number of bins")
    width = number("width of each time series bin")
    if None in (low, span, chan_bw, bins, width):
        raise ValueError(f"{path}: .inf is missing a required field")
    return {
        "low_mhz": low,
        "span_mhz": span,
        "chan_bw_mhz": chan_bw,
        "bins": int(round(bins)),
        "bin_width_s": width,
        "centre_mhz": low + span / 2.0,
        "nchans": int(round(span / chan_bw)),
        "source": fields.get("object being observed", "unknown"),
        "telescope": fields.get("telescope used", "?"),
        "instrument": fields.get("instrument used", "?"),
        "mjd": number("epoch of observation"),
    }


def lag1(values: np.ndarray) -> float:
    if values.size < 3:
        return float("nan")
    v = values.astype(np.float64)
    if v.std() == 0:
        return float("nan")
    return float(np.corrcoef(v[:-1], v[1:])[0, 1])


def hypotheses(raw: np.ndarray, info: dict) -> dict:
    """Which shapes the payload size allows, with the evidence for each."""
    total_time = info["bins"] * info["bin_width_s"]
    out: dict = {"n_floats": int(raw.size), "total_time_s": total_time}

    as_series = raw.size == info["bins"]
    out["timeseries"] = {
        "fits": bool(as_series),
        "lag1": lag1(raw) if as_series else None,
    }

    nch = info["nchans"]
    as_waterfall = raw.size % nch == 0
    out["waterfall"] = {"fits": bool(as_waterfall), "nchans": nch}
    if as_waterfall:
        nsamp = raw.size // nch
        wf = raw.reshape(nsamp, nch)
        col = wf.mean(axis=0)
        out["waterfall"].update(
            n_samples=int(nsamp),
            sample_time_s=total_time / nsamp if nsamp else None,
            per_channel_lag1=float(np.nanmedian([
                lag1(wf[:, c]) for c in range(0, nch, max(1, nch // 16))
            ])),
            channel_spread=float(col.std() / col.mean()) if col.mean() else None,
        )
    return out


def load(dat_path: str, inf_path: str, shape: str | None = None):
    """Return (2-D float32 array, info) with columns as frequency channels."""
    info = read_inf(inf_path)
    raw = np.fromfile(dat_path, dtype=FLOAT)
    if raw.size == 0:
        raise ValueError(f"{dat_path}: no floats")
    if not np.all(np.isfinite(raw)):
        raise ValueError(f"{dat_path}: payload contains non-finite values")

    total_time = info["bins"] * info["bin_width_s"]
    choice = shape
    if choice is None:
        if raw.size == info["bins"]:
            choice = "timeseries"
        elif raw.size % info["nchans"] == 0:
            choice = "waterfall"
        else:
            raise ValueError(
                f"{dat_path}: {raw.size} floats is neither {info['bins']} bins "
                f"(time series) nor a multiple of {info['nchans']} channels"
            )

    if choice == "timeseries":
        if raw.size != info["bins"]:
            raise ValueError(
                f"{dat_path}: asked for a time series but holds {raw.size} "
                f"samples, not {info['bins']}"
            )
        data = raw.reshape(-1, 1)
        tsamp = info["bin_width_s"]
    else:
        if raw.size % info["nchans"]:
            raise ValueError(
                f"{dat_path}: {raw.size} floats is not a multiple of "
                f"{info['nchans']} channels"
            )
        data = raw.reshape(-1, info["nchans"])
        tsamp = total_time / data.shape[0]

    return data, {
        **info,
        "shape": choice,
        "tsamp_s": tsamp,
        "nints": int(data.shape[0]),
        "ncols": int(data.shape[1]),
    }


def convert(dat_path: str, inf_path: str, out_path: str,
            shape: str | None = None) -> dict:
    data, info = load(dat_path, inf_path, shape)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    write_filterbank(
        out_path,
        data.astype(np.float32),
        fch1=info["low_mhz"],
        foff=info["chan_bw_mhz"],
        tsamp=info["tsamp_s"],
        source=f"{info['source']} {info['instrument']}".strip()[:32],
        nbits=32,
    )

    # read it back with this project's parser and compare every sample
    fb = Filterbank.open(out_path)
    back = fb.read()
    exact = back.shape == data.shape and bool(np.array_equal(back, data.astype(np.float32)))
    return {
        "out": out_path,
        "nints": int(data.shape[0]),
        "ncols": int(data.shape[1]),
        "tsamp_s": info["tsamp_s"],
        "total_time_s": fb.total_time,
        "frequencies_mhz": (float(np.min(fb.frequencies())), float(np.max(fb.frequencies()))),
        "bytes": os.path.getsize(out_path),
        "round_trip_exact": exact,
        "max_abs_diff": float(np.max(np.abs(back - data))) if back.shape == data.shape else None,
    }


def print_inspection(dat_path: str, inf_path: str) -> dict:
    info = read_inf(inf_path)
    raw = np.fromfile(dat_path, dtype=FLOAT)
    h = hypotheses(raw, info)

    print(f"   {os.path.basename(dat_path)}: {h['n_floats']:,} float32 samples")
    print(f"   .inf: {info['low_mhz']:.3f} MHz + {info['span_mhz']:.1f} MHz, "
          f"{info['nchans']} channels of {info['chan_bw_mhz'] * 1e3:.4f} kHz")
    print(f"   .inf: {info['bins']:,} bins at {info['bin_width_s'] * 1e6:.4f} us "
          f"= {h['total_time_s']:.2f} s")
    print(f"   payload value range {raw.min():.4g} .. {raw.max():.4g}, "
          f"median {np.median(raw):.4g}")

    ts = h["timeseries"]
    print(f"   time series?  floats == bins: {ts['fits']}"
          + (f", lag-1 corr {ts['lag1']:+.4f}" if ts["fits"] else ""))
    wf = h["waterfall"]
    if wf["fits"]:
        print(f"   waterfall?    {wf['n_samples']:,} rows x {wf['nchans']} channels, "
              f"sample time {wf['sample_time_s']:.6g} s")
        print(f"                 median per-channel lag-1 corr "
              f"{wf['per_channel_lag1']:+.4f}, channel-mean spread "
              f"{wf['channel_spread']:.4f}")
    else:
        print(f"   waterfall?    no: {h['n_floats']} is not a multiple of "
              f"{wf['nchans']} channels")
    return h


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dat", required=True)
    parser.add_argument("--inf", required=True)
    parser.add_argument("--out")
    parser.add_argument("--shape", choices=("timeseries", "waterfall"))
    parser.add_argument("--inspect", action="store_true")
    args = parser.parse_args()

    if args.inspect or not args.out:
        h = print_inspection(args.dat, args.inf)
        if not args.out:
            print()
            print("   nothing written: pass --out to convert")
            return 0
        print()

    result = convert(args.dat, args.inf, args.out, args.shape)
    print(f"   wrote {result['out']}  ({result['bytes']:,} bytes)")
    print(f"   {result['nints']:,} integrations x {result['ncols']} channels, "
          f"{result['tsamp_s']:.6g} s each = {result['total_time_s']:.2f} s")
    print(f"   frequencies {result['frequencies_mhz'][0]:.6f} - "
          f"{result['frequencies_mhz'][1]:.6f} MHz")
    print(f"   round trip exact: {result['round_trip_exact']} "
          f"(max |diff| {result['max_abs_diff']:.3g})")
    return 0 if result["round_trip_exact"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
