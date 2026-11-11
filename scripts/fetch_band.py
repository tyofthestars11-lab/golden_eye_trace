"""Pull a band from an open archive, verify what it actually is, and refuse to
call it a second band unless it is far enough from the first.

The pre-registration needs two bands at least a decade apart. Nothing in this
project enforced that, so a plausible-looking download could quietly become "the
second band" while sitting at 1520 MHz, 0.00 decades from the 1405 MHz first
band. This script checks the file's own header, computes the separation, and
writes a ledger entry that either passes the gate or states the deficit.

    python3 scripts/fetch_band.py --list
    python3 scripts/fetch_band.py --check lwa-archive lofar-hba
    python3 scripts/fetch_band.py --get zenodo-vela --max-mb 64
    python3 scripts/fetch_band.py --gate

A SIGPROC file is read with this project's own parser; a FITS pulsar file is
read far enough to get OBSFREQ and OBSBW, because the one .sf file that
downloads turns out to be FITS and a .sf name is not evidence of a format.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import urllib.error
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bands import (  # noqa: E402
    MIN_DECADES,
    REFERENCE_MHZ,
    SOURCES,
    centre,
    decades,
    qualifies,
)
from seti_search.filterbank import Filterbank  # noqa: E402

LEDGER_DIR = os.path.join(ROOT, "bands")
DOWNLOAD_DIR = os.path.join(ROOT, "data", "bands")
USER_AGENT = "golden_eye_trace/1.0 (band fetch; python urllib)"
TIMEOUT = 30.0


def source(source_id: str) -> dict:
    for src in SOURCES:
        if src["id"] == source_id:
            return src
    raise SystemExit(
        f"unknown source {source_id!r}; try --list\n"
        + "\n".join(f"  {s['id']}" for s in SOURCES)
    )


# -- inspection ------------------------------------------------------------


def fits_band(path: str) -> tuple[float, float]:
    """(centre MHz, bandwidth MHz) from a FITS pulsar header."""
    cards: dict[str, str] = {}
    with open(path, "rb") as fh:
        while True:
            block = fh.read(2880)
            if not block:
                break
            text = block.decode("ascii", errors="replace")
            for i in range(0, len(text) - 79, 80):
                card = text[i:i + 80]
                key = card[:8].strip()
                if not key:
                    continue
                value = card[10:].split("/")[0].strip() if card[8:10] == "= " else ""
                cards.setdefault(key, value)
            if "END" in cards:
                break
    for key in ("OBSFREQ", "FREQ", "RESTFREQ"):
        if key in cards:
            try:
                center = float(cards[key])
            except ValueError:
                continue
            span = 0.0
            for span_key in ("OBSBW", "BW"):
                if span_key in cards:
                    try:
                        span = abs(float(cards[span_key]))
                    except ValueError:
                        span = 0.0
                    break
            return center, span
    raise ValueError("FITS header has no OBSFREQ/FREQ card")


def presto_inf(path: str) -> dict:
    """Band, channels and scan length from a PRESTO ``.inf`` sidecar.

    PRESTO writes the metadata as text next to the binary, so the band of an
    LWA DRX file is readable without touching the 15 MB payload. The sample
    interval here is the *search*-mode one; the folded series uses a different
    bin width, which is why the two are read separately.
    """
    text = open(path, encoding="utf-8", errors="replace").read()
    fields: dict[str, str] = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        fields[key.strip().lower()] = value.strip()

    def pick(*prefixes) -> float | None:
        for prefix in prefixes:
            for key, value in fields.items():
                if key.startswith(prefix):
                    try:
                        return float(value)
                    except ValueError:
                        return None
        return None

    low = pick("central freq of low channel")
    span = pick("total bandwidth")
    chan_bw = pick("channel bandwidth")
    bins = pick("number of bins")
    width = pick("width of each time series bin")
    if None in (low, span, chan_bw, bins, width) or width <= 0:
        raise ValueError(f"{path}: .inf lacks a required field")
    total_time = bins * width
    return {
        "fmt": "presto-inf",
        "center_mhz": low + span / 2,
        "span_mhz": span,
        "channel_khz": chan_bw * 1e3,
        "nchans": int(round(span / chan_bw)),
        "total_time_s": total_time,
        "drift_resolution_hz_s": chan_bw * 1e6 / total_time,
        "source_name": fields.get("object being observed"),
    }


def inspect(path: str) -> dict:
    """Identify a downloaded file from its own bytes and report its band."""
    info: dict = {"path": os.path.relpath(path, ROOT),
                  "bytes": os.path.getsize(path)}
    if path.endswith(".inf"):
        info.update(presto_inf(path))
        return _add_band_verdict(info)
    try:
        fb = Filterbank.open(path)
    except (ValueError, OSError):
        center, span = fits_band(path)
        info.update(
            fmt="fits-pulsar",
            center_mhz=center,
            span_mhz=span,
            channel_khz=None,
            total_time_s=None,
            drift_resolution_hz_s=None,
            source_name=None,
        )
    else:
        freqs = fb.frequencies()
        channel_hz = abs(fb.foff) * 1e6
        total_time = fb.total_time
        info.update(
            fmt="sigproc",
            center_mhz=float(0.5 * (freqs.min() + freqs.max())),
            span_mhz=float(freqs.max() - freqs.min()),
            channel_khz=channel_hz / 1e3,
            total_time_s=float(total_time) if total_time == total_time else None,
            drift_resolution_hz_s=(
                channel_hz / total_time if total_time and total_time > 0 else None
            ),
            source_name=fb.source,
            nchans=fb.nchans,
            nints=fb.nints,
        )
    return _add_band_verdict(info)


def _add_band_verdict(info: dict) -> dict:
    info["decades_from_first"] = decades(REFERENCE_MHZ, info["center_mhz"])
    info["qualifies"] = bool(info["decades_from_first"] >= MIN_DECADES)
    return info


# -- network --------------------------------------------------------------


def probe(url: str) -> dict:
    if url.startswith("local:"):
        path = os.path.join(ROOT, url.split(":", 1)[1])
        return {"status": 200 if os.path.exists(path) else 404,
                "bytes": os.path.getsize(path) if os.path.exists(path) else 0,
                "type": "local file", "error": None}
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
            return {
                "status": int(response.status),
                "bytes": int(response.headers.get("Content-Length") or 0),
                "type": response.headers.get("Content-Type", ""),
                "error": None,
            }
    except urllib.error.HTTPError as exc:
        return {"status": int(exc.code), "bytes": 0, "type": "", "error": str(exc)}
    except Exception as exc:  # connection refused, DNS, TLS, timeout
        return {"status": 0, "bytes": 0, "type": "", "error": f"{type(exc).__name__}: {exc}"}


def download(url: str, dest: str, max_bytes: int) -> tuple[int, str]:
    """Stream to disk, refusing anything over max_bytes; return (bytes, sha256)."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    digest = hashlib.sha256()
    written = 0
    partial = dest + ".part"
    with urllib.request.urlopen(request, timeout=TIMEOUT) as response, open(partial, "wb") as fh:
        while True:
            chunk = response.read(1 << 20)
            if not chunk:
                break
            written += len(chunk)
            if written > max_bytes:
                fh.close()
                os.unlink(partial)
                raise SystemExit(
                    f"refusing {url}: larger than --max-mb {max_bytes / 2**20:.0f}"
                )
            digest.update(chunk)
            fh.write(chunk)
    os.replace(partial, dest)
    return written, digest.hexdigest()


# -- commands -------------------------------------------------------------


def cmd_list() -> int:
    print(f"{'id':<15}{'telescope':<32}{'centre MHz':>11}{'decades':>9}{'ok':>5}")
    for src in SOURCES:
        c = centre(src)
        d = decades(REFERENCE_MHZ, c)
        print(f"{src['id']:<15}{src['telescope'][:31]:<32}{c:>11.1f}{d:>9.2f}"
              f"{'yes' if d >= MIN_DECADES else 'no':>5}")
    print()
    print("decades are measured from the STARE2 first band at "
          f"{REFERENCE_MHZ:.0f} MHz; a decade is a factor 10")
    return 0


def cmd_check(ids: list[str]) -> int:
    for src in (SOURCES if not ids else [source(i) for i in ids]):
        result = probe(src["url"])
        ok = result["status"] == 200
        detail = result["error"] or f"{result['type']} {result['bytes']} bytes"
        print(f"{src['id']:<15}{result['status']:>4}  "
              f"{'reachable' if ok else 'no':<10}{detail[:60]}")
        print(f"{'':<15}     {src['url']}")
    return 0


def cmd_get(source_id: str, max_mb: float, with_payload: bool) -> int:
    src = source(source_id)
    url = src["url"]
    os.makedirs(LEDGER_DIR, exist_ok=True)
    os.makedirs(DOWNLOAD_DIR, exist_ok=True)

    if url.startswith("local:"):
        info = inspect(os.path.join(ROOT, url.split(":", 1)[1]))
        info.update(id=source_id, sha256=None)
    elif src.get("meta_url"):
        # metadata-only source: the band comes from the sidecar, the payload
        # stays on the server unless asked for
        name = os.path.basename(src["meta_url"])
        meta_path = os.path.join(DOWNLOAD_DIR, f"{source_id}__{name}")
        meta_size, meta_sha = download(src["meta_url"], meta_path, 1 << 20)
        info = inspect(meta_path)
        info.update(id=source_id, meta_sha256=meta_sha, meta_bytes=meta_size,
                    payload_url=url, payload_readable=False)
        if with_payload:
            payload_name = os.path.basename(url.split("?")[0])
            payload_path = os.path.join(DOWNLOAD_DIR, f"{source_id}__{payload_name}")
            size, sha = download(url, payload_path, int(max_mb * 2**20))
            info.update(payload_bytes=size, payload_sha256=sha,
                        payload_path=os.path.relpath(payload_path, ROOT))
    else:
        name = os.path.basename(url.split("?")[0]) or f"{source_id}.bin"
        dest = os.path.join(DOWNLOAD_DIR, f"{source_id}__{name}")
        size, sha = download(url, dest, int(max_mb * 2**20))
        info = inspect(dest)
        info.update(id=source_id, sha256=sha, bytes=size)

    ledger = os.path.join(LEDGER_DIR, f"{source_id}.json")
    with open(ledger, "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(f"{source_id}: {info['fmt']}, {info['bytes']:,} bytes")
    print(f"  centre {info['center_mhz']:.1f} MHz, span {info['span_mhz']:.1f} MHz")
    if info.get("drift_resolution_hz_s"):
        print(f"  drift resolution {info['drift_resolution_hz_s']:.3g} Hz/s")
    print(f"  {info['decades_from_first']:.2f} decades from {REFERENCE_MHZ:.0f} MHz: "
          f"{'qualifies' if info['qualifies'] else 'too close'}")
    if not info["qualifies"]:
        need = REFERENCE_MHZ / 10 ** MIN_DECADES
        print(f"  a decade needs {need:.0f} MHz or below, or "
              f"{REFERENCE_MHZ * 10 ** MIN_DECADES:.0f} MHz or above")
    if info.get("payload_url"):
        state = "fetched" if info.get("payload_bytes") else "still on the server"
        print(f"  payload: {os.path.basename(info['payload_url'])} ({state}), "
              "PRESTO .dat")
        print("  the waterfall needs dat2fil (or equivalent) before the search can")
        print("  read it; --with-payload downloads it now")
    print(f"  ledger {os.path.relpath(ledger, ROOT)}")
    return 0 if info["qualifies"] else 1


def cmd_gate() -> int:
    entries = []
    if os.path.isdir(LEDGER_DIR):
        for name in sorted(os.listdir(LEDGER_DIR)):
            if not name.endswith(".json"):
                continue
            with open(os.path.join(LEDGER_DIR, name), encoding="utf-8") as fh:
                entries.append(json.load(fh))
    if not entries:
        print("no ledger entries: run --get on a source first")
        return 1

    print(f"{'id':<15}{'centre MHz':>11}{'decades':>9}{'ok':>5}")
    for entry in entries:
        print(f"{entry['id']:<15}{entry['center_mhz']:>11.1f}"
              f"{entry['decades_from_first']:>9.2f}"
              f"{'yes' if entry['qualifies'] else 'no':>5}")

    good = [e for e in entries if e["qualifies"]]
    print()
    if good:
        print(f"{len(good)} entry passes the decade requirement: "
              + ", ".join(e["id"] for e in good))
        print("the two-band test can now be attempted with scripts/phi_twoband.py")
        return 0
    best = max(entries, key=lambda e: e["decades_from_first"])
    print(f"nothing clears {MIN_DECADES:g} decade; the closest is {best['id']} at "
          f"{best['decades_from_first']:.2f}.")
    print("the pre-registered test stays NOT MEASURED until one does.")
    return 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="show the registry")
    parser.add_argument("--check", nargs="*", metavar="ID",
                        help="probe source URLs over HTTP")
    parser.add_argument("--get", metavar="ID", help="download and verify one source")
    parser.add_argument("--with-payload", action="store_true",
                        help="also fetch the data file behind a metadata sidecar")
    parser.add_argument("--max-mb", type=float, default=256.0)
    parser.add_argument("--gate", action="store_true", help="test the ledger")
    args = parser.parse_args()

    if args.list:
        return cmd_list()
    if args.check is not None:
        return cmd_check(args.check)
    if args.get:
        return cmd_get(args.get, args.max_mb, args.with_payload)
    if args.gate:
        return cmd_gate()
    print(__doc__.strip().splitlines()[0])
    print("  --list, --check, --get ID, --gate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
