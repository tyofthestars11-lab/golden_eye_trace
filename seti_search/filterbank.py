"""SIGPROC filterbank (.fil) I/O and a synthetic file generator.

Only the subset a narrowband technosignature search needs: header parsing, a
memory-mapped data view, a writer, and an injector for recovery tests.

Layout (SIGPROC v3.x), all little endian:

    HEADER_START
    int32   header_size        optional; 3093750 in SIGPROC output, absent in some files
    key/value pairs, each key a length-prefixed ASCII string:
      source_name(str) telescope_id(i) machine_id(i) data_type(i)
      barycentric(i) pulsarcentric(i) az_start(d) za_start(d) src_raj(d)
      src_dej(d) tstart(d) tsamp(d) nbits(i) fch1(d) foff(d) nchans(i) nifs(i)
    HEADER_END
    nints records of nchans * nifs * nbits / 8 bytes

Real files vary in two ways that are handled here: the HEADER_START/HEADER_END
markers are sometimes written length-prefixed and sometimes raw, and a few
telescopes add extra keys. Unknown keys cannot be typed, so they raise a
message naming the key and byte offset rather than silently misaligning.

Sample widths: 8 and 16 bit are read as unsigned integers, 32 bit as IEEE float
(what Breakthrough Listen filterbank products use), 64 bit as double.
"""

from __future__ import annotations

import os
import struct
from dataclasses import dataclass, field, replace

import numpy as np

HEADER_MAGIC = 3093750

INT_KEYS = {
    "header_size",
    "telescope_id",
    "machine_id",
    "data_type",
    "barycentric",
    "pulsarcentric",
    "nbits",
    "nchans",
    "nifs",
    "nbeams",
    "ibeam",
    "srclen",
    # seen in STARE2 candidate files
    "start_sample",
    "end_sample",
    "cand_location",
    "MJD_start",
    "MJD_hour",
    "MJD_minute",
}

FLOAT_KEYS = {
    "az_start",
    "za_start",
    "src_raj",
    "src_dej",
    "tstart",
    "tsamp",
    "fch1",
    "foff",
    "refdm",
    "period",
    "fchannel",
    "epoch",
    "MJD_second",
}

STRING_KEYS = {
    "source_name",
    "srcname",
    "pol_type",
    "telescope",
}

DTYPES = {8: np.uint8, 16: np.uint16, 32: np.float32, 64: np.float64}


class _Cursor:
    """Little-endian reader over an in-memory header buffer."""

    def __init__(self, buf: bytes):
        self.buf = buf
        self.pos = 0

    def take(self, n: int) -> bytes:
        if self.pos + n > len(self.buf):
            raise ValueError(
                f"filterbank header is truncated at byte {self.pos} "
                f"(need {n} more bytes, have {len(self.buf) - self.pos})"
            )
        chunk = self.buf[self.pos:self.pos + n]
        self.pos += n
        return chunk

    def int32(self) -> int:
        return struct.unpack("<i", self.take(4))[0]

    def float64(self) -> float:
        return struct.unpack("<d", self.take(8))[0]

    def string(self) -> str:
        n = self.int32()
        if not 0 <= n <= 8192:
            raise ValueError(
                f"implausible string length {n} at byte {self.pos - 4}: "
                "this does not look like a SIGPROC filterbank header"
            )
        return self.take(n).decode("ascii", "replace")


def parse_header(buf: bytes) -> tuple[dict, int]:
    """Return (header dict, byte offset where the data records start).

    Two header styles are in the wild and both are accepted:

    * key/value pairs, where every field is a key string followed by a value:
      SIGPROC's own `filterbank` writer, and the STARE2 candidate files;
    * the fixed positional sequence used by GUPPI and Parkes products - source
      name first with no key, then telescope_id, machine_id, ... nifs - which is
      the layout blimpy's reader assumes.

    The keyed form is tried first, so a genuinely unknown key is reported as an
    error instead of being misread as a source name.
    """
    prefix, body_start = _split_prefix(buf)

    keyed = _Cursor(buf)
    keyed.pos = body_start
    try:
        parsed, offset = _keyed_fields(keyed)
        return {**prefix, **parsed}, offset
    except ValueError as keyed_error:
        positional = _Cursor(buf)
        positional.pos = body_start
        try:
            parsed, offset = _positional_fields(positional)
        except ValueError:
            raise keyed_error from None
        return {**prefix, **parsed}, offset


def _split_prefix(buf: bytes) -> tuple[dict, int]:
    cur = _Cursor(buf)
    if buf[:4] == b"HEAD":
        cur.take(len(b"HEADER_START"))
    else:
        marker = cur.string()
        if not marker.startswith("HEADER_START"):
            raise ValueError(
                f"not a SIGPROC filterbank file: expected HEADER_START, got {marker!r}"
            )

    prefix: dict = {}
    # SIGPROC sometimes writes the header length here; other writers (the STARE2
    # files, for instance) go straight into key/value pairs. A key length prefix
    # is small and a real header length is not, so the cases cannot be confused.
    if len(buf) >= cur.pos + 4:
        candidate = struct.unpack("<i", buf[cur.pos:cur.pos + 4])[0]
        if candidate == HEADER_MAGIC or candidate == 0 or candidate > 4096:
            prefix["header_size"] = cur.int32()
    return prefix, cur.pos


def _keyed_fields(cur: "_Cursor") -> tuple[dict, int]:
    header: dict = {}
    while True:
        if cur.buf[cur.pos:cur.pos + len(b"HEADER_END")] == b"HEADER_END":
            cur.take(len(b"HEADER_END"))
            break
        key = cur.string()
        if key == "HEADER_END":
            break
        if key in STRING_KEYS:
            value = cur.string()
        elif key in INT_KEYS:
            value = cur.int32()
        elif key in FLOAT_KEYS:
            value = cur.float64()
        else:
            raise ValueError(
                f"unknown filterbank header key {key!r} at byte {cur.pos}; "
                f"known keys: {sorted(STRING_KEYS | INT_KEYS | FLOAT_KEYS)}"
            )
        header[key] = value
    return header, cur.pos


def _positional_fields(cur: "_Cursor") -> tuple[dict, int]:
    """Fixed-order header: source name first, with no key name."""
    header: dict = {"source_name": cur.string()}
    header["telescope_id"] = cur.int32()
    header["machine_id"] = cur.int32()
    header["data_type"] = cur.int32()
    header["barycentric"] = cur.int32()
    header["pulsarcentric"] = cur.int32()
    for key in ("az_start", "za_start", "src_raj", "src_dej", "tstart", "tsamp"):
        header[key] = cur.float64()
    header["nbits"] = cur.int32()
    header["fch1"] = cur.float64()
    header["foff"] = cur.float64()
    header["nchans"] = cur.int32()
    header["nifs"] = cur.int32()

    # a positional header that decoded to nonsense must not be accepted just
    # because the bytes happened to line up
    if header["nchans"] <= 0 or header["nbits"] not in DTYPES:
        raise ValueError("positional header did not validate")

    extra, offset = _keyed_fields(cur)
    header.update(extra)
    return header, offset


@dataclass
class Filterbank:
    """A read-only view of a filterbank file."""

    path: str
    header: dict
    nints: int
    nchans: int
    nifs: int
    nbits: int
    tsamp: float
    fch1: float
    foff: float
    data_offset: int
    file_size: int

    @classmethod
    def open(cls, path) -> "Filterbank":
        path = str(path)
        size = os.path.getsize(path)
        with open(path, "rb") as fh:
            buf = fh.read(min(size, 1 << 20))
        header, offset = parse_header(buf)

        missing = [k for k in ("nchans", "nbits", "fch1", "foff") if k not in header]
        if missing:
            raise ValueError(f"filterbank header is missing {missing}")

        nchans = int(header["nchans"])
        nbits = int(header["nbits"])
        if nbits not in DTYPES:
            raise ValueError(f"unsupported nbits={nbits}; expected one of {sorted(DTYPES)}")
        nifs = int(header.get("nifs", 1) or 1)
        record = nchans * nifs * nbits // 8
        if record <= 0:
            raise ValueError(f"bad record size (nchans={nchans}, nifs={nifs}, nbits={nbits})")
        nints = (size - offset) // record
        if nints <= 0:
            raise ValueError(
                f"no complete integration in {path}: header ends at byte {offset}, "
                f"file is {size} bytes, each integration is {record} bytes"
            )
        return cls(
            path=path,
            header=header,
            nints=int(nints),
            nchans=nchans,
            nifs=nifs,
            nbits=nbits,
            tsamp=float(header.get("tsamp", float("nan"))),
            fch1=float(header["fch1"]),
            foff=float(header["foff"]),
            data_offset=int(offset),
            file_size=int(size),
        )

    # -- properties -------------------------------------------------------
    @property
    def source(self) -> str:
        return str(self.header.get("source_name", "unknown")).strip() or "unknown"

    @property
    def tstart(self) -> float:
        return float(self.header.get("tstart", 0.0))

    @property
    def total_time(self) -> float:
        return self.nints * self.tsamp

    @property
    def bandwidth(self) -> float:
        return abs(self.foff) * self.nchans

    def frequencies(self) -> np.ndarray:
        """Channel centre frequencies in MHz, in stored channel order.

        foff is negative for most telescopes, so index 0 is usually the
        highest frequency in the band.
        """
        return self.fch1 + self.foff * np.arange(self.nchans, dtype=np.float64)

    def describe(self) -> str:
        if np.isfinite(self.tsamp) and self.tsamp > 0:
            span = f"{self.nints:,}  ({self.total_time:g} s at {self.tsamp:g} s)"
        else:
            span = f"{self.nints:,}  (tsamp missing from header; pass --tsamp)"
        lines = [
            f"file            {self.path}",
            f"source          {self.source}",
            f"telescope_id    {self.header.get('telescope_id', '?')}",
            f"channels        {self.nchans:,}",
            f"integrations    {span}",
            f"bits/sample     {self.nbits}",
            f"polarisations    {self.nifs}",
            f"channel width   {abs(self.foff) * 1e3:.6g} kHz  (foff={self.foff:+.6g} MHz)",
            f"bandwidth       {self.bandwidth * 1e3:.6g} kHz",
        ]
        if self.nchans > 1:
            freqs = self.frequencies()
            lines.append(f"frequency range {freqs.min():.6f} - {freqs.max():.6f} MHz")
        lines.append(f"data starts     byte {self.data_offset:,} of {self.file_size:,}")
        return "\n".join(lines)

    # -- data access ------------------------------------------------------
    def read(self, start: int = 0, stop: int | None = None,
             chan_start: int = 0, chan_stop: int | None = None) -> np.ndarray:
        """Read a (n_ints, n_channels) float32 block; polarisations are averaged."""
        stop = self.nints if stop is None else min(int(stop), self.nints)
        chan_stop = self.nchans if chan_stop is None else min(int(chan_stop), self.nchans)
        if not 0 <= start < stop:
            raise ValueError(f"bad integration range {start}:{stop}")
        if not 0 <= chan_start < chan_stop:
            raise ValueError(f"bad channel range {chan_start}:{chan_stop}")

        dtype = np.dtype(DTYPES[self.nbits]).newbyteorder("<")
        view = np.memmap(
            self.path,
            dtype=dtype,
            mode="r",
            offset=self.data_offset,
            shape=(self.nints, self.nifs, self.nchans),
        )
        block = np.asarray(view[start:stop, :, chan_start:chan_stop], dtype=np.float32)
        if self.nifs > 1:
            block = block.mean(axis=1)
        else:
            block = block[:, 0, :]
        return np.ascontiguousarray(block, dtype=np.float32)


def _str(value: str) -> bytes:
    raw = value.encode("ascii")
    return struct.pack("<i", len(raw)) + raw


def write_filterbank(
    path,
    data: np.ndarray,
    *,
    fch1: float,
    foff: float,
    tsamp: float,
    source: str = "SYNTHETIC",
    nbits: int = 32,
    tstart: float = 0.0,
    telescope_id: int = 6,
    machine_id: int = 0,
) -> int:
    """Write a (n_ints, n_channels) array as a SIGPROC filterbank file."""
    data = np.asarray(data)
    if data.ndim != 2:
        raise ValueError(f"expected a 2-D array, got shape {data.shape}")
    if nbits not in DTYPES:
        raise ValueError(f"unsupported nbits={nbits}")
    nints, nchans = data.shape

    header = bytearray()
    header += _str("HEADER_START")
    header += _str(source)
    header += struct.pack("<i", telescope_id)
    header += struct.pack("<i", machine_id)
    header += struct.pack("<i", 1)  # data_type: 1 = filterbank
    header += struct.pack("<i", 0)  # barycentric
    header += struct.pack("<i", 0)  # pulsarcentric
    header += struct.pack("<d", 0.0)  # az_start
    header += struct.pack("<d", 0.0)  # za_start
    header += struct.pack("<d", 0.0)  # src_raj
    header += struct.pack("<d", 0.0)  # src_dej
    header += struct.pack("<d", float(tstart))
    header += struct.pack("<d", float(tsamp))
    header += struct.pack("<i", int(nbits))
    header += struct.pack("<d", float(fch1))
    header += struct.pack("<d", float(foff))
    header += struct.pack("<i", int(nchans))
    header += struct.pack("<i", 1)  # nifs
    header += _str("HEADER_END")

    payload = data.astype(np.dtype(DTYPES[nbits]).newbyteorder("<"), copy=False).tobytes()
    with open(path, "wb") as fh:
        fh.write(bytes(header))
        fh.write(payload)
    return len(header)


@dataclass
class InjectedSignal:
    """A synthetic narrowband signal, in physical units."""

    freq_mhz: float
    drift_hz_s: float = 0.0
    amplitude: float = 1.0
    width_channels: float = 1.5

    def as_dict(self) -> dict:
        return {
            "freq_mhz": self.freq_mhz,
            "drift_hz_s": self.drift_hz_s,
            "amplitude": self.amplitude,
            "width_channels": self.width_channels,
        }


@dataclass
class RFILine:
    """A stationary carrier: constant in time and frequency."""

    channel: int
    amplitude: float = 40.0


def simulate(
    path,
    *,
    nchans: int = 8192,
    nints: int = 200,
    tsamp: float = 1.0,
    fch1: float = 1420.4,
    foff: float = -2.86e-6,
    signals: tuple[InjectedSignal, ...] = (),
    rfi: tuple[RFILine, ...] = (),
    bandpass_ripple: float = 0.2,
    glitch: bool = True,
    seed: int = 0,
    nbits: int = 32,
    sky: tuple[float, float] | None = None,
    observer=None,
    mjd: float | None = None,
) -> dict:
    """Generate a filterbank file with known signals buried in noise.

    Noise is exponential (total-power statistics: mean 1, std 1) so that
    dividing by the bandpass leaves unit-variance fluctuations, which is what
    the search assumes. A smooth bandpass gain is applied so bandpass control
    is exercised, and an optional broadband glitch row is injected so row
    excision is exercised.

    Pass sky=(ra_deg, dec_deg) together with observer and mjd to make each
    signal's drift_hz_s the source's *barycentric* drift: the file then records
    the topocentric drift a telescope at that site would measure at that epoch,
    which is what an ON/OFF cadence test needs.
    """
    if (sky is None) != (observer is None):
        raise ValueError("sky and observer must be given together")
    if sky is not None and mjd is None:
        raise ValueError("sky and observer need mjd to place the observation")
    # mjd on its own is fine: it just records the epoch in the header, which is
    # what a run needs in order to link scans
    rng = np.random.default_rng(seed)
    data = rng.exponential(size=(nints, nchans)).astype(np.float32)

    phase = np.linspace(0.0, 3.0 * np.pi, nchans)
    gain = (1.0 + bandpass_ripple * np.sin(phase) + 0.5 * bandpass_ripple * np.cos(2.7 * phase))
    gain = gain.astype(np.float32)
    data *= gain[None, :]

    freqs = fch1 + foff * np.arange(nchans, dtype=np.float64)
    times = np.arange(nints, dtype=np.float64) * tsamp

    if sky is not None and observer is not None and mjd is not None:
        from .cadence import SkyPosition, topocentric_drift

        position = SkyPosition(float(sky[0]), float(sky[1]))
        signals = tuple(
            replace(
                sig,
                drift_hz_s=topocentric_drift(
                    sig.drift_hz_s, sig.freq_mhz * 1e6, observer, mjd, position
                ),
            )
            for sig in signals
        )

    for sig in signals:
        centre = sig.freq_mhz + sig.drift_hz_s * 1e-6 * times  # MHz, one per integration
        sigma_mhz = max(abs(sig.width_channels * foff), 1e-9)
        profile = np.exp(-0.5 * ((freqs[None, :] - centre[:, None]) / sigma_mhz) ** 2)
        data += (sig.amplitude * profile).astype(np.float32)

    for line in rfi:
        data[:, int(line.channel)] += float(line.amplitude)

    if glitch and nints > 10:
        rows = slice(nints // 3, nints // 3 + 2)
        data[rows, :] *= 5.0

    header_bytes = write_filterbank(
        path, data, fch1=fch1, foff=foff, tsamp=tsamp,
        source="SYNTHETIC", nbits=nbits, telescope_id=6,
        tstart=0.0 if mjd is None else float(mjd),
    )
    return {
        "path": str(path),
        "nchans": nchans,
        "nints": nints,
        "tsamp": tsamp,
        "fch1": fch1,
        "foff": foff,
        "nbits": nbits,
        "header_bytes": header_bytes,
        "signals": [s.as_dict() for s in signals],
        "rfi": [{"channel": r.channel, "amplitude": r.amplitude} for r in rfi],
        "channel_freqs_mhz": (float(freqs[0]), float(freqs[-1])),
    }
