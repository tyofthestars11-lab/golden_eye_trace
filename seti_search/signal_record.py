"""A typed record for one candidate, with the six layers kept separate.

Every field answers one of the questions the framework refuses to merge: what was
measured, how significant it is when the noise is measured locally, what the
transform did, and what we are allowed to conclude. Nothing here computes a
disposition; `disposition` is supplied by the pipeline that earned it, so a record
cannot claim more than the run that produced it.

The record is frozen: a candidate's provenance is not something to edit in place
after seeing a later result. Revise by building a new record with a new
transform_version.
"""

from __future__ import annotations

import datetime as dt
import math
from dataclasses import asdict, dataclass, field

# The vocabulary is closed on purpose. A disposition outside this set is a bug in
# the caller, not a new category, and is rejected rather than stored.
DISPOSITIONS = (
    "CONFIRMED TECHNOSIGNATURE CANDIDATE",
    "NONLINEAR-COHERENT BUT UNRESOLVED",
    "ASTROPHYSICAL",
    "TERRESTRIAL RFI",
    "INSTRUMENTAL ALIAS",
    "PIPELINE ARTIFACT",
    "UNRESOLVED",
    "NULL RESULT",
)

ON_SOURCE = ("on-source", "off-source", "untestable")


def channel_rate(drift_hz_per_s: float, channel_hz: float, tsamp_s: float) -> float:
    """Drift in channels per sample, which is the quantity a dedoppler can see.

    A rate beyond 0.5 folds: 0.6 cycles per sample is recorded as -0.4, and the
    sign flips with the fold. That is why the fold is applied here rather than
    left to a reader to remember.
    """
    if channel_hz <= 0 or tsamp_s <= 0:
        raise ValueError("channel width and sample time must be positive")
    # signed on purpose: folding 0.6 gives -0.4, not +0.4. Taking the magnitude
    # first throws away the sign the fold is supposed to flip, which is the exact
    # bug an earlier injection test carried.
    nu = float(drift_hz_per_s) * tsamp_s / channel_hz
    if nu > 0.5:
        return nu - 1.0
    if nu < -0.5:
        return nu + 1.0
    return nu


@dataclass(frozen=True)
class SignalRecord:
    source_file: str
    observation_date: str
    instrument: str
    frequency_hz: float
    drift_hz_per_s: float
    snr_global: float
    snr_at_cell: float
    local_scale: float
    resolution_hz_per_s: float
    channel_rate: float
    duration_s: float
    on_source_status: str
    transform_version: str
    disposition: str
    metadata: dict = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        missing = [name for name in ("source_file", "instrument",
                                     "transform_version")
                   if not str(getattr(self, name)).strip()]
        if missing:
            raise ValueError(f"empty provenance field(s): {', '.join(missing)}")
        date = str(self.observation_date)
        if len(date) != 8 or not date.isdigit():
            raise ValueError(f"observation_date must be YYYYMMDD, got {date!r}")
        try:
            dt.date(int(date[:4]), int(date[4:6]), int(date[6:]))
        except ValueError as exc:
            raise ValueError(f"observation_date is not a real date: {date}") from exc
        for name in ("frequency_hz", "resolution_hz_per_s", "duration_s",
                     "local_scale"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive, got {value}")
        for name in ("drift_hz_per_s", "snr_global", "snr_at_cell"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite, got {value}")
        if self.on_source_status not in ON_SOURCE:
            raise ValueError(f"on_source_status must be one of {ON_SOURCE}, "
                             f"got {self.on_source_status!r}")
        if self.disposition not in DISPOSITIONS:
            raise ValueError(f"disposition must be one of {DISPOSITIONS}, "
                             f"got {self.disposition!r}")
        if self.snr_global <= 0 and self.snr_at_cell <= 0:
            raise ValueError("a record with no positive significance anywhere is "
                             "not a candidate; keep it as a null-result summary")
        if not isinstance(self.metadata, dict):
            raise ValueError("metadata must be a mapping")

    def as_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "SignalRecord":
        known = {f for f in cls.__dataclass_fields__}  # type: ignore[attr-defined]
        unknown = set(data) - known
        if unknown:
            raise ValueError(f"unknown field(s): {sorted(unknown)}")
        return cls(**data)

    def significance_ratio(self) -> float:
        """How much the local noise reference moved this candidate."""
        if self.snr_global == 0:
            return float("nan")
        return self.snr_at_cell / self.snr_global
