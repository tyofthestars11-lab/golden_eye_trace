"""The diamond NV-centre lane: a local field witness that can only detract.

A nitrogen-vacancy ensemble in diamond measures magnetic field, and a radio wave
carries one: for a plane wave ``E = c B``, so 1 pT is 0.3 mV/m. That is why an
NV sensor is worth integrating as an **independent local witness** at the
telescope site: if a radio transient is accompanied by a magnetic transient of
the same frequency inside the same second, the simplest explanation is a local
emitter, and the candidate detracts.

The asymmetry is the whole design, and it is the same asymmetry the phi layer
carries:

* A **positive** local-field record detracts. It is evidence *against* an
  astrophysical origin, and it is allowed to lower a disposition.
* A **null** local-field record proves nothing at all. The NV sensitivity is
  worse than a dipole antenna's, its sensing volume is millimetres, and its
  coupling to a wave arriving from a given direction depends on the NV axis. No
  transient in the sensor is compatible with a real local emitter that is
  polarised away from the axis, or that couples to the antenna but not the
  diamond.
* Therefore the lane may **never** support a disposition, and
  ``guard_support`` raises rather than returning a flag.

Reach, stated rather than assumed: the NV ground-state splitting is 2.87 GHz, and
a bias field tunes it by 2.8 MHz/G, so sensing at 50 MHz needs about 100 mT. The
dynamical-decoupling sequences that reach fT-to-pT sensitivity work at kHz to a
few MHz (their upper limit is the Rabi frequency); above that the practical route
is RF-dressed CW-ODMR, whose demonstrated bandwidth is a few MHz. The LWA band
13.4-86.9 MHz is therefore reachable, but only in slices and only with a
substantial bias field. None of that is hidden in a constant here: the
``zeeman_bias_t`` and ``band_covered`` functions compute it per request.
"""

from __future__ import annotations

import datetime as dt
import math
from dataclasses import dataclass

__all__ = [
    "NVAsEvidence",
    "NVSensorRecord",
    "SENSING_MODES",
    "VERDICTS",
    "band_covered",
    "field_from_magnetic",
    "guard_support",
    "magnetic_from_field",
    "significance",
    "verdict",
    "zeeman_bias_t",
]

#: NV ground-state zero-field splitting, and the Zeeman tuning rate.
ZFS_HZ = 2.87e9
ZEEMAN_HZ_PER_T = 2.8e10  # 2.8 MHz/G
C_M_PER_S = 299_792_458.0

#: Closed vocabulary of acquisition schemes. Each has a different usable band:
#: the dynamical-decoupling family (echo, CPMG) is narrowband and reaches pT
#: sensitivity but is limited above by the Rabi frequency; CW-ODMR is broadband
#: and far less sensitive.
SENSING_MODES = (
    "cw-odmr",           # broadband, uT-nT/sqrt(Hz), continuous
    "ramsey",            # 1e3-1e5 Hz, uT-nT/sqrt(Hz)
    "hahn-echo",         # 1e3-1e6 Hz, nT-pT/sqrt(Hz)
    "cpmg",              # 1e4-1e7 Hz, nT-pT/sqrt(Hz)
    "rf-dressed",        # tuned to the dressed resonance, few MHz bandwidth
    "quantum-heterodyne",
)

VERDICTS = (
    "LOCAL TRANSIENT (DETRACTS)",
    "NO LOCAL TRANSIENT (NO INFORMATION)",
    "NOT TESTABLE",
)


class NVAsEvidence(Exception):
    """Raised when an NV record is asked to support a disposition."""


def zeeman_bias_t(frequency_hz: float) -> float:
    """Bias field that brings the NV resonance down to ``frequency_hz``.

    ``f = ZFS - gamma B`` for the transition that moves downwards with field, so
    ``B = (ZFS - f) / gamma``. A frequency above the zero-field splitting is not
    reachable with a positive bias and raises rather than returning a negative
    field, which would look like a sensor that exists.
    """
    if frequency_hz <= 0:
        raise ValueError("a frequency must be positive")
    if frequency_hz >= ZFS_HZ:
        raise ValueError(
            f"{frequency_hz:.3e} Hz is at or above the 2.87 GHz zero-field "
            "splitting: a positive bias field cannot reach it")
    return (ZFS_HZ - frequency_hz) / ZEEMAN_HZ_PER_T


def field_from_magnetic(b_t: float) -> float:
    """Electric field of a plane wave with this magnetic amplitude, in V/m."""
    return C_M_PER_S * b_t


def magnetic_from_field(e_v_m: float) -> float:
    """Magnetic amplitude of a plane wave with this electric amplitude, in T."""
    return e_v_m / C_M_PER_S


def band_covered(band_low_hz: float, band_high_hz: float,
                 sensor_low_hz: float, sensor_high_hz: float) -> float:
    """Fraction of ``[band_low, band_high]`` the sensor's band covers, 0 to 1."""
    for value in (band_low_hz, band_high_hz, sensor_low_hz, sensor_high_hz):
        if value < 0:
            raise ValueError("frequencies cannot be negative")
    if band_high_hz <= band_low_hz:
        raise ValueError("the band is empty or reversed")
    overlap = (min(band_high_hz, sensor_high_hz)
               - max(band_low_hz, sensor_low_hz))
    return max(0.0, overlap) / (band_high_hz - band_low_hz)


def significance(b_amplitude_t: float, sensitivity_t_per_sqrt_hz: float,
                 bandwidth_hz: float) -> float:
    """Signal-to-noise of a coherent amplitude in a stated bandwidth.

    ``SNR = B / (eta sqrt(delta f))``: the sensitivity is per root hertz, so a
    wider bandwidth costs signal-to-noise. Quoting an NV sensitivity without its
    bandwidth is the same error as quoting a sigma without its null.
    """
    if sensitivity_t_per_sqrt_hz <= 0 or bandwidth_hz <= 0:
        raise ValueError("sensitivity and bandwidth must be positive")
    if b_amplitude_t < 0:
        raise ValueError("an amplitude cannot be negative")
    return b_amplitude_t / (sensitivity_t_per_sqrt_hz * math.sqrt(bandwidth_hz))


@dataclass(frozen=True)
class NVSensorRecord:
    """One NV acquisition: what it could see, when, and from where.

    ``sensitivity_t_per_sqrt_hz`` and ``bandwidth_hz`` are both required, because
    either alone is meaningless. ``bias_field_t`` is required for the same
    reason: an NV ensemble at zero field cannot sense the LWA band at all, and a
    record that omits it cannot be checked.
    """

    sensor_id: str
    instrument: str
    sensing_mode: str
    band_low_hz: float
    band_high_hz: float
    sensitivity_t_per_sqrt_hz: float
    bandwidth_hz: float
    bias_field_t: float
    start_time_utc: dt.datetime
    time_uncertainty_s: float
    axis_tilt_deg: float
    data_release: str
    public_source_url: str

    def __post_init__(self) -> None:
        for name in ("sensor_id", "instrument", "data_release",
                     "public_source_url"):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"an NV record needs {name}")
        if self.sensing_mode not in SENSING_MODES:
            raise ValueError(f"sensing_mode {self.sensing_mode!r} not one of "
                             f"{SENSING_MODES}")
        if self.band_high_hz <= self.band_low_hz:
            raise ValueError("the sensor band is empty or reversed")
        if self.sensitivity_t_per_sqrt_hz <= 0 or self.bandwidth_hz <= 0:
            raise ValueError("sensitivity and bandwidth must be positive")
        if self.bias_field_t < 0:
            raise ValueError("a bias field cannot be negative")
        if self.time_uncertainty_s < 0:
            raise ValueError("a time uncertainty cannot be negative")
        if not (0.0 <= self.axis_tilt_deg <= 90.0):
            raise ValueError("axis_tilt_deg is an angle in [0, 90]")
        stamp = self.start_time_utc
        if not isinstance(stamp, dt.datetime):
            raise ValueError("start_time_utc must be a datetime")
        if stamp.tzinfo is None or stamp.utcoffset() is None:
            raise ValueError("start_time_utc must be timezone-aware: a local "
                             "wall clock cannot be compared with UTC")

    @property
    def field_sensitivity_v_m(self) -> float:
        """Sensitivity expressed as an electric-field amplitude per root hertz."""
        return field_from_magnetic(self.sensitivity_t_per_sqrt_hz)

    def covers(self, frequency_hz: float) -> bool:
        return self.band_low_hz <= frequency_hz <= self.band_high_hz

    def to_dict(self) -> dict:
        return {
            "sensor_id": self.sensor_id,
            "instrument": self.instrument,
            "sensing_mode": self.sensing_mode,
            "band_low_hz": self.band_low_hz,
            "band_high_hz": self.band_high_hz,
            "sensitivity_t_per_sqrt_hz": self.sensitivity_t_per_sqrt_hz,
            "sensitivity_v_m_per_sqrt_hz": self.field_sensitivity_v_m,
            "bandwidth_hz": self.bandwidth_hz,
            "bias_field_t": self.bias_field_t,
            "start_time_utc": self.start_time_utc.isoformat(),
            "time_uncertainty_s": self.time_uncertainty_s,
            "axis_tilt_deg": self.axis_tilt_deg,
            "data_release": self.data_release,
            "public_source_url": self.public_source_url,
        }


def verdict(sensor: NVSensorRecord | None, *, frequency_hz: float,
            radio_time_utc: dt.datetime | None,
            radio_duration_s: float = 0.0,
            b_amplitude_t: float | None = None,
            sigma: float = 3.0) -> dict:
    """Grade the local-field lane for one radio record.

    Three outcomes, and only the first is evidence:

    ``LOCAL TRANSIENT (DETRACTS)`` - the sensor covers the frequency, overlaps the
    window, and saw an amplitude above threshold. It may lower a disposition.
    ``NO LOCAL TRANSIENT (NO INFORMATION)`` - the sensor looked and saw nothing.
    This is compatible with a real local emitter, so it says nothing.
    ``NOT TESTABLE`` - there is no sensor, or none covering that frequency or
    time. Also says nothing, and says why.
    """
    if sensor is None:
        return {"verdict": "NOT TESTABLE",
                "reasons": ["no NV sensor record was supplied"]}
    if radio_time_utc is None:
        return {"verdict": "NOT TESTABLE",
                "reasons": ["the radio record has no UTC timestamp, so no local "
                            "field can be compared with it"]}
    if not sensor.covers(frequency_hz):
        return {"verdict": "NOT TESTABLE",
                "reasons": [f"the sensor covers "
                            f"{sensor.band_low_hz / 1e6:.1f}-"
                            f"{sensor.band_high_hz / 1e6:.1f} MHz and the record "
                            f"is at {frequency_hz / 1e6:.4f} MHz"]}
    dt_s = abs((sensor.start_time_utc - radio_time_utc).total_seconds())
    window_s = sigma * (sensor.time_uncertainty_s + radio_duration_s)
    if dt_s > max(window_s, sensor.time_uncertainty_s):
        return {"verdict": "NOT TESTABLE",
                "reasons": [f"the sensor acquisition is {dt_s:.0f} s from the "
                            "record, outside any window the timing supports"]}
    if b_amplitude_t is None:
        return {"verdict": "NO LOCAL TRANSIENT (NO INFORMATION)",
                "reasons": ["no amplitude was reported for the window; a null "
                            "here is compatible with a real local emitter, so it "
                            "is not evidence"],
                "window_s": window_s}
    snr = significance(b_amplitude_t, sensor.sensitivity_t_per_sqrt_hz,
                       sensor.bandwidth_hz)
    if snr < sigma:
        return {"verdict": "NO LOCAL TRANSIENT (NO INFORMATION)",
                "reasons": [f"{snr:.1f} sigma in {sensor.bandwidth_hz:.0f} Hz "
                            "is below threshold; the sensor is mm-scale and "
                            "polarisation-dependent, so a null is not evidence"],
                "snr": snr, "window_s": window_s}
    return {"verdict": "LOCAL TRANSIENT (DETRACTS)",
            "reasons": [f"{snr:.1f} sigma at {frequency_hz / 1e6:.4f} MHz in the "
                        "same window: the simplest explanation is a local "
                        "emitter"],
            "snr": snr, "window_s": window_s,
            "b_t": b_amplitude_t, "e_v_m": field_from_magnetic(b_amplitude_t)}


def guard_support(record, *, claimed: str | None = None) -> dict:
    """Refuse to let a local-field record support anything.

    Two refusals, both aimed at the same mistake: treating the absence of a
    local transient as evidence of an astrophysical origin. That inference needs
    a sensor that couples to every direction and polarisation at full
    sensitivity, which a millimetre of diamond is not.
    """
    detail = f"claimed {claimed!r}" if claimed else f"record {record}"
    raise NVAsEvidence(
        "a local-field record may only detract: absence of a transient in an NV "
        "sensor is compatible with a real local emitter that is polarised away "
        "from the NV axis or couples to the antenna and not to the diamond "
        f"({detail})")
