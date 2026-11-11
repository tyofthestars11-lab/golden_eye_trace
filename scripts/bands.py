"""Open archives that could supply a second band, and the arithmetic that
decides whether one is far enough away to move the lever arm.

The pre-registered test in ``prereg/phi_2band.json`` measures
``p = d ln|drift| / d ln f`` and separates two statements: a flat rung rate
(``p = 1``, 4.7849 rungs per decade) from the pitch read as a slope
(``p = tan 17.0324 deg = 0.3063``, 1.4659 rungs per decade). Over a span of D
decades the two predictions differ by 3.3190 * D rungs per decade, so a slope
measured to better than 1.1065 * D rungs per decade can tell them apart at 3
sigma. That is the whole reason the registration asks for a decade: at D = 1 the
requirement is 1.1 rungs/decade, which is a loose measurement; at D = 0.2 it is
0.22 and no archive here delivers it.

Statuses below were checked from this workspace on 2026-11-11 and say what the
HTTP layer did, not what a search engine claims.
"""

from __future__ import annotations

import math

#: Frequency the first band sits at, in MHz (STARE2 OVRO, 1280-1530 MHz).
REFERENCE_MHZ = 1405.0

#: The pre-registration's evidence requirement.
MIN_DECADES = 1.0

#: ln(10)/ln(phi): rungs spanned by one decade, by definition.
RUNGS_PER_DECADE = math.log(10.0) / math.log((1.0 + math.sqrt(5.0)) / 2.0)

#: tan(17.0324 deg), the pitch read as a log-log slope.
PITCH_EXPONENT = math.tan(math.radians(17.0324))

PITCH_RUNGS_PER_DECADE = RUNGS_PER_DECADE * PITCH_EXPONENT

#: How far apart the two predictions are, in rungs per decade.
SEPARATION = RUNGS_PER_DECADE - PITCH_RUNGS_PER_DECADE

#: field-by-field: id, telescope, band, url, status seen here, format, note
SOURCES = [
    {
        "id": "stare2-ovro",
        "telescope": "STARE2 OVRO",
        "low_mhz": 1280.2,
        "high_mhz": 1530.1,
        "url": "local:data/stare2_ovro.fil",
        "status": "local file, 3 candidates at 3 sigma",
        "format": "sigproc",
        "auth": "none",
        "note": "the first band; drift resolution 3.2e5 Hz/s, no resolved drifts",
    },
    {
        "id": "lwa-archive",
        "telescope": "LWA1 (UNM)",
        "low_mhz": 10.0,
        "high_mhz": 88.0,
        "url": "https://lda10g.alliance.unm.edu/",
        "status": "200 directory listing",
        "format": "tbn/tbw time series",
        "auth": "none",
        "note": "only open source that clears a decade (10-88 MHz); needs "
                "channelising to a filterbank",
    },
    {
        "id": "lwa-drx",
        "telescope": "LWA1 (UNM), DRX beam 3",
        "low_mhz": 32.2,
        "high_mhz": 51.8,
        "url": "https://lda10g.alliance.unm.edu/PulsarArchive/B0329+54/56797/"
               "drx_56797_B0329+54_b3t1_DM0.00_topo.dat",
        "meta_url": "https://lda10g.alliance.unm.edu/PulsarArchive/B0329+54/56797/"
                    "drx_56797_B0329+54_b3t1_DM0.00_topo.inf",
        "status": "200, 16138240-byte .dat + 1.2 kB .inf",
        "format": "presto-dat",
        "auth": "none",
        "note": "the .inf gives 4096 channels x 4.785 kHz over 19.6 MHz and an "
                "843 s scan, i.e. 5.7 Hz/s drift resolution: 55000x finer than "
                "the first band, and 1.52 decades away",
    },
    {
        "id": "lofar-hba",
        "telescope": "LOFAR HBA",
        "low_mhz": 110.0,
        "high_mhz": 190.0,
        "url": "https://lta.lofar.eu/",
        "status": "no connection from here (000)",
        "format": "hdf5",
        "auth": "LTA account required",
        "note": "1.07 decades at 120 MHz; needs login, so not scriptable here",
    },
    {
        "id": "mwa",
        "telescope": "Murchison Widefield Array",
        "low_mhz": 72.0,
        "high_mhz": 300.0,
        "url": "https://data.mwatelescope.org/",
        "status": "no connection from here (000)",
        "format": "raw visibility",
        "auth": "Pawsey account required",
        "note": "1.29 decades at 72 MHz; not open",
    },
    {
        "id": "bl-opendata",
        "telescope": "Breakthrough Listen (GBT/Parkes)",
        "low_mhz": 1100.0,
        "high_mhz": 11200.0,
        "url": "https://seti.berkeley.edu/opendata",
        "status": "200 frontend; /api/* 404 from here",
        "format": "sigproc / hdf5",
        "auth": "none",
        "note": "X band at 11.2 GHz is 0.90 decades: the closest any open "
                "SETI archive comes, and still short of a decade",
    },
    {
        "id": "zenodo-vela",
        "telescope": "Parkes",
        "low_mhz": 1360.0,
        "high_mhz": 1680.0,
        "url": "https://zenodo.org/api/records/22848961/files/S53006_1.sf/content",
        "status": "200, 377280 bytes, header OBSFREQ 1520 OBSBW 320",
        "format": "fits (despite the .sf name)",
        "auth": "none",
        "note": "fetches cleanly and is 0.00 decades from the first band: the "
                "concrete case the gate exists to reject",
    },
    {
        "id": "fast",
        "telescope": "FAST",
        "low_mhz": 1000.0,
        "high_mhz": 1500.0,
        "url": "https://fast.bao.ac.cn/",
        "status": "200 landing page",
        "format": "psrfits",
        "auth": "account required",
        "note": "same band as the first: useless for the lever arm",
    },
]


def decades(a_mhz: float, b_mhz: float) -> float:
    """Separation of two frequencies in decades."""
    if a_mhz <= 0 or b_mhz <= 0:
        raise ValueError("frequencies must be positive")
    return abs(math.log10(b_mhz / a_mhz))


def qualifies(center_mhz: float, reference_mhz: float = REFERENCE_MHZ,
              min_decades: float = MIN_DECADES) -> bool:
    return decades(reference_mhz, center_mhz) >= min_decades


def slope_sigma_ceiling(span_decades: float, sigmas: float = 3.0) -> float:
    """Largest slope uncertainty, in rungs/decade, that still separates them."""
    return SEPARATION * span_decades / sigmas


def power_table(spans=(0.1, 0.3, 0.6, 0.9, 1.0, 1.5)) -> list[tuple[float, float]]:
    return [(span, slope_sigma_ceiling(span)) for span in spans]


def scan_time_for(drift_hz_s: float, channel_hz: float) -> float:
    """Contiguous scan a drift needs to move one full channel width.

    Drift resolution is channel width over integration time, so resolving
    small drifts needs either a long scan or narrow channels - and the two are
    not interchangeable when the channels are fixed by the instrument's FFT.
    """
    if drift_hz_s <= 0 or channel_hz <= 0:
        raise ValueError("drift and channel width must be positive")
    return channel_hz / drift_hz_s


def channel_width_for(drift_hz_s: float, scan_s: float) -> float:
    """Channel width that would resolve a given drift in a given scan."""
    if drift_hz_s <= 0 or scan_s <= 0:
        raise ValueError("drift and scan length must be positive")
    return drift_hz_s * scan_s


def centre(src: dict) -> float:
    if src.get("center_mhz"):
        return float(src["center_mhz"])
    return 0.5 * (float(src["low_mhz"]) + float(src["high_mhz"]))


def report() -> str:
    lines = [
        f"reference band {REFERENCE_MHZ:.0f} MHz; a second band needs "
        f"{MIN_DECADES:g} decade = a factor {10 ** MIN_DECADES:g}",
        "",
        f"{'id':<15}{'centre MHz':>11}{'decades':>9}{'ok':>5}  status",
    ]
    for src in SOURCES:
        c = centre(src)
        d = decades(REFERENCE_MHZ, c)
        lines.append(
            f"{src['id']:<15}{c:>11.1f}{d:>9.2f}{'yes' if d >= MIN_DECADES else 'no':>5}"
            f"  {src['status'][:52]}"
        )
    lines += [
        "",
        f"the two statements differ by {SEPARATION:.4f} rungs per decade, so a",
        "slope known to this precision separates them at 3 sigma:",
        f"{'span (decades)':>14}{'sigma ceiling':>16}",
    ]
    for span, ceiling in power_table():
        lines.append(f"{span:>14.2f}{ceiling:>16.4f}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
