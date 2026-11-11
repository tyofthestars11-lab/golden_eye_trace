"""A narrowband drift-rate search for technosignatures in SIGPROC filterbank data.

Quick start:

    python -m seti_search simulate data/demo.fil
    python -m seti_search search data/demo.fil
    python -m seti_search info data/demo.fil
"""

from .cadence import (
    SITES,
    Observer,
    SkyPosition,
    barycentric_drift,
    doppler_factor,
    predict_observation,
    topocentric_drift,
)
from .filterbank import Filterbank, InjectedSignal, RFILine, simulate, write_filterbank
from .flux import (
    COMPLETENESS_MARGIN,
    FluxEstimate,
    boxcar_scale,
    combine_fluxes,
    flux_density_jy,
    limit_jy,
    noise_per_sample_jy,
    run_limit,
    sefd_from_gain,
    sefd_from_tsys,
)
from .rungs import RungAddress, address, ladder_step, repeated_scales, rung, rung_uniformity
from .run import Event, Scan, link_run
from .search import Candidate, SearchConfig, search_band, search_file
from .verify import CadenceResult, Verdict, verify_cadence

__version__ = "0.4.0"

__all__ = [
    "Candidate",
    "CadenceResult",
    "Event",
    "RungAddress",
    "Scan",
    "Filterbank",
    "FluxEstimate",
    "InjectedSignal",
    "Observer",
    "RFILine",
    "SITES",
    "SearchConfig",
    "SkyPosition",
    "Verdict",
    "address",
    "barycentric_drift",
    "doppler_factor",
    "ladder_step",
    "link_run",
    "predict_observation",
    "repeated_scales",
    "rung",
    "rung_uniformity",
    "search_band",
    "search_file",
    "simulate",
    "topocentric_drift",
    "verify_cadence",
    "write_filterbank",
]
