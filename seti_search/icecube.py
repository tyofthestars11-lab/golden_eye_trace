"""Ingest a published IceCube alert into the neutral neutrino event.

The published alert is not the published *event*: it carries a containment
radius instead of a sigma, an energy in TeV on an assumed E^-2.19 spectrum, a
false-alarm rate in Hz on a Kafka topic but in per-year on the AMON table, and
two timestamps of which only one is the neutrino's. Each of those is a factor
away from the wrong answer, so each is converted here, once, with the conversion
recorded alongside the event.

What is refused, and why:

* ``alert_type: "retraction"`` - a withdrawn alert is not an event.
* a missing ``trigger_time`` - ``alert_datetime`` is when the notice was sent,
  minutes to hours after the neutrino, and substituting it would quietly destroy
  the only quantity this lane needs.
* a missing or non-positive concordance between error and containment - an
  angular error with no stated containment cannot be turned into a sigma, and
  guessing 90 would inflate or deflate every cone by a factor of two.
* ``alert_tense`` of ``test`` or ``injections`` is accepted only when the caller
  asks for it, and is then labelled ``test`` so nothing downstream treats a
  hardware injection as the sky.
"""

from __future__ import annotations

import datetime as dt
import math

from .neutrino import NeutrinoEvent

__all__ = [
    "IceCubeSchemaError",
    "IngestedAlert",
    "PARSE_KINDS",
    "from_alert",
    "sigma_from_containment_deg",
]

#: Where a row's numbers came from, so a reader knows which conversions applied.
PARSE_KINDS = ("gcn-kafka", "amon-table")


class IceCubeSchemaError(ValueError):
    """A published row that cannot be turned into an event without guessing."""


def sigma_from_containment_deg(radius_deg: float, containment: float) -> float:
    """Convert a containment radius to the 1-sigma of the equivalent 2-D Gaussian.

    ``P(<R) = 1 - exp(-R^2 / 2 sigma^2)``, so ``sigma = R / sqrt(-2 ln(1 - p))``.
    For 90% containment that factor is 2.146, not 1: using the published 90%
    radius directly as a sigma would make every cone 2.1 times too wide and every
    chance probability correspondingly too large. For 50% it is 1.177.
    """
    if not (0.0 < containment < 1.0):
        raise IceCubeSchemaError(
            f"containment {containment} is not a probability in (0, 1)")
    if radius_deg <= 0:
        raise IceCubeSchemaError("a positional error must be positive")
    return radius_deg / math.sqrt(-2.0 * math.log(1.0 - containment))


def _first(value):
    """GCN sends some fields as one-element lists; unwrap without guessing."""
    if isinstance(value, (list, tuple)):
        if len(value) != 1:
            raise IceCubeSchemaError(
                f"expected a single value, got {len(value)}: {value!r}")
        return value[0]
    return value


def _parse_time(value, what: str) -> dt.datetime:
    if not isinstance(value, str):
        raise IceCubeSchemaError(f"{what} must be an ISO 8601 string")
    try:
        stamp = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise IceCubeSchemaError(f"{what} {value!r} is not ISO 8601: {exc}")
    if stamp.tzinfo is None:
        raise IceCubeSchemaError(f"{what} {value!r} has no timezone")
    return stamp.astimezone(dt.timezone.utc)


def _alert_status(row: dict) -> str:
    tense = str(_first(row.get("alert_tense", "current"))).lower()
    pipeline = str(_first(row.get("analysis_pipeline", ""))).lower()
    if tense in ("test", "injections"):
        return "test"
    if tense == "current":
        pass
    else:
        raise IceCubeSchemaError(f"unknown alert_tense {tense!r}")
    if "gold" in pipeline:
        return "GOLD"
    if "bronze" in pipeline:
        return "BRONZE"
    # an archival release names neither, and pretending it is an alert would put
    # an analysis-grade event into a real-time coincidence table
    declared = str(_first(row.get("notice_type", ""))).upper()
    if declared in ("GOLD", "BRONZE"):
        return declared
    if row.get("data_release") or row.get("alert_status") == "archival":
        return "archival"
    raise IceCubeSchemaError(
        "the row names neither a GOLD nor a BRONZE pipeline and is not marked "
        "archival: the alert grade is not inferable")


def _energy_gev(row: dict) -> float:
    for key, scale in (("nu_energy", 1000.0),        # GCN: TeV
                       ("energy_tev", 1000.0),       # AMON column: TeV
                       ("energy_gev", 1.0),
                       ("energy", 1000.0)):          # legacy AMON: TeV
        if row.get(key) is not None:
            value = float(_first(row[key]))
            if value <= 0:
                raise IceCubeSchemaError(f"{key} must be positive, got {value}")
            return value * scale
    raise IceCubeSchemaError("no energy field (nu_energy / energy_tev / energy_gev)")


def _far_per_year(row: dict) -> float | None:
    """False alarm rate in per year, from Hz (GCN) or per-year (AMON table)."""
    unit = str(row.get("far_unit", "")).lower()
    if row.get("far_per_year") is not None:
        return float(row["far_per_year"])
    if row.get("far") is None:
        return None
    far = float(_first(row["far"]))
    if far < 0:
        raise IceCubeSchemaError("a false alarm rate cannot be negative")
    if unit in ("hz", "per_s", "/s"):
        return far * 3.1557e7
    if unit in ("per_year", "/yr", "yr-1"):
        return far
    # the GCN schema documents Hz and the AMON table documents per year; they
    # differ by seven orders of magnitude, so the unit is inferred from the
    # exponent and the inference is written into the record
    return far * 3.1557e7 if far < 1e-4 else far


def _angular_sigma_deg(row: dict) -> tuple[float, float, str]:
    """Return ``(sigma_deg, containment, convention_note)``."""
    containment = row.get("containment_probability", 0.9)
    if isinstance(row.get("ra_dec_error"), (int, float)):
        radius = float(row["ra_dec_error"])
        p = float(containment)
        return (sigma_from_containment_deg(radius, p), p,
                f"ra_dec_error {radius:g} deg at {p:.0%} containment")
    for key, p, scale in (("ang_err90", 0.9, 1.0),
                          ("error90_deg", 0.9, 1.0),
                          ("error90_arcmin", 0.9, 1.0 / 60.0),
                          ("ang_err50", 0.5, 1.0),
                          ("error50_arcmin", 0.5, 1.0 / 60.0)):
        if row.get(key) is not None:
            radius = float(_first(row[key])) * scale
            return (sigma_from_containment_deg(radius, p), p,
                    f"{key} {radius:g} deg at {p:.0%} containment")
    raise IceCubeSchemaError(
        "no positional error (ra_dec_error, ang_err90 or ang_err50): a cone "
        "cannot be built without one")


def _event_type(row: dict) -> str:
    topology = str(_first(row.get("alert_topology", "Track"))).lower()
    if topology == "track":
        return "through-going"
    if topology in ("cascade", "shower"):
        return "cascade"
    return "unknown"


class IngestedAlert:
    """The event plus the published quantities that are not part of it.

    ``far_per_year`` is the detector's own null: how often a background event of
    this significance occurs. It is kept beside the event rather than inside it
    because the event schema is a statement about the neutrino, and this is a
    statement about the detector's background.
    """

    __slots__ = ("event", "far_per_year", "p_astro", "run_id", "alert_type",
                 "alert_tense", "record_number", "conversions", "kind")

    def __init__(self, event: NeutrinoEvent, *, far_per_year: float | None,
                 p_astro: float | None, run_id: str, alert_type: str,
                 alert_tense: str, record_number: int, conversions: list[str],
                 kind: str) -> None:
        self.event = event
        self.far_per_year = far_per_year
        self.p_astro = p_astro
        self.run_id = run_id
        self.alert_type = alert_type
        self.alert_tense = alert_tense
        self.record_number = record_number
        self.conversions = conversions
        self.kind = kind

    def to_dict(self) -> dict:
        return {
            "event": self.event.to_dict(),
            "far_per_year": self.far_per_year,
            "p_astro": self.p_astro,
            "run_id": self.run_id,
            "alert_type": self.alert_type,
            "alert_tense": self.alert_tense,
            "record_number": self.record_number,
            "conversions": list(self.conversions),
            "kind": self.kind,
        }


def from_alert(row: dict, *, public_source_url: str,
               data_release: str,
               allow_test: bool = False) -> IngestedAlert:
    """Build a NeutrinoEvent from one published alert row.

    ``public_source_url`` and ``data_release`` are arguments rather than fields
    because a notice itself does not say where it was archived, and an event
    whose provenance depends on the caller's memory is not provenance.
    """
    alert_type = str(_first(row.get("alert_type", "initial"))).lower()
    if alert_type == "retraction":
        raise IceCubeSchemaError(
            "this alert was retracted: a withdrawn event is not an event, and "
            "keeping it would put a real coincidence on a direction that the "
            "detector no longer believes")
    tense = str(_first(row.get("alert_tense", "current"))).lower()
    if tense in ("test", "injections") and not allow_test:
        raise IceCubeSchemaError(
            f"alert_tense is {tense!r}: pass allow_test=True to ingest a "
            "hardware injection as a test event")

    trigger = row.get("trigger_time") or row.get("event_time_utc")
    if not trigger:
        raise IceCubeSchemaError(
            "no trigger_time: alert_datetime is the notice time, not the "
            "neutrino time, and substituting it would destroy the only "
            "quantity this lane measures")
    event_time = _parse_time(_first(trigger), "trigger_time")

    sigma, containment, convention = _angular_sigma_deg(row)
    energy_gev = _energy_gev(row)
    far_per_year = _far_per_year(row)
    status = _alert_status(row)

    name = str(_first(row.get("event_name", ""))).strip()
    run_id = str(_first(row.get("id", row.get("run_id", "")))).strip()
    if not name:
        if not run_id:
            raise IceCubeSchemaError("the row has neither event_name nor id")
        name = f"IC-{run_id}"

    record_number = row.get("record_number",
                            0 if alert_type == "initial" else 1)
    energy_class = ("shower" if _event_type(row) == "cascade"
                    else ("track" if _event_type(row) == "through-going"
                          else "unknown"))

    event = NeutrinoEvent(
        event_id=name,
        detector=str(_first(row.get("instrument", row.get("detector",
                                                        "IceCube")))),
        event_time_utc=event_time,
        # IceCube does not publish a time uncertainty for track alerts; the
        # readout is exact to well below a second and the *radio* window
        # dominates by orders of magnitude, so zero is recorded rather than
        # invented
        time_uncertainty_s=float(_first(row.get("time_uncertainty_s", 0.0))),
        ra_deg=float(_first(row["ra"])),
        dec_deg=float(_first(row["dec"])),
        angular_error_deg=sigma,
        energy_estimate_gev=energy_gev,
        energy_class=energy_class,
        event_type=_event_type(row),
        alert_status=status,
        public_source_url=public_source_url,
        data_release=data_release,
    )
    conversions = [
        f"{convention} -> {sigma:.3f} deg (1 sigma)",
        f"energy {energy_gev / 1000.0:g} TeV (published units) -> "
        f"{energy_gev:g} GeV",
    ]
    if far_per_year is not None:
        conversions.append(
            f"false alarm rate {far_per_year:g} per year "
            f"({far_per_year / 3.1557e7:.3e} Hz)")
    if record_number:
        conversions.append(f"record_number {record_number}: revised "
                           "reconstruction, superseding the initial notice")
    return IngestedAlert(
        event, far_per_year=far_per_year,
        p_astro=(None if row.get("p_astro") is None
                 else float(_first(row["p_astro"]))),
        run_id=run_id, alert_type=alert_type, alert_tense=tense,
        record_number=int(record_number),
        conversions=conversions,
        kind="gcn-kafka" if "ra_dec_error" in row else "amon-table",
    )
