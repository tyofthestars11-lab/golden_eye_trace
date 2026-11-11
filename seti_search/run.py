"""Link candidates across a whole observing run.

One cadence can confirm a candidate. A run - several ON scans of the same target
with OFF scans interleaved - can do much better, because a real source has to
keep the *same barycentric drift* every time it is observed, however the
topocentric drift wanders as Earth turns. Two ON scans six hours apart see the
same source at drifts that differ by a tenth of a hertz per second or more, so
a set of hits that agrees only after the barycentric correction is a much
stronger statement than a set that agrees as measured.

Each hit is therefore referred to two run-fixed quantities:

* its barycentric drift, which belongs to the source;
* its barycentric frequency at a common reference epoch, which also belongs to
  the source - the raw frequency does not, because a drifting source's own
  frequency changes across the run.

Events are built by clustering ON hits on those two numbers, requiring hits from
several distinct ON scans. An event is then tested against every OFF scan, at
the position and drift the prediction gives for that scan's epoch *and* at the
topocentric position of the event's own hits (a carrier does not move). Any OFF
match means the signal survived pointing away, so the event does not count.

Without a site, a sky position and epochs, the run still links, but on
topocentric quantities: only scans close enough together in time will link, and
nothing corrects for the wandering drift. The result says which mode it used.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .cadence import (
    Observer,
    SkyPosition,
    barycentric_drift,
    doppler_factor,
    topocentric_drift,
    topocentric_frequency,
)
from .flux import FluxEstimate, combine_fluxes, flux_density_jy
from .rungs import RungAddress, address
from .search import Candidate

__all__ = [
    "Event",
    "HitRef",
    "OffMatch",
    "RunResult",
    "Scan",
    "link_run",
    "scan_sensitivity",
]

SURVIVED = "survived"
SEEN_IN_OFF = "seen in off"


@dataclass
class Scan:
    """One searched file, with what it takes to place it in the run."""

    path: str
    kind: str  # "on" or "off"
    mjd: float | None
    hits: list[Candidate]
    meta: dict
    # system equivalent flux density in Jy; without it the run reports noise
    # units and no physical flux
    sefd_jy: float | None = None


@dataclass
class HitRef:
    """A candidate plus the run-fixed quantities derived from it."""

    scan: str
    kind: str
    mjd: float | None
    candidate: Candidate
    freq_hz: float  # barycentric frequency referred to the run epoch, else topocentric
    drift_hz_s: float  # barycentric drift, else topocentric
    topocentric_freq_mhz: float
    sensitivity: float = 1.0
    sefd_jy: float | None = None
    df_hz: float = 0.0
    tsamp_s: float = 0.0
    boxcar: int = 1

    @property
    def amplitude(self) -> float:
        """SNR per unit sensitivity: the source strength this scan implies."""
        return self.candidate.snr / self.sensitivity if self.sensitivity > 0 else 0.0

    def as_dict(self) -> dict:
        return {
            "scan": self.scan,
            "kind": self.kind,
            "mjd": self.mjd,
            "freq_mhz": round(self.candidate.freq_mhz, 9),
            "drift_hz_s": round(self.candidate.drift_hz_s, 6),
            "snr": round(self.candidate.snr, 2),
        }


@dataclass
class OffMatch:
    """An OFF-scan hit that sits where the event predicts.

    ``window_hz`` is the full width of the frequency window the hit was found
    in, which is what makes the match's chance computable afterwards: the OFF
    test accepts any drift inside that window, so the cells it searched are the
    window's channels times the scan's drift trials.
    """

    scan: str
    hit: Candidate
    window_hz: float


@dataclass
class Event:
    """A candidate that recurred across several ON scans."""

    freq_hz: float  # barycentric at the reference epoch (or topocentric)
    drift_hz_s: float
    score: float
    snr: float
    # the same event in the telescope's own frame: something that does not move
    # is linked here rather than in the barycentric frame, and then rejected
    topo_freq_hz: float = 0.0
    topo_drift_hz_s: float = 0.0
    members: list[HitRef] = field(default_factory=list)
    off_matches: list[OffMatch] = field(default_factory=list)
    freq_scatter_mhz: float = 0.0
    drift_scatter_hz_s: float = 0.0
    status: str = SURVIVED
    # sensitivity-aware combination, see _recompute
    amplitude: float = 0.0
    amplitude_scatter: float = 0.0
    total_sensitivity: float = 0.0
    expected_snr: dict[str, float] = field(default_factory=dict)
    missed_scans: list[str] = field(default_factory=list)
    flux: FluxEstimate | None = None
    rung_freq: RungAddress | None = None
    rung_drift: RungAddress | None = None
    rung_snr: RungAddress | None = None

    @property
    def n_on(self) -> int:
        return len(self.members)

    @property
    def n_scans(self) -> int:
        return len({m.scan for m in self.members})

    @property
    def freq_mhz(self) -> float:
        return self.freq_hz / 1e6

    @property
    def survived(self) -> bool:
        return self.status == SURVIVED


@dataclass
class RunResult:
    events: list[Event] = field(default_factory=list)
    predicted: bool = False
    notes: list[str] = field(default_factory=list)

    @property
    def survivors(self) -> list[Event]:
        return [e for e in self.events if e.survived]


def link_run(
    scans: list[Scan],
    *,
    observer: Observer | None = None,
    sky: SkyPosition | None = None,
    drift_tolerance_hz_s: float = 0.05,
    freq_tolerance_hz: float | None = None,
    min_scans: int = 2,
    min_ons: int | None = None,
    line_width_channels: float = 1.0,
) -> RunResult:
    """Cluster ON hits across scans and test the clusters against the OFF scans.

    Args:
        scans: searched files, in any order, each with its kind and epoch.
        observer, sky: needed together with scan epochs for the barycentric
            conversion. Without them hits are linked as measured.
        drift_tolerance_hz_s: how far two hits' barycentric drifts may differ
            and still be the same source.
        freq_tolerance_hz: how far their referred frequencies may differ.
            Default is four channels of the widest scan.
        min_scans: how many distinct ON scans an event must appear in.
        min_ons: alias for ``min_scans`` kept for readability at the call site.
    """
    if min_ons is not None:
        min_scans = min_ons
    on_scans = [s for s in scans if s.kind == "on"]
    off_scans = [s for s in scans if s.kind == "off"]
    if not on_scans:
        raise ValueError("a run needs at least one ON scan")

    result = RunResult()
    can_convert = (
        observer is not None
        and sky is not None
        and all(s.mjd is not None for s in on_scans)
    )
    result.predicted = bool(can_convert)
    if not can_convert:
        missing = []
        if observer is None:
            missing.append("site latitude/longitude")
        if sky is None:
            missing.append("RA/Dec")
        if any(s.mjd is None for s in on_scans):
            missing.append("an ON epoch")
        result.notes.append(
            "barycentric conversion unavailable (missing " + ", ".join(missing) + "); "
            "linking on topocentric drift, which only works while the drift has not"
            " moved much between scans"
        )

    epochs = [s.mjd for s in on_scans if s.mjd is not None]
    epoch_ref = min(epochs) if epochs else None
    widest_df_hz = max(abs(float(s.meta["df_mhz"])) for s in scans) * 1e6
    # Referring every hit to one epoch means a drift error of e Hz/s moves it by
    # e * elapsed, so the window has to grow with the run. Two sources closer
    # together than that, with drifts within tolerance, would merge - which is
    # why short runs link more sharply than long ones.
    max_elapsed_s = 86400.0 * max((abs(e - epoch_ref) for e in epochs), default=0.0)
    freq_window_hz = (
        float(freq_tolerance_hz)
        if freq_tolerance_hz
        else max(4.0 * widest_df_hz, 0.02 * max_elapsed_s)
    )

    events: list[Event] = []
    for scan in on_scans:
        for hit in sorted(scan.hits, key=lambda h: -h.snr):
            freq_hz, drift = _hit_key(hit, scan, observer, sky, epoch_ref, can_convert)
            target = None
            for event in events:
                if scan.path in {m.scan for m in event.members}:
                    continue
                if _links(
                    event,
                    freq_hz,
                    drift,
                    hit.freq_mhz * 1e6,
                    hit.drift_hz_s,
                    drift_tolerance_hz_s,
                    freq_window_hz,
                ):
                    target = event
                    break
            member = HitRef(
                scan=scan.path,
                kind=scan.kind,
                mjd=scan.mjd,
                candidate=hit,
                freq_hz=freq_hz,
                drift_hz_s=drift,
                topocentric_freq_mhz=hit.freq_mhz,
                sensitivity=scan_sensitivity(scan.meta),
                sefd_jy=scan.sefd_jy,
                df_hz=abs(float(scan.meta["df_mhz"])) * 1e6,
                tsamp_s=float(scan.meta["tsamp"]),
                boxcar=int(scan.meta["boxcar"]),
            )
            if target is None:
                events.append(
                    Event(
                        freq_hz=freq_hz,
                        drift_hz_s=drift,
                        score=0.0,
                        snr=hit.snr,
                        topo_freq_hz=hit.freq_mhz * 1e6,
                        topo_drift_hz_s=hit.drift_hz_s,
                        members=[member],
                    )
                )
            else:
                target.members.append(member)
                _recompute(target, line_width_channels)

    survivors = [event for event in events if event.n_scans >= min_scans]
    for event in survivors:
        _check_completeness(event, on_scans)
        event.rung_freq = address(event.freq_hz, "Hz")
        # rung is undefined at zero, and a zero drift is the one value where the
        # address would say nothing at all
        event.rung_drift = (
            address(abs(event.drift_hz_s), "Hz/s") if event.drift_hz_s else None
        )
        event.rung_snr = address(event.snr, "sigma")
        _test_off_scans(
            event,
            off_scans,
            observer,
            sky,
            can_convert,
            freq_tolerance_hz=freq_window_hz,
        )

    survivors.sort(key=lambda e: (-e.survived, -e.score))
    result.events = survivors
    if can_convert:
        result.notes.append(
            "hits were referred to a common epoch before linking, so drift a"
            " telescope measures moving across the run does not break a link"
        )
    return result


def _hit_key(hit, scan, observer, sky, epoch_ref, can_convert) -> tuple[float, float]:
    if not can_convert or scan.mjd is None:
        return hit.freq_mhz * 1e6, hit.drift_hz_s
    factor = doppler_factor(observer, scan.mjd, sky)
    bary_freq_hz = hit.freq_mhz * 1e6 / factor
    drift = barycentric_drift(hit.drift_hz_s, hit.freq_mhz * 1e6, observer, scan.mjd, sky)
    # a drifting source's own frequency changes with time: refer it to one epoch
    elapsed = (float(scan.mjd) - float(epoch_ref)) * 86400.0
    return bary_freq_hz - drift * elapsed, drift


def _links(
    event: Event,
    freq_hz: float,
    drift: float,
    topo_freq_hz: float,
    topo_drift_hz_s: float,
    drift_tolerance_hz_s: float,
    freq_window_hz: float,
) -> bool:
    """Two hits are the same event if they agree in either frame.

    Barycentric agreement is the source test: it survives the drift wandering as
    Earth turns. Topocentric agreement is the interference test: a carrier sits
    at the same frequency and drift in every scan no matter the epoch, so it
    links too - and is then thrown out by the OFF scans rather than vanishing
    from the report without a verdict.
    """
    if abs(drift - event.drift_hz_s) <= drift_tolerance_hz_s and abs(
        freq_hz - event.freq_hz
    ) <= freq_window_hz:
        return True
    return abs(topo_drift_hz_s - event.topo_drift_hz_s) <= drift_tolerance_hz_s and abs(
        topo_freq_hz - event.topo_freq_hz
    ) <= freq_window_hz


def scan_sensitivity(meta: dict) -> float:
    """How strongly a scan measures a fixed intrinsic signal, in sqrt(N*W).

    The detection statistic is a dedoppler sum over N integrations, boxcar
    smoothed over W channels, normalised by the robust noise level. For a
    signal that adds a fixed amount to each sample, SNR = A * sqrt(N * W): the
    signal grows as N*W and the noise as its square root. So this is the factor
    that turns a measured SNR back into the source amplitude, and its square is
    the inverse variance of that estimate.

    It deliberately does not involve the detection threshold. The threshold
    decides which hits exist, never how strong one is: two hits at 12 sigma are
    equally strong whether the cut was 10 or 20. What the threshold does govern
    is whether an *absence* is informative, which is what _check_completeness
    uses it for.

    Assumes comparable channel widths. A scan with much wider channels puts
    more continuum flux behind a narrowband signal, and then SNR is no longer
    proportional to sqrt(N*W) alone; for such a run, combine scans of matching
    resolution or normalise the files first.
    """
    n_ints = float(meta.get("n_ints", 1) or 1)
    boxcar = float(meta.get("boxcar", 1) or 1)
    return float(np.sqrt(max(n_ints, 1.0) * max(boxcar, 1.0)))


def spectrum_flux(
    snr: float,
    meta: dict,
    sefd_jy: float | None,
    line_width_channels: float = 1.0,
) -> tuple[float, float] | None:
    """Flux density and its uncertainty for one hit in one file, in Jy.

    The single-file form of what an event does across a run: the same noise
    model, applied to a spectrum whose sensitivity is fixed by the file's own
    integration count, so a candidate found in one scan can still be quoted on a
    physical scale. None when the caller supplied no SEFD.
    """
    if sefd_jy is None:
        return None
    sensitivity = scan_sensitivity(meta)
    df_hz = abs(float(meta.get("df_mhz", 0.0))) * 1e6
    tsamp_s = float(meta.get("tsamp", 0.0) or 0.0)
    if sensitivity <= 0 or df_hz <= 0 or tsamp_s <= 0:
        return None
    return flux_density_jy(
        amplitude=float(snr) / sensitivity,
        sensitivity=sensitivity,
        sefd_jy=sefd_jy,
        df_hz=df_hz,
        tsamp_s=tsamp_s,
        boxcar=int(meta.get("boxcar", 1) or 1),
        line_width_channels=line_width_channels,
    )


def _recompute(event: Event, line_width_channels: float = 1.0) -> None:
    snr = np.array([m.candidate.snr for m in event.members], dtype=np.float64)
    sens = np.array([m.sensitivity for m in event.members], dtype=np.float64)
    freqs = np.array([m.freq_hz for m in event.members], dtype=np.float64)
    drifts = np.array([m.drift_hz_s for m in event.members], dtype=np.float64)
    topo_freqs = np.array(
        [m.candidate.freq_mhz * 1e6 for m in event.members], dtype=np.float64
    )
    topo_drifts = np.array([m.candidate.drift_hz_s for m in event.members], dtype=np.float64)

    # inverse-variance combination: each scan's amplitude estimate is
    # snr/sensitivity with variance 1/sensitivity^2, so the best estimate is
    # the sensitivity-weighted mean and its significance is that mean times the
    # total sensitivity of the scans that saw it. With equal scans this reduces
    # to sum of squares, the old behaviour.
    total = float(np.sqrt(np.sum(sens**2)))
    event.total_sensitivity = total
    if total > 0:
        event.score = float(np.sum(snr * sens) / total)
        event.amplitude = event.score / total
    else:
        event.score = float(snr.max())
        event.amplitude = 0.0
    if sens.size > 1:
        amplitudes = snr / np.where(sens > 0, sens, 1.0)
        event.amplitude_scatter = float(np.std(amplitudes))
    else:
        event.amplitude_scatter = 0.0

    # Physical flux, when the scans carry a SEFD. Each scan contributes its own
    # estimate from its own noise level and its own receiver, combined by
    # inverse variance - that is what makes a run across two telescopes add up
    # to one number instead of a weighted guess.
    estimates = []
    for member in event.members:
        if member.sefd_jy is None:
            continue
        flux, error = flux_density_jy(
            amplitude=member.amplitude,
            sensitivity=member.sensitivity,
            sefd_jy=member.sefd_jy,
            df_hz=member.df_hz,
            tsamp_s=member.tsamp_s,
            boxcar=member.boxcar,
            line_width_channels=line_width_channels,
        )
        estimates.append((flux, error, member.scan))
    event.flux = combine_fluxes(estimates, line_width_channels=line_width_channels)

    weights = snr
    event.freq_hz = float(np.average(freqs, weights=weights))
    event.drift_hz_s = float(np.average(drifts, weights=weights))
    event.topo_freq_hz = float(np.average(topo_freqs, weights=weights))
    event.topo_drift_hz_s = float(np.average(topo_drifts, weights=weights))
    event.snr = float(snr.max())
    event.freq_scatter_mhz = float(np.std(freqs) / 1e6) if freqs.size > 1 else 0.0
    event.drift_scatter_hz_s = float(np.std(drifts)) if drifts.size > 1 else 0.0


def _raw_spread(event: Event) -> float:
    raw = [m.candidate.drift_hz_s for m in event.members]
    return float(max(raw) - min(raw)) if raw else 0.0


def _check_completeness(event: Event, on_scans: list[Scan]) -> None:
    """Did any ON scan fail to see an event it was deep enough to see?

    This is where the threshold earns its place. From the combined amplitude,
    the SNR an absent scan *should* have shown is predictable; if that beats the
    scan's own threshold, its silence is a real discrepancy, not a shallow scan
    doing its best. Reported rather than scored: a source that fades is a
    physical possibility, and turning a fade into a lower score would hide it.
    """
    if not event.members or event.amplitude <= 0:
        return
    present = {m.scan for m in event.members}
    for scan in on_scans:
        if scan.path in present:
            continue
        expected = event.amplitude * scan_sensitivity(scan.meta)
        event.expected_snr[scan.path] = expected
        threshold = float(scan.meta.get("snr_threshold", 10.0))
        if expected >= threshold:
            event.missed_scans.append(scan.path)


def _test_off_scans(event, off_scans, observer, sky, can_convert, *, freq_tolerance_hz):
    """Test each ON member against each OFF scan.

    Per member, not per event: a member is a real observation of a real
    direction, and the two hypotheses (the source where the prediction puts it,
    interference where the telescope's own frame puts it) are each defined by
    that observation's epoch.
    """
    ref = _epoch_ref(event)
    for scan in off_scans:
        df_hz = abs(float(scan.meta["df_mhz"])) * 1e6
        static_window_mhz = max(3.0 * df_hz, 3.0) / 1e6
        for member in event.members:
            anchors = [(member.candidate.freq_mhz, static_window_mhz)]
            if can_convert and scan.mjd is not None and member.mjd is not None:
                elapsed = (float(scan.mjd) - ref) * 86400.0
                predicted_topo = topocentric_frequency(
                    member.freq_hz + member.drift_hz_s * elapsed,
                    observer,
                    scan.mjd,
                    sky,
                )
                window_mhz = max(freq_tolerance_hz, 4.0 * df_hz) / 1e6
                anchors.append((predicted_topo / 1e6, window_mhz))
            for hit in scan.hits:
                matched = [
                    width
                    for anchor, width in anchors
                    if abs(hit.freq_mhz - anchor) <= width
                ]
                if matched:
                    event.off_matches.append(
                        OffMatch(
                            scan=scan.path,
                            hit=hit,
                            window_hz=2.0 * max(matched) * 1e6,
                        )
                    )
                    event.status = SEEN_IN_OFF
                    return


def _epoch_ref(event: Event) -> float:
    epochs = [m.mjd for m in event.members if m.mjd is not None]
    return min(epochs) if epochs else 0.0
