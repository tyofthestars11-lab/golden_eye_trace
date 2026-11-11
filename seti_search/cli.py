"""Command line interface: `info`, `search`, `simulate`, `verify`, `link`."""

from __future__ import annotations

import argparse
import csv
import json
import sys

import numpy as np

from .cadence import SITES, Observer, SkyPosition
from .filterbank import Filterbank, InjectedSignal, RFILine, simulate
from .flux import COMPLETENESS_MARGIN, limit_jy, run_limit, sefd_from_gain
from .rungs import repeated_scales, rung_uniformity
from .run import Scan, link_run, scan_sensitivity, spectrum_flux
from .trials import (
    drift_trials,
    expected_false_alarms,
    log10_tail,
    sigma_for_band_false_alarm,
    tail_probability,
)
from .search import SearchConfig, search_file
from .verify import verify_cadence

CSV_FIELDS = [
    "freq_mhz",
    "drift_hz_s",
    "snr",
    "channel",
    "source",
    "tstart",
    "tsamp",
    "n_ints",
    "total_time_s",
    "boxcar",
    "snr_threshold",
    "flux_jy",
    "flux_error_jy",
]

DISCLAIMER = (
    "A candidate is not a detection: a real signal must reappear in an ON scan and "
    "vanish in an OFF scan at a drift consistent with a fixed sky position."
)

VERIFY_FIELDS = [
    "freq_mhz",
    "drift_hz_s",
    "snr",
    "channel",
    "status",
    "reason",
    "predicted_freq_mhz",
    "predicted_drift_hz_s",
    "off_freq_mhz",
    "off_drift_hz_s",
    "off_snr",
    "flux_jy",
    "flux_error_jy",
    "threshold_sigma",
    "trial_cells",
    "snr_tail_log10",
    "expected_noise_cells",
]


def _verdict_table(
    verdicts,
    on_meta=None,
    sefd_jy=None,
    line_width_channels: float = 1.0,
    sigma: float | None = None,
    n_cells: int = 0,
) -> str:
    show_flux = on_meta is not None and sefd_jy is not None
    header = (
        f"{'#':>3}  {'frequency (MHz)':>16}  {'drift (Hz/s)':>12}  {'SNR':>6}"
    )
    if show_flux:
        header += f"  {'flux (Jy)':>12}"
    header += "  verdict"
    lines = [header, "-" * len(header)]
    for index, verdict in enumerate(verdicts, 1):
        hit = verdict.on
        row = (
            f"{index:>3}  {hit.freq_mhz:>16.6f}  {hit.drift_hz_s:>+12.4f}"
            f"  {hit.snr:>6.1f}"
        )
        if show_flux:
            flux = spectrum_flux(hit.snr, on_meta, sefd_jy, line_width_channels)
            row += f"  {'-':>12}" if flux is None else f"  {flux[0]:>12.4g}"
        lines.append(row + f"  {verdict.status}")

    for index, verdict in enumerate(verdicts, 1):
        detail = [f"  #{index} {verdict.status}: {verdict.reason}"]
        if show_flux:
            flux = spectrum_flux(verdict.on.snr, on_meta, sefd_jy, line_width_channels)
            if flux is not None:
                detail.append(
                    f"      flux {flux[0]:.4g} +/- {flux[1]:.2g} Jy in ON",
                )
        if verdict.predicted_freq_mhz is not None:
            predicted = f"      predicted in OFF: {verdict.predicted_freq_mhz:.6f} MHz"
            if verdict.predicted_drift_hz_s is not None:
                predicted += f", {verdict.predicted_drift_hz_s:+.4f} Hz/s"
                if verdict.drift_discrimination_hz_s is not None:
                    predicted += (
                        f"  (drift moves {verdict.drift_discrimination_hz_s:.4f} Hz/s"
                        " between scans)"
                    )
            detail.append(predicted)
        if verdict.off_hit is not None:
            detail.append(
                f"      OFF hit: {verdict.off_hit.freq_mhz:.6f} MHz @"
                f" {verdict.off_hit.drift_hz_s:+.4f} Hz/s, {verdict.off_hit.snr:.1f} sigma"
            )
        if not verdict.confirmed and sigma is not None and n_cells > 0:
            detail.append(_rejection_budget(verdict.on, sigma, n_cells))
        lines.append("\n".join(detail))
    return "\n".join(lines)


def _write_verdicts(
    path: str,
    verdicts,
    on_meta=None,
    sefd_jy=None,
    line_width_channels: float = 1.0,
    sigma: float | None = None,
    n_cells: int = 0,
) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=VERIFY_FIELDS)
        writer.writeheader()
        for verdict in verdicts:
            off = verdict.off_hit
            flux = (
                None
                if on_meta is None
                else spectrum_flux(
                    verdict.on.snr, on_meta, sefd_jy, line_width_channels
                )
            )
            writer.writerow(
                {
                    "freq_mhz": f"{verdict.on.freq_mhz:.9f}",
                    "drift_hz_s": f"{verdict.on.drift_hz_s:.6f}",
                    "snr": f"{verdict.on.snr:.2f}",
                    "channel": verdict.on.channel,
                    "status": verdict.status,
                    "reason": verdict.reason,
                    "predicted_freq_mhz": (
                        "" if verdict.predicted_freq_mhz is None
                        else f"{verdict.predicted_freq_mhz:.9f}"
                    ),
                    "predicted_drift_hz_s": (
                        "" if verdict.predicted_drift_hz_s is None
                        else f"{verdict.predicted_drift_hz_s:.6f}"
                    ),
                    "off_freq_mhz": "" if off is None else f"{off.freq_mhz:.9f}",
                    "off_drift_hz_s": "" if off is None else f"{off.drift_hz_s:.6f}",
                    "off_snr": "" if off is None else f"{off.snr:.2f}",
                    "flux_jy": "" if flux is None else f"{flux[0]:.6g}",
                    "flux_error_jy": "" if flux is None else f"{flux[1]:.4g}",
                    "threshold_sigma": "" if sigma is None else f"{sigma:.6g}",
                    "trial_cells": "" if not n_cells else f"{int(n_cells)}",
                    "snr_tail_log10": f"{log10_tail(verdict.on.snr):.4g}",
                    "expected_noise_cells": (
                        ""
                        if not n_cells
                        else f"{tail_probability(verdict.on.snr) * n_cells:.6g}"
                    ),
                }
            )


def _observation_mjd(fb) -> float | None:
    """Start time of a file in MJD, if its header says one."""
    header = fb.header
    if "MJD_start" in header:
        seconds = (
            float(header.get("MJD_hour", 0)) * 3600.0
            + float(header.get("MJD_minute", 0)) * 60.0
            + float(header.get("MJD_second", 0))
        )
        return float(header["MJD_start"]) + seconds / 86400.0
    tstart = float(header.get("tstart", 0.0) or 0.0)
    if 30000.0 < tstart < 80000.0:
        return tstart  # SIGPROC writes MJD here
    if tstart > 1e8:
        return tstart / 86400.0 + 40587.0  # unix seconds
    return None


def _sexagesimal(value: float, *, right_ascension: bool) -> float:
    """HHMMSS.SS / DDMMSS.SS as written in SIGPROC headers -> degrees."""
    sign = -1.0 if value < 0 else 1.0
    magnitude = abs(float(value))
    hours = int(magnitude // 10000)
    minutes = int((magnitude % 10000) // 100)
    seconds = magnitude % 100
    sexagesimal = hours + minutes / 60.0 + seconds / 3600.0
    return sign * sexagesimal * (15.0 if right_ascension else 1.0)


def _header_sky(fb):
    header = fb.header
    if "src_raj" not in header or "src_dej" not in header:
        return None
    ra_raw = float(header["src_raj"])
    dec_raw = float(header["src_dej"])
    if ra_raw == 0.0 and dec_raw == 0.0:
        return None
    dec = _sexagesimal(dec_raw, right_ascension=False)
    if not -90.0 <= dec <= 90.0:
        return None
    return SkyPosition(_sexagesimal(ra_raw, right_ascension=True) % 360.0, dec)


def _resolve_sky(args, fb=None):
    if args.ra is not None and args.dec is not None:
        return SkyPosition(float(args.ra), float(args.dec))
    if fb is not None and args.ra is None and args.dec is None:
        return _header_sky(fb)
    return None


def _resolve_observer(args):
    if args.lat is not None and args.lon is not None:
        return Observer(float(args.lat), float(args.lon))
    if args.site:
        key = args.site.strip().lower()
        if key not in SITES:
            raise SystemExit(
                f"unknown --site {args.site!r}; known sites: {', '.join(sorted(SITES))}"
            )
        return SITES[key]
    return None


def _table(hits, meta=None, sefd_jy=None, line_width_channels: float = 1.0) -> str:
    show_flux = meta is not None and sefd_jy is not None
    header = (
        f"{'#':>3}  {'frequency (MHz)':>16}  {'drift (Hz/s)':>13}"
        f"  {'SNR':>7}  {'channel':>9}"
    )
    if show_flux:
        header += f"  {'flux (Jy)':>12}  {'+/-':>8}"
    lines = [header, "-" * len(header)]
    for index, hit in enumerate(hits, 1):
        row = (
            f"{index:>3}  {hit.freq_mhz:>16.6f}  {hit.drift_hz_s:>+13.4f}"
            f"  {hit.snr:>7.1f}  {hit.channel:>9d}"
        )
        if show_flux:
            flux = spectrum_flux(hit.snr, meta, sefd_jy, line_width_channels)
            if flux is None:
                row += f"  {'-':>12}  {'-':>8}"
            else:
                row += f"  {flux[0]:>12.4g}  {flux[1]:>8.2g}"
        lines.append(row)
    return "\n".join(lines)


def _write_csv(
    path: str,
    hits,
    meta: dict,
    sefd_jy: float | None = None,
    line_width_channels: float = 1.0,
) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for hit in hits:
            flux = spectrum_flux(hit.snr, meta, sefd_jy, line_width_channels)
            writer.writerow(
                {
                    "freq_mhz": f"{hit.freq_mhz:.9f}",
                    "drift_hz_s": f"{hit.drift_hz_s:.6f}",
                    "snr": f"{hit.snr:.2f}",
                    "channel": hit.channel,
                    "source": meta["source"],
                    "tstart": f"{meta['tstart']:.6f}",
                    "tsamp": meta["tsamp"],
                    "n_ints": meta["n_ints"],
                    "total_time_s": f"{meta['total_time']:.3f}",
                    "boxcar": meta["boxcar"],
                    "snr_threshold": meta["snr_threshold"],
                    "flux_jy": "" if flux is None else f"{flux[0]:.6g}",
                    "flux_error_jy": "" if flux is None else f"{flux[1]:.4g}",
                }
            )


def cmd_info(args) -> int:
    fb = Filterbank.open(args.file)
    print(fb.describe())
    print("\nheader")
    for key, value in fb.header.items():
        print(f"  {key:<16} {value}")
    return 0


def cmd_search(args) -> int:
    fb = Filterbank.open(args.file)
    if args.false_alarm is not None:
        _apply_false_alarm(args, _nominal_trials(fb, args), quiet=args.quiet)

    config = SearchConfig(
        boxcar=args.boxcar,
        max_drift_hz_s=args.max_drift,
        drift_step_hz_s=args.drift_step,
        snr_threshold=args.snr,
        max_hits=args.max_hits,
        rfi_sigma=args.rfi_sigma,
        carrier_sigma=args.carrier_sigma,
    )
    if not args.quiet:
        tsamp_shown = args.tsamp if args.tsamp else fb.tsamp
        print(f"searching {args.file}")
        print(
            f"  source {fb.source}, {fb.nchans:,} channels"
            f" x {fb.nints:,} integrations @ {tsamp_shown:g} s"
        )
        print(
            f"  drift range +/-{args.max_drift:g} Hz/s, boxcar {args.boxcar},"
            f" detection threshold {args.snr:g} sigma"
        )

    def progress(c0: int, c1: int, total: int, found: int) -> None:
        print(f"  ... {c1:,}/{total:,} channels, {found} candidate(s)", file=sys.stderr)

    hits, meta = search_file(
        args.file,
        config,
        chan_chunk=args.chan_chunk,
        max_ints=(args.max_ints or None),
        progress=None if args.quiet else progress,
        tsamp=args.tsamp,
    )

    if not args.quiet:
        blanked = meta["blanked_channels"]
        print(
            f"  blanked {blanked:,}/{meta['nchans']:,} channels"
            f" and {meta['blanked_rows']}/{meta['n_ints']} integrations as RFI"
        )
        calibration = float(meta.get("snr_calibration", 1.0))
        if abs(calibration - 1.0) > 0.02:
            print(
                f"  noise calibration: this file's dedoppler map is {calibration:.3f}x"
                f" wider than the robust estimate, so the sigmas below are nominal"
                f" divided by {calibration:.3f}"
            )
        tails = meta.get("snr_tail") or {}
        if tails:
            parts = []
            for sigma in sorted(tails):
                entry = tails[sigma]
                if entry["expected"] > 0.0:
                    parts.append(
                        f"{entry['observed']:,} above {sigma:g} sigma"
                        f" ({entry['observed'] / entry['expected']:.3g}x the Gaussian count)"
                    )
                else:
                    parts.append(
                        f"{entry['observed']:,} above {sigma:g} sigma (none expected)"
                    )
            print("  map tails: " + "; ".join(parts))
        if meta["truncated"]:
            print(
                f"  note: searched the first {meta['n_ints']:,} of"
                f" {meta['n_ints_in_file']:,} integrations (use --max-ints 0 for all)"
            )
        print()

    sefd = _scan_sefd(
        args.file,
        _resolve_quantity(args.sefd),
        _resolve_quantity(args.tsys),
        _resolve_quantity(args.gain),
    )
    if not hits:
        print(f"no candidates at {args.snr:g} sigma or above")
    else:
        print(_table(hits, meta, sefd, args.line_width))

    if sefd is not None:
        fifty, ninety = _limit_values(args, meta, sefd)
        lo, hi = _band_edges(meta)
        searched = 100.0 * (1.0 - meta["blanked_channels"] / max(meta["nchans"], 1))
        print(
            f"\nlimit: any {args.line_width:g}-channel line brighter than"
            f" {fifty:.3g} Jy in {lo:.6f}-{hi:.6f} MHz at |drift| <="
            f" {args.max_drift:g} Hz/s would have been found"
            f" ({ninety:.3g} Jy for 95% completeness), over the {searched:.1f}%"
            " of channels searched"
        )
        print(_false_alarm_note(args.snr, _band_trials(meta, args)[1]))

    out = args.out or f"{args.file}.hits.csv"
    _write_csv(out, hits, meta, sefd, args.line_width)
    print(f"\nwrote {len(hits)} candidate(s) to {out}")
    if sefd is not None:
        limits_out = _limits_path(out)
        _write_limits(limits_out, [_limit_row(args, args.file, meta, sefd)])
        print(f"wrote 1 limit row to {limits_out}")
    if hits:
        print(DISCLAIMER)
    return 0


def cmd_verify(args) -> int:
    on_fb = Filterbank.open(args.on)
    off_fb = Filterbank.open(args.off)
    # the threshold has to hold for both scans, so the wider of the two budgets
    # is the one it was set against
    n_cells = max(_nominal_trials(on_fb, args), _nominal_trials(off_fb, args))
    if args.false_alarm is not None:
        _apply_false_alarm(args, n_cells, quiet=args.quiet)

    config = SearchConfig(
        boxcar=args.boxcar,
        max_drift_hz_s=args.max_drift,
        drift_step_hz_s=args.drift_step,
        snr_threshold=args.snr,
        max_hits=args.max_hits,
        rfi_sigma=args.rfi_sigma,
        carrier_sigma=args.carrier_sigma,
    )

    def progress(label):
        if args.quiet:
            return None

        def report(c0: int, c1: int, total: int, found: int) -> None:
            print(
                f"  {label} {c1:,}/{total:,} channels, {found} candidate(s)",
                file=sys.stderr,
            )

        return report

    on_hits, on_meta = search_file(
        args.on,
        config,
        chan_chunk=args.chan_chunk,
        max_ints=(args.max_ints or None),
        tsamp=args.tsamp,
        progress=progress("ON "),
    )
    off_hits, off_meta = search_file(
        args.off,
        config,
        chan_chunk=args.chan_chunk,
        max_ints=(args.max_ints or None),
        tsamp=args.tsamp,
        progress=progress("OFF"),
    )

    observer = _resolve_observer(args)
    sky = _resolve_sky(args, on_fb)
    on_mjd = args.on_mjd if args.on_mjd is not None else _observation_mjd(on_fb)
    off_mjd = args.off_mjd if args.off_mjd is not None else _observation_mjd(off_fb)

    result = verify_cadence(
        on_hits,
        off_hits,
        on_meta=on_meta,
        off_meta=off_meta,
        observer=observer,
        sky=sky,
        on_mjd=on_mjd,
        off_mjd=off_mjd,
        drift_tolerance_hz_s=args.drift_tol,
        freq_tolerance_hz=args.freq_tol,
    )

    if not args.quiet:
        print(f"  {len(on_hits)} candidate(s) in ON, {len(off_hits)} in OFF at {args.snr:g} sigma")
        if observer is not None and sky is not None and on_mjd is not None and off_mjd is not None:
            print(f"  site {observer}, target RA {sky.ra_deg:.3f} Dec {sky.dec_deg:+.3f}")
            print(
                f"  epochs {on_mjd:.6f} / {off_mjd:.6f} MJD"
                f"  (cadence {abs(off_mjd - on_mjd) * 1440.0:.1f} min)"
            )
        print()

    on_sefd = _scan_sefd(
        args.on,
        _resolve_quantity(args.sefd),
        _resolve_quantity(args.tsys),
        _resolve_quantity(args.gain),
    )
    if not result.verdicts:
        print("nothing to verify: the ON scan produced no candidates")
    else:
        print(
            _verdict_table(
                result.verdicts,
                on_meta,
                on_sefd,
                args.line_width,
                sigma=args.snr,
                n_cells=n_cells,
            )
        )
        print(f"\n{len(result.confirmed)} confirmed, {len(result.rejected)} rejected")
        print(_false_alarm_note(args.snr, n_cells))
    for note in result.notes:
        print(f"note: {note}")
    if result.confirmed:
        print(
            "\nA confirmed candidate behaved like a sky source in this cadence."
            " That is not a detection: repeat the cadence, and check the drift"
            " against the position."
        )

    out = args.out or f"{args.on}.verified.csv"
    _write_verdicts(
        out,
        result.verdicts,
        on_meta,
        on_sefd,
        args.line_width,
        sigma=args.snr,
        n_cells=n_cells,
    )
    print(f"\nwrote {len(result.verdicts)} verdict(s) to {out}")
    return 0


LIMIT_FIELDS = [
    "kind",
    "path",
    "threshold_sigma",
    "line_width_channels",
    "max_drift_hz_s",
    "freq_min_mhz",
    "freq_max_mhz",
    "searched_fraction",
    "sefd_jy",
    "n_ints",
    "boxcar",
    "df_hz",
    "tsamp_s",
    "limit_jy",
    "limit_jy_95",
    "trials_per_channel",
    "total_trials",
    "false_alarm_log10",
    "expected_false_alarms",
    "scans",
    "min_scans",
    "set_by",
]

LINK_FIELDS = [
    "freq_mhz",
    "bary_freq_mhz",
    "drift_hz_s",
    "snr",
    "score",
    "n_on",
    "n_scans",
    "status",
    "scans",
    "drift_scatter_hz_s",
    "amplitude",
    "amplitude_scatter",
    "missed",
    "flux_jy",
    "flux_error_jy",
    "threshold_sigma",
    "trial_cells",
    "chance_log10",
    "off_window_hz",
    "off_cells",
    "off_chance_log10",
    "rung_freq",
    "rung_freq_nearest",
    "rung_freq_drift",
    "rung_drift",
    "rung_drift_nearest",
    "rung_drift_drift",
]


def _parse_scan_spec(text: str) -> tuple[str, float | None]:
    """PATH, or PATH@MJD to override the epoch the header claims."""
    if "@" in text:
        path, _, mjd = text.rpartition("@")
        try:
            return path, float(mjd)
        except ValueError:
            raise SystemExit(f"bad scan {text!r}: expected PATH or PATH@MJD")
    return text, None


def _rung_report(events) -> str:
    lines = [
        "",
        "golden-angle rungs: rung(x) = ln x / ln phi. These can only detract from a",
        "candidate, never support one - the fractional part of a rung is uniform for",
        "any smooth population, so a small drift is not a finding by itself.",
    ]
    for index, event in enumerate(events, 1):
        parts = []
        if event.rung_freq is not None:
            parts.append(
                f"freq {event.rung_freq.rung:.4f} (near {event.rung_freq.nearest},"
                f" drift {event.rung_freq.drift:.4f})"
            )
        if event.rung_drift is not None:
            parts.append(
                f"drift {event.rung_drift.rung:.4f} (near {event.rung_drift.nearest},"
                f" drift {event.rung_drift.drift:.4f})"
            )
        else:
            parts.append("drift rung undefined at zero")
        if event.rung_snr is not None:
            parts.append(
                f"snr {event.rung_snr.rung:.4f} (near {event.rung_snr.nearest},"
                f" drift {event.rung_snr.drift:.4f})"
            )
        lines.append(f"  #{index}  " + "   ".join(parts))

    freqs = [event.freq_hz for event in events]
    uniform = rung_uniformity(freqs)
    if uniform["n"] >= 4:
        verdict = (
            "uniform"
            if uniform["p_value"] > 0.01
            else "NOT uniform - the population looks quantized by the pipeline"
        )
        lines.append(
            f"  frequencies: {uniform['n']} values, {verdict}"
            f" (KS p = {uniform['p_value']:.3g})"
        )
    else:
        lines.append(
            f"  frequencies: {uniform['n']} value(s), too few for a uniformity test"
        )

    for label, values in (
        ("frequencies", freqs),
        ("drifts", [abs(e.drift_hz_s) for e in events if e.drift_hz_s]),
    ):
        found = repeated_scales(values, min_pairs=3)
        if found:
            for repetition in found:
                lines.append(
                    f"  {label}: RECURRING SCALE - {repetition} - one scaling"
                    " operation, so an artifact of the processing"
                )
        else:
            lines.append(f"  {label}: no repeated scale factor across independent pairs")
    return "\n".join(lines)


def _write_links(path: str, events, budgets=None) -> None:
    budgets = budgets or []
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LINK_FIELDS)
        writer.writeheader()
        for index, event in enumerate(events):
            budget = budgets[index] if index < len(budgets) else None
            off = None if budget is None else budget["off"]
            strongest = max(event.members, key=lambda m: m.candidate.snr)
            writer.writerow(
                {
                    "freq_mhz": f"{strongest.candidate.freq_mhz:.9f}",
                    "bary_freq_mhz": f"{event.freq_mhz:.9f}",
                    "drift_hz_s": f"{event.drift_hz_s:.6f}",
                    "snr": f"{event.snr:.2f}",
                    "score": f"{event.score:.2f}",
                    "n_on": event.n_on,
                    "n_scans": event.n_scans,
                    "status": event.status,
                    "scans": ";".join(m.scan for m in event.members),
                    "drift_scatter_hz_s": f"{event.drift_scatter_hz_s:.6f}",
                    "amplitude": f"{event.amplitude:.4f}",
                    "amplitude_scatter": f"{event.amplitude_scatter:.4f}",
                    "missed": ";".join(event.missed_scans),
                    "flux_jy": (
                        "" if event.flux is None else f"{event.flux.jy:.6g}"
                    ),
                    "flux_error_jy": (
                        "" if event.flux is None else f"{event.flux.uncertainty_jy:.4g}"
                    ),
                    "threshold_sigma": (
                        "" if budget is None else f"{budget['sigma']:.6g}"
                    ),
                    "trial_cells": "" if budget is None else f"{budget['cells_total']}",
                    "chance_log10": (
                        "" if budget is None else f"{budget['chance_log10']:.4g}"
                    ),
                    "off_window_hz": (
                        "" if off is None else f"{off['window_hz']:.4g}"
                    ),
                    "off_cells": "" if off is None else f"{off['cells']}",
                    "off_chance_log10": (
                        "" if off is None else f"{off['chance_log10']:.4g}"
                    ),
                    "rung_freq": (
                        "" if event.rung_freq is None else f"{event.rung_freq.rung:.6f}"
                    ),
                    "rung_freq_nearest": (
                        "" if event.rung_freq is None else event.rung_freq.nearest
                    ),
                    "rung_freq_drift": (
                        "" if event.rung_freq is None else f"{event.rung_freq.drift:.6f}"
                    ),
                    "rung_drift": (
                        "" if event.rung_drift is None else f"{event.rung_drift.rung:.6f}"
                    ),
                    "rung_drift_nearest": (
                        "" if event.rung_drift is None else event.rung_drift.nearest
                    ),
                    "rung_drift_drift": (
                        "" if event.rung_drift is None else f"{event.rung_drift.drift:.6f}"
                    ),
                }
            )


def _resolve_quantity(specs) -> tuple[float | None, dict[str, float]]:
    """Parse repeatable VALUE or PATH=VALUE flags into a default and overrides."""
    default: float | None = None
    by_path: dict[str, float] = {}
    for spec in specs or []:
        if "=" in spec:
            path, _, value = spec.partition("=")
            try:
                by_path[path.strip()] = float(value)
            except ValueError:
                raise SystemExit(f"bad value in {spec!r}: expected PATH=NUMBER or NUMBER")
        else:
            try:
                default = float(spec)
            except ValueError:
                raise SystemExit(f"bad value {spec!r}: expected PATH=NUMBER or NUMBER")
    return default, by_path


def _scan_sefd(path, sefd, tsys, gain) -> float | None:
    """SEFD in Jy for one scan: explicit value first, then Tsys divided by gain."""
    default, by_path = sefd
    if path in by_path:
        return by_path[path]
    if default is not None:
        return default
    t = tsys[1].get(path, tsys[0])
    g = gain[1].get(path, gain[0])
    if t is not None and g:
        return sefd_from_gain(t, g)
    return None


def cmd_link(args) -> int:
    if not args.on:
        raise SystemExit("give at least one --on scan")

    entries = [(spec, "on") for spec in args.on]
    entries += [(spec, "off") for spec in args.off]

    if args.false_alarm is not None:
        # the threshold has to hold for every scan searched, so it comes from the
        # file testing the most cells rather than from whichever ran first
        worst = 0
        for spec, _kind in entries:
            path, _mjd = _parse_scan_spec(spec)
            worst = max(worst, _nominal_trials(Filterbank.open(path), args))
        _apply_false_alarm(args, worst, quiet=args.quiet)

    config = SearchConfig(
        boxcar=args.boxcar,
        max_drift_hz_s=args.max_drift,
        drift_step_hz_s=args.drift_step,
        snr_threshold=args.snr,
        max_hits=args.max_hits,
        rfi_sigma=args.rfi_sigma,
        carrier_sigma=args.carrier_sigma,
    )

    sefd_flags = _resolve_quantity(args.sefd)
    tsys_flags = _resolve_quantity(args.tsys)
    gain_flags = _resolve_quantity(args.gain)
    calibrated = 0

    scans: list[Scan] = []
    first_fb = None
    for spec, kind in entries:
        path, mjd = _parse_scan_spec(spec)
        fb = Filterbank.open(path)
        first_fb = first_fb or fb
        hits, meta = search_file(
            path,
            config,
            chan_chunk=args.chan_chunk,
            max_ints=(args.max_ints or None),
            tsamp=args.tsamp,
        )
        resolved = mjd if mjd is not None else _observation_mjd(fb)
        sefd = _scan_sefd(path, sefd_flags, tsys_flags, gain_flags)
        calibrated += sefd is not None
        scans.append(
            Scan(path=path, kind=kind, mjd=resolved, hits=hits, meta=meta, sefd_jy=sefd)
        )
        if not args.quiet:
            epoch = "no epoch" if resolved is None else f"{resolved:.6f} MJD"
            print(
                f"  {kind.upper():3} {path}  {len(hits)} candidate(s)  {epoch}"
                f"  sensitivity {scan_sensitivity(meta):.1f}"
                + ("" if sefd is None else f"  SEFD {sefd:.3g} Jy")
            )

    by_path = {scan.path: scan for scan in scans}
    observer = _resolve_observer(args)
    sky = _resolve_sky(args, first_fb)
    result = link_run(
        scans,
        observer=observer,
        sky=sky,
        drift_tolerance_hz_s=args.drift_tol,
        freq_tolerance_hz=args.freq_tol,
        min_scans=args.min_scans,
        line_width_channels=args.line_width,
    )
    if not args.quiet and 0 < calibrated < len(scans):
        print(
            f"  note: {len(scans) - calibrated} scan(s) had no SEFD;"
            " they contribute to the score but not to the flux density"
        )

    if not args.quiet:
        n_on = sum(1 for s in scans if s.kind == "on")
        print(
            f"\nrun: {n_on} ON, {len(scans) - n_on} OFF, "
            f"{sum(len(s.hits) for s in scans)} candidate(s)"
        )
        if result.predicted:
            print(f"  site {observer}, target RA {sky.ra_deg:.3f} Dec {sky.dec_deg:+.3f}")
        print()

    budgets = [_event_budget(event, by_path, args) for event in result.events]
    if not result.events:
        print(f"no event appears in {args.min_scans} or more ON scans")
    else:
        header = (
            f"  {'#':>3}  {'frequency (MHz)':>16}  {'drift (Hz/s)':>12}"
            f"  {'score':>7}  {'log10 p':>9}  {'ON':>3}  {'scans':>5}  status"
        )
        print(header)
        print("-" * len(header))
        for index, event in enumerate(result.events, 1):
            print(
                f"  {index:>3}  {event.freq_mhz:>16.6f}  {event.drift_hz_s:>+12.4f}"
                f"  {event.score:>7.1f}  {budgets[index - 1]['chance_log10']:>9.1f}"
                f"  {event.n_on:>3d}  {event.n_scans:>5d}  {event.status}"
            )

        for index, event in enumerate(result.events, 1):
            raws = [m.candidate.drift_hz_s for m in event.members]
            detail = [
                f"  #{index} {event.status}: the same barycentric drift in"
                f" {event.n_scans} ON scan(s)",
                f"      as measured: {'  '.join(f'{r:+.4f}' for r in raws)}"
                f"   (spread {max(raws) - min(raws):.4f} Hz/s)",
                f"      frequency scatter {event.freq_scatter_mhz * 1e3:.3f} kHz,"
                f" drift scatter {event.drift_scatter_hz_s:.4f} Hz/s",
                f"      amplitude {event.amplitude:.3f} sigma per unit sensitivity"
                f" (scatter {event.amplitude_scatter:.3f}),"
                f" combined score {event.score:.1f}",
                _event_budget_note(budgets[index - 1], event.n_scans),
            ]
            if event.flux is not None:
                detail.append(
                    f"      flux {event.flux.jy:.4g} +/- {event.flux.uncertainty_jy:.2g} Jy"
                    f" over {event.flux.n_scans} calibrated scan(s), assuming"
                    f" {args.line_width:g}-channel lines"
                )
            for off in event.off_matches:
                detail.append(
                    f"      OFF hit in {off.scan}: {off.hit.freq_mhz:.6f} MHz @"
                    f" {off.hit.drift_hz_s:+.4f} Hz/s, {off.hit.snr:.1f} sigma"
                )
                off_budget = budgets[index - 1]["off"]
                if off_budget is not None:
                    detail.append(_off_budget_note(off_budget, args.snr))
            for missed in event.missed_scans:
                detail.append(
                    f"      MISSED in {missed}: at this amplitude it should have"
                    f" shown {event.expected_snr[missed]:.1f} sigma there"
                )
            print("\n".join(detail))

        print(_rung_report(result.events))

    if calibrated:
        lines = [
            "",
            "limits: what the searched band is ruled out to, at |drift| <="
            f" {args.max_drift:g} Hz/s with {args.line_width:g}-channel lines",
            "  (50% / 95% completeness)",
        ]
        values: dict[str, tuple[float, float]] = {}
        rows = []
        worst_trials = 0
        worst_row = None
        for scan in scans:
            if scan.kind != "on" or scan.sefd_jy is None:
                continue
            fifty, ninety = _limit_values(args, scan.meta, scan.sefd_jy)
            values[scan.path] = (fifty, ninety)
            rows.append(_limit_row(args, scan.path, scan.meta, scan.sefd_jy))
            lines.append(f"  {scan.path}  {fifty:.3g} / {ninety:.3g} Jy")
            total = int(rows[-1]["total_trials"])
            if total > worst_trials:
                worst_trials, worst_row = total, rows[-1]
        best = run_limit({p: v[0] for p, v in values.items()}, args.min_scans)
        if best is not None:
            if result.survivors:
                lines.append(
                    f"  run  {best[0]:.3g} Jy - every other frequency brighter than"
                    f" this is ruled out across {args.min_scans} of"
                    f" {len(values)} ON scans"
                )
            else:
                lines.append(
                    f"  run  {best[0]:.3g} Jy - brighter than this is found in"
                    f" {args.min_scans} of {len(values)} ON scans, so this null"
                    " result rules it out"
                )
            run_row = {field: "" for field in LIMIT_FIELDS}
            run_row.update(
                kind="run",
                threshold_sigma=f"{args.snr:g}",
                line_width_channels=f"{args.line_width:g}",
                max_drift_hz_s=f"{args.max_drift:g}",
                limit_jy=f"{best[0]:.6g}",
                limit_jy_95=f"{values[best[1]][1]:.6g}",
                # the threshold had to hold for the scan testing the most cells,
                # so that scan's budget is the run's, not the shallower one's
                trials_per_channel=worst_row["trials_per_channel"],
                total_trials=worst_row["total_trials"],
                false_alarm_log10=worst_row["false_alarm_log10"],
                expected_false_alarms=worst_row["expected_false_alarms"],
                scans=str(len(values)),
                min_scans=str(args.min_scans),
                set_by=best[1],
            )
            rows.append(run_row)
        print("\n".join(lines))
        print(_false_alarm_note(args.snr, worst_trials))
        limits_out = _limits_path(args.out or "run.linked.csv")
        _write_limits(limits_out, rows)
        print(f"\nwrote {len(rows)} limit row(s) to {limits_out}")

    for note in result.notes:
        print(f"note: {note}")
    if result.survivors:
        print(
            "\nA surviving event is a candidate that behaved like a source across the"
            " run. It is still not a detection."
        )

    out = args.out or "run.linked.csv"
    _write_links(out, result.events, budgets)
    print(f"\nwrote {len(result.events)} event(s) to {out}")
    return 0


def _band_trials(meta: dict, args) -> tuple[int, int]:
    """Drift trials per channel, and total trial cells, for one scanned file.

    The cell count uses the whole band rather than only the surviving channels:
    the threshold has to be set before blanking is known, so the number it was
    set from travels with the result instead of being recomputed afterwards.
    """
    per_channel = drift_trials(
        n_ints=int(meta["n_ints"]),
        tsamp_s=float(meta["tsamp"]),
        df_mhz=float(meta["df_mhz"]),
        boxcar=int(meta["boxcar"]),
        max_drift_hz_s=args.max_drift,
        drift_step_hz_s=args.drift_step,
    )
    return per_channel, per_channel * int(max(meta["nchans"], 1))


def _nominal_trials(fb, args) -> int:
    """Trial cells in a whole file, from the header alone - before any searching."""
    n_used = int(fb.nints)
    if args.max_ints:
        n_used = min(n_used, int(args.max_ints))
    per_channel = drift_trials(
        n_ints=n_used,
        tsamp_s=float(args.tsamp or fb.tsamp),
        df_mhz=abs(float(fb.foff)),
        boxcar=args.boxcar,
        max_drift_hz_s=args.max_drift,
        drift_step_hz_s=args.drift_step,
    )
    return per_channel * int(max(fb.nchans, 1))


def _apply_false_alarm(args, total_trials: int, quiet: bool = False) -> None:
    """Turn a requested band false-alarm rate into the sigma actually searched.

    The derived threshold replaces --snr rather than competing with it, and the
    replacement is announced: a threshold nobody knows was changed is worse than
    no threshold at all.
    """
    sigma = sigma_for_band_false_alarm(args.false_alarm, total_trials)
    # quiet mode still hears this: it changes which search ran, so it is not
    # chatter. It goes to stderr so a script reading stdout sees only results.
    stream = sys.stderr if quiet else sys.stdout
    print(
        f"  threshold {sigma:.2f} sigma for a {args.false_alarm:g} chance of"
        f" any false alarm over {total_trials:,} trial cells",
        file=stream,
    )
    if abs(float(args.snr) - sigma) > 1e-9:
        print(
            f"  note: --snr {args.snr:g} is not used when --false-alarm is given",
            file=stream,
        )
    args.snr = sigma


def _false_alarm_note(sigma: float, n_trials: int) -> str:
    """What the threshold means in false-alarm terms, for the limits block."""
    p = tail_probability(sigma)
    expected = p * n_trials
    if p > 0:
        text = f"p = {p:.3g} per trial"
    else:
        text = f"p below double precision (log10 p = {log10_tail(sigma):.1f})"
    if expected > 0:
        text += f", {expected:.3g} expected over {n_trials:,} trial cells"
    else:
        text += f", none expected over {n_trials:,} trial cells"
    return f"  false alarm at this threshold: {text}"


def _rejection_budget(hit, sigma: float, n_cells: int) -> str:
    """A rejected candidate in the terms its threshold was set in.

    A rejection says the OFF scan saw the same transmitter. How much that proves
    depends on whether the ON hit was a real signal or a noise peak that cleared
    the threshold, so quote the noise count at both - the threshold it was tested
    against, and the significance it actually reached.
    """

    def count(value: float) -> str:
        return "none" if value == 0.0 else f"{value:.3g}"

    at_threshold = tail_probability(sigma) * n_cells
    at_hit = tail_probability(hit.snr) * n_cells
    lines = [
        f"      false-alarm budget: {count(at_threshold)} noise peaks expected over"
        f" {n_cells:,} cells at the {sigma:.2f} sigma threshold,"
        f" {count(at_hit)} at this {hit.snr:.1f} sigma hit"
    ]
    if at_hit < 1.0:
        lines.append(
            "      so the rejection reads as interference: noise does not reach"
            " this strength here, and the OFF detection is the same transmitter"
        )
    else:
        lines.append(
            "      so the rejection is not decisive: a band this size expects noise"
            " this strong, and a chance meeting with an OFF fluctuation would look"
            " the same"
        )
    return "\n".join(lines)


def _event_budget(event, by_path: dict, args) -> dict:
    """What noise would have to do to fake this event in every scan it appears.

    Each member is a peak in one scan's search, so the chance that noise alone
    puts a peak that strong somewhere in that band is cells * tail(snr), and
    independent scans multiply. The product ignores the frequency and drift
    agreement the link also requires, so it overstates a chance coincidence
    rather than understating it.
    """
    sigma = float(args.snr)
    cells_total = 0
    per_scan_cells = 0
    expected_per_scan = 0.0
    chance_log10 = 0.0
    for member in event.members:
        scan = by_path.get(member.scan)
        if scan is None:
            continue
        cells = _band_trials(scan.meta, args)[1]
        cells_total += cells
        expected = expected_false_alarms(sigma, cells)
        if expected > expected_per_scan:
            expected_per_scan, per_scan_cells = expected, cells
        chance_log10 += log10_tail(member.candidate.snr) + float(
            np.log10(max(cells, 1))
        )
    return {
        "sigma": sigma,
        "cells_total": cells_total,
        "per_scan_cells": per_scan_cells,
        "expected_per_scan": expected_per_scan,
        "chance_log10": chance_log10,
        "off": _off_budget(event, by_path, args),
    }


def _off_budget(event, by_path: dict, args) -> dict | None:
    """How surprising the OFF match is, in the same units as the threshold.

    The OFF test looks for a hit inside a frequency window, at any drift, so the
    cells it searched are the window's channels times that scan's drift trials.
    A rejection is only worth as much as the improbability of a noise peak
    landing in that window: quote that, not the bare fact that a peak was found.
    """
    if not event.off_matches:
        return None
    match = event.off_matches[0]
    scan = by_path.get(match.scan)
    if scan is None:
        return None
    df_hz = abs(float(scan.meta["df_mhz"])) * 1e6
    per_channel = _band_trials(scan.meta, args)[0]
    channels = max(1, int(round(match.window_hz / df_hz))) if df_hz > 0 else 1
    cells = channels * per_channel
    return {
        "scan": match.scan,
        "snr": float(match.hit.snr),
        "window_hz": float(match.window_hz),
        "channels": channels,
        "cells": cells,
        "expected": expected_false_alarms(float(args.snr), cells),
        "chance_log10": log10_tail(match.hit.snr)
        + float(np.log10(max(cells, 1))),
    }


def _off_budget_note(off: dict, sigma: float) -> str:
    """Why an OFF match does, or does not, settle the event."""
    verdict = (
        "but a window this wide makes peaks that strong by chance, so the"
        " rejection is not decisive"
        if off["chance_log10"] >= 0
        else "so the OFF peak is not chance noise, and pointing away did not remove it"
    )
    return "\n".join(
        [
            f"      OFF test searched {off['channels']:,} channel(s) x"
            f" {off['cells'] // off['channels']} drift trials in a"
            f" {off['window_hz'] / 1e3:.3g} kHz window ({off['cells']:,} cells)",
            f"      {off['expected']:.3g} noise peaks expected there at the"
            f" {sigma:.2f} sigma threshold; a noise peak reaching {off['snr']:.1f}"
            f" sigma has probability 10^{off['chance_log10']:.1f} per scan,"
            f" {verdict}",
        ]
    )


def _event_budget_note(budget: dict, n_scans: int) -> str:
    """The same budget the threshold was set in, applied to a whole event."""
    return "\n".join(
        [
            f"      false-alarm budget: {budget['expected_per_scan']:.3g} noise peaks"
            f" expected per scan at the {budget['sigma']:.2f} sigma threshold"
            f" ({budget['per_scan_cells']:,} cells)",
            f"      noise alone matching this event in all {n_scans} ON scans:"
            f" 10^{budget['chance_log10']:.1f}, and that ignores the frequency and"
            " drift agreement, so it overstates the chance",
        ]
    )


def _add_threshold_arguments(parser) -> None:
    """The false-alarm flag shared by `search`, `verify` and `link`."""
    parser.add_argument(
        "--false-alarm", type=float, default=None, metavar="P",
        help="set the threshold from a per-band chance P of any false alarm"
        " instead of --snr; P counts every channel-drift cell searched",
    )


def _limits_path(out: str) -> str:
    """Where the limits record goes, beside the results it belongs to."""
    if out.endswith(".csv"):
        return out[: -len(".csv")] + ".limits.csv"
    return out + ".limits.csv"


def _limit_row(args, label: str, meta: dict, sefd_jy: float) -> dict:
    """One limits-CSV row, carrying every number the limit was computed from.

    A limit nobody can recompute is a claim, so the threshold, line width,
    drift range and noise terms travel with it rather than living in a log file.
    """
    fifty, ninety = _limit_values(args, meta, sefd_jy)
    lo, hi = _band_edges(meta)
    searched = 1.0 - meta["blanked_channels"] / max(meta["nchans"], 1)
    per_channel, total_trials = _band_trials(meta, args)
    return {
        "kind": "file",
        "path": label,
        "threshold_sigma": f"{args.snr:g}",
        "line_width_channels": f"{args.line_width:g}",
        "max_drift_hz_s": f"{args.max_drift:g}",
        "freq_min_mhz": f"{lo:.6f}",
        "freq_max_mhz": f"{hi:.6f}",
        "searched_fraction": f"{searched:.4g}",
        "sefd_jy": f"{sefd_jy:g}",
        "n_ints": meta.get("n_ints", ""),
        "boxcar": meta.get("boxcar", ""),
        "df_hz": f"{abs(float(meta['df_mhz'])) * 1e6:.6g}",
        "tsamp_s": f"{float(meta['tsamp']):.6g}",
        "limit_jy": f"{fifty:.6g}",
        "limit_jy_95": f"{ninety:.6g}",
        "trials_per_channel": per_channel,
        "total_trials": total_trials,
        "false_alarm_log10": f"{log10_tail(args.snr):.4f}",
        "expected_false_alarms": f"{expected_false_alarms(args.snr, total_trials):.6g}",
        "scans": "",
        "min_scans": "",
        "set_by": "",
    }


def _write_limits(path: str, rows) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LIMIT_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _limit_values(args, meta: dict, sefd_jy: float) -> tuple[float, float]:
    """The 50% and 95% completeness flux limits for one file, in Jy."""
    common = dict(
        sensitivity=scan_sensitivity(meta),
        sefd_jy=sefd_jy,
        df_hz=abs(float(meta["df_mhz"])) * 1e6,
        tsamp_s=float(meta["tsamp"]),
        boxcar=int(meta["boxcar"]),
        line_width_channels=args.line_width,
    )
    return (
        limit_jy(threshold=args.snr, margin_sigma=0.0, **common),
        limit_jy(threshold=args.snr, margin_sigma=COMPLETENESS_MARGIN, **common),
    )


def _band_edges(meta: dict) -> tuple[float, float]:
    fch1 = float(meta.get("fch1", 0.0) or 0.0)
    foff = float(meta.get("foff", 0.0) or 0.0)
    nchans = int(meta.get("nchans", 1) or 1)
    edge = fch1 + foff * (nchans - 1)
    return (min(fch1, edge), max(fch1, edge))


def _add_flux_arguments(parser) -> None:
    """The flux-calibration flags shared by `search`, `verify` and `link`."""
    parser.add_argument(
        "--sefd", action="append", default=[], metavar="[PATH=]JY",
        help="system equivalent flux density in Jy: one value for every scan,"
        " or PATH=JY per file (repeatable)",
    )
    parser.add_argument(
        "--tsys", action="append", default=[], metavar="[PATH=]K",
        help="system temperature in K, used with --gain when --sefd is absent",
    )
    parser.add_argument(
        "--gain", action="append", default=[], metavar="[PATH=]K/JY",
        help="antenna gain in K/Jy, used with --tsys when --sefd is absent",
    )
    parser.add_argument(
        "--line-width", type=float, default=1.0, metavar="CHANNELS",
        help="assumed line width for the flux conversion (default 1, unresolved);"
        " wider lines carry no boxcar penalty",
    )


def _parse_signal(text: str) -> InjectedSignal:
    parts = [p.strip() for p in text.split(",")]
    if len(parts) not in (3, 4):
        raise argparse.ArgumentTypeError(
            "expected FREQ_MHZ,DRIFT_HZ_S,AMPLITUDE[,WIDTH_CHANNELS]"
        )
    values = [float(p) for p in parts]
    return InjectedSignal(
        freq_mhz=values[0],
        drift_hz_s=values[1],
        amplitude=values[2],
        width_channels=values[3] if len(values) == 4 else 1.5,
    )


def cmd_simulate(args) -> int:
    if args.signal:
        signals = [_parse_signal(text) for text in args.signal]
    elif args.no_signals:
        signals = []
    else:
        freqs = args.fch1 + args.foff * np.arange(args.nchans, dtype=np.float64)
        picks = [
            (args.nchans // 8, 0.5, 1.0),
            (5 * args.nchans // 8, -1.0, 0.9),
            (7 * args.nchans // 8, 0.05, 1.2),
        ]
        signals = [InjectedSignal(float(freqs[c]), drift, amp) for c, drift, amp in picks]

    rfi = ()
    if not args.no_rfi:
        if args.rfi:
            chan, amp = [float(p) for p in args.rfi.split(",")]
            rfi = (RFILine(channel=int(chan), amplitude=amp),)
        else:
            rfi = (RFILine(channel=3 * args.nchans // 8, amplitude=40.0),)

    observer = _resolve_observer(args)
    sky = _resolve_sky(args)
    if (args.ra is None) != (args.dec is None):
        raise SystemExit("give both --ra and --dec, or neither")
    if sky is not None and (observer is None or args.mjd is None):
        raise SystemExit(
            "--ra/--dec also need a site (--lat/--lon or --site) and --mjd, "
            "so the barycentric drift can be turned into what the telescope sees"
        )

    truth = simulate(
        args.out,
        nchans=args.nchans,
        nints=args.nints,
        tsamp=args.tsamp,
        fch1=args.fch1,
        foff=args.foff,
        signals=tuple(signals),
        rfi=rfi,
        glitch=not args.no_glitch,
        seed=args.seed,
        sky=None if sky is None else (sky.ra_deg, sky.dec_deg),
        observer=observer,
        mjd=args.mjd,
    )

    print(f"wrote {args.out}")
    print(
        f"  {args.nchans:,} channels @ {args.foff:+.6g} MHz,"
        f" {args.nints:,} integrations @ {args.tsamp:g} s"
    )
    print("  injected signals")
    for index, sig in enumerate(truth["signals"], 1):
        print(
            f"    {index}  {sig['freq_mhz']:.6f} MHz  {sig['drift_hz_s']:+.4f} Hz/s"
            f"  amplitude {sig['amplitude']:g}"
        )
    for line in truth["rfi"]:
        print(f"  rfi: channel {line['channel']} amplitude {line['amplitude']:g}")

    if args.truth:
        with open(args.truth, "w") as handle:
            json.dump(truth, handle, indent=2)
        print(f"  truth written to {args.truth}")

    print(f"\nnext: python -m seti_search search {args.out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="seti_search",
        description="Narrowband drift-rate (dedoppler) search for technosignatures in SIGPROC filterbank files.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    info = sub.add_parser("info", help="print a filterbank header and layout")
    info.add_argument("file")
    info.set_defaults(func=cmd_info)

    search = sub.add_parser("search", help="search a .fil file for drifting narrowband signals")
    search.add_argument("file")
    search.add_argument("-o", "--out", help="CSV output path (default: FILE.hits.csv)")
    search.add_argument("--snr", type=float, default=10.0, help="detection threshold in sigma")
    search.add_argument("--max-drift", type=float, default=4.0, help="maximum |drift| in Hz/s")
    search.add_argument(
        "--drift-step", type=float, default=None,
        help="drift trial spacing in Hz/s (default: half a boxcar of smear)",
    )
    search.add_argument("--boxcar", type=int, default=3, help="frequency smoothing, in channels")
    search.add_argument("--max-hits", type=int, default=50)
    search.add_argument("--rfi-sigma", type=float, default=6.0)
    search.add_argument(
        "--carrier-sigma", type=float, default=20.0,
        help="blank channels whose median power exceeds this many sigma (RFI)",
    )
    search.add_argument("--chan-chunk", type=int, default=1 << 18)
    search.add_argument(
        "--max-ints", type=int, default=4096,
        help="integrations to use, 0 for every integration in the file",
    )
    search.add_argument(
        "--tsamp", type=float, default=None,
        help="override tsamp in seconds, for headers that omit it",
    )
    _add_threshold_arguments(search)
    _add_flux_arguments(search)
    search.add_argument("-q", "--quiet", action="store_true")
    search.set_defaults(func=cmd_search)

    verify = sub.add_parser("verify", help="test ON candidates against an OFF scan")
    verify.add_argument("on", help="ON scan, pointed at the target")
    verify.add_argument("off", help="OFF scan, pointed away")
    verify.add_argument("-o", "--out", help="CSV output path (default: ON.verified.csv)")
    verify.add_argument("--snr", type=float, default=10.0)
    verify.add_argument("--max-drift", type=float, default=4.0)
    verify.add_argument("--drift-step", type=float, default=None)
    verify.add_argument("--boxcar", type=int, default=3)
    verify.add_argument("--max-hits", type=int, default=50)
    verify.add_argument("--rfi-sigma", type=float, default=6.0)
    verify.add_argument("--carrier-sigma", type=float, default=20.0)
    verify.add_argument("--chan-chunk", type=int, default=1 << 18)
    verify.add_argument("--max-ints", type=int, default=4096)
    verify.add_argument("--tsamp", type=float, default=None)
    verify.add_argument("--ra", type=float, default=None, help="target RA, degrees (J2000)")
    verify.add_argument("--dec", type=float, default=None, help="target Dec, degrees (J2000)")
    verify.add_argument("--lat", type=float, default=None, help="site latitude, degrees north")
    verify.add_argument("--lon", type=float, default=None, help="site longitude, degrees east")
    verify.add_argument("--site", help=f"named site: {', '.join(sorted(SITES))}")
    verify.add_argument("--on-mjd", type=float, default=None, help="ON start time (MJD)")
    verify.add_argument("--off-mjd", type=float, default=None, help="OFF start time (MJD)")
    verify.add_argument(
        "--drift-tol", type=float, default=0.05,
        help="how far an OFF hit's drift may sit from the prediction (Hz/s)",
    )
    verify.add_argument(
        "--freq-tol", type=float, default=None,
        help="frequency window in Hz (default: scales with the cadence)",
    )
    _add_threshold_arguments(verify)
    _add_flux_arguments(verify)
    verify.add_argument("-q", "--quiet", action="store_true")
    verify.set_defaults(func=cmd_verify)

    link = sub.add_parser(
        "link", help="link candidates across several ON and OFF scans of one target"
    )
    link.add_argument(
        "--on", action="append", default=[], metavar="PATH[@MJD]",
        help="an ON scan of the target; repeat for each",
    )
    link.add_argument(
        "--off", action="append", default=[], metavar="PATH[@MJD]",
        help="an OFF scan, pointed away; repeat for each",
    )
    link.add_argument("-o", "--out", help="CSV output path (default: run.linked.csv)")
    link.add_argument("--snr", type=float, default=10.0)
    link.add_argument("--max-drift", type=float, default=4.0)
    link.add_argument("--drift-step", type=float, default=None)
    link.add_argument("--boxcar", type=int, default=3)
    link.add_argument("--max-hits", type=int, default=50)
    link.add_argument("--rfi-sigma", type=float, default=6.0)
    link.add_argument("--carrier-sigma", type=float, default=20.0)
    link.add_argument("--chan-chunk", type=int, default=1 << 18)
    link.add_argument("--max-ints", type=int, default=4096)
    link.add_argument("--tsamp", type=float, default=None)
    link.add_argument("--ra", type=float, default=None, help="target RA, degrees (J2000)")
    link.add_argument("--dec", type=float, default=None, help="target Dec, degrees (J2000)")
    link.add_argument("--lat", type=float, default=None, help="site latitude, degrees north")
    link.add_argument("--lon", type=float, default=None, help="site longitude, degrees east")
    link.add_argument("--site", help=f"named site: {', '.join(sorted(SITES))}")
    link.add_argument(
        "--min-scans", type=int, default=2,
        help="how many distinct ON scans an event must appear in",
    )
    link.add_argument(
        "--drift-tol", type=float, default=0.05,
        help="how far barycentric drifts may differ and still link (Hz/s)",
    )
    link.add_argument(
        "--freq-tol", type=float, default=None,
        help="frequency window in Hz (default: grows with the run length)",
    )
    _add_threshold_arguments(link)
    _add_flux_arguments(link)
    link.add_argument("-q", "--quiet", action="store_true")
    link.set_defaults(func=cmd_link)

    sim = sub.add_parser("simulate", help="write a filterbank file with known signals")
    sim.add_argument("out")
    sim.add_argument("--nchans", type=int, default=8192)
    sim.add_argument("--nints", type=int, default=200)
    sim.add_argument("--tsamp", type=float, default=1.0)
    sim.add_argument("--fch1", type=float, default=1420.4)
    sim.add_argument("--foff", type=float, default=-2.86e-6)
    sim.add_argument("--signal", action="append", help="FREQ_MHZ,DRIFT_HZ_S,AMP[,WIDTH]")
    sim.add_argument("--rfi", help="CHANNEL,AMPLITUDE for a stationary carrier")
    sim.add_argument("--no-rfi", action="store_true")
    sim.add_argument(
        "--no-signals", action="store_true", help="write noise only, e.g. an OFF scan"
    )
    sim.add_argument(
        "--mjd", type=float, default=None, help="start time (MJD), recorded in the header"
    )
    sim.add_argument(
        "--ra", type=float, default=None,
        help="target RA in degrees; --signal drifts then mean barycentric",
    )
    sim.add_argument("--dec", type=float, default=None, help="target Dec in degrees")
    sim.add_argument("--lat", type=float, default=None, help="site latitude, degrees north")
    sim.add_argument("--lon", type=float, default=None, help="site longitude, degrees east")
    sim.add_argument("--site", help=f"named site: {', '.join(sorted(SITES))}")
    sim.add_argument("--no-glitch", action="store_true")
    sim.add_argument("--seed", type=int, default=0)
    sim.add_argument("--truth", help="write the injected truth to this JSON file")
    sim.set_defaults(func=cmd_simulate)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
