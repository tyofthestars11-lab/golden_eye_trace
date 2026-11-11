"""Execute the phi read on real candidates, with its own null attached.

The directive names three coordinates to read: the rung address
r(x) = ln(x)/ln(phi), the golden pitch mu = 17.0324 deg, and the golden angle
137.5078 deg. This script computes all three for every candidate in a real
file, and then computes the same three for noise draws and for injected
signals that contain no phi at all.

A read that separates a broadcast from a natural event has to return a
different answer for those groups. Printing the group pass rates next to the
candidates' values is the whole point: if every group passes at the same
rate, the read is a description of the coordinates rather than a test.

    python3 scripts/phi_read.py --fil data/stare2_ovro.fil
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import numpy as np  # noqa: E402

from phi_native import FOFF, MU_DEG, SLOPE, TSAMP_S  # noqa: E402
from seti_search.filterbank import (  # noqa: E402
    Filterbank,
    InjectedSignal,
    simulate,
)
from seti_search.rungs import LN_PHI, PHI, rung  # noqa: E402
from seti_search.search import SearchConfig, search_file  # noqa: E402

GOLDEN_ANGLE_DEG = 360.0 / (PHI * PHI)
TOL = 0.05


def _phases(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return (np.log(values) / LN_PHI) % 1.0


def _dist_to_rung(phases: np.ndarray) -> np.ndarray:
    phases = np.asarray(phases, dtype=float)
    return np.abs(phases - np.round(phases))


def read_candidates(path: str):
    print(f"READ 0. THE CANDIDATES IN {os.path.basename(path)}")
    cfg = SearchConfig(boxcar=3, max_drift_hz_s=4.0, snr_threshold=3.0, max_hits=200)
    hits, _meta = search_file(path, cfg, chan_chunk=2048, max_ints=None)
    fb = Filterbank.open(path)
    freqs = fb.frequencies()
    low, high = float(freqs.min()), float(freqs.max())
    print(f"   {len(hits)} hits at {cfg.snr_threshold:g} sigma in "
          f"{low:.3f}-{high:.3f} MHz")
    print(f"   {fb.nchans} channels, {fb.nints} rows, {fb.total_time:g} s total, "
          f"channel {abs(fb.foff) * 1e3:.4f} kHz")
    print("   ranked by SNR, which is the only ordering the search has:")
    print(f"   {'MHz':>14}{'SNR':>8}{'drift Hz/s':>13}")
    for hit in sorted(hits, key=lambda h: -h.snr)[:8]:
        drift = getattr(hit, "drift_hz_s", 0.0)
        print(f"   {hit.freq_mhz:>14.6f}{hit.snr:>8.1f}{drift:>13.3g}")
    print()
    return hits, fb, low, high


def read_rung(hits, fb, low: float, high: float) -> None:
    print("READ 1. THE RUNG ADDRESS  r(x) = ln(x) / ln(phi)")
    freqs = np.array([h.freq_mhz for h in hits], dtype=float)
    phases = _phases(freqs)
    dist = _dist_to_rung(phases)
    print(f"   one rung is {low * (PHI - 1.0):.1f} MHz wide at this "
          f"frequency; the band spans {rung(high) - rung(low):.4f} rungs")
    print(f"   {'MHz':>14}{'rung':>10}{'phase':>9}{'dist':>8}")
    for f_mhz, ph, d in zip(freqs, phases, dist):
        print(f"   {f_mhz:>14.6f}{rung(float(f_mhz)):>10.4f}{ph:>9.4f}{d:>8.4f}")

    rng = np.random.default_rng(23)
    draws = rng.uniform(low, high, 40000)
    chance = float(np.mean(_dist_to_rung(_phases(draws)) <= TOL))
    inside = int(np.sum(dist <= TOL))
    print(f"   {inside} of {len(hits)} are within {TOL} of an integer rung")
    print(f"   40000 uniform frequencies in this same band: {chance:.3f} pass")
    print("   so the read returns that value for noise at the same rate. With")
    print(f"   {rung(high) - rung(low):.3f} rungs spread across the band every "
          "channel shares one address.")

    # relabelling test: the same physics, expressed in Hz instead of MHz
    shift = math.log(1e6) / LN_PHI
    shifted = (_phases(draws) + shift) % 1.0
    chance_hz = float(np.mean(_dist_to_rung(shifted) <= TOL))
    flips = int(np.sum((_dist_to_rung(phases) <= TOL)
                       != (_dist_to_rung((phases + shift) % 1.0) <= TOL)))
    print(f"   writing the same frequencies in Hz shifts every rung by "
          f"{shift:.4f}")
    print(f"   and moves the chance rate to {chance_hz:.3f}; {flips} of "
          f"{len(hits)} candidates change verdict.")
    print()


def read_pitch(hits, fb) -> None:
    print("READ 2. THE GOLDEN PITCH: WHICH BASELINE MAKES A DRIFT THE PITCH?")
    delta_hz = abs(fb.foff) * 1e6
    resolution = delta_hz / fb.total_time
    print(f"   tan(mu) = {SLOPE:.6f} at mu = {MU_DEG:.4f} deg, so the pitch is")
    print("   drift = tan(mu) * f / T, and any drift is the pitch at some T:")
    print(f"   {'MHz':>14}{'drift Hz/s':>13}{'T* (s)':>14}{'T* (yr)':>12}")
    baselines = []
    for hit in sorted(hits, key=lambda h: -h.snr)[:8]:
        f_hz = hit.freq_mhz * 1e6
        drift = float(getattr(hit, "drift_hz_s", 0.0))
        if drift == 0.0:
            print(f"   {hit.freq_mhz:>14.6f}{drift:>13.3f}{'-':>14}{'inf':>12}")
            continue
        t_s = SLOPE * f_hz / abs(drift)
        baselines.append(t_s)
        print(f"   {hit.freq_mhz:>14.6f}{drift:>13.3f}{t_s:>14.3e}"
              f"{t_s / 3.156e7:>12.1f}")
    if len(baselines) > 1:
        print(f"   the baselines these candidates require span a factor of "
              f"{max(baselines) / min(baselines):.1f};")
        print("   one pitch cannot put them all on the lattice at once.")
    elif not baselines:
        print("   every peak sits at zero drift: with one drift trial this file")
        print("   cannot place any candidate on a pitch baseline at all.")
    print(f"   the drift resolution of this file is df / T = {resolution:.3e} Hz/s,")
    print("   against drifts of order 0.1 Hz/s, so the drifts are not measured:")
    print("   every T* above is unconstrained, and the read names a baseline")
    print("   rather than a source.")
    window_t = SLOPE * 1.4206e9 / 4.0
    print(f"   the whole +/-4 Hz/s window is the pitch at baselines of "
          f"{window_t / 3.156e7:.2f} yr and longer,")
    print("   so 100% of the drifts this search can see are pitch-consistent.")
    print()


def _galactic_lon(ra_deg: float, dec_deg: float) -> float:
    a_ngp, d_ngp, l_ncp = 192.85948, 27.12825, 122.93192
    a, d = math.radians(ra_deg), math.radians(dec_deg)
    an, dn = math.radians(a_ngp), math.radians(d_ngp)
    y = math.cos(d) * math.sin(a - an)
    x = (math.sin(d) * math.cos(dn)
         - math.cos(d) * math.sin(dn) * math.cos(a - an))
    return (l_ncp - math.degrees(math.atan2(y, x))) % 360.0


def _ecliptic_lon(ra_deg: float, dec_deg: float, obl_deg: float = 23.4392911) -> float:
    a, d = math.radians(ra_deg), math.radians(dec_deg)
    e = math.radians(obl_deg)
    lam = math.atan2(math.sin(a) * math.cos(e) + math.tan(d) * math.sin(e),
                     math.cos(a))
    return math.degrees(lam) % 360.0


def _offset(lon_deg: float) -> float:
    return (lon_deg - GOLDEN_ANGLE_DEG + 180.0) % 360.0 - 180.0


def _decode_raj(value: float):
    if not value:
        return None
    if value > 24.0001:
        hh = int(value / 10000.0)
        mm = int(value / 100.0) % 100
        ss = value - 100.0 * int(value / 100.0)
        return 15.0 * (hh + mm / 60.0 + ss / 3600.0)
    return 15.0 * float(value)


def _decode_dej(value: float):
    if not value:
        return None
    if abs(value) > 90.0001:
        sign = -1.0 if value < 0 else 1.0
        v = abs(value)
        dd = int(v / 10000.0)
        mm = int(v / 100.0) % 100
        ss = v - 100.0 * int(v / 100.0)
        return sign * (dd + mm / 60.0 + ss / 3600.0)
    return float(value)


def read_angle(fb, ra_deg, dec_deg) -> None:
    print("READ 3. THE GOLDEN ANGLE ADDRESS  theta = 137.5078 deg")
    raw_ra = fb.header.get("src_raj", 0.0)
    raw_dej = fb.header.get("src_dej", 0.0)
    if ra_deg is None:
        ra_deg = _decode_raj(float(raw_ra or 0.0))
        dec_deg = _decode_dej(float(raw_dej or 0.0))
    print(f"   header carries src_raj={raw_ra} src_dej={raw_dej}"
          f" (source {fb.source})")
    demonstration = False
    if ra_deg is None or dec_deg is None:
        ra_deg, dec_deg = 266.4168, -29.0078
        demonstration = True
        print("   no position is recorded for this candidate, so the frame test is")
        print("   run on a reference direction instead (Sgr A*, 266.4168 -29.0078):")
    rows = [("equatorial RA", ra_deg),
            ("galactic l", _galactic_lon(ra_deg, dec_deg)),
            ("ecliptic lambda", _ecliptic_lon(ra_deg, dec_deg))]
    if not demonstration:
        print(f"   position {ra_deg:.4f} deg, {dec_deg:.4f} deg read in three frames:")
    else:
        print(f"   reference {ra_deg:.4f} deg, {dec_deg:.4f} deg read in three frames:")
    print(f"   {'frame':>16}{'longitude':>12}{'offset':>10}")
    for name, lon in rows:
        print(f"   {name:>16}{lon:>12.4f}{_offset(lon):>10.4f}")
    offs = [_offset(lon) for _, lon in rows]
    print(f"   the offsets differ by up to "
          f"{max(abs((a - b + 180.0) % 360.0 - 180.0) for a in offs for b in offs):.1f} deg")
    print("   for one fixed direction, so membership depends on the frame.")
    print(f"   chance of landing within 5 deg of {GOLDEN_ANGLE_DEG:.1f} in one frame")
    print(f"   is {10.0 / 360.0:.3f}; with three frames it is "
          f"{1.0 - (1.0 - 10.0 / 360.0) ** 3:.3f}. Rotating the frame by any angle")
    print("   moves the offset by that angle, so the criterion can be met by")
    print("   choosing coordinates rather than by choosing a source.")
    print()


def read_injections(count: int = 24) -> None:
    print("READ 4. THE SAME READ ON SIGNALS WITH NO PHI IN THEM")
    rng = np.random.default_rng(11)
    centres = rng.uniform(1300.0, 2300.0, count)
    passed = 0
    found = 0
    near = 0
    t_stars: list[float] = []
    with tempfile.TemporaryDirectory() as folder:
        for i, f_mhz in enumerate(centres):
            path = os.path.join(folder, f"inj{i:02d}.fil")
            simulate(path, nchans=1024, nints=48, tsamp=TSAMP_S,
                     fch1=float(f_mhz) + 0.0015, foff=FOFF,
                     signals=(InjectedSignal(float(f_mhz), drift_hz_s=0.05,
                                             amplitude=12.0, width_channels=0.4),),
                     rfi=(), glitch=False, seed=500 + i)
            hit = None
            for boxcar in (1, 3):
                cfg = SearchConfig(boxcar=boxcar, max_drift_hz_s=4.0,
                                   snr_threshold=3.0, max_hits=20)
                hits, _meta = search_file(path, cfg, chan_chunk=1024, max_ints=None)
                if hits:
                    hit = max(hits, key=lambda h: h.snr)
                    break
            if hit is None:
                continue
            found += 1
            if abs(float(hit.freq_mhz) - float(f_mhz)) < 0.01:
                near += 1
            phase = float(rung(float(hit.freq_mhz)) % 1.0)
            if abs(phase - round(phase)) <= TOL:
                passed += 1
            drift = abs(float(getattr(hit, "drift_hz_s", 0.0)))
            t_stars.append(SLOPE * float(hit.freq_mhz) * 1e6 / drift
                           if drift else math.inf)
    expected = found * 2.0 * TOL
    spread = math.sqrt(max(found, 1) * 2.0 * TOL * (1.0 - 2.0 * TOL))
    finite = sum(1 for t in t_stars if math.isfinite(t))
    width_hz = 1024 * abs(FOFF) * 1e6
    print(f"   {found} of {count} plain narrowband beacons were recovered, "
          f"{near} of them")
    print("   at the injected frequency; none was constructed with any phi in it.")
    print(f"   {passed} land within {TOL} of a rung, expected {expected:.1f} "
          f"+/- {spread:.1f} from a coin.")
    print(f"   each file is {width_hz:.0f} Hz wide, or "
          f"{width_hz / 1420.4e6 / LN_PHI:.1e} rungs, so the centres from 1300")
    print("   to 2300 MHz are what varies the address, and those phases are")
    print("   uniform by construction.")
    print(f"   {finite} of {found} have a finite pitch baseline, spanning seconds")
    print("   to millennia, so the pitch accepts every one of them. A read that")
    print("   accepts a signal with no phi in it is measuring the ruler, not the")
    print("   source.")
    print()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fil", required=True, help="filterbank file to read")
    parser.add_argument("--ra", type=float, default=None,
                        help="right ascension in degrees (overrides the header)")
    parser.add_argument("--dec", type=float, default=None,
                        help="declination in degrees (overrides the header)")
    parser.add_argument("--injections", type=int, default=24,
                        help="how many phi-free beacons to read")
    args = parser.parse_args()

    hits, fb, low, high = read_candidates(args.fil)
    if hits:
        read_rung(hits, fb, low, high)
        read_pitch(hits, fb)
    else:
        print("READ 1. no candidates above threshold, nothing to read")
        print()
    read_angle(fb, args.ra, args.dec)
    read_injections(args.injections)
    print("SUMMARY")
    print("   rung: every frequency in the band shares one address, so the read")
    print("         returns the noise value to a signal and to noise alike")
    print("   pitch: every drift is the pitch at some baseline; the read names a")
    print("          baseline, and this file's drifts are unmeasured anyway")
    print("   angle: one position gives three different offsets in three frames")
    print("   beacon: phi-free injections pass at the coin's rate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
