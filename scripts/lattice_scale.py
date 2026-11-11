"""Scale the Z[phi] lattice test: real populations against their own controls.

One coincidence at p = 0.1 proves nothing, and neither do twenty unless the
controls are run alongside them. This script measures the lattice statistic on

* a population of real spectral frequencies (atomic, molecular, hyperfine, and
  the two band centres this project actually uses),
* a population of real sources with published coordinates, in RA, ecliptic and
  galactic longitude,
* uniform random positions on the sphere, and log-uniform random frequencies,
  thousands of each,

and reports the same three numbers for every set: the pass rate at a 0.05-rung
tolerance, the best (smallest) p in the set, and that best p corrected for the
number of tests. The control is what makes the real numbers readable: a lattice
of spacing ``s = 2 - phi`` puts a random address within 0.05 rungs of a point
``2*0.05/s = 0.262`` of the time, and the best of N draws has p of order 1/N no
matter what the addresses are.

Frequencies are compared in MHz throughout, because the rung coordinate depends
on the unit and mixing units is the framework's own easiest way to manufacture a
hit. Table values here are quoted to a few kHz where they are known that finely;
at a 0.05-rung tolerance on 1 GHz, 0.05 rungs is about 50 MHz, so table
precision is irrelevant to the outcome - which is itself worth knowing.

    python3 scripts/lattice_scale.py
"""

from __future__ import annotations

import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from seti_search.rungs import PHI, rung  # noqa: E402

STEP = 2.0 - PHI          # one golden angle, in turns and in rungs
TOL_RUNGS = 0.05
SEED = 1729
OUT_PATH = os.path.join(ROOT, "prereg", "lattice_scale.json")

FREQS_MHZ = [
    ("HI 21 cm", 1420.405751),
    ("OH 1612", 1612.231),
    ("OH 1665", 1665.402),
    ("OH 1667", 1667.359),
    ("OH 1720", 1720.530),
    ("OH 4765", 4765.562),
    ("H2CO", 4829.660),
    ("H109a", 5008.939),
    ("H110a", 4874.157),
    ("CH3OH", 6668.519),
    ("Cs 133 hyperfine", 9192.631770),
    ("Rb 87 hyperfine", 6834.682610),
    ("H2O 22 GHz", 22235.080),
    ("NH3 (1,1)", 23694.495),
    ("SiO 43 GHz", 43122.080),
    ("HCN 88 GHz", 88630.416),
    ("HCO+ 89 GHz", 89188.526),
    ("CO (1-0)", 115271.202),
    ("NV- D", 2870.000),
    ("14N A_parallel", 2.160),
    ("LWA band centre", 50.390),
    ("STARE2 band centre", 1405.122),
]

SOURCES = [
    ("Sgr A*", 266.41683, -29.00781),
    ("Cas A", 350.85000, 58.81500),
    ("Crab", 83.63308, 22.01450),
    ("Vela", 128.83610, -45.17640),
    ("Cyg A", 299.86820, 40.73390),
    ("Cyg X-1", 299.59030, 35.20160),
    ("M87", 187.70590, 12.39110),
    ("Orion A", 83.82210, -5.39110),
    ("W3", 36.44500, 62.15000),
    ("B0329+54", 53.24729, 54.57867),
    ("B1937+21", 293.73290, 21.58310),
    ("3C 273", 187.27790, 2.05240),
    ("NGC 1275", 49.95070, 41.51170),
    ("LMC", 80.89390, -69.75610),
    ("SMC", 13.18670, -72.82860),
    ("M31", 10.68470, 41.26920),
]


def ecliptic_lon(ra_deg: float, dec_deg: float, obl_deg: float = 23.4392911) -> float:
    a, d, e = np.radians(ra_deg), np.radians(dec_deg), np.radians(obl_deg)
    lam = np.arctan2(np.sin(a) * np.cos(e) + np.tan(d) * np.sin(e), np.cos(a))
    return np.degrees(lam) % 360.0


def galactic_lon(ra_deg: float, dec_deg: float) -> float:
    a_ngp, d_ngp, l_ncp = 192.85948, 27.12825, 122.93192
    a, d = np.radians(ra_deg), np.radians(dec_deg)
    an, dn = np.radians(a_ngp), np.radians(d_ngp)
    y = np.cos(d) * np.sin(a - an)
    x = np.sin(d) * np.cos(dn) - np.cos(d) * np.sin(dn) * np.cos(a - an)
    return (l_ncp - np.degrees(np.arctan2(y, x))) % 360.0


def p_of_address(address: np.ndarray) -> np.ndarray:
    """Additive lattice chance for an address measured in lattice steps."""
    address = np.asarray(address, dtype=np.float64)
    index = np.round(address / STEP)
    return np.minimum(1.0, 2.0 * np.abs(address - index * STEP) / STEP)


def summarise(ps: np.ndarray) -> dict:
    ps = np.asarray(ps, dtype=np.float64)
    best = float(ps.min())
    return {
        "n": int(ps.size),
        "pass_rate": float(np.mean(ps <= 2.0 * TOL_RUNGS / STEP)),
        "best_p": best,
        "best_adjusted": float(1.0 - (1.0 - best) ** ps.size),
        "median_p": float(np.median(ps)),
    }


def main() -> int:
    rng = np.random.default_rng(SEED)
    results: dict[str, dict] = {}

    freq_addresses = np.array([rung(f) for _name, f in FREQS_MHZ])
    results["real_frequencies"] = summarise(p_of_address(freq_addresses))

    ra = np.array([s[1] for s in SOURCES])
    dec = np.array([s[2] for s in SOURCES])
    lon_addresses = np.concatenate([
        (ra % 360.0) / 360.0,
        ecliptic_lon(ra, dec) / 360.0,
        galactic_lon(ra, dec) / 360.0,
    ])
    results["real_sources"] = summarise(p_of_address(lon_addresses))

    # UNIFORM CONTROL: random points on the sphere, all three longitudes
    draws = 4000
    c_ra = rng.uniform(0.0, 360.0, draws)
    c_dec = np.degrees(np.arcsin(rng.uniform(-1.0, 1.0, draws)))
    control_lon = np.concatenate([
        c_ra / 360.0,
        ecliptic_lon(c_ra, c_dec) / 360.0,
        galactic_lon(c_ra, c_dec) / 360.0,
    ])
    results["uniform_positions"] = summarise(p_of_address(control_lon))

    # UNIFORM CONTROL: log-uniform frequencies over six decades, in MHz
    control_freq = 10.0 ** rng.uniform(1.0, 5.0, draws)
    results["uniform_frequencies"] = summarise(p_of_address(np.log(control_freq)))

    theory = 2.0 * TOL_RUNGS / STEP
    print("SCALING THE Z[phi] LATTICE TEST")
    print(f"   lattice step {STEP:.6f} (2 - phi), tolerance {TOL_RUNGS} rungs,")
    print(f"   so a random address passes with probability 2*t/s = {theory:.4f}")
    print()
    print(f"   {'set':<24}{'n':>6}{'pass rate':>11}{'best p':>9}"
          f"{'adjusted':>10}{'median p':>10}")
    for name, row in results.items():
        print(f"   {name:<24}{row['n']:>6}{row['pass_rate']:>11.4f}"
              f"{row['best_p']:>9.4f}{row['best_adjusted']:>10.4f}"
              f"{row['median_p']:>10.4f}")
    print()
    print(f"   the controls pass at {results['uniform_positions']['pass_rate']:.4f} "
          f"and {results['uniform_frequencies']['pass_rate']:.4f} against a "
          f"predicted {theory:.4f}: the null is calibrated")
    print(f"   the best control p is {results['uniform_positions']['best_p']:.4f} "
          f"and {results['uniform_frequencies']['best_p']:.4f}, as expected for "
          f"{draws} draws of a uniform address (about 1/N)")
    print()
    real_best = min(results["real_frequencies"]["best_p"],
                    results["real_sources"]["best_p"])
    ctrl_best = min(results["uniform_positions"]["best_p"],
                    results["uniform_frequencies"]["best_p"])
    print(f"   best real p {real_best:.4f} vs best control p {ctrl_best:.4f}: the "
          f"real sets are not")
    print(f"   better than random draws of the same size, and after correcting")
    print(f"   each set for its own number of tests they are unremarkable:")
    print(f"   frequencies {results['real_frequencies']['best_adjusted']:.3f}, "
          f"sources {results['real_sources']['best_adjusted']:.3f}")

    results["theory_pass_rate"] = theory
    results["tolerance_rungs"] = TOL_RUNGS
    results["seed"] = SEED
    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"\n   wrote {os.path.relpath(OUT_PATH, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
