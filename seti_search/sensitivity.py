"""What the neutrino lane could detect, stated as numbers instead of hope.

A lane that reports ``UNINFORMATIVE`` has to say what it would take to be
informative, or the report is an excuse. Everything here answers one of three
questions:

* **Significance.** How often does a random event land in the box? That is the
  background, and it grows with the solid angle and with the time window.
* **Sensitivity.** How much flux does it take to put one event in the box? That
  falls with effective area and with duration, so a longer window is *better* for
  a flux limit and *worse* for a coincidence.
* **Reach.** Given a source that emits a total energy in neutrinos, how far away
  could this lane still see it? That is the number worth quoting to an observer,
  because it converts the instrument into a distance.

The two uses pull in opposite directions and the code keeps them apart rather
than blending them into a single "sensitivity": ``background_events`` for the
first, ``fluence_limit_cm2`` for the second, ``horizon_mpc`` for the third.
Every function takes its assumptions as arguments, because the honest output of
this module is a curve with the assumptions attached, not a single number.

Units are fixed and checked: metres and seconds for the instrument, cm^-2 for
fluence, erg for energy, Mpc for distance. A mixed-unit answer here would be
wrong by 10^4 and look plausible, which is the failure this module is written to
avoid.
"""

from __future__ import annotations

import math

from .neutrino import solid_angle_sr

__all__ = [
    "background_events",
    "chance_probability",
    "energy_fluence_limit_erg_cm2",
    "fluence_limit_cm2",
    "horizon_mpc",
    "power",
    "required_events",
    "signal_events",
    "window_for_chance_probability",
]

SECONDS_PER_YEAR = 3.1557e7
CM2_PER_M2 = 1.0e4
ERG_PER_GEV = 1.602176634e-3
CM_PER_MPC = 3.0856775814913673e24


def background_events(rate_per_year: float, half_angle_deg: float,
                      window_s: float) -> float:
    """Expected accidental events in a cone of this size and a window this long.

    The rate is the detector's whole-sky alert rate, so the fraction of the sky
    inside the cone is ``Omega / 4 pi``. A rate of zero means a catalog with no
    background model, which is not the same as a background of zero - callers
    pass ``None`` upstream when the rate is unknown.
    """
    if rate_per_year < 0:
        raise ValueError("a rate cannot be negative")
    if window_s < 0:
        raise ValueError("a window cannot be negative")
    omega = solid_angle_sr(half_angle_deg)
    return rate_per_year * (omega / (4.0 * math.pi)) * (window_s / SECONDS_PER_YEAR)


def chance_probability(mu_background: float) -> float:
    """``P(>= 1) = 1 - exp(-mu)`` for a Poisson background."""
    if mu_background < 0:
        raise ValueError("mu cannot be negative")
    return 1.0 - math.exp(-mu_background)


def window_for_chance_probability(rate_per_year: float, half_angle_deg: float,
                                  p_max: float) -> float:
    """Longest window whose accidental-coincidence probability stays under ``p_max``.

    This is the timing precision the radio side has to supply: past this window
    the box is crowded with background and a coincidence means nothing.
    """
    if not (0.0 < p_max < 1.0):
        raise ValueError("p_max is a probability in (0, 1)")
    if rate_per_year <= 0:
        raise ValueError("a rate limit needs a positive rate")
    omega = solid_angle_sr(half_angle_deg)
    mu_max = -math.log(1.0 - p_max)
    return mu_max * 4.0 * math.pi * SECONDS_PER_YEAR / (rate_per_year * omega)


def required_events(detection_power: float) -> float:
    """Events needed for a given probability of seeing at least one."""
    if not (0.0 < detection_power < 1.0):
        raise ValueError("power is a probability in (0, 1)")
    return -math.log(1.0 - detection_power)


def power(mu_signal: float) -> float:
    """Probability of at least one signal event, given the expected count."""
    if mu_signal < 0:
        raise ValueError("mu cannot be negative")
    return 1.0 - math.exp(-mu_signal)


def signal_events(flux_cm2_s: float, effective_area_m2: float,
                  duration_s: float) -> float:
    """Expected detected neutrinos from a steady flux.

    ``N = flux * area * time``, with the area in m^2 converted to cm^2 once. A
    flux quoted per m^2 would be 10^4 times smaller and would look like a better
    limit, so the conversion happens here rather than at the call sites.
    """
    if flux_cm2_s < 0:
        raise ValueError("a flux cannot be negative")
    if effective_area_m2 <= 0 or duration_s <= 0:
        raise ValueError("area and duration must be positive")
    return flux_cm2_s * effective_area_m2 * CM2_PER_M2 * duration_s


def fluence_limit_cm2(detection_power: float, effective_area_m2: float,
                      duration_s: float) -> float:
    """Number fluence (neutrinos per cm^2) for the stated power, zero background.

    The background is set to zero on purpose: this is the sensitivity question,
    and adding a background it does not depend on would make the number look more
    conservative without being so.
    """
    if effective_area_m2 <= 0 or duration_s <= 0:
        # a zero area divides to infinity and looks like a perfect instrument
        raise ValueError("area and duration must be positive")
    needed = required_events(detection_power)
    return needed / (effective_area_m2 * CM2_PER_M2 * duration_s)


def energy_fluence_limit_erg_cm2(detection_power: float,
                                 effective_area_m2: float, duration_s: float,
                                 pivot_gev: float) -> float:
    """Energy fluence, treating each detected neutrino as the pivot energy.

    That is a convention, not a measurement: a real E^-2 source deposits a
    distribution of energies, and the effective area is a function of energy, so
    the number is an order-of-magnitude statement and is labelled as one.
    """
    if pivot_gev <= 0:
        raise ValueError("pivot energy must be positive")
    return (fluence_limit_cm2(detection_power, effective_area_m2, duration_s)
            * pivot_gev * ERG_PER_GEV)


def horizon_mpc(energy_erg: float, energy_fluence_limit_erg_cm2: float) -> float:
    """Distance at which a source emitting ``energy_erg`` sits at the limit.

    Isotropic emission: ``F = E / (4 pi d^2)``, so
    ``d = sqrt(E / (4 pi F))``. A jet beamed at us would reach further by the
    beaming factor, which is not included, so this is the conservative end.
    """
    if energy_erg <= 0 or energy_fluence_limit_erg_cm2 <= 0:
        raise ValueError("energy and fluence limit must be positive")
    d_cm = math.sqrt(energy_erg / (4.0 * math.pi * energy_fluence_limit_erg_cm2))
    return d_cm / CM_PER_MPC
