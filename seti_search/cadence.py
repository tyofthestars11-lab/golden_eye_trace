"""Where a fixed sky position ought to drift.

A narrowband signal from a fixed point on the sky has a *topocentric* drift set
by the observer's motion: Earth's orbit around the Sun (about 30 km/s), Earth's
rotation (up to 465 m/s at the equator) and any acceleration the transmitter
itself has. Between an ON scan and an OFF scan taken minutes apart, the rotation
term changes, so the same source shows a slightly different drift in each file.
Predicting that difference is what lets an ON/OFF pair test a candidate instead
of only counting it: drift is the one quantity a celestial signal cannot fake.

Model
-----
For a source at unit vector ``n`` the observed frequency is

    nu_obs(t) = nu_bary(t) * f(t),      f(t) = 1 - v(t).n / c

to first order, with ``v`` the observer's barycentric velocity. Differentiating,

    drift_obs(t) = f(t) * drift_bary - (nu_obs(t) / f(t)) * (a(t).n) / c

with ``a = dv/dt``. Given a measured topocentric drift at one epoch this inverts
to the barycentric drift, which is what stays constant for a source, and from
there predicts the topocentric drift at any other epoch.

Accuracy
--------
The solar position uses the standard low-precision series (about 0.01 degrees),
which is ample here: a drift prediction needs only the *change* in projected
acceleration, and the series' own error changes very little over minutes. The
lunar term is not modelled; it contributes up to about 12 m/s of velocity and
3e-5 m/s^2 of acceleration, i.e. roughly 1e-4 Hz/s of drift at 1.4 GHz, far
below any practical drift trial spacing. Second-order Doppler terms are ~1e-8
and ignored.

Absolute velocity errors of order 10 m/s matter for comparing *absolute*
frequencies across epochs only through their variation, which is tiny over a
cadence; the frequency ratio used here inherits that, not the absolute error.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "Observer",
    "SkyPosition",
    "SITES",
    "barycentric_acceleration",
    "barycentric_drift",
    "barycentric_velocity",
    "doppler_factor",
    "gmst_deg",
    "predict_observation",
    "topocentric_drift",
    "topocentric_frequency",
]

C_LIGHT = 299792458.0
AU_METRES = 1.495978707e11
EARTH_RADIUS_M = 6.371008e6
SECONDS_PER_DAY = 86400.0


@dataclass(frozen=True)
class Observer:
    """A telescope site. Longitude is east-positive, as in most maps."""

    latitude_deg: float
    longitude_deg: float
    elevation_m: float = 0.0

    def __str__(self) -> str:
        return f"{self.latitude_deg:+.4f}, {self.longitude_deg:+.4f}"


# A few sites that appear in public narrowband data, so a file can be checked
# without looking up coordinates.
SITES = {
    "gbt": Observer(38.4331, -79.8398, 806.0),
    "parkes": Observer(-32.9983, 148.2635, 395.0),
    "fast": Observer(25.6529, 106.8564, 1110.0),
    "meerkat": Observer(-30.7215, 21.4103, 1038.0),
    "arecibo": Observer(18.3442, -66.7528, 497.0),
    "ovro": Observer(37.2340, -118.2820, 1222.0),
    "greenbank": Observer(38.4331, -79.8398, 806.0),
    "effelsberg": Observer(50.5247, 6.8836, 369.0),
    "nancay": Observer(47.3803, 2.1972, 150.0),
}


@dataclass(frozen=True)
class SkyPosition:
    """Equatorial coordinates of the target, J2000."""

    ra_deg: float
    dec_deg: float

    def unit_vector(self) -> np.ndarray:
        ra = np.radians(self.ra_deg)
        dec = np.radians(self.dec_deg)
        return np.array(
            [np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]
        )


def gmst_deg(mjd: float) -> float:
    """Greenwich mean sidereal time in degrees (IAU 1982 series).

    Good to well under an arcsecond over the era of modern radio data; the
    sub-arcsecond terms are irrelevant for a drift of order 0.01 Hz/s.
    """
    jd = float(mjd) + 2400000.5
    days = jd - 2451545.0
    centuries = days / 36525.0
    theta = (
        280.46061837
        + 360.98564736629 * days
        + 0.000387933 * centuries**2
        - centuries**3 / 38710000.0
    )
    return float(theta % 360.0)


def _sun_geocentric(mjd: float) -> np.ndarray:
    """Apparent geocentric position of the Sun in AU, equatorial frame."""
    days = float(mjd) - 51544.5  # from J2000.0
    mean_longitude = np.radians(280.460 + 0.9856474 * days)
    mean_anomaly = np.radians(357.528 + 0.9856003 * days)
    ecliptic_longitude = (
        mean_longitude
        + np.radians(1.915) * np.sin(mean_anomaly)
        + np.radians(0.020) * np.sin(2.0 * mean_anomaly)
    )
    radius = 1.00014 - 0.01671 * np.cos(mean_anomaly) - 0.00014 * np.cos(
        2.0 * mean_anomaly
    )
    obliquity = np.radians(23.4392911 - 3.563e-7 * days)
    return radius * np.array(
        [
            np.cos(ecliptic_longitude),
            np.cos(obliquity) * np.sin(ecliptic_longitude),
            np.sin(obliquity) * np.sin(ecliptic_longitude),
        ]
    )


def _earth_rotation_position(observer: Observer, mjd: float) -> np.ndarray:
    """Position of the site relative to Earth's centre, equatorial frame (m)."""
    theta = np.radians(gmst_deg(mjd) + observer.longitude_deg)
    latitude = np.radians(observer.latitude_deg)
    radius = EARTH_RADIUS_M + observer.elevation_m
    return radius * np.array(
        [
            np.cos(latitude) * np.cos(theta),
            np.cos(latitude) * np.sin(theta),
            np.sin(latitude),
        ]
    )


def _barycentric_position(observer: Observer, mjd: float) -> np.ndarray:
    """Observer's position relative to the solar system barycentre, metres.

    Only the Sun's pull dominates at the accuracy this module claims; the
    barycentre itself is within about 0.005 AU of the Sun, which changes the
    orbital velocity by well under a percent and the acceleration by less.
    """
    return -_sun_geocentric(mjd) * AU_METRES + _earth_rotation_position(observer, mjd)


def barycentric_velocity(observer: Observer, mjd: float, dt_s: float = 1.0) -> np.ndarray:
    """Observer's barycentric velocity in m/s, by central difference."""
    step = dt_s / SECONDS_PER_DAY
    return (
        _barycentric_position(observer, mjd + step)
        - _barycentric_position(observer, mjd - step)
    ) / (2.0 * dt_s)


def barycentric_acceleration(observer: Observer, mjd: float, dt_s: float = 300.0):
    """Observer's barycentric acceleration in m/s^2, by second difference.

    A 300 s baseline keeps the truncation error around 1e-10 m/s^2 while the
    double-precision roundoff, roughly ``1e-5 m / (300 s)^2``, stays smaller
    still - both are thousands of times below the 0.03 m/s^2 rotation term.
    """
    step = dt_s / SECONDS_PER_DAY
    return (
        _barycentric_position(observer, mjd + step)
        - 2.0 * _barycentric_position(observer, mjd)
        + _barycentric_position(observer, mjd - step)
    ) / dt_s**2


def doppler_factor(observer: Observer, mjd: float, sky: SkyPosition) -> float:
    """First-order factor mapping barycentric to topocentric frequency."""
    return 1.0 - float(np.dot(barycentric_velocity(observer, mjd), sky.unit_vector())) / C_LIGHT


def topocentric_frequency(
    barycentric_hz: float, observer: Observer, mjd: float, sky: SkyPosition
) -> float:
    """Frequency a telescope sees for a source at this barycentric frequency."""
    return barycentric_hz * doppler_factor(observer, mjd, sky)


def topocentric_drift(
    barycentric_drift_hz_s: float,
    freq_hz: float,
    observer: Observer,
    mjd: float,
    sky: SkyPosition,
) -> float:
    """Drift a source with this barycentric drift shows at this site and time."""
    factor = doppler_factor(observer, mjd, sky)
    acceleration = float(np.dot(barycentric_acceleration(observer, mjd), sky.unit_vector()))
    return factor * barycentric_drift_hz_s - (freq_hz / factor) * acceleration / C_LIGHT


def barycentric_drift(
    observed_drift_hz_s: float,
    freq_hz: float,
    observer: Observer,
    mjd: float,
    sky: SkyPosition,
) -> float:
    """Invert :func:`topocentric_drift`: what a measured drift implies about the source."""
    factor = doppler_factor(observer, mjd, sky)
    acceleration = float(np.dot(barycentric_acceleration(observer, mjd), sky.unit_vector()))
    return (observed_drift_hz_s + (freq_hz / factor) * acceleration / C_LIGHT) / factor


def predict_observation(
    freq_hz: float,
    drift_hz_s: float,
    observer: Observer,
    mjd_from: float,
    mjd_to: float,
    sky: SkyPosition,
) -> tuple[float, float]:
    """Carry one observation forward: (predicted frequency, predicted drift).

    ``freq_hz`` and ``drift_hz_s`` are the topocentric values measured at
    ``mjd_from``; the return values are what the same source should show at
    ``mjd_to``. The source's own drift accumulates across the gap, so a long
    gap magnifies any error in the measured drift - the caller has to allow for
    that when choosing a frequency window.
    """
    factor_from = doppler_factor(observer, mjd_from, sky)
    factor_to = doppler_factor(observer, mjd_to, sky)
    intrinsic = barycentric_drift(drift_hz_s, freq_hz, observer, mjd_from, sky)
    gap_s = (float(mjd_to) - float(mjd_from)) * SECONDS_PER_DAY
    predicted_hz = (freq_hz / factor_from + intrinsic * gap_s) * factor_to
    predicted_drift = topocentric_drift(intrinsic, predicted_hz, observer, mjd_to, sky)
    return predicted_hz, predicted_drift
