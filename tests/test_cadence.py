"""Does the drift prediction reproduce real celestial drift?

If these tests pass, the ON/OFF check is testing physics and not just arithmetic:
the ephemeris has to put Earth where it belongs, the rotation term has to have
the right size, and the forward and inverse drift models have to be consistent.
"""

import unittest

import numpy as np

from seti_search.cadence import (
    SITES,
    Observer,
    SkyPosition,
    barycentric_acceleration,
    barycentric_drift,
    barycentric_velocity,
    doppler_factor,
    gmst_deg,
    predict_observation,
    topocentric_drift,
)

MJD = 58967.23  # 2020-04-28, the STARE2 candidate epoch


class EphemerisTests(unittest.TestCase):
    def test_gmst_at_j2000(self):
        # 2000-01-01 12:00 UT is 18h 41m 50.5s of sidereal time
        self.assertAlmostEqual(gmst_deg(51544.5), 280.46061837, delta=1e-4)

    def test_gmst_advances_by_one_sidereal_day(self):
        # 360.98564736 degrees per solar day
        self.assertAlmostEqual(
            (gmst_deg(MJD + 1.0) - gmst_deg(MJD)) % 360.0, 0.98564736, delta=1e-3
        )

    def test_orbital_speed_is_earths(self):
        # at the pole there is no rotation term, so this is the orbit alone
        speed = float(np.linalg.norm(barycentric_velocity(Observer(90.0, 0.0), MJD)))
        self.assertGreater(speed / 1e3, 28.5)
        self.assertLess(speed / 1e3, 31.0)

    def test_orbital_speed_peaks_at_perihelion(self):
        speeds = [
            float(np.linalg.norm(barycentric_velocity(Observer(90.0, 0.0), MJD + 91.3 * k)))
            for k in range(4)
        ]
        # eccentricity 0.0167 gives about 3.4 per cent between perihelion and
        # aphelion; a wrong obliquity or a wrong sign would not reproduce it
        self.assertGreater(max(speeds) / min(speeds), 1.02)
        self.assertLess(max(speeds) / min(speeds), 1.08)

    def test_rotation_acceleration_is_centripetal(self):
        equator = barycentric_acceleration(Observer(0.0, 0.0), MJD)
        pole = barycentric_acceleration(Observer(90.0, 0.0), MJD)
        # omega^2 R = 0.0337 m/s^2 at the equator, absent at the pole
        self.assertAlmostEqual(np.linalg.norm(equator - pole), 0.0337, delta=0.002)

    def test_doppler_factor_is_within_the_expected_range(self):
        sky = SkyPosition(75.0, 20.0)
        factor = doppler_factor(SITES["ovro"], MJD, sky)
        self.assertLess(abs(factor - 1.0), 2e-4)  # 30 km/s is 1e-4 of c


class DriftModelTests(unittest.TestCase):
    def test_round_trip(self):
        observer, sky = SITES["parkes"], SkyPosition(300.0, -45.0)
        freq = 1420.4e6
        for barycentric in (-1.5, -0.05, 0.0, 0.3, 2.0):
            observed = topocentric_drift(barycentric, freq, observer, MJD, sky)
            self.assertAlmostEqual(
                barycentric_drift(observed, freq, observer, MJD, sky), barycentric, places=9
            )

    def test_prediction_to_the_same_epoch_is_the_identity(self):
        observer, sky = SITES["ovro"], SkyPosition(75.0, 20.0)
        freq = 1420.4e6
        drift = topocentric_drift(0.2, freq, observer, MJD, sky)
        pred_freq, pred_drift = predict_observation(freq, drift, observer, MJD, MJD, sky)
        self.assertAlmostEqual(pred_freq, freq, places=3)
        self.assertAlmostEqual(pred_drift, drift, places=6)

    def test_prediction_moves_the_drift_by_a_physical_amount(self):
        observer, sky = SITES["ovro"], SkyPosition(75.0, 20.0)
        freq = 1420.4e6
        start = topocentric_drift(0.0, freq, observer, MJD, sky)
        _, later = predict_observation(freq, start, observer, MJD, MJD + 1.0 / 24.0, sky)
        change = abs(later - start)
        # an hour of rotation moves the projected acceleration, but not the sky
        self.assertGreater(change, 1e-4)
        self.assertLess(change, 0.3)

    def test_a_fixed_source_cannot_drift_arbitrarily_fast(self):
        freq = 1420.4e6
        worst = 0.0
        for dec in (-80.0, -40.0, 0.0, 40.0, 80.0):
            for ra in range(0, 360, 30):
                for hour in (0.0, 0.25, 0.5, 0.75):
                    drift = topocentric_drift(
                        0.0, freq, SITES["gbt"], MJD + hour, SkyPosition(float(ra), dec)
                    )
                    worst = max(worst, abs(drift))
        # nu * a_max / c with a_max about 0.04 m/s^2
        self.assertLess(worst, 0.35)
        self.assertGreater(worst, 0.02)


if __name__ == "__main__":
    unittest.main()
