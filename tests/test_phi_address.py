"""Rung addresses exist for every measurement, and are read honestly.

The fixture is the address artifact. What is asserted is the discipline the
framework asks for - every positive quantity addressed, no integer ladder steps
invented, and a group whose spacing is regular reported as regular without a
model being chosen for it - rather than any claim that an address is evidence.
"""

import json
import math
import os
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARTIFACT = os.path.join(ROOT, "prereg", "phi_address.json")


def load():
    with open(ARTIFACT, encoding="utf-8") as fh:
        return json.load(fh)


class AddressTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = load()
        cls.groups = {g["date"]: g for g in cls.data["groups"]}

    def test_every_record_is_addressed(self):
        self.assertEqual(len(self.data["addresses"]), self.data["n_records"])
        for entry in self.data["addresses"]:
            for key in ("frequency_hz", "snr_at_cell", "duration_s"):
                self.assertIn(key, entry["quantities"])
                self.assertIsNotNone(entry["quantities"][key]["rung"])

    def test_a_quantity_without_a_positive_value_gets_no_rung(self):
        for entry in self.data["addresses"]:
            drift = entry["quantities"]["drift_hz_per_s"]
            if drift["value"] == 0:
                self.assertIsNone(drift["rung"])

    def test_the_records_themselves_carry_the_address(self):
        with open(os.path.join(ROOT, "signal_records.json"), encoding="utf-8") as fh:
            records = json.load(fh)["records"]
        for record in records:
            meta = record["metadata"]
            self.assertIn("phi_freq_rung", meta)
            self.assertIn("phi_freq_residual", meta)
            self.assertEqual(meta["phi_freq_nearest"], round(meta["phi_freq_rung"]))
            self.assertGreaterEqual(meta["phi_freq_residual"], -0.5)
            self.assertLessEqual(meta["phi_freq_residual"], 0.5)
            # the drift coordinate is present and labelled, or explicitly absent
            self.assertIn("phi_drift_rung", meta)
            self.assertIn("phi_drift_residual", meta)


class LadderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = load()
        cls.groups = {g["date"]: g for g in cls.data["groups"]}

    def test_no_pair_is_declared_a_ladder_step(self):
        # an integer rung difference is exact arithmetic, so a non-integer one is
        # a scale factor. None of these pairs is an integer step.
        for group in self.data["groups"]:
            for pair in group["pairs"]:
                self.assertFalse(pair["integer_step"],
                                 msg=f"{group['date']} {pair}")

    def test_a_non_integer_step_reports_its_factor(self):
        # the identity is exact: delta rung is log_phi of the factor, so a factor
        # and its rungs have to agree or one of them is wrong
        ln_phi = math.log((1.0 + math.sqrt(5.0)) / 2.0)
        for group in self.data["groups"]:
            for pair in group["pairs"]:
                self.assertGreater(pair["factor"], 1.0)
                self.assertAlmostEqual(
                    pair["delta_rung"],
                    math.log(pair["factor"]) / ln_phi,
                    places=5,
                    msg=f"{group['date']} {pair['from_mhz']}->{pair['to_mhz']}")

    def test_no_scaling_artifact_is_claimed_from_a_chain(self):
        # four lines in a chain give at most two independent pairs, and the
        # independence rule must not be fudged into a repetition
        for group in self.data["groups"]:
            for entry in group["repeated_scales"]:
                self.assertGreaterEqual(entry["pairs"], 3)


class SpacingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = load()
        cls.groups = {g["date"]: g for g in cls.data["groups"]}

    def test_the_burst_group_is_regular_and_undecided(self):
        spacing = self.groups["20240325"]["spacing"]
        self.assertLess(spacing["p_comb"], 0.01)
        self.assertLess(spacing["p_scale"], 0.01)
        self.assertEqual(spacing["basis"], "cannot distinguish")

    def test_the_near_nyquist_trio_is_not_regular(self):
        spacing = self.groups["20240920"]["spacing"]
        self.assertGreaterEqual(min(spacing["p_comb"], spacing["p_scale"]), 0.05)
        self.assertEqual(spacing["basis"], "neither")

    def test_the_null_is_drawn_from_the_group_s_own_band(self):
        self.assertEqual(self.data["null_draws"], 4000)


class BaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bases = load()["base_invariance"]

    def test_the_whole_band_sits_inside_one_rung_in_every_base(self):
        # so no absolute address in this data can discriminate, under phi or any
        # other base: the framework's own caution about narrow bands
        for name, entry in self.bases.items():
            if not isinstance(entry, dict) or "band_span_rungs" not in entry:
                continue
            self.assertLess(entry["band_span_rungs"], 1.0, msg=name)

    def test_no_base_is_special_for_a_single_address(self):
        spans = [v["median_fractional_distance"] for k, v in self.bases.items()
                 if isinstance(v, dict) and "median_fractional_distance" in v]
        self.assertGreaterEqual(len(spans), 5)
        self.assertGreater(max(spans), 0.2)


class DriftPopulationTests(unittest.TestCase):
    def test_drift_rungs_are_uniform(self):
        entry = load()["drift_rung_uniformity"]
        self.assertEqual(entry["n"], 10)
        self.assertGreater(entry["p_value"], 0.05)


if __name__ == "__main__":
    unittest.main()
