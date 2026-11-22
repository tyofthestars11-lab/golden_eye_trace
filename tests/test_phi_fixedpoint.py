"""phi^2 = phi + 1 as a recurrence: the rate, the depth, and the rejected window."""

from __future__ import annotations

import csv
import json
import math
import os
import sys
import unittest
from fractions import Fraction

from seti_search.phi_fixedpoint import (
    CONTRACTION,
    PHI,
    STEPS_PER_DECADE,
    agreement_depth,
    continued_fraction,
    convergents,
    derivative,
    effective_depth,
    excess,
    fibonacci,
    is_phi_shaped,
    orbit,
    step,
)

LN_PHI = math.log(PHI)


def rung_width_hz(tolerance_rungs: float, freq_hz: float) -> float:
    """Physical width spanned by a rung tolerance, used by the fixture below."""
    return (PHI ** tolerance_rungs - PHI ** (-tolerance_rungs)) * freq_hz


class RecurrenceTests(unittest.TestCase):
    def test_phi_is_the_fixed_point(self):
        self.assertAlmostEqual(step(PHI), PHI, places=15)
        self.assertAlmostEqual(PHI ** 2, PHI + 1.0, places=15)

    def test_the_derivative_at_the_fixed_point_is_the_contraction(self):
        self.assertAlmostEqual(derivative(PHI), -CONTRACTION, places=15)
        self.assertAlmostEqual(CONTRACTION, 1.0 / PHI ** 2, places=16)

    def test_every_start_contracts_at_the_same_rate(self):
        for start in (0.5, 1.0, 2.0, 6.485, 100.0):
            ob = orbit(start, turns=30)
            self.assertAlmostEqual(ob.contraction, CONTRACTION, places=4)

    def test_the_residual_shrinks_two_turns_at_a_time(self):
        # T has a negative derivative at phi, so the orbit straddles the fixed
        # point: residuals alternate but every second one is smaller
        res = orbit(0.5, turns=20).residuals
        for i in range(1, len(res) - 2):
            self.assertLess(res[i + 2], res[i])

    def test_turns_per_decade_match_the_rate(self):
        self.assertAlmostEqual(STEPS_PER_DECADE,
                               1.0 / abs(math.log10(CONTRACTION)), places=12)
        self.assertAlmostEqual(STEPS_PER_DECADE, 2.3925, places=4)

    def test_non_positive_inputs_are_refused(self):
        with self.assertRaises(ValueError):
            step(0.0)
        with self.assertRaises(ValueError):
            derivative(-1.0)
        with self.assertRaises(ValueError):
            orbit(0.0)


class LadderTests(unittest.TestCase):
    def test_phi_reads_add_one_at_every_level(self):
        part = continued_fraction(PHI, depth=12)
        self.assertEqual(part[:2], [1, 1])
        self.assertTrue(is_phi_shaped(PHI))
        self.assertTrue(is_phi_shaped(1.0 / PHI))
        self.assertTrue(is_phi_shaped(PHI ** 2))

    def test_the_channel_grid_says_add_one_exactly_once(self):
        # (k+1)/k = 1 + 1/k, so [1; k]: the +1 is carried after one level
        part = continued_fraction(Fraction(318, 317))
        self.assertEqual(part, [1, 317])
        self.assertEqual(agreement_depth(part), 0)
        self.assertEqual(excess(part), 316)
        self.assertFalse(is_phi_shaped(318 / 317))

    def test_a_convergent_matches_until_its_terminating_two(self):
        part = continued_fraction(Fraction(13, 8))
        self.assertEqual(part, [1, 1, 1, 1, 2])
        self.assertEqual(agreement_depth(part), 3)
        self.assertEqual(excess(part), 1)

    def test_depth_is_monotone_in_closeness_to_phi(self):
        depths = [agreement_depth(continued_fraction(PHI * (1 + j)))
                  for j in (1e-1, 1e-2, 1e-3, 1e-4, 1e-6)]
        self.assertEqual(depths, sorted(depths))
        self.assertGreater(depths[-1], depths[0] + 5)

    def test_effective_depth_is_continuous_and_can_be_negative(self):
        self.assertTrue(math.isinf(effective_depth(PHI)))
        self.assertAlmostEqual(effective_depth(PHI * (1 + 1e-4)), 9.57, places=1)
        self.assertAlmostEqual(effective_depth(PHI * (1 + 1e-2)), 4.78, places=1)
        self.assertLess(effective_depth(869 / 134), 0.0)

    def test_bad_depth_is_refused(self):
        with self.assertRaises(ValueError):
            continued_fraction(0.0)
        with self.assertRaises(ValueError):
            continued_fraction(2.0, depth=0)


class ConvergentDecayTests(unittest.TestCase):
    def test_the_ladder_is_the_fibonacci_ladder(self):
        conv = convergents(PHI, depth=8)
        self.assertEqual([(p, q) for p, q, _ in conv[:5]],
                         [(1, 1), (2, 1), (3, 2), (5, 3), (8, 5)])
        fibs = {fibonacci(n) for n in range(1, 12)}
        for p, q, _ in conv:
            self.assertIn(p, fibs)
            self.assertIn(q, fibs)

    def test_each_level_costs_one_over_phi_squared(self):
        conv = convergents(PHI, depth=14)
        ratios = [conv[i][2] / conv[i - 1][2] for i in range(1, len(conv))]
        for r in ratios[-6:]:
            self.assertAlmostEqual(r, CONTRACTION, places=3)


class RejectedToleranceTest(unittest.TestCase):
    """The ±0.01 rung window, kept as a fixture for why it was withdrawn."""

    def test_a_rung_window_is_not_scale_invariant(self):
        low = rung_width_hz(0.01, 13.4e6)
        high = rung_width_hz(0.01, 86.9e6)
        self.assertAlmostEqual(high / low, 86.9 / 13.4, places=9)
        self.assertGreater(high / low, 6.0)

    def test_the_tolerated_width_is_not_a_constant_of_the_instrument(self):
        widths = [rung_width_hz(0.01, f) for f in (13.4e6, 47e6, 86.9e6)]
        self.assertAlmostEqual(widths[0] / 1e6, 0.1290, places=3)
        self.assertAlmostEqual(widths[-1] / 1e6, 0.8363, places=3)
        # a fixed rung tolerance therefore tests a different physical width
        # at every frequency, which is why the threshold was withdrawn
        self.assertGreater(max(widths) / min(widths), 6.0)

    def test_the_rate_itself_is_scale_free(self):
        # the rate, unlike a window, has no unit: it is the same number at
        # every frequency because it is the derivative at the fixed point
        for f in (13.4e6, 47e6, 86.9e6):
            self.assertAlmostEqual(math.exp(math.log(PHI ** 2) * 1.0),
                                   PHI ** 2, places=12)
            self.assertAlmostEqual(derivative(PHI), -1.0 / PHI ** 2, places=16)


class CampaignRatioTests(unittest.TestCase):
    """The depth reading applied to the ratios the campaign actually measured."""

    @classmethod
    def setUpClass(cls):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        table = os.path.join(root, "results", "comb_vs_scale.csv")
        if not os.path.exists(table):
            raise unittest.SkipTest("run scripts/phi_model_discrimination.py first")
        sys.path.insert(0, os.path.join(root, "scripts"))
        import phi_plus_one
        cls.table = table
        cls.script = phi_plus_one
        cls.payload = phi_plus_one.build()
        with open(table, newline="", encoding="utf-8") as fh:
            cls.rows = list(csv.DictReader(fh))

    def test_every_row_is_read_from_the_table(self):
        camp = self.payload["campaign"]
        self.assertEqual(len(camp), len(self.rows))
        for row, entry in zip(self.rows, camp):
            self.assertEqual(entry["date"], row["date"])
            self.assertEqual(entry["source"], row["source"])
            self.assertEqual(entry["verdict"], row["verdict"])
            # the fitted interval is carried as printed, not re-derived
            self.assertEqual(entry["ratios"][0]["printed"], row["scale_q"])
            self.assertEqual(entry["ratios"][0]["value"], float(row["scale_q"]))

    def test_the_span_ratio_is_two_printed_intervals(self):
        for row, entry in zip(self.rows, self.payload["campaign"]):
            span = Fraction(row["span_mhz"])
            step = Fraction(row["comb_df_mhz"])
            self.assertEqual(entry["ratios"][1]["value"], float(span / step))

    def test_no_measured_ratio_reaches_phi(self):
        flat = [r for e in self.payload["campaign"] for r in e["ratios"]]
        self.assertTrue(flat)
        for r in flat:
            self.assertFalse(r["unbounded"], msg=r["label"])
            self.assertLessEqual(r["depth"], 1, msg=r["label"])
        self.assertTrue(any(r["effective_depth"] < 0 for r in flat))

    def test_a_single_level_of_depth_is_cheap(self):
        # depth 1 is granted by any ratio whose fractional part exceeds 1/2, so
        # one level says nothing about phi without the continuous reading
        flat = [r for e in self.payload["campaign"] for r in e["ratios"]]
        one = [r for r in flat if r["depth"] == 1]
        self.assertEqual(len(one), 2)
        for r in one:
            self.assertGreater(r["value"] - math.floor(r["value"]), 0.5)
            self.assertGreater(r["value"], 2.5)

    def test_an_empty_recursion_looks_like_depth_zero(self):
        # 20240920's span is exactly two fitted comb steps, [2]: depth 0 and
        # excess 0 with no fractional level at all, the ambiguity that keeps
        # `unbounded` beside every row
        entry = [e for e in self.payload["campaign"] if e["date"] == "20240920"]
        self.assertEqual(len(entry), 1)
        ratio = entry[0]["ratios"][1]
        self.assertTrue(ratio["empty_recursion"])
        self.assertEqual((ratio["depth"], ratio["excess"]), (0, 0))
        self.assertFalse(ratio["unbounded"])

    def test_masked_rows_stay_unclaimed(self):
        camp = self.payload["campaign"]
        masked = [e for e in camp if e["verdict"] == "MASK-MANUFACTURED SIGNAL"]
        self.assertEqual(len(masked), 2)
        self.assertTrue(all(e["diagnostic_only"] for e in masked))
        self.assertEqual(self.payload["campaign_summary"]["rows"], len(camp))

    def test_the_additive_ratio_is_1_plus_step_over_f_low(self):
        # the comb's own adjacent-line ratio at each group's lowest line, built
        # from the fit's step and the fixture's frequency, exactly
        with open(os.path.join(os.path.dirname(self.table),
                               "expanded_phi_groups.json"), encoding="utf-8") as fh:
            groups = json.load(fh)["groups"]
        index = {(g["source"], g["date"], int(g["n_lines"])): g
                 for g in groups}
        self.assertEqual(len(index), len(self.rows))
        for row, entry in zip(self.rows, self.payload["campaign"]):
            group = index[(row["source"], row["date"], int(row["n_lines"]))]
            f_low = Fraction(str(group["freqs_mhz"][0]))
            self.assertEqual(float(f_low), min(group["freqs_mhz"]))
            expected = 1 + Fraction(row["comb_df_mhz"]) / f_low
            self.assertEqual(entry["additive"]["value"], float(expected))
            self.assertEqual(entry["lowest_line_mhz"], float(f_low))

    def test_the_additive_ratio_is_diagnostic_only(self):
        # predicted by the fit, not measured from the lines: it must not enter
        # the count of measurements, and it must not reach phi
        add = [e["additive"] for e in self.payload["campaign"]]
        self.assertEqual(len(add), len(self.payload["campaign"]))
        measured = sum(len(e["ratios"]) for e in self.payload["campaign"])
        self.assertEqual(self.payload["campaign_summary"]["ratios"], measured)
        for r in add:
            self.assertTrue(r["diagnostic_only"])
            self.assertFalse(r["unbounded"], msg=r["label"])
            self.assertEqual(r["depth"], 0, msg=r["label"])
        summary = self.payload["campaign_summary"]["additive"]
        self.assertTrue(summary["diagnostic_only"])
        self.assertEqual(summary["ratios"], len(add))
        self.assertEqual(sum(summary["depth_counts"].values()), len(add))

    def _groups(self):
        with open(os.path.join(os.path.dirname(self.table),
                               "expanded_phi_groups.json"),
                  encoding="utf-8") as fh:
            groups = json.load(fh)["groups"]
        return {(g["source"], g["date"], int(g["n_lines"])): g
                for g in groups}

    def test_every_adjacent_pair_is_read(self):
        # the additive model's ratio is carried to every pair in the group, not
        # only to the first line where the two models coincide
        index = self._groups()
        for row, entry in zip(self.rows, self.payload["campaign"]):
            group = index[(row["source"], row["date"], int(row["n_lines"]))]
            freqs = [Fraction(str(f)) for f in group["freqs_mhz"]]
            step = Fraction(row["comb_df_mhz"])
            pairs = entry["additive_pairs"]
            self.assertEqual(len(pairs), len(freqs) - 1)
            for i, pair in enumerate(pairs):
                self.assertEqual(pair["value"], float(1 + step / freqs[i]))
                self.assertEqual(pair["line_low_mhz"], float(freqs[i]))
                self.assertEqual(pair["line_high_mhz"], float(freqs[i + 1]))
                self.assertTrue(pair["diagnostic_only"])

    def test_the_pairs_fall_and_cross_the_fitted_interval_once(self):
        # 1 + d/f falls as f rises, so each group brackets q at exactly one
        # frequency: the models touch there and separate away from it
        for entry in self.payload["campaign"]:
            pairs = entry["additive_pairs"]
            values = [p["value"] for p in pairs]
            self.assertEqual(values, sorted(values, reverse=True))
            deltas = [p["signed_delta_vs_scale_q"] for p in pairs]
            self.assertGreater(deltas[0], 0.0)
            self.assertLess(deltas[-1], 0.0)
            self.assertEqual(
                sum(1 for a, b in zip(deltas, deltas[1:]) if a * b < 0), 1)
            span = entry["additive_span"]
            self.assertTrue(span["crossing_inside_span"])
            self.assertEqual(span["gap_sign_changes"], 1)
            # the first line sees only the near side of the crossing
            self.assertGreater(span["drift_over_first_gap"], 1.5)

    def test_the_compounded_value_is_the_product_of_the_pairs(self):
        index = self._groups()
        for row, entry in zip(self.rows, self.payload["campaign"]):
            group = index[(row["source"], row["date"], int(row["n_lines"]))]
            freqs = [Fraction(str(f)) for f in group["freqs_mhz"]]
            step = Fraction(row["comb_df_mhz"])
            product = Fraction(1)
            for f in freqs[:-1]:
                product *= 1 + step / f
            span = entry["additive_span"]
            self.assertEqual(span["value"], float(product))
            self.assertTrue(span["diagnostic_only"])
            self.assertEqual(span["pairs"], len(freqs) - 1)
            scale_cumulative = Fraction(row["scale_q"]) ** (len(freqs) - 1)
            self.assertEqual(span["reference_scale"]["value"],
                             float(scale_cumulative))
            self.assertEqual(span["relative_gap_to_scale"],
                             float(product / scale_cumulative - 1))
            self.assertTrue(span["reference_scale"]["diagnostic_only"])

    def test_the_pair_and_span_readings_stay_out_of_the_counts(self):
        camp = self.payload["campaign"]
        summary = self.payload["campaign_summary"]["additive"]
        pairs = sum(len(e["additive_pairs"]) for e in camp)
        self.assertEqual(summary["pairs"], pairs)
        self.assertEqual(sum(summary["pair_depth_counts"].values()), pairs)
        self.assertEqual(sum(summary["span_depth_counts"].values()), len(camp))
        self.assertEqual(len(summary["span_relative_gap_to_scale"]), len(camp))
        self.assertTrue(all(c["inside_span"] for c in summary["crossings"]))
        # predicted, never measured: the measured count is unchanged
        measured = sum(len(e["ratios"]) for e in camp)
        self.assertEqual(self.payload["campaign_summary"]["ratios"], measured)
        self.assertLess(summary["max_pair_gap_to_scale_q"], 1e-2)

    def test_the_first_line_entry_is_the_first_pair(self):
        # one code path: the additive entry is additive_pairs[0], so the lowest
        # line cannot be computed twice and disagree with itself
        for entry in self.payload["campaign"]:
            first = entry["additive_pairs"][0]
            add = entry["additive"]
            self.assertEqual(add["pair"], 1)
            for key in ("value", "depth", "excess", "effective_depth",
                        "levels", "unbounded", "printed"):
                self.assertEqual(add[key], first[key], msg=key)
            self.assertEqual(add["lowest_line_mhz"], first["line_low_mhz"])
            self.assertEqual(add["gap_to_scale_q"], first["gap_to_scale_q"])
            self.assertTrue(add["diagnostic_only"])

    def test_the_dead_crossing_field_is_gone(self):
        # crossing_after_pair was written and nothing read it
        for entry in self.payload["campaign"]:
            span = entry["additive_span"]
            self.assertNotIn("crossing_after_pair", span)
            self.assertIn("scale_crossing_mhz", span)

    def test_the_measured_pair_ratio_is_scored_against_both_fits(self):
        # f_high/f_low measured on the lines, beside both predictions: the
        # divergence is scored on the data, not only between the two fits
        index = self._groups()
        scored = 0
        for row, entry in zip(self.rows, self.payload["campaign"]):
            group = index[(row["source"], row["date"], int(row["n_lines"]))]
            freqs = [Fraction(str(f)) for f in group["freqs_mhz"]]
            step = Fraction(row["comb_df_mhz"])
            q = Fraction(row["scale_q"])
            for i, m in enumerate(entry["pair_residuals"]):
                measured = freqs[i + 1] / freqs[i]
                additive = 1 + step / freqs[i]
                self.assertEqual(m["measured"], float(measured))
                self.assertEqual(m["relative_vs_scale"],
                                 float(measured / q - 1))
                self.assertEqual(m["relative_vs_additive"],
                                 float(measured / additive - 1))
                self.assertEqual(
                    m["closer"],
                    "additive" if abs(measured / additive - 1)
                    < abs(measured / q - 1) else "scale")
                self.assertTrue(m["diagnostic_only"])
                scored += 1
        summary = self.payload["campaign_summary"]["line_fit"]
        self.assertEqual(summary["pairs"], scored)
        self.assertEqual(sum(summary["closer_counts"].values()), scored)
        self.assertGreaterEqual(summary["max_abs_relative_vs_scale"],
                                summary["mean_abs_relative_vs_scale"])

    def test_the_span_residual_is_measured_from_the_lines(self):
        index = self._groups()
        for row, entry in zip(self.rows, self.payload["campaign"]):
            group = index[(row["source"], row["date"], int(row["n_lines"]))]
            freqs = [Fraction(str(f)) for f in group["freqs_mhz"]]
            span = entry["additive_span"]
            measured = freqs[-1] / freqs[0]
            scale_cumulative = Fraction(row["scale_q"]) ** (len(freqs) - 1)
            product = Fraction(1)
            for f in freqs[:-1]:
                product *= 1 + Fraction(row["comb_df_mhz"]) / f
            self.assertEqual(span["measured_span_ratio"], float(measured))
            self.assertEqual(span["relative_vs_scale"],
                             float(measured / scale_cumulative - 1))
            self.assertEqual(span["relative_vs_additive"],
                             float(measured / product - 1))
            self.assertEqual(
                span["closer"],
                "additive" if abs(measured / product - 1)
                < abs(measured / scale_cumulative - 1) else "scale")
            self.assertEqual(
                span["measured_between_predictions"],
                min(product, scale_cumulative) <= measured
                <= max(product, scale_cumulative))

    def test_the_scatter_is_read_from_the_campaign_table(self):
        # the rms columns are ingest facts, held as printed, never re-derived
        for row, entry in zip(self.rows, self.payload["campaign"]):
            noise = entry["noise"]
            self.assertEqual(noise["as_printed"],
                             [row["rms_comb_channels"],
                              row["rms_scale_channels"]])
            self.assertEqual(noise["rms_comb_channels"],
                             float(Fraction(row["rms_comb_channels"])))
            self.assertEqual(noise["rms_scale_channels"],
                             float(Fraction(row["rms_scale_channels"])))
            self.assertGreater(noise["rms_comb_channels"], 0.0)
            self.assertGreater(noise["rms_scale_channels"], 0.0)

    def test_each_pair_miss_is_a_line_offset_over_its_own_scatter(self):
        # a ratio miss displaces the pair's upper line, so the miss in
        # channels is relative * predicted high line, and one rms divides it
        index = self._groups()
        channel_hz = self.script.CHANNEL_HZ
        for row, entry in zip(self.rows, self.payload["campaign"]):
            group = index[(row["source"], row["date"], int(row["n_lines"]))]
            freqs = [Fraction(str(f)) for f in group["freqs_mhz"]]
            step = Fraction(row["comb_df_mhz"])
            q = Fraction(row["scale_q"])
            rms_comb = float(Fraction(row["rms_comb_channels"]))
            rms_scale = float(Fraction(row["rms_scale_channels"]))
            for i, m in enumerate(entry["pair_residuals"]):
                measured = freqs[i + 1] / freqs[i]
                additive = 1 + step / freqs[i]
                offset_add = (float((freqs[i] + step)
                                    * (measured / additive - 1))
                              * 1e6 / channel_hz)
                offset_scale = (float(q * freqs[i] * (measured / q - 1))
                                * 1e6 / channel_hz)
                self.assertAlmostEqual(m["offset_channels_vs_additive"],
                                       offset_add, places=9)
                self.assertAlmostEqual(m["offset_channels_vs_scale"],
                                       offset_scale, places=9)
                self.assertAlmostEqual(m["over_comb_rms"],
                                       offset_add / rms_comb, places=9)
                self.assertAlmostEqual(m["over_scale_rms"],
                                       offset_scale / rms_scale, places=9)

    def test_the_span_miss_is_a_high_line_offset_over_the_scatter(self):
        index = self._groups()
        channel_hz = self.script.CHANNEL_HZ
        for row, entry in zip(self.rows, self.payload["campaign"]):
            group = index[(row["source"], row["date"], int(row["n_lines"]))]
            freqs = [Fraction(str(f)) for f in group["freqs_mhz"]]
            step = Fraction(row["comb_df_mhz"])
            product = Fraction(1)
            for f in freqs[:-1]:
                product *= 1 + step / f
            measured = freqs[-1] / freqs[0]
            span = entry["additive_span"]
            rms_comb = float(Fraction(row["rms_comb_channels"]))
            rms_scale = float(Fraction(row["rms_scale_channels"]))
            offset_add = (float(product * freqs[0] * (measured / product - 1))
                          * 1e6 / channel_hz)
            self.assertAlmostEqual(span["offset_channels_vs_additive"],
                                   offset_add, places=9)
            self.assertAlmostEqual(span["over_comb_rms"],
                                   offset_add / rms_comb, places=9)
            self.assertAlmostEqual(
                span["over_scale_rms"],
                float(span["offset_channels_vs_scale"]) / rms_scale,
                places=9)

    def test_the_two_fits_gap_is_read_in_scatter_units(self):
        channel_hz = self.script.CHANNEL_HZ
        for row, entry in zip(self.rows, self.payload["campaign"]):
            rms_comb = float(Fraction(row["rms_comb_channels"]))
            rms_scale = float(Fraction(row["rms_scale_channels"]))
            for m in entry["pair_residuals"]:
                gap = (m["predicted_high_mhz"]["additive"]
                       - m["predicted_high_mhz"]["scale"]) * 1e6 / channel_hz
                self.assertAlmostEqual(m["fit_gap_channels"], gap, places=9)
                self.assertAlmostEqual(m["fit_gap_over_comb_rms"],
                                       gap / rms_comb, places=9)
                self.assertAlmostEqual(m["fit_gap_over_scale_rms"],
                                       gap / rms_scale, places=9)
            span = entry["additive_span"]
            gap = (span["predicted_high_mhz"]["additive"]
                   - span["predicted_high_mhz"]["scale"]) * 1e6 / channel_hz
            self.assertAlmostEqual(span["fit_gap_channels"], gap, places=9)
            self.assertAlmostEqual(span["fit_gap_over_comb_rms"],
                                   gap / rms_comb, places=9)
            self.assertAlmostEqual(span["fit_gap_over_scale_rms"],
                                   gap / rms_scale, places=9)

    def test_the_noise_scaling_stays_out_of_the_model_selection_counts(self):
        summary = self.payload["campaign_summary"]
        line_fit = summary["line_fit"]
        noise = line_fit["noise_scale"]
        self.assertTrue(noise["diagnostic_only"])
        self.assertEqual(summary["ratios"],
                         sum(len(e["ratios"])
                             for e in self.payload["campaign"]))
        self.assertEqual(summary["additive"]["pairs"],
                         sum(len(e["additive_pairs"])
                             for e in self.payload["campaign"]))
        self.assertEqual(noise["pairs"], line_fit["pairs"])
        for key in ("aic", "bic", "delta_aic", "null_p", "p_value",
                    "verdict"):
            self.assertNotIn(key, noise)
        self.assertLessEqual(noise["pairs_within_comb_rms"], noise["pairs"])
        self.assertLessEqual(noise["pairs_within_scale_rms"], noise["pairs"])
        self.assertEqual(len(noise["rms_comb_channels"]),
                         len(self.payload["campaign"]))

    def test_the_tables_carry_the_miss_in_scatter_units(self):
        root = os.path.dirname(os.path.dirname(self.table))
        with open(os.path.join(root, "PHI_PLUS_ONE.md"),
                  encoding="utf-8") as fh:
            text = fh.read()
        headers = [line for line in text.splitlines()
                   if line.startswith("| group |")]
        pair_header = [h for h in headers if "| pair |" in h]
        self.assertEqual(len(pair_header), 1)
        self.assertIn("miss / comb rms", pair_header[0])
        self.assertIn("miss / scale rms", pair_header[0])
        span_header = [h for h in headers if "| pairs |" in h]
        self.assertEqual(len(span_header), 1)
        self.assertIn("miss / comb rms", span_header[0])
        self.assertIn("miss / scale rms", span_header[0])

    def test_the_tables_carry_no_excess_column(self):
        # the pair and span tables were trimmed to depth and the continuous
        # reading; excess still travels in the JSON but is not tabulated
        root = os.path.dirname(os.path.dirname(self.table))
        with open(os.path.join(root, "PHI_PLUS_ONE.md"), encoding="utf-8") as fh:
            text = fh.read()
        headers = [line for line in text.splitlines()
                   if line.startswith("| group |")]
        pair_header = [h for h in headers if "| pair |" in h]
        self.assertEqual(len(pair_header), 1)
        self.assertNotIn("excess", pair_header[0])
        self.assertIn("measured", pair_header[0])
        span_header = [h for h in headers if "| pairs |" in h]
        self.assertEqual(len(span_header), 1)
        self.assertNotIn("excess", span_header[0])
        self.assertIn("measured", span_header[0])


if __name__ == "__main__":
    unittest.main()
