"""The omni layer refuses a disposition its gates do not support."""

import unittest

from seti_search.omni import (
    LAYERS,
    REQUIRED,
    STATUSES,
    Audit,
    EvidenceError,
    Gate,
    PhiAsEvidence,
    audit,
)
from seti_search.signal_record import DISPOSITIONS


def gate(gate_id, status="PASS", layer="measurement", evidence="a number"):
    return Gate(id=gate_id, layer=layer, status=status, evidence=evidence)


CONFIRMING = (
    gate("detected"),
    gate("window"),
    gate("repeat", layer="confirmation"),
    gate("off_source", layer="confirmation"),
    gate("independent_telescope", layer="confirmation"),
    gate("barycentric", layer="confirmation"),
)


class VocabularyTests(unittest.TestCase):
    def test_every_disposition_has_a_gate_list(self):
        self.assertEqual(set(REQUIRED), set(DISPOSITIONS))

    def test_only_the_negative_dispositions_are_free(self):
        free = {d for d, req in REQUIRED.items() if not req}
        self.assertEqual(free, {"NULL RESULT"})

    def test_a_confirmation_requires_the_gates_a_single_beam_cannot_supply(self):
        required = REQUIRED["CONFIRMED TECHNOSIGNATURE CANDIDATE"]
        for gate_id in ("off_source", "independent_telescope", "barycentric",
                        "repeat"):
            self.assertIn(gate_id, required)

    def test_no_disposition_requires_a_phi_gate(self):
        for disposition, required in REQUIRED.items():
            for gate_id in required:
                self.assertNotIn("phi", gate_id.lower(),
                                 msg=f"{disposition} requires {gate_id}")

    def test_the_layers_are_ordered_from_measurement_to_confirmation(self):
        self.assertEqual(LAYERS[0], "measurement")
        self.assertEqual(LAYERS[-1], "confirmation")
        self.assertIn("phi", LAYERS)


class GateTests(unittest.TestCase):
    def test_a_pass_without_evidence_is_rejected(self):
        with self.assertRaises(EvidenceError) as caught:
            gate("detected", evidence="   ")
        self.assertIn("no evidence", str(caught.exception))

    def test_an_unknown_status_is_rejected(self):
        with self.assertRaises(EvidenceError):
            gate("detected", status="MAYBE")

    def test_an_unknown_layer_is_rejected(self):
        with self.assertRaises(EvidenceError):
            gate("detected", layer="vibes")

    def test_the_statuses_are_closed(self):
        self.assertEqual(STATUSES, ("PASS", "FAIL", "UNTESTABLE", "NOT RUN"))
        for status in STATUSES:
            gate("detected", status=status)


class AuditTests(unittest.TestCase):
    def test_a_supported_claim_stands(self):
        result = audit("CONFIRMED TECHNOSIGNATURE CANDIDATE", CONFIRMING)
        self.assertFalse(result.downgraded)
        self.assertFalse(result.escalation)
        self.assertEqual(result.deficits, ())

    def test_an_untestable_gate_blocks_a_confirmation(self):
        gates = [
            gate("detected"),
            gate("window"),
            gate("repeat", layer="confirmation"),
            gate("off_source", status="UNTESTABLE", layer="confirmation",
                 evidence="one disk-integrated beam, no off-source scan"),
            gate("independent_telescope", status="UNTESTABLE",
                 layer="confirmation", evidence="one telescope"),
            gate("barycentric", status="UNTESTABLE", layer="confirmation",
                 evidence="no WCS, no sky position"),
        ]
        result = audit("CONFIRMED TECHNOSIGNATURE CANDIDATE", gates)
        self.assertTrue(result.downgraded)
        self.assertTrue(result.escalation)
        self.assertEqual(result.permitted, "UNRESOLVED")
        self.assertEqual(len(result.deficits), 3)
        self.assertTrue(any("off_source: UNTESTABLE" in d for d in result.deficits))
        self.assertTrue(any("no off-source scan" in d for d in result.deficits))

    def test_a_missing_gate_is_named_as_not_run(self):
        result = audit("CONFIRMED TECHNOSIGNATURE CANDIDATE", [gate("detected")])
        self.assertTrue(any("NOT RUN" in d for d in result.deficits))
        self.assertIn("window", " ".join(result.deficits))

    def test_a_failed_detection_leaves_only_a_null_result(self):
        result = audit("UNRESOLVED", [gate("detected", status="FAIL",
                                           evidence="5.16 sigma vs 5.5")])
        self.assertEqual(result.permitted, "NULL RESULT")
        self.assertTrue(result.downgraded)
        self.assertFalse(result.escalation)

    def test_an_unresolved_candidate_keeps_its_name_when_detected(self):
        result = audit("UNRESOLVED", [gate("detected")])
        self.assertFalse(result.downgraded)

    def test_a_null_result_needs_nothing(self):
        result = audit("NULL RESULT", [])
        self.assertFalse(result.downgraded)

    def test_a_nonlinear_claim_needs_the_control(self):
        gates = [gate("detected"),
                 gate("nonlinear_beats_control", status="FAIL",
                      layer="timedomain",
                      evidence="gain 0.57 over linear, warp control wins")]
        result = audit("NONLINEAR-COHERENT BUT UNRESOLVED", gates)
        self.assertEqual(result.permitted, "UNRESOLVED")
        self.assertIn("warp control", " ".join(result.deficits))

    def test_an_astrophysical_claim_needs_its_mechanism(self):
        gates = [gate("detected")]
        result = audit("ASTROPHYSICAL", gates)
        self.assertTrue(result.downgraded)
        self.assertIn("astrophysical_explanation", " ".join(result.deficits))

    def test_an_alias_claim_needs_the_alias_signature(self):
        result = audit("INSTRUMENTAL ALIAS", [gate("detected")])
        self.assertEqual(result.permitted, "UNRESOLVED")
        self.assertIn("alias_signature", " ".join(result.deficits))

    def test_a_duplicate_gate_id_is_a_bug(self):
        with self.assertRaises(EvidenceError) as caught:
            audit("UNRESOLVED", [gate("detected"), gate("detected")])
        self.assertIn("share the id", str(caught.exception))

    def test_an_unknown_disposition_is_refused(self):
        with self.assertRaises(EvidenceError):
            audit("PROBABLY ALIENS", [gate("detected")])


class PhiTests(unittest.TestCase):
    def test_phi_cannot_be_cited_as_support(self):
        gates = list(CONFIRMING) + [
            gate("phi_lock", layer="phi", evidence="lock distance 0.001"),
        ]
        with self.assertRaises(PhiAsEvidence) as caught:
            audit("CONFIRMED TECHNOSIGNATURE CANDIDATE", gates,
                  supporting=("phi_lock",))
        self.assertIn("only ever detract", str(caught.exception))

    def test_a_phi_gate_may_exist_and_change_nothing(self):
        locked = [gate("detected"),
                  gate("phi_lock", layer="phi",
                       evidence="lock distance 0.0001 at rung 3")]
        unlocked = [gate("detected"),
                    gate("phi_lock", status="FAIL", layer="phi",
                         evidence="lock distance 0.42 at rung 3")]
        self.assertEqual(audit("UNRESOLVED", locked).permitted,
                         audit("UNRESOLVED", unlocked).permitted)

    def test_citing_an_absent_gate_is_refused(self):
        with self.assertRaises(EvidenceError):
            audit("UNRESOLVED", [gate("detected")], supporting=("phi_lock",))


class ReportTests(unittest.TestCase):
    def test_the_record_carries_the_deficits_and_the_gates(self):
        result = audit("ASTROPHYSICAL", [gate("detected")])
        payload = result.as_dict()
        self.assertEqual(payload["claimed"], "ASTROPHYSICAL")
        self.assertEqual(payload["permitted"], "UNRESOLVED")
        self.assertTrue(payload["downgraded"])
        self.assertFalse(payload["escalation"])
        self.assertTrue(payload["deficits"])
        self.assertEqual(payload["gates"][0]["id"], "detected")

    def test_gates_are_looked_up_by_id(self):
        result = audit("UNRESOLVED", [gate("detected")])
        self.assertEqual(result.status("detected"), "PASS")
        self.assertEqual(result.status("nonsense"), "NOT RUN")
        self.assertIsNone(result.gate("nonsense"))

    def test_the_audit_is_frozen(self):
        result = audit("UNRESOLVED", [gate("detected")])
        self.assertIsInstance(result, Audit)
        with self.assertRaises(Exception):
            result.permitted = "ASTROPHYSICAL"


class PipelineTests(unittest.TestCase):
    """The pipeline's exit code, against a doctored copy of the records."""

    @staticmethod
    def pipeline():
        import importlib.util
        import os

        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        spec = importlib.util.spec_from_file_location(
            "omni_pipeline", os.path.join(root, "scripts", "omni_pipeline.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def doctored(self, module, disposition):
        real_load = module.load

        def load(name):
            payload = real_load(name)
            if name != "signal_records.json":
                return payload
            records = [dict(r) for r in payload["records"]]
            # doctor one record whatever it currently is: the fixture must not
            # depend on which disposition the data happens to carry this week
            records[0] = {**records[0], "disposition": disposition}
            return {**payload, "records": records}

        module.load = load
        module.write_outputs = lambda *args, **kwargs: None
        return module

    def test_a_claim_every_gate_supports_still_runs_clean(self):
        module = self.doctored(self.pipeline(), "ASTROPHYSICAL")
        self.assertEqual(module.run(strict=False), 0)

    def test_an_unsupported_confirmation_fails_the_run(self):
        module = self.doctored(self.pipeline(),
                               "CONFIRMED TECHNOSIGNATURE CANDIDATE")
        self.assertEqual(module.run(strict=False), 1)

    def test_a_lowered_midlevel_claim_is_a_report_not_a_failure(self):
        module = self.doctored(self.pipeline(), "ASTROPHYSICAL")
        self.assertEqual(module.run(strict=False), 0)
        self.assertEqual(module.run(strict=True), 1)


if __name__ == "__main__":
    unittest.main()
