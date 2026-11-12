"""The omni layer: one disposition per candidate, audited against its gates.

The framework's six layers each produce a verdict, and the failure this module
exists to prevent is those verdicts being quietly merged into a stronger one.
A "CONFIRMED TECHNOSIGNATURE CANDIDATE" is not a stronger "ASTROPHYSICAL": it
needs an off-source scan, a second telescope and a repeat, none of which any
single-beam product can supply. So the pipeline does not derive a disposition; it
checks the one it was handed and reports the strongest disposition the gates
actually permit.

Three rules are enforced here rather than trusted:

1. **A claim needs its gates.** Every disposition names the gates it requires,
   and a missing or non-PASS gate is a deficit, not a rounding error.
2. **A pass needs evidence.** A gate with an empty evidence string is rejected:
   otherwise "PASS" becomes an assertion with no number behind it.
3. **phi may not support anything.** A gate in the ``phi`` layer is refused as
   support, and no disposition requires one. phi can be recorded, ranked and
   reported, and it may detract; it can never be the reason.

When a claim outruns its gates, ``audit`` falls back rather than fails: a
classification whose own evidence is missing drops to ``UNRESOLVED`` if the
detection itself still holds, and to ``NULL RESULT`` if it does not. It never
falls sideways into a different classification, because "we do not know what this
is" is the honest answer and "it is probably RFI" is a different claim.

An escalation to ``CONFIRMED TECHNOSIGNATURE CANDIDATE`` that survives no gates
is the one case the caller must not paper over, so ``Audit.escalation`` marks it
and the pipeline exits non-zero.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from .signal_record import DISPOSITIONS

__all__ = [
    "Audit",
    "EvidenceError",
    "Gate",
    "LAYERS",
    "PhiAsEvidence",
    "REQUIRED",
    "STATUSES",
    "audit",
]

#: The layers, in the order they are allowed to speak. Measurement first, the
#: interpretation last: a disposition is a statement about the last layer, so it
#: cannot be formed before the earlier ones have run.
LAYERS = (
    "measurement",
    "calibration",
    "timedomain",
    "phi",
    "interpretation",
    "confirmation",
)

STATUSES = ("PASS", "FAIL", "UNTESTABLE", "NOT RUN")

#: The layer that may only ever detract from a candidate.
DETRACTING_LAYERS = ("phi",)

#: What each disposition requires. Empty means the claim asserts nothing beyond
#: "this is what the run produced", which is why only the negative dispositions
#: are free.
REQUIRED = {
    "CONFIRMED TECHNOSIGNATURE CANDIDATE": (
        "detected",
        "window",
        "repeat",
        "off_source",
        "independent_telescope",
        "barycentric",
    ),
    "NONLINEAR-COHERENT BUT UNRESOLVED": ("detected", "nonlinear_beats_control"),
    "ASTROPHYSICAL": ("detected", "astrophysical_explanation"),
    "TERRESTRIAL RFI": ("detected",),
    "INSTRUMENTAL ALIAS": ("detected", "alias_signature"),
    "PIPELINE ARTIFACT": ("detected",),
    "UNRESOLVED": ("detected",),
    "NULL RESULT": (),
}

_UNKNOWN = set(REQUIRED) - set(DISPOSITIONS)
if _UNKNOWN:  # pragma: no cover - import-time guard, not a runtime path
    raise RuntimeError(f"REQUIRED names dispositions outside the vocabulary: {_UNKNOWN}")
if set(REQUIRED) != set(DISPOSITIONS):  # pragma: no cover
    raise RuntimeError("every disposition needs a gate list, including the empty one")


class EvidenceError(ValueError):
    """A gate set that cannot be audited at all: a bug, not a weak result."""


class PhiAsEvidence(EvidenceError):
    """phi was asked to support a claim, which is the one thing it may not do."""


@dataclass(frozen=True)
class Gate:
    """One layer's verdict on one candidate, with the number that produced it."""

    id: str
    layer: str
    status: str
    evidence: str

    def __post_init__(self) -> None:
        if not str(self.id).strip():
            raise EvidenceError("a gate needs an id")
        if self.layer not in LAYERS:
            raise EvidenceError(f"gate {self.id!r} has layer {self.layer!r}, "
                                f"not one of {LAYERS}")
        if self.status not in STATUSES:
            raise EvidenceError(f"gate {self.id!r} has status {self.status!r}, "
                                f"not one of {STATUSES}")
        if not str(self.evidence).strip():
            raise EvidenceError(f"gate {self.id!r} is {self.status} with no "
                                "evidence; a verdict without a number behind it "
                                "is an assertion")
        if self.status != "PASS" and self.layer in DETRACTING_LAYERS:
            # a failing phi gate is the normal case and needs no special care
            pass


@dataclass(frozen=True)
class Audit:
    """What was claimed, what the gates permit, and why they differ."""

    claimed: str
    permitted: str
    deficits: tuple[str, ...]
    gates: tuple[Gate, ...]

    @property
    def downgraded(self) -> bool:
        return self.permitted != self.claimed

    @property
    def escalation(self) -> bool:
        """A confirmation claimed and not permitted: the caller must report it."""
        return (self.claimed == "CONFIRMED TECHNOSIGNATURE CANDIDATE"
                and self.downgraded)

    def status(self, gate_id: str) -> str:
        for gate in self.gates:
            if gate.id == gate_id:
                return gate.status
        return "NOT RUN"

    def gate(self, gate_id: str) -> Gate | None:
        for gate in self.gates:
            if gate.id == gate_id:
                return gate
        return None

    def as_dict(self) -> dict:
        return {
            "claimed": self.claimed,
            "permitted": self.permitted,
            "downgraded": self.downgraded,
            "escalation": self.escalation,
            "deficits": list(self.deficits),
            "gates": [{"id": g.id, "layer": g.layer, "status": g.status,
                       "evidence": g.evidence} for g in self.gates],
        }


def _index(gates: Sequence[Gate]) -> dict[str, Gate]:
    seen: dict[str, Gate] = {}
    for gate in gates:
        if gate.id in seen:
            raise EvidenceError(f"two gates share the id {gate.id!r}; the second "
                                "would silently replace the first")
        seen[gate.id] = gate
    return seen


def audit(claimed: str, gates: Iterable[Gate],
          supporting: Iterable[str] = ()) -> Audit:
    """Check one claimed disposition against the gates that produced it.

    ``supporting`` names the gates the caller believes carried the claim. It is
    optional because the requirement list already decides that; it exists so a
    caller that cites phi gets an error instead of a silent upgrade.
    """
    if claimed not in DISPOSITIONS:
        raise EvidenceError(f"{claimed!r} is not a disposition; the vocabulary is "
                            f"closed: {DISPOSITIONS}")
    gates = tuple(gates)
    by_id = _index(gates)

    for gate_id in supporting:
        gate = by_id.get(gate_id)
        if gate is None:
            raise EvidenceError(f"claim cites {gate_id!r}, which no gate reported")
        if gate.layer in DETRACTING_LAYERS:
            raise PhiAsEvidence(
                f"claim cites the {gate.layer} gate {gate_id!r} as support; a "
                "golden-ratio feature may only ever detract from a candidate"
            )

    required = REQUIRED[claimed]
    for gate_id in required:
        gate = by_id.get(gate_id)
        if gate is not None and gate.layer in DETRACTING_LAYERS:  # pragma: no cover
            raise PhiAsEvidence(f"{claimed} requires {gate_id!r}, a {gate.layer} gate")

    deficits: list[str] = []
    for gate_id in required:
        gate = by_id.get(gate_id)
        if gate is None:
            deficits.append(f"{gate_id}: NOT RUN (no run reported this gate)")
        elif gate.status != "PASS":
            deficits.append(f"{gate_id}: {gate.status} ({gate.evidence})")

    if not deficits:
        return Audit(claimed=claimed, permitted=claimed, deficits=(),
                     gates=gates)

    detected = by_id.get("detected")
    if detected is not None and detected.status == "PASS":
        permitted = "UNRESOLVED"
    else:
        permitted = "NULL RESULT"
    return Audit(claimed=claimed, permitted=permitted,
                 deficits=tuple(deficits), gates=gates)
