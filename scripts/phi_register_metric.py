"""Register the Z[phi] metric as revision 3 of the two-band registration.

Revision 2 replaced the drift-blind band. This one fixes the *measurement* of
lattice agreement, before any drift slope is computed, because a metric chosen
after seeing which one is kinder is not a test:

* the lattice is Z[phi]: addresses are multiples of ``phi**-2 = 2 - phi``, one
  golden angle, in turns for sky addresses and in rungs for frequencies;
* the distance is additive in rungs, ``|r - n*(2-phi)|``, because the rung
  coordinate is already logarithmic and distances there are differences;
* the null is ``2*d/s`` with ``s = 2 - phi``;
* the multiplicative alternative ``|log_phi(r/(n*s))|`` is recorded with its own
  null ``2*d*ln(phi)*r/s``, since without that factor it inflates significance by
  the address scale - measured, the two agree at p = 0.213 against 0.214;
* multiple comparisons are corrected as ``1 - (1-p)**N`` over the number of
  quantities tested, and the calibration run in ``prereg/lattice_scale.json`` is
  attached as evidence that the null is right.

    python3 scripts/lattice_scale.py && python3 scripts/phi_register_metric.py
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bands import PITCH_EXPONENT  # noqa: E402
from seti_search.rungs import PHI  # noqa: E402

R2_PATH = os.path.join(ROOT, "prereg", "phi_2band_r2.json")
SCALE_PATH = os.path.join(ROOT, "prereg", "lattice_scale.json")
OUT_PATH = os.path.join(ROOT, "prereg", "phi_2band_r3.json")


def sha256(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main() -> int:
    r2 = json.load(open(R2_PATH, encoding="utf-8"))
    scale = json.load(open(SCALE_PATH, encoding="utf-8")) if os.path.exists(SCALE_PATH) else None

    step = 2.0 - PHI
    revision = dict(r2)
    revision.update({
        "revision": 3,
        "supersedes": {
            "path": "prereg/phi_2band_r2.json",
            "sha256": sha256(R2_PATH),
            "also_supersedes": r2.get("supersedes"),
        },
        "date": "2026-11-11",
        "metric": {
            "lattice": "Z[phi]: addresses are multiples of phi**-2 = 2 - phi",
            "angle_step_turns": step,
            "angle_step_exact": "2 - phi",
            "rung_step": step,
            "frequency_ratio_per_step": float(PHI ** step),
            "distance": "additive in rungs: |r - n*(2 - phi)|; turns for sky",
            "null": "2*d/(2-phi)",
            "multiplicative_alternative": {
                "distance": "|log_phi(r / (n*(2-phi)))|",
                "null": "2*d*ln(phi)*r/(2-phi)",
                "measured_comparison": "0.213 additive against 0.214 multiplicative "
                                       "for the band separation; the apparent 0.30x "
                                       "improvement in distance is exactly cancelled",
            },
            "look_elsewhere": "1 - (1-p)**N over the quantities tested",
            "why_registered_now": "the metric cannot be chosen after seeing which "
                                  "reading favours the framework",
            "calibration": scale,
            "predictions_restated": {
                "flat_drift_ratio_high_over_low": r2["golden_axis"][
                    "predicted_drift_ratio_flat"],
                "axis_drift_ratio_high_over_low": r2["golden_axis"][
                    "predicted_drift_ratio_axis"],
            },
        },
    })

    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(revision, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print("REVISION 3: THE METRIC IS NOW REGISTERED")
    print(f"   supersedes r2 ({revision['supersedes']['sha256'][:16]}...)")
    print(f"   lattice step {step:.9f} = 2 - phi, ratio per step {PHI ** step:.6f}")
    print("   distance additive in rungs; null 2d/(2-phi); look-elsewhere "
          "1-(1-p)**N")
    if scale:
        print(f"   calibration: control pass {scale['uniform_positions']['pass_rate']:.4f} "
              f"(theory {scale['theory_pass_rate']:.4f})")
        print(f"   best real p {min(scale['real_frequencies']['best_p'], scale['real_sources']['best_p']):.4f} "
              f"vs best control p {min(scale['uniform_positions']['best_p'], scale['uniform_frequencies']['best_p']):.4f}")
    else:
        print("   WARNING: no calibration file; run scripts/lattice_scale.py first")
    print(f"   wrote {os.path.relpath(OUT_PATH, ROOT)}")
    print(f"   sha256 {sha256(OUT_PATH)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
