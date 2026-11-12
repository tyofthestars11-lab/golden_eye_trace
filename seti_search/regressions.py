"""Regression fixtures for methods that were tested and rejected.

Each entry is a record of a run that produced a wrong answer, why it was wrong,
and what replaced it. The point is not sentiment: a rejected method that is
deleted will be reinvented, and the reinvented version will look plausible for
the same reason the first one did.
"""

from __future__ import annotations

REGRESSIONS = (
    {
        "id": "narrow-window-truncation",
        "date": "2026-11-12",
        "method": "32-channel search windows, local noise reference",
        "result": "0 of 10 candidates survive, including a 954 sigma burst group",
        "verdict": "INVALID",
        "why": (
            "The 20240325 lanes drift about 0.105 channels per sample, which is "
            "105 channels across the 1000-sample integration. A 32-channel window "
            "holds a third of that trail, so the dedoppler sum is truncated "
            "mid-flight and the significance collapses. The 954 sigma line was "
            "not lost because the noise reference was local; it was lost because "
            "the window cut the signal in half."
        ),
        "tell": "A 954 sigma line reporting 0 sigma is a statement about the window.",
        "superseded_by": (
            "per-channel local scale on the survey's own 96-channel window, "
            "plus a pre-flight guard (seti_search/stabilizer.py) that refuses a "
            "window which cannot hold the trail"
        ),
        "superseded_result": "9 of 10 survive",
        "fixture": "tests/test_window_truncation.py",
        "artifact": "prereg/regression_narrow_window.json",
    },
    {
        "id": "per-trial-calibration",
        "date": "2026-11-12",
        "method": "calibrating each drift trial of the map by its own spread",
        "result": (
            "a 32-channel window reported 174 sigma for a lane the 96-channel "
            "window reported at 135"
        ),
        "verdict": "INVALID",
        "why": (
            "A trial whose accumulated shift leaves the band sums fewer "
            "integrations, so its spread really is smaller - and dividing each "
            "trial by its own spread normalises away exactly the sensitivity loss "
            "truncation causes. The truncated search then outscored the complete "
            "one, which inverts the quantity being measured."
        ),
        "tell": (
            "The narrow window scoring higher than the wide window on the same "
            "synthetic lane."
        ),
        "superseded_by": "one spread for the whole map, in calibrate_map()",
        "superseded_result": "narrow below 0.8x the wide, as the fixture requires",
        "fixture": "tests/test_window_truncation.py",
        "artifact": "prereg/regression_narrow_window.json",
    },
    {
        "id": "chi2-lock-k0",
        "date": "2026-11-11",
        "method": "golden-ratio lock distance scored over all integer k",
        "result": "p = 0.0265 from a 172 Hz/s line",
        "verdict": "INVALID",
        "why": (
            "k = 0 is the stationary lattice point. Including it lets a near-zero "
            "drift score an arbitrarily small lock distance, because distance to "
            "k = 0 shrinks with the rate itself. The registered criterion excludes "
            "k = 0: a stationary line is not a phase advance."
        ),
        "tell": "A 0.0004-cycle drift scoring closer to the lattice than 0.05.",
        "superseded_by": "k = max(1, round(nu/theta)), recorded in prereg/phi_lock.json",
        "superseded_result": "p = 1.0000",
        "fixture": "scripts/phi_echo.py, scripts/phi_lock_prereg.py",
        "artifact": "prereg/phi_lock.json",
    },
    {
        "id": "alias-sign",
        "date": "2026-11-12",
        "method": "aliased injection recovery, expected rate +0.4 for an injected 0.6",
        "result": "0 of 2 aliased injections recovered",
        "verdict": "INVALID",
        "why": (
            "The fold changes the sign as well as the magnitude: 0.6 cycles per "
            "sample is recorded as -0.4, not +0.4. Comparing against +0.4 missed "
            "a signal that had been found."
        ),
        "tell": "A recovered alias at -157 kHz/s scored against an expectation of +157.",
        "superseded_by": "fold() in scripts/ovro_confirm.py, verified exactly in scripts/nonlinear_time.py",
        "superseded_result": "fold exact at all six probed rates",
        "fixture": "scripts/nonlinear_time.py:alias_probe",
        "artifact": "nonlinear_candidates.json",
    },
    {
        "id": "real-data-injection",
        "date": "2026-11-12",
        "method": "completeness measured by injecting into the real solar spectra",
        "result": "0 of 8 recovered",
        "verdict": "INVALID",
        "why": (
            "The solar spectra carry bandpass structure that the excision and "
            "normalisation stages act on, so the injection measured those stages "
            "rather than the detection pipeline. Completeness is a property of the "
            "search, so it has to be measured on data the search is responsible for."
        ),
        "tell": "Every rate failing equally, including a stationary line at 14 sigma elsewhere.",
        "superseded_by": "synthetic unit-variance exponential noise",
        "superseded_result": "6 of 8 recovered, 0 of 3 false positives",
        "fixture": "scripts/ovro_confirm.py:inject_synthetic",
        "artifact": "prereg/ovro_confirm.json",
    },
    {
        "id": "shuffled-null-for-amplitude",
        "date": "2026-11-12",
        "method": "shuffled-time nulls as the primary control for nonlinear coherence",
        "result": "null percentile 1.000 for every candidate and coordinate",
        "verdict": "INVALID as a primary control",
        "why": (
            "A permutation keeps the sample values, so an amplitude statistic "
            "cannot move: the null sits on the observation rather than below it. "
            "Measured worst-case deviation across twelve candidates is 1.8e-15."
        ),
        "tell": "Every candidate passing every null at the 100th percentile.",
        "superseded_by": "best-of-five randomised monotone warps, best-of-five against best-of-five",
        "superseded_result": "0 of 12 beat the warp-gain control",
        "fixture": "scripts/nonlinear_time.py",
        "artifact": "nonlinear_candidates.json",
    },
    {
        "id": "centroid-chunks",
        "date": "2026-11-12",
        "method": "line centroid fitted in eight time chunks to measure drift",
        "result": "the two days came out swapped: 47.304 for the 47.242 line and 47.232 for the 47.314",
        "verdict": "INVALID",
        "why": (
            "At roughly 3 sigma per chunk the centroid is driven by whichever "
            "channel is brightest in the local window, not by the candidate. The "
            "recovered centroids were two different features."
        ),
        "tell": "Each day recovering the other day's frequency.",
        "superseded_by": "full-integration dedoppler plus a within-day window scan",
        "superseded_result": "neither line persists; verdict NOT A REPEAT",
        "fixture": "scripts/rfi_pair_test.py",
        "artifact": "prereg/rfi_pair.json",
    },
    {
        "id": "band-wide-noise-reference",
        "date": "2026-11-12",
        "method": "significance against one MAD across the whole 9.66 MHz window",
        "result": "significances wrong by factors from 0.03 to 3.0",
        "verdict": "SUPERSEDED, not invalid",
        "why": (
            "The statistic is not wrong, it is mis-scaled: a line is measured "
            "against a band average that its own neighbourhood may beat. The local "
            "scale spans 0.30 to 31.18 of the band-wide value across the ten "
            "candidates, so a fixed sigma threshold was not a constant false-alarm "
            "rate across the band."
        ),
        "tell": "A 6.49 sigma line and an 8.66 sigma line whose neighbourhoods differ by 2x.",
        "superseded_by": "per-channel local scale",
        "superseded_result": "9 of 10 survive; the survivor set changes by one",
        "fixture": "scripts/local_noise_survey.py",
        "artifact": "prereg/local_noise.json",
    },
    {
        "id": "phi-as-evidence",
        "date": "2026-11-11",
        "method": "golden-ratio phase as a detection criterion",
        "result": "no lock survives any control",
        "verdict": "REJECTED as evidence; retained as a preregistered feature",
        "why": (
            "Fractional rungs of a smooth population are uniform, so a small drift "
            "is not a finding. The population test gives adjusted p of 0.455 and "
            "0.935; the echo harness orders theta above an arbitrary rate but its "
            "own null overlaps. phi may only ever detract from a candidate."
        ),
        "tell": "A lock distance small enough to matter appearing in a random draw at the same rate.",
        "superseded_by": "feature only, with a randomized-phase control beside it",
        "superseded_result": "p = 1.0000",
        "fixture": "tests/test_rungs.py",
        "artifact": "prereg/phi_lock.json",
    },
    {
        "id": "window-dependent-visibility",
        "date": "2026-11-12",
        "method": "reading a candidate's sigma, or its existence, from one window width",
        "result": (
            "the 20240325 four-line group scores 954 sigma in a 96-channel window "
            "and 14 sigma in the 731-channel band, same data, same threshold, same "
            "1000 integrations; at full band the wide-band carrier test blanks 208 "
            "of 731 channels and the group's channels are among them, so the "
            "expanded search does not find it at all"
        ),
        "verdict": "INVALID as a verdict either way",
        "why": (
            "Two gates depend on the window's extent, not on the sky. calibrate_map "
            "measures one spread for the whole map, and the spread of a 73.5 MHz band "
            "full of low-frequency interference is 1.965 where the 96-channel window's "
            "is 0.377 - a factor 5.2. Separately, the carrier test compares each "
            "channel against the band it sits in, so widening the band changes which "
            "channels are blanked. A candidate that disappears when the band widens "
            "has been masked, not refuted, and one that brightens has not been "
            "confirmed."
        ),
        "tell": "The same cell scoring 954 sigma and 14 sigma in two runs of the same pipeline.",
        "superseded_by": (
            "the preregistered 2x2 (96 or 731 channels) x (1000 or all integrations) "
            "with the window guard printed beside every result, and "
            "scripts/band_calibration_check.py recording both gates per candidate"
        ),
        "superseded_result": (
            "variant B holds the drift window at 100% and finds no group in "
            "45.4-54.5 MHz; the group is reported as masked, and the label for it "
            "stays INDETERMINATE DUE TO LIMITED SPAN rather than becoming a null"
        ),
        "fixture": "scripts/band_calibration_check.py",
        "artifact": "prereg/band_calibration.json",
    },
    {
        "id": "theory-scaled-row-test",
        "date": "2026-11-12",
        "method": "row excision at sigma * 1.25/sqrt(n_f) * typical_std",
        "result": (
            "505 of 1000 integrations blanked on 20240325 and 228 on 20240920; the "
            "four-line group reads 954 sigma through that mask and 1.44 sigma with no "
            "mask at all, at the same cells and the same drift"
        ),
        "verdict": "INVALID: the mask manufactures the candidates",
        "why": (
            "1.25/sqrt(n_f) is the standard error of a median only if the channels are "
            "independent. In a solar dynamic spectrum one source spans the band, so the "
            "channel gains are correlated, the row medians scatter far more, and a "
            "two-sided 6 sigma threshold flags ordinary integrations. It removes both "
            "bright and faint rows, so the surviving subset is selected on the same "
            "quantity the dedoppler then measures, and the sum over that subset is a "
            "much quieter statistic. The same data, the same drift and the same "
            "calibration read 1.44 sigma without the mask."
        ),
        "tell": "A sigma that shrinks to 1.5 when nothing but the excision threshold changes.",
        "superseded_by": (
            "row_noise=\"empirical\", which scales the same test by "
            "1.4826 * MAD(row medians); scripts/light_path.py traces every stage with "
            "all three masks and results/light_path_ledger.csv carries the numbers"
        ),
        "superseded_result": (
            "0 of 1000 rows blanked on 20240325 and 1 on 20240920; no cell in the 2x2 "
            "reaches 1.5 sigma, so the group is not a candidate"
        ),
        "fixture": "scripts/light_path.py",
        "artifact": "results/light_path_stages.json",
    },
)
