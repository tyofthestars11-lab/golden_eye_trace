# The omni pipeline

**9 of 10 records keep the disposition their run gave them. 0 claim a confirmation, and 0 of those would be unsupported.**

Every other script here answers one question. This one orders the questions - measurement, calibration, nonlinear time, phi, interpretation, confirmation - and then checks the claimed disposition against the gates underneath it. It cannot raise a disposition; it can leave one alone or lower it.

## The window guard

The survey's window is 96 channels by 1000 samples, and the stabilizer asks whether that window can hold the trail a candidate's drift lays down. At this geometry the honest drift range is 58.9 kHz/s, while the dedoppler evaluates trials out to the fold at 196.4 kHz/s: a rate between those two is searched, found, and measured through a window that loses part of it.

| date | MHz | cycles/sample | trail | held | amplitude lost | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| 20240207 | 53.7833 | +0.0177 | 18 | 100.0% | 0.0% | HOLDS |
| 20240322 | 47.2416 | -0.0004 | 0 | 100.0% | 0.0% | HOLDS |
| 20240324 | 47.3143 | -0.0013 | 1 | 100.0% | 0.0% | HOLDS |
| 20240325 | 45.9843 | +0.1049 | 105 | 91.5% | 4.3% | HOLDS |
| 20240325 | 46.9828 | +0.1049 | 105 | 91.5% | 4.3% | HOLDS |
| 20240325 | 47.9380 | +0.1049 | 105 | 91.6% | 4.3% | HOLDS |
| 20240325 | 48.9937 | +0.1048 | 105 | 91.6% | 4.3% | HOLDS |
| 20240920 | 45.4206 | +0.4965 | 497 | 19.3% | 56.0% | TRUNCATES |
| 20240920 | 46.5890 | +0.4962 | 496 | 19.3% | 56.0% | TRUNCATES |
| 20240920 | 47.5358 | +0.4962 | 496 | 19.3% | 56.0% | TRUNCATES |

The 20240325 group holds at 91.5% of its trail, so its sigma carries a 4.3% truncation footnote and no more. The 20240920 trio sits at 0.4965 cycles per sample, holds 19% of its trail, and is the set this geometry is least able to measure - which is another reason its strength says nothing about its origin.

The measured cost of that loss, as drift rises past the window's limit, is in DRIFT_SENSITIVITY.md: the guard's sqrt(held) form overstates the loss, so a window it passes keeps at least the significance it claims.

## The gate matrix

`.` pass, `F` fail, `U` untestable with this product, `-` not run.

| date | MHz | window | detected | alias_signature | noise_reference | nonlinear_beats_control | phi_lock | astrophysical_explanation | off_source | independent_telescope | barycentric | repeat | claimed | permitted |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20240207 | 53.7833 | . | . | F | . | F | F | F | U | U | U | F | UNRESOLVED | UNRESOLVED |
| 20240322 | 47.2416 | . | F | F | . | F | F | F | U | U | U | F | UNRESOLVED | NULL RESULT |
| 20240324 | 47.3143 | . | . | F | . | F | F | F | U | U | U | F | UNRESOLVED | UNRESOLVED |
| 20240325 | 45.9843 | . | . | F | . | F | F | . | U | U | U | F | ASTROPHYSICAL | ASTROPHYSICAL |
| 20240325 | 46.9828 | . | . | F | . | F | F | . | U | U | U | F | ASTROPHYSICAL | ASTROPHYSICAL |
| 20240325 | 47.9380 | . | . | F | . | F | F | . | U | U | U | F | ASTROPHYSICAL | ASTROPHYSICAL |
| 20240325 | 48.9937 | . | . | F | . | F | F | . | U | U | U | F | ASTROPHYSICAL | ASTROPHYSICAL |
| 20240920 | 45.4206 | F | . | . | . | F | F | F | U | U | U | F | INSTRUMENTAL ALIAS | INSTRUMENTAL ALIAS |
| 20240920 | 46.5890 | F | . | . | . | F | F | F | U | U | U | F | INSTRUMENTAL ALIAS | INSTRUMENTAL ALIAS |
| 20240920 | 47.5358 | F | . | . | . | F | F | F | U | U | U | F | INSTRUMENTAL ALIAS | INSTRUMENTAL ALIAS |

## Where the claims were lowered

**20240322 47.2416 MHz**: UNRESOLVED -> NULL RESULT.
- detected: FAIL (5.16 sigma at the reported cell against a 5.5 sigma threshold)


## phi

The golden-ratio gate is computed for every record and is not allowed to support anything: it is not in any disposition's requirement list, and citing it raises an error. No record's disposition depends on whether its rung is close to an integer.

## What would raise a disposition

Nothing here can be repaired by more analysis of this data. A confirmation needs an off-source scan, a second telescope and a barycentric correction, and the product has one disk-integrated beam, one instrument and no world coordinate system. Those three gates are `UNTESTABLE` rather than `FAIL`, which is the difference between a method that found nothing and a method that was never run.
