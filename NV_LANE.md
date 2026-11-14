# The diamond NV-centre lane

Framework: φ² = φ + 1. Attribution: Tyree Jones (tyofthestarz), prior art from December 2025.

A nitrogen-vacancy ensemble in diamond measures magnetic field, and a
radio wave carries one: for a plane wave E = cB, so 1 pT is 0.3 mV/m. A
millimetre of diamond at the telescope is therefore a **local field
witness** - it can state that a magnetic transient happened at the site
at the same frequency and second as a radio candidate.

## The asymmetry, which is the whole point

| outcome | what it means | allowed to |
| --- | --- | --- |
| a local transient in band and in window | the simplest explanation is a
local emitter | **detract** |
| no local transient | compatible with a real emitter polarised away from
the NV axis, or coupled to the antenna and not the diamond | nothing |
| no sensor, no band, no window | the question was not asked | nothing |

`guard_support` raises `NVAsEvidence`, so the lane cannot be cited as
support even by accident, and the OMNI gate `local_field` is required by
no disposition. This is the same rule the phi layer carries, for the same
reason: a check that can only ever remove belief must never add it.

## Feasibility: can a diamond see 13.4-86.9 MHz at all?

The NV ground-state splitting is 2.87 GHz and tunes at 2.8 MHz/G, so the
band is reachable only with a substantial bias field:

| frequency | bias field needed |
| --- | --- |
| 13.40 MHz | 102.0 mT |
| 45.98 MHz | 100.9 mT |
| 47.24 MHz | 100.8 mT |
| 53.78 MHz | 100.6 mT |
| 86.90 MHz | 99.4 mT |

Sensitivity per technique, with the bandwidth that makes the number
mean something (threshold = sensitivity x sqrt(bandwidth)):

| mode | sensitivity | as E-field | bandwidth | LWA band covered | threshold | SNR for 1 V/m |
| --- | --- | --- | --- | --- | --- | --- |
| cw-odmr | 1.0e-09 T/rtHz | 299.79 mV/m/rtHz | 1.0e+03 Hz | 0% | 9480.27 mV/m | 0.1 |
| ramsey | 5.0e-10 T/rtHz | 149.90 mV/m/rtHz | 1.0e+04 Hz | 0% | 14989.62 mV/m | 0.1 |
| hahn-echo | 2.1e-13 T/rtHz | 0.06 mV/m/rtHz | 1.0e+05 Hz | 0% | 19.91 mV/m | 50.2 |
| cpmg | 2.3e-10 T/rtHz | 68.95 mV/m/rtHz | 1.0e+05 Hz | 0% | 21804.62 mV/m | 0.0 |
| rf-dressed | 1.0e-09 T/rtHz | 299.79 mV/m/rtHz | 7.5e+06 Hz | 100% | 821015.46 mV/m | 0.0 |
| quantum-heterodyne | 1.0e-12 T/rtHz | 0.30 mV/m/rtHz | 1.0e+06 Hz | 100% | 299.79 mV/m | 3.3 |

Read the two columns together and the conclusion is narrower than the
sensitivity figures suggest, because band and sensitivity are mutually
exclusive with present hardware. The sequences that reach pT sensitivity
(echo, CPMG) stop at a few MHz and cover **none** of the survey band; the
technique that does reach the band (RF-dressed CW-ODMR) is nT-sensitive
over 7.5 MHz, so its threshold is ~0.8 V/m. A diamond can therefore see
**strong local RFI** across 13.4-86.9 MHz, or pT-level fields in narrow
slices below 10 MHz, but not both. The 20 mV/m row with 0% coverage is
the trap: the impressive number belongs to a band the survey never
observes.

## This survey's records

| date | MHz | verdict | reason |
| --- | --- | --- | --- |
| 20240207 | 53.7833 | NOT TESTABLE | no NV sensor record was supplied |
| 20240322 | 47.2416 | NOT TESTABLE | no NV sensor record was supplied |
| 20240324 | 47.3143 | NOT TESTABLE | no NV sensor record was supplied |
| 20240325 | 45.9843 | NOT TESTABLE | no NV sensor record was supplied |
| 20240325 | 46.9828 | NOT TESTABLE | no NV sensor record was supplied |
| 20240325 | 47.9380 | NOT TESTABLE | no NV sensor record was supplied |
| 20240325 | 48.9937 | NOT TESTABLE | no NV sensor record was supplied |
| 20240920 | 45.4206 | NOT TESTABLE | no NV sensor record was supplied |
| 20240920 | 46.5890 | NOT TESTABLE | no NV sensor record was supplied |
| 20240920 | 47.5358 | NOT TESTABLE | no NV sensor record was supplied |

**No NV sensor data exists for these observations**, so every record is
`NOT TESTABLE` and none is detracted. That is the current state and not a
result: the lane needs a sensor at the site, running during the
observation, with a bias field and a calibration - three things this
product does not have on disk.

## What integration would require

1. A sensor record per night: band, sensitivity, bandwidth, bias field,
   axis orientation, and a UTC start with an uncertainty. The schema
   refuses a record missing any of them, because each is a way for a null
   to look like evidence.
2. A calibration against a known local emitter, so the sensitivity is a
   measurement and not a datasheet number.
3. A stated coupling model: the NV axis responds to one field component,
   so an unpolarised local emitter can be invisible to it. Without that
   model a null cannot even be interpreted as a weak null.

Until then this lane subtracts nothing and adds nothing, which is the
correct behaviour for an instrument that does not exist yet.

