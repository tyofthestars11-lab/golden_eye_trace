# The neutrino lane

Framework: φ² = φ + 1. Attribution: Tyree Jones (tyofthestarz), prior art from December 2025.

A neutrino is a second messenger: it leaves the source without being
absorbed the way radio is, and it arrives with no dispersion to speak of,
so a coincidence in time and angle is a coincidence between two
independent channels. That is what makes it worth asking about. It is
also why the answer has to be stated carefully, because a neutrino
detector is a poor instrument for ruling things out:

| limitation | consequence |
| --- | --- |
| energy threshold (~TeV and up) | a source that is radio-bright but
neutrino-quiet sits below it and is invisible |
| angular resolution (0.5-15 deg) | a one-degree cone contains many
square degrees of sky, so chance coincidences are normal |
| sky coverage (up-going vs down-going) | half the sky is seen through
the Earth, at a different effective area |
| duty cycle and livetime | an event can be missed while the detector is
not taking data |
| flux | a single source may simply not produce a detectable event in
the exposure |

Every one of those is a reason that **an absent neutrino is not evidence
that a radio feature is artificial**, and this module refuses to represent
it otherwise: `Coincidence.absence_is_informative` is a field fixed at
`False`, and an empty query must name the limit that makes it empty.

## The rule that outranks the rest

A neutrino coincidence cannot restore an invalid radio candidate. If the
radio pipeline found that a signal was manufactured by its own mask, then
the radio event does not exist, and a real particle that happens to point
there is a coincidence with an artifact. The run below calls the guard on
every invalid record and records 7 refusals, so
the rule is exercised rather than promised.

## The records

| date | MHz | disposition | status | chance p | limits |
| --- | --- | --- | --- | --- | --- |
| 20240207 | 53.7833 | UNRESOLVED | UNINFORMATIVE | no null available | the radio record has no pointing: a single-beam product with no WCS cannot be crossmatched in angle |
| 20240322 | 47.2416 | UNRESOLVED | UNINFORMATIVE | no null available | the radio record has no pointing: a single-beam product with no WCS cannot be crossmatched in angle |
| 20240324 | 47.3143 | UNRESOLVED | UNINFORMATIVE | no null available | the radio record has no pointing: a single-beam product with no WCS cannot be crossmatched in angle |
| 20240325 | 45.9843 | MASK-MANUFACTURED SIGNAL | EXCLUDED_FROM_CROSSMATCH | no null available | the radio event was manufactured by the excision mask; a neutrino coincidence is not evidence about an artifact, so the lane does not run on it |
| 20240325 | 46.9828 | MASK-MANUFACTURED SIGNAL | EXCLUDED_FROM_CROSSMATCH | no null available | the radio event was manufactured by the excision mask; a neutrino coincidence is not evidence about an artifact, so the lane does not run on it |
| 20240325 | 47.9380 | MASK-MANUFACTURED SIGNAL | EXCLUDED_FROM_CROSSMATCH | no null available | the radio event was manufactured by the excision mask; a neutrino coincidence is not evidence about an artifact, so the lane does not run on it |
| 20240325 | 48.9937 | MASK-MANUFACTURED SIGNAL | EXCLUDED_FROM_CROSSMATCH | no null available | the radio event was manufactured by the excision mask; a neutrino coincidence is not evidence about an artifact, so the lane does not run on it |
| 20240920 | 45.4206 | MASK-MANUFACTURED SIGNAL | EXCLUDED_FROM_CROSSMATCH | no null available | the radio event was manufactured by the excision mask; a neutrino coincidence is not evidence about an artifact, so the lane does not run on it |
| 20240920 | 46.5890 | MASK-MANUFACTURED SIGNAL | EXCLUDED_FROM_CROSSMATCH | no null available | the radio event was manufactured by the excision mask; a neutrino coincidence is not evidence about an artifact, so the lane does not run on it |
| 20240920 | 47.5358 | MASK-MANUFACTURED SIGNAL | EXCLUDED_FROM_CROSSMATCH | no null available | the radio event was manufactured by the excision mask; a neutrino coincidence is not evidence about an artifact, so the lane does not run on it |

## The radio-time gate

A coincidence needs a clock and a cone. Each record is graded against
both, and a coincidence may not be reported for one that is not
testable - `enforce` raises instead of writing the row.

| date | MHz | radio_time_status | multimessenger_status | window |
| --- | --- | --- | --- | --- |
| 20240207 | 53.7833 | MISSING | NOT TESTABLE | - |
| 20240322 | 47.2416 | MISSING | NOT TESTABLE | - |
| 20240324 | 47.3143 | MISSING | NOT TESTABLE | - |
| 20240325 | 45.9843 | MISSING | NOT TESTABLE | - |
| 20240325 | 46.9828 | MISSING | NOT TESTABLE | - |
| 20240325 | 47.9380 | MISSING | NOT TESTABLE | - |
| 20240325 | 48.9937 | MISSING | NOT TESTABLE | - |
| 20240920 | 45.4206 | MISSING | NOT TESTABLE | - |
| 20240920 | 46.5890 | MISSING | NOT TESTABLE | - |
| 20240920 | 47.5358 | MISSING | NOT TESTABLE | - |

**0 of 10 records are testable.** The verdict is `MISSING` rather than `PARTIAL` because a
date is not a clock, and a 24 h window is not free: it carries ~4.6e-4
chance probability in a 4.68 deg cone, about 300x the 256 s integration
the survey actually used.

## What the answer currently is

**UNINFORMATIVE, and the reason is the product, not the sky.** These
records have no world coordinate system and no UTC timestamp: they are a
single disk-integrated beam, named by date. There is therefore no position
to put a cone around and no time to draw a window on, and the lane reports
that instead of reporting a null. This is the same limitation that makes
the off-source and barycentric gates untestable, and it is not repaired by
having a catalog.

What would make it informative: a beam-formed or targeted follow-up with a
pointing, a UTC start time, and a localisation of a degree or better. Then
the same function decides, and the number it reports is a chance
probability against the stated exposure.

## Worked example, hypothetical geometry

- catalog: synthetic (data/neutrino_events_synthetic.json), used only for this example
- radio position: RA 100.00 deg, Dec -20.00 deg (**invented to sit inside the synthetic cone**)
- radio time: 2024-03-22T08:15:00+00:00, +-256 s
- exposure: 100 events over 3.156e+07 s and the whole sky, i.e. 2.521e-07 per sr per s (order-of-magnitude alert rate, stated for the worked example only)

Nearest event: `SYN-0001` at 0.48 deg and 4 s, inside 4.69 deg and 768 s.

Chance probability: **8.13e-06**. Read it as the probability that one random event from this exposure lands in that
box, which is what makes it small: the box is small in solid angle,
not because the event is remarkable.

The example is not a candidate and is not written to the record list. It
demonstrates the arithmetic, and it shows the asymmetry: a coincidence is
an invitation to point a telescope, never a disposition.

