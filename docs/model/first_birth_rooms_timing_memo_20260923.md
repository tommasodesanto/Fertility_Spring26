# First-birth rooms: timing correction and provisional 0.6-room decision

Date: 2026-09-23. Author decision: record **0.60 rooms** (two-year window +3/+4
versus −3/−2, corrected dates) as the provisional first-birth housing response,
replacing the 0.770 placeholder. Not yet a calibration-contract change: the
frozen target, weight and fingerprint are untouched until the horizon is set on
model grounds and the target contract is renamed.

## What was wrong

In the merged panel `PSID/PSIDSHELF_MOBILITY.dta`, `ACTUALROOMS_` on the row
for interview year \(t\) holds the answer given at the **next** interview
(\(t+1\) through 1997, \(t+2\) after). Every other variable on the row is
contemporaneous: on 512 persons, the wave-specific family ID and reported age
match their own wave 100% and the next wave 0% among changers; rooms match the
next wave 100% and their own wave 0%. The full-panel census (September 12) has
811,486/811,486 matches after shifting. Rooms is an add-on merged from
`mobility_long_withadd.dta`, not part of the official PSID-SHELF build; the
construction code was never found. Ownership (`HOMEOWN`) is built by the
official pipeline with year-correct sources and is **not** affected. The moving
variables from the same add-on file are unverified; do not use them on either
clock without a check.

Consequence: the original "+3" coefficient (0.770) read rooms four to five years
after the birth against rooms in the year before or the year of the birth.

## Designs run tonight (Torch 18388849 smoke, 18388850 array, all pass)

Author's original Sun–Abraham specification throughout: all adults 18+,
unweighted, ID and year fixed effects, age and education covariates, ID
clustering, last first-birth cohort as `control_cohort`, people without a
recorded first birth retained, original room codes. Only the room date
assignment and the design change. Estimator
`code/data/psid_followup_mar2026/audit_rooms_window_designs.do`; results in
`code/data/psid_followup_mar2026/output/first_birth_correction_review/window_designs/`.

| Design, corrected dates | Rows / persons | Horizon | Rooms (SE) | Original dates |
|---|---|---|---|---|
| Two-year windows, baseline −3/−2 | 315,737 / 34,191 | +1/+2 | 0.58 (0.12) | 0.58 |
| | | **+3/+4** | **0.60 (0.13)** | 0.79 |
| | | +5/+6 | 1.07 (0.13) | 1.02 |
| Annual, cohorts observed at −2 only | 285,155 / 30,430 | +2 | 0.58 (0.18) | 0.58 |
| | | +3 | 0.66 (0.18) | 0.40 |
| | | +4 | 0.56 (0.18) | 0.80 |
| | | +6 | 1.13 (0.19) | 1.11 |

Exclusions: the window design drops first births 1968–1971 and 1978 (no
interview two or more years before the birth in the common sample); the
−2-only design additionally drops even birth years 2000–2018, which after 1997
are observed only at odd event times. From cohort 1999 on, no cohort observes
both −2 and +3, so an annual "−2 to +3" contrast is an annual-era object.

Findings: the corrected +3 point in the full annual fit (0.40) was a
normalization artifact from cohorts lacking a −2 interview; with those removed
the curve is smooth (0.40, 0.58, 0.66, 0.56, 0.95, 1.13 at +1…+6). The
short-run response clusters at 0.6 rooms; the completed response is 1.0–1.1
rooms. The original 0.77 cannot be recovered at any clean horizon: it blended
true +4 and true +5 readings.

## If 0.6 is not satisfactory: prepared follow-ups, none launched

1. Horizon on model grounds: if one model period after birth maps to a shorter
   or longer calendar span, read the window curve above (+1/+2 = 0.58,
   +5/+6 = 1.07); no new regression needed.
2. Rebuild rooms from the 41 raw year-specific variables via the saved
   crosswalk (recovers 1969 and 1976 answers, ~9k rows; restores 1977–1978
   references).
3. Controlled sequence toward the August household specification (0.720):
   add −6 dummy; explicit confirmed-childless controls; household unit (women
   ref/spouse, one per household-year, single-FU dwellings); biological birth
   history and room-code cleaning; IW weights. One change per fit, six fits,
   about one Torch array. Isolates which of the six differences moves the
   number and where the household pre-birth trend comes from.
4. Same window design with a −2/−1 baseline (four-year span for every cohort,
   anticipation year inside the reference) as a sensitivity row.

## Not changed

Frozen target 0.720246 and its weight; slides; paper text. Bundle for the
run: `/tmp/psid_rooms_window_20260923/task_root` (private; not in Git).
