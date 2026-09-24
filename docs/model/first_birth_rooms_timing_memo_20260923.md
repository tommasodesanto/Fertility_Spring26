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

## Overnight follow-ups (September 23–24): sequence to the household specification

Torch jobs 18391624/18391626, 18395474/18395475, 18396745/18396746; sixteen
fits, all pass, receipts and covariances verified. Extended frozen sample
`/tmp/psid_rooms_sequence_20260923/analysis_sample.dta` (rows and both room
columns asserted identical to the reference sample). Estimator
`code/data/psid_followup_mar2026/audit_rooms_sequence.do`; results and figure
in `code/data/psid_followup_mar2026/output/first_birth_correction_review/sequence_designs/`.
All arms: corrected room dates, two-year windows, headline = +3/+4 window
versus the −3/−2 baseline window, Sun–Abraham as in the original code.

### Finding 1: non-room codes must be cleaned everywhere

The PSID codes "don't know / not answered" as 9 (through 1984), 99 (1985–1993)
and 98/99 (1994 on). The original code, the September 12 timing comparison and
tonight's first pass left them as room counts. In the full common sample there
are 4,559 such rows (1.3%), 912 inside the −3 to +4 event window; each is a
~90-room outlier against a mean of six. With them set to missing:

| Arm | +3/+4 rooms (SE) | N |
|---|---|---|
| S0 original spec, corrected dates, codes retained | 0.602 (0.126) | 315,737 |
| **S0c same, codes cleaned** | **0.726 (0.033)** | 311,453 |
| B21c baseline −2/−1, window +2/+3, codes cleaned | 0.553 (0.025) | 317,403 |

The standard error falls four-fold because the outliers dominated the residual
variance. Cleaned full-sample path (S0c): −0.13 at −7/−6, −0.14 at −5/−4, 0 at
−3/−2, 0.18 at −1/0, 0.49 at +1/+2, 0.73 at +3/+4, 0.97 at +5/+6, 1.10, 1.24,
1.35 at ≥+11. Note the mild pre-birth slope (about 0.07 rooms per year over the
seven years before the baseline), now statistically visible.

**Recommendation for the morning:** replace the provisional 0.60 with
**0.73 (SE 0.03)**. This is the same design and sample as the 0.60; the only
change is treating the codebook's non-answers as missing, which is not a
modelling choice. The 0.60 should not be used.

### Finding 2: what separates the full-sample number from the August household number

One change per fit, codes cleaned throughout, corrected dates, −3/−2 baseline:

| Step | Change added | +3/+4 rooms (SE) | N | Never-treated rows |
|---|---|---|---|---|
| S0c | original specification | 0.73 (0.03) | 311,453 | 100,889 |
| S1c | comparison group = confirmed-childless only (unknown histories dropped) | 0.73 (0.03) | 285,205 | 74,641 |
| S2c | women who are current reference person or spouse | 0.89 (0.06) | 118,928 | 22,693 |
| S3c | single-family-unit dwellings, one woman per household-year | 0.90 (0.06) | 101,852 | 20,049 |
| S4ac | first birth from full biological-child history | 0.93 (0.06) | 101,160 | 20,049 |
| S4c | first-birth-after-first-observation rule (as in the August code) | 0.88 (0.07) | 60,178 | 20,049 |
| S6 | PSID longitudinal weights | 0.93 (0.08) | 41,419 | 14,796 |

The S6 endpoint is the August specification in window form and matches its
annual numbers (+3 = 1.20, +4 = 0.83 relative to −2; −1 = 0.48). So the
0.73-versus-0.93 gap is the unit (+0.16: women who already head a household
or are spouses, versus all adults including adult children living with
parents) and the survey weights (+0.06). The comparison group, the household
de-duplication and the birth-history definition each move it by 0.03 or less.
The August entry rule removes 40% of the household sample and lowers the
estimate by 0.06. With codes retained the same sequence looked very different
(S4 = 1.34 with SE 0.24), which was the outlier contamination, not design.

The household-unit arms (S2c onward) show a larger pre-birth rise (−0.33 at
−5/−4 versus −0.14 in the full sample): women who are already heads or spouses
are moving into larger dwellings in the years before the first birth. This is
the household-formation selection discussed above; it is a population
definition, not an error, and it should be stated when either number is used.

### Not done

The annual −2-only design and the September 12 timing pair were not rerun with
codes cleaned; if the annual curve is needed for the write-up, that is one
Torch array. Rooms were not rebuilt from the raw yearly variables (would
recover 1969 and 1976 answers and the 1977–1978 references). Moving variables
remain unverified.
