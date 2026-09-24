# First-birth housing-space event study, version 2 (final empirical version)

Date: 2026-09-24. Author decisions: household unit, PSID weights, full
biological birth history, entry at first adult observation, −3/−2 baseline
window as headline with −2/−1 as sensitivity. Horizon for the model target is
a separate, later decision; the full window path is reported for that purpose.

## Specification

Sun–Abraham interacted-cohort event study (`eventstudyinteract`) in two-year
interview windows, person and survey-year fixed effects, age and education
covariates, standard errors clustered by person. Outcome: rooms in the
dwelling, taken from the official year-specific PSID variables (41 variables,
crosswalk in `../first_birth_correction_review/all_wave_variable_crosswalk.csv`)
merged to the interview year; codes 9 (through 1984), 99 (1985–1993) and 98/99
(1994 on) set to missing; zero retained. Event: first biological child from the
full RELCHI1–20 history. Unit: one woman, current reference person or spouse,
per single-family-unit household-year, reference person first. Weights: PSID
longitudinal weight IW. Controls: women whose full history reports zero
children. Entry: first birth no earlier than the woman's first observation as a
current adult in the PSID. Support: every treated cohort must be observed in
the baseline window. Windows: ≤−8, −7/−6, −5/−4, [−3/−2], −1/0, +1/+2, +3/+4,
+5/+6, +7/+8, +9/+10, ≥+11.

Preparation `../../prepare_first_birth_rooms_v2.do` (rebuilt rooms verified
against the shifted shelf column: 811,486 overlapping values, 0 mismatches;
52,317 values only the official variables supply; 12,426 non-room codes among
current rows set to missing; 1,029 zeros retained). Estimator
`../../sa_rooms_first_birth_v2.do`; launchers `code/cluster/run_rooms_v2*.sh`;
collector `../../collect_rooms_v2.py`. Torch jobs 18424705/18424708 and
18425245/18425246, all pass; covariance symmetric and positive definite;
receipts reproduce the exported coefficients.

## Headline

| Design | Window | Rooms (SE) | Rows | Persons | Treated / control women |
|---|---|---|---|---|---|
| **Household, baseline −3/−2** | **+3/+4** | **1.03 (0.07)** | 63,338 | 5,309 | 3,655 / 1,654 |
| Household, baseline −2/−1 (sensitivity) | +2/+3 | 0.76 (0.06) | 65,053 | 5,381 | 3,727 / 1,654 |
| All current adults, unweighted (robustness) | +3/+4 | 0.82 (0.03) | 300,684 | 28,663 | 12,564 / last cohort |

Full paths with standard errors are in `window_path.csv`; figure
`first_birth_rooms_v2.png/pdf`.

Household path relative to −3/−2: −0.50 (≤−8), −0.50 (−7/−6), −0.34 (−5/−4),
0, +0.41 (−1/0), +0.82 (+1/+2), +1.03 (+3/+4), +1.17 (+5/+6), +1.20, +1.25,
+1.16 (≥+11). Baseline mean 4.70 rooms (weighted). Excluded by data coverage:
first births 1968–1970 and 1985–1986 (no interview in the baseline window in
this sample).

Reading: relative to two to three years before the birth, a household adds
0.4 rooms by the birth year, 0.8 rooms in the two years after, 1.0 rooms at
three to four years, and about 1.2 rooms from five years on. The rise in the
years before the baseline (from −0.5 at seven years out) is household
formation: these are women who already head a household, and they are moving
into larger dwellings as the birth approaches. Under the −2/−1 baseline the
anticipation year is inside the reference and the +2/+3 window reads 0.76.

## Reconciliation with every earlier number

| Object | Rooms | Why it differs |
|---|---|---|
| May slide, +3 coefficient | 0.77 | rooms dated one interview late; single-year +3 with an unpinned reference for post-1998 cohorts; codes as counts |
| August household target, +3 minus −1 | 0.72 | annual clock, −2 omitted, unpinned reference for 16 cohorts; −1 start removes anticipation |
| Sept 23 window, all adults, codes as counts | 0.60 | 4,559 code rows as ~90-room outliers |
| Sept 24 window, all adults, codes cleaned (S0c) | 0.73 | includes 11,003 non-current person-years and the original block's entry rule |
| **Final, all adults (A)** | **0.82 (0.03)** | current adults only, rebuilt rooms, entry at first adult observation |
| Sept 24 sequence endpoint (S6, August rules) | 0.93 | August entry rule (−0.05), shifted rooms, last-cohort control |
| Final household with August entry rule (Hentry) | 0.97 (0.07) | 48,225 rows |
| Final household with shifted rooms (Hshift) | 1.00 (0.07) | rebuilt rooms add 0.02 |
| Final household with last-cohort control (Hctrl) | 1.02 (0.07) | designation of the control cohort is immaterial |
| **Final household (H)** | **1.03 (0.07)** | |

The gap between the household and all-adult designs (0.20) is the unit: women
who already head a household versus every adult, including adult children in
their parents' dwelling, plus survey weights (Sept 24 sequence: unit +0.16,
weights +0.06, everything else ≤0.03). The residual between S0c (0.73) and
Ashift (0.83), which share the shifted rooms, is the non-current rows and the
entry-rule definition; it was not split further because it concerns the
robustness design only.

## Open items for the author

1. Horizon of the model target: read from the household path once the model's
   observation rule is fixed. Not chosen here.
2. Interpretation of the pre-baseline rise in the household design (household
   formation before the birth) in the write-up.
3. Frozen calibration target 0.720246 and its weight are unchanged; a new
   target contract name is in `target_receipt.csv` for when the horizon is set.

## Reproduce

```bash
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -bq do /tmp/psid_rooms_v2_20260924/prepare_first_birth_rooms_v2.do   # writes analysis_sample.dta (private)
sbatch code/cluster/run_rooms_v2.sh toy <task_root>; sbatch --array=0-2 --dependency=afterok:<smoke> code/cluster/run_rooms_v2.sh full <task_root>
code/model/.venv/bin/python3 code/data/psid_followup_mar2026/collect_rooms_v2.py
```
