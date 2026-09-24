# Opus 5.5: write the empirical evidence section of the JMP mock

Snapshot: September 24, 2026, end of day. Author: Tommaso De Santo. This
merges the September 24 ChatGPT handoff with the Fable session that produced
the final PSID estimates. Later explicit author decisions take precedence.

Project root: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.
Paths below are relative to it. Read the saved sources; this document is a map
with the verified numbers, not a substitute for the receipts.

## 1. Deliverable and ownership

Write into the EXISTING `latex/JMP_DS_mock/sections/02_empirical_evidence.tex`
(43 lines, uncommitted edits present: read it first, keep what is useful, and
replace obsolete claims deliberately; its line 17–18 sentence "0.72 rooms,
woman-clustered" is superseded). Update the EXISTING, currently unwired
`latex/JMP_DS_mock/data_appendix.tex` (label `app:empirical-data`, already
referenced from the section) and add its `\input` to `JMP_DS_mock.tex` beside
the other inputs. Do not create another appendix, suggestions folder, or
manuscript. Use existing figures and tables. Use the existing bibliography.

Do not edit `latex/JMP_DS_draft/` (author-owned), the slides, the model, data,
regression scripts, or `sections/04_quantification.tex`. No new estimation,
cluster jobs, calibrations, or policy exercises. Do not change any calibration
target by writing a revised estimate into prose; the frozen target 0.720246 and
its weight are unchanged and its replacement is an author decision. Preserve
unrelated dirty files.

Complete one grounded writing pass, one claim-to-source check, one typesetting
pass (compile the mock twice with build products outside its source directory,
deliver `output/pdf/JMP_DS_mock.pdf`, inspect changed pages), and one compact
receipt. If a decision is missing, finish the independent text and list the
exact unresolved choice in the receipt. Then stop.

## 2. Read first

Repository startup (`memory/AGENT_MEMORY.md`, `memory/daily/2026-09-24.md`,
`CALIBRATION_STATUS.md` September 24 entries), then:

- `docs/style/econ_writing_style_guide.md` in full.
- `latex/README.md` for the three document roles.
- `latex/JMP_DS_draft/sections/03_model.tex` (read only, for voice and notation).
- `latex/JMP_DS_mock/sections/02_empirical_evidence.tex` and
  `sections/04_quantification.tex` (the latter read only; its first-birth row
  still says 0.600 and is stale; flag in the receipt, do not edit).
- `docs/model/ACTIVE_DECISION_LEDGER.md` for recorded author decisions.
- `docs/model/first_birth_rooms_timing_memo_20260923.md`: the full chronology
  of the PSID correction. Its final three sections are the current state.

## 3. The PSID evidence: what was established on September 23–24

Everything below is verified against saved receipts. Use it to write; cite
the packages, not this file.

### 3a. The data correction (belongs in the data appendix, two or three sentences, no workflow narrative)

The merged PSID panel used before September 2026 stored three dwelling items
(rooms, moved since last interview, reason for move) on the row one interview
EARLIER than the interview at which they were reported (one year while the PSID
was annual through 1997, two years after). Ownership, age and the family
identifier were correctly dated. Verified on the full panel: 811,486 of
811,486 room values match the next wave's official variable, none match their
own wave among people whose answer changed. The final estimates rebuild all
four items from the official year-specific PSID family variables merged to the
interview year (`code/data/psid_followup_mar2026/psid_family_item_crosswalk.csv`
and `output/first_birth_correction_review/all_wave_variable_crosswalk.csv`).
Non-room codes (9 through 1984; 99 in 1985–1993; 98 and 99 from 1994) are set
to missing; a reported zero is retained as a valid shared-room answer. In the
reason-for-move item, code 9 is DK/NA before 2019 and "homeless" in 2019 and is
treated as missing throughout. The "moved since" question refers to the spring
of the prior survey year through 2001 and to January 1 of the prior year from
2003, so it is an interval outcome. Consequences for the old figures: the May
"+3" coefficient of 0.77 measured rooms four to five years after the birth
against the year before or the year of the birth. State this once, plainly.

### 3b. Estimator (identical in every selected fit; verify in each `fit_receipt.csv` and `estimation.log`)

Sun–Abraham interacted-cohort event study (`eventstudyinteract`), cohort = year
of first biological birth from the full 20-slot PSID birth history; event time
grouped into two-year interview windows so that every cohort is observed in the
reference window despite biennial interviews after 1997 (from the 1999 birth
cohort on, no cohort is observed at both single years −2 and +3, which made the
old annual single-year coefficients depend on an arbitrary normalization for
half the cohorts). Omitted window −3/−2 in the headline; the last pre-birth
year is deliberately outside the reference so anticipatory moves count as
response. Windows: ≤−8, −7/−6, −5/−4, [−3/−2], −1/0, +1/+2, +3/+4, +5/+6,
+7/+8, +9/+10, ≥+11. Person and survey-year fixed effects, age and years-of-
education dummies, standard errors clustered by person, PSID longitudinal
individual weight IW (which restricts to sample persons). Control cohort:
adults whose full relationship history reports zero children ("confirmed
childless"); unknown histories dropped. Treated persons enter only if the first
birth is no earlier than their first observation as a current adult in the
PSID. Every treated cohort must have observations in the reference window;
cohorts without are excluded and listed in the receipts. Standard errors come
from the estimator's interaction-weighted variance including cohort-share
estimation (`e(V_iw)`), from which contrasts are computed. Each outcome uses its
own complete-case sample. Missing observations outside a window are not zeros.

Do not describe this as generic two-way fixed effects. Do not use "parity".

### 3c. Two well-defined household populations (never call them the same sample)

- **Household design (H).** One woman who is the household's reference person
  or spouse per single-family-unit household-year, reference person first.
  Author's stated unit choice on September 24 (avoids counting one dwelling
  through several adults). Its rows exist only while the woman holds that
  status, so its far pre-birth points come from women who headed a household
  six or more years before the birth.
- **Status-at-baseline design (A2h).** All current adults, men and women, who
  were reference person or spouse in the −3/−2 window, a characteristic fixed
  before the birth; all rows of those persons. Same weights, history, controls,
  entry rule, codes. Includes both partners of a sample couple when both are
  sample persons, so its person-clustered standard errors are slightly
  optimistic; say so once if it is used.
- The complement (A2n), adults NOT heading a household at baseline, and the
  pooled all-adult fit (A2) exist and explain the household-formation margin.

### 3d. Rooms: `code/data/psid_followup_mar2026/output/sa_rooms_first_birth_v2/` (README, `window_path.csv`, `summary.csv`, `target_receipt.csv`, figure `first_birth_rooms_v2.png/pdf`)

Rooms relative to the −3/−2 window (SE), full paths:

| Window | H (women heads/spouses) | A2h (heads/spouses at baseline) |
|---|---|---|
| ≤−8 | −0.50 (0.12) | +0.15 (0.08) |
| −7/−6 | −0.50 (0.10) | +0.04 (0.06) |
| −5/−4 | −0.34 (0.06) | −0.03 (0.04) |
| −3/−2 | 0 | 0 |
| −1/0 | +0.41 (0.05) | +0.61 (0.04) |
| +1/+2 | +0.82 (0.06) | +1.19 (0.04) |
| **+3/+4** | **+1.03 (0.07)** | **+1.47 (0.05)** |
| +5/+6 | +1.17 (0.08) | +1.62 (0.06) |
| +7/+8 | +1.20 (0.09) | +1.64 (0.07) |
| +9/+10 | +1.25 (0.10) | +1.83 (0.08) |
| ≥+11 | +1.16 (0.13) | +1.85 (0.09) |

Samples: H 63,338 household-years, 5,309 women, 3,655 treated, 1,654 controls,
weighted baseline mean 4.70 rooms, cohorts 1968–1970 and 1985–1986 excluded
for lack of a reference window. A2h 117,853 person-years, 9,310 persons, 3,302
treated, 6,008 controls, no cohort excluded.

Other saved rooms fits (same package): H with −2/−1 baseline, +2/+3 window
0.76 (0.06), N 65,053; H with −5/−4 baseline, +3/+4 = 1.30 (0.10); with −7/−6
baseline 1.43 (0.13) and a flat path further back (+0.22, +0.14, +0.04); the
increment from −3/−2 to +3/+4 is 1.03 / 1.00 / 0.94 across the three baselines.
Pooled all adults A2 0.97 (0.05), N 187,873; A2n −0.38 (0.08), N 87,239, path
−0.83 at the birth recovering to zero by +7/+8 (they leave a parent's larger
dwelling). Reconciliation fits: rebuilt versus shifted rooms ±0.02; control-
cohort designation 0.00; the August "first birth after first head/spouse
observation" rule −0.05 with 24% fewer rows.

Interpretation the evidence supports: relative to six or more years before the
first birth the path is flat; adjustment starts about five years before the
birth and is monotone through it; a household that already exists two to three
years before the birth adds about 0.6 rooms by the birth, 1.2 in the next two
years, 1.5 at three to four years and 1.6 to 1.8 thereafter; the household
design's pre-baseline rise from −0.5 is a consequence of requiring head/spouse
status row by row (a row-selection effect), not a background trend, since
fixing status at baseline gives a flat pre-period. Adults not yet heading a
household lose rooms at the birth. Present H and A2h as two comparisons in one
table with sample, window and outcome columns; do not merge them into one
"effect" column.

### 3e. Ownership and moving: `code/data/psid_followup_mar2026/output/sa_first_birth_outcomes_v3/` (README, `window_path.csv`, `summary.csv`, figure `first_birth_outcomes_v3.png/pdf`)

+3/+4 versus −3/−2, H / A2h, with baseline levels and SEs:

| Outcome | Baseline (H / A2h) | −1/0 | +1/+2 | +3/+4 | +5/+6 |
|---|---|---|---|---|---|
| Owns the dwelling (share) | 0.41 / 0.43 | +0.11 (0.01) / +0.16 (0.01) | +0.16 (0.02) / +0.25 (0.01) | +0.16 (0.02) / +0.27 (0.01) | +0.14 (0.02) / +0.26 (0.01) |
| Moved since last interview (share) | 0.58 / 0.54 | −0.04 (0.02) / −0.08 (0.01) | −0.12 (0.02) / −0.17 (0.01) | −0.16 (0.02) / −0.21 (0.01) | −0.16 (0.02) / −0.21 (0.01) |
| Moved for more space, share of all household-years | 0.070 / 0.064 | +0.014 (0.009) / +0.015 (0.006) | +0.034 (0.009) / +0.021 (0.006) | +0.020 (0.010) / +0.010 (0.006) | +0.007 (0.010) / +0.009 (0.006) |
| Moved for neighbourhood, share of all | 0.032 / 0.030 | −0.004 (0.006) / −0.009 (0.004) | −0.005 (0.006) / −0.008 (0.004) | −0.010 (0.005) / −0.009 (0.004) | +0.002 (0.006) / −0.003 (0.004) |
| Moved for more space, share of MOVES | 0.122 / 0.123 | +0.054 (0.017) / +0.050 (0.013) | +0.140 (0.021) / +0.114 (0.017) | +0.133 (0.023) / +0.142 (0.020) | +0.108 (0.026) / +0.147 (0.022) |
| Moved for neighbourhood, share of MOVES | 0.054 / 0.058 | −0.017 (0.012) / −0.022 (0.009) | −0.011 (0.014) / +0.005 (0.011) | −0.017 (0.014) / −0.002 (0.013) | +0.022 (0.018) / +0.012 (0.015) |

Sample sizes in `summary.csv` (own H 61,229; moved H 63,817; reason-of-all H
62,813; movers-only H 19,456; A2h 111,272 / 118,875 / 116,697 / 27,462). A2h
pre-periods for ownership are flat (0.00, 0.00, −0.02). Moving rises INTO the
baseline window (A2h −0.13, −0.10, −0.03 before it) and falls after the birth.

The economic distinction to state plainly: the unconditional space-move share
rises little because total moving falls by 16 to 21 points after the birth,
while among households that do move the share moving for more space doubles
from 12% to about 25% from the birth window on and stays elevated for a decade.
Conditioning on movers selects on an outcome the birth changes; it describes
the composition of moves, not a causal effect on space-motivated moving. The
neighbourhood outcome does not respond in either version and serves as a
placebo. Separate outcome paths do not prove a single observed joint event.
Ownership: +16 to +27 points on a 42% base, concentrated in the birth window
and the two years after.

### 3f. Reconciliation with earlier numbers (one footnote or one appendix sentence at most)

May slide +3 coefficient 0.77: rooms dated one interview late, single-year
coefficient with an unpinned reference for post-1998 cohorts, non-room codes
as counts. August household target 0.72 (+3 minus −1, annual clock): unpinned
reference for 16 cohorts. September 23 window estimate 0.60: 4,559 code rows
as ~90-room outliers; 0.73 with codes cleaned, all adults, old choices. The
final numbers above supersede all of these. Do not quote 0.60, 0.72, 0.73, 0.77
or 0.93 as current estimates.

### 3g. Open author decisions (receipt, not prose)

1. Headline unit: H (author's stated unit choice) versus A2h (status fixed
   before treatment, flat pre-period, cleaner analogue of a model household).
   Present both; in prose lead with H as the author's chosen unit and describe
   A2h as the comparison that fixes household position before the birth.
2. Baseline: −3/−2 (headline) versus an earlier window that counts the
   anticipatory adjustment (1.30 / 1.43 for H). Report the −3/−2 numbers;
   mention that the total adjustment from the flat period is larger.
3. Horizon for the model target: not an empirical-section matter.
4. Second-birth results: none re-estimated on this specification; do not
   import old second-birth numbers.

## 4. Other evidence blocks (from the ChatGPT handoff; not re-verified on September 24, so check each source before quoting)

### Housing stock and tenure (descriptive motivation)

`docs/model/evidence_tenure_segmentation_20260918.md`,
`code/data/ahs_supply_snapshot/README.md`, and
`code/data/ahs_supply_snapshot/output_ahs_family_unit_menu_national/AHS_2023_FAMILY_UNIT_MENU.md`
with its tables and figure packet. Verify denominators (share of owners in
large homes is not the share of large homes rented). Occupied stock, not a
supply elasticity. Keep rooms and bedrooms distinct. Do not mix metro,
national or period definitions. The mock's current sentences (three-bedroom
shares, 9.2% of four-plus-bedroom homes rented) must be re-checked against the
packet before being kept.

### ACS matched pseudo-panel (complementary, descriptive)

`code/empirical/acs/kleven_pseudo/README.md`,
`first_birth_housing_run_report_20260920_job18079576.md`,
`first_birth_sensitivity_bundle_receipt_18080591.md`,
`second_birth_housing_run_report_20260920_job18080900.md`,
`second_birth_proxy_design.md`, and the linked receipts; national continuation
under `output/national_acs_comparison/national_continuation_20260921b/`
(`continuation_receipt.json`, selected `*_receipt.json`, contrast and covariance
files). Distinguish completed fits from unavailable ones yourself; regional
curves are not national estimates; constructed pre-birth donors are not
repeated observations of the same women; roster-based birth order is not a
biological history. Use only what adds to the PSID story.

### Fertility instruments (separate identification exercise)

`output/acs_fertility_iv/data_appendix/SOURCE_AUDIT.md`,
`output/acs_fertility_iv/national_128g_results/review_summary.json` and
`national_18row_table.csv`, `output/acs_fertility_iv/samesex_diagnosis/report.md`,
`latex/JMP_DS_draft/sections/appendix_acs_fertility_iv.tex` (read only, adapt
substance if useful), and for PSID
`code/data/psid_followup_mar2026/output/iv_housing_reaudit_20260809/README.md`
(supersedes March interpretations) and
`code/empirical/acs/kleven_pseudo/acs_twins_samesex_housing_contract.md`.
Twins versus first-two-child sex composition have different realization dates,
samples, first stages and weak-instrument uncertainty; do not restrict the
sex-composition sample to exactly two children; the same-sex housing pattern
raises exclusion concerns without formally rejecting exclusion. Roster proxies
do not establish biological motherhood or exact twin births. Detailed
diagnostics go to an appendix; do not frame the section around negative tests.

## 5. Paper question and narrative

How housing costs and access to family-sized housing interact with fertility,
tenure and housing allocation over the lifecycle. Organize the evidence around
three facts, letting the verified numbers set the claims: (1) how housing size
is distributed across tenure; (2) how rooms, ownership and mobility change
around the first birth, including adjustment before the birth and the
household-formation margin; (3) what the cross-sectional and instrument
exercises add and what variation they identify. Facts motivate model
ingredients; they do not prove that housing constraints cause low fertility.
Model counterfactuals stay out of this section. Do not resurrect the old
spatial or center–periphery evidence.

## 6. Writing and presentation

First person, the author's voice (read the draft's model section), the economic
fact first, then the comparison, then the evidence; restrained exposition in
the manner of Menzio/Fernández and Boar–Gorea–Midrigan (consult passages, do not
copy). Compact subsections with the mock's run-in subheaders. No job IDs, arm
names (H, A2h, S0c and the like), workflow narrative, calibration receipts or
repeated caveat lists in the paper; name the designs in words. State each
economically material limitation once where it matters. Three decimals at most;
meaningful units, confidence intervals or standard errors, and sample sizes. A
table identifies sample, window and outcome; it never hides incomparable
estimates in one column. Do not use "parity". Cite data and estimators from
verified references (Sun and Abraham 2021 for the estimator; PSID documentation
for items), not from memory. Use existing figures: `first_birth_rooms_v2.pdf`
and `first_birth_outcomes_v3.pdf` are the current ones; if a figure is used,
copy or reference it from its output folder without regenerating results.

## 7. Completion check

Record a claim-to-source map outside the manuscript (every number to its CSV
or receipt). Compile twice, inspect changed pages, deliver
`output/pdf/JMP_DS_mock.pdf`. Report files changed, narrative chosen,
unresolved decisions (section 3g), cross-section discrepancies (quantification
row 0.600; frozen target 0.720246), and any claim you could not verify. Stop
and return control to Tommaso.
