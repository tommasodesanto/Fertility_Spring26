# PSID Follow-Up (March 2026)

This folder contains standalone scripts for:

1. Replicating the single Sun-Abraham `own_f_c_y_all` plot used in slides.
2. Running first-pass wealth regressions centered on ownership transitions around first birth.
3. Running a twins-IV and same-sex-IV implementation on ownership outcomes.
4. Running a second-born-gender design explicitly tied to mixed-sex preference.

No existing project files were overwritten.

## Files

- `sa_replication_own_only.do`
- `wealth_regressions_v1.do`
- `twins_and_gender_iv_v1.do`
- `secondborn_gender_design_v1.do`
- `compile_no_location_family_space_packet.py`
- `future_parent_dp_moments_v1.do`
- `build_future_parent_dp_moments_v1.py`
- `audit_intergen_bequest_family_size_targets.R`
- `build_intergen_bequest_retention_targets.R`
- `output/` (generated artifacts and logs)

## Exact Replication Target

The target plot is the PSID homeownership event-study around first birth:

- `own_f_c_y_all.png` used in:
  - `latex/Slides_20Jan.tex`
  - `latex/Slides_21Jan_model_first.tex`

The replication script keeps the same core logic as:

- `/Users/tommasodesanto/Desktop/Projects/Fertility/Codes/code_per tommi_addingcontrolsandfixingthings_modifiedbyludo.do`

## Wealth Regression Design (v1)

Design choice: use first-birth event time and focus on **pre-birth renters** (`t=-2` ownership = renter), years `1984-2019` where wealth variables are populated.

### Outcomes

- `own_post3`: became owner in years `[0,3]` after first birth
- `moved_to_own_post3`: transition rent->own that occurs with a move in `[0,3]`
- `rooms_change_post3`: mean rooms in `[0,3]` minus rooms at `t=-2`

### Key regressors

- `asinh(networth2r_pre)` (net worth excluding home, pre-birth)
- wealth components (`asinh(wlthsavetotr_pre)`, `asinh(wlthfundtotr_pre)`, `asinh(wlthodebtotr_pre)`)
- controls: age, education FE, first-child-year FE, female indicator

### Why this design

It isolates the tenure threshold margin that matters for the model: whether pre-birth wealth allows renters to cross into ownership when children arrive.

## Twins/Gender IV Design (implemented)

`twins_and_gender_iv_v1.do` implements two IV designs on pre-birth renters:

1. `twin_firstbirth` instrument for `two_plus` children (twins surprise at first birth).
2. `same_sex_first2` instrument for `three_plus` children (mixed-sex preference logic).

Outcomes:

- `own_post3`: owner at least once in years `0..3` after first birth.
- `moved_to_own_post3`: rent->own transition with a move in years `0..3`.

Core outputs:

- `output/twins_gender_iv_v1/iv_first_stage_v1.tex`
- `output/twins_gender_iv_v1/iv_vs_ols_own_post3_v1.tex`
- `output/twins_gender_iv_v1/iv_vs_ols_moved_to_own_post3_v1.tex`
- `output/twins_gender_iv_v1/iv_design_summary_v1.csv`

## Second-Born Gender Design (implemented)

`secondborn_gender_design_v1.do` explicitly tests the “want two kids of different sex” margin using:

- `mixed_first2 = 1(child1_sex != child2_sex)` as the design variable/instrument for `three_plus`.
- split first-stage checks for `second_girl` conditional on first child sex.

Core outputs:

- `output/secondborn_gender_v1/secondborn_first_stage_v1.tex`
- `output/secondborn_gender_v1/secondborn_own_post3_ols_iv_v1.tex`
- `output/secondborn_gender_v1/secondborn_moved_to_own_post3_ols_iv_v1.tex`
- `output/secondborn_gender_v1/secondborn_design_summary_v1.csv`
- `output/secondborn_gender_v1/secondborn_design_teststats_v1.csv`

## Execution Notes

Run scripts from StataMP in batch mode, for example:

```bash
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -b do /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/sa_replication_own_only.do
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -b do /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/wealth_regressions_v1.do
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -b do /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/twins_and_gender_iv_v1.do
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -b do /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/secondborn_gender_design_v1.do
```

## Run Status (March 7, 2026)

- `wealth_regressions_v1.do`: completed successfully.
  - Log: `output/wealth_v1/wealth_regressions_v1.log`
  - Regression table: `output/wealth_v1/wealth_transition_regressions_v1.tex`
  - Event sample: `output/wealth_v1/wealth_event_sample_pre_renters_v1.dta`
  - Tercile summary: `output/wealth_v1/wealth_tercile_outcomes_v1.csv`

- `sa_replication_own_only.do`: completed successfully.
  - Log: `output/sa_replication/sa_replication_own_only.log`
  - Graph: `output/sa_replication/own_f_c_y_all_repl.png`
  - Estimates: `output/sa_replication/own_f_c_y_all_repl_estimates.dta`
  - Repl-vs-original merge: `output/sa_replication/own_f_c_y_all_repl_vs_original.dta`

- `twins_and_gender_iv_v1.do`: completed successfully.
  - Log: `output/twins_gender_iv_v1/twins_and_gender_iv_v1.log`
  - Main sample: `output/twins_gender_iv_v1/analysis_sample_iv_v1.dta`
  - Parent instruments: `output/twins_gender_iv_v1/parent_instruments_v1.dta`

- `secondborn_gender_design_v1.do`: completed successfully.
  - Log: `output/secondborn_gender_v1/secondborn_gender_v1.log`
  - First-stage diagnostics: `output/secondborn_gender_v1/secondborn_first_stage_v1.tex`
  - Test-stat export: `output/secondborn_gender_v1/secondborn_design_teststats_v1.csv`

## Future-Parent Down-Payment Moments (June 30, 2026)

`future_parent_dp_moments_v1.do` builds an ID-level PSID sample of respondents
observed as renters and still childless at ages 25-30, then splits pre-birth
liquid wealth and down-payment distance by eventual parenthood horizons. The
Python companion `build_future_parent_dp_moments_v1.py` writes:

- `output/future_parent_dp_moments_v1/future_parent_dp_entry_sample_v1.csv`
- `output/future_parent_dp_moments_v1/future_parent_dp_moments_v1.csv`
- `output/future_parent_dp_moments_v1/future_parent_dp_comparison_v1.csv`
- `output/future_parent_dp_moments_v1/future_parent_dp_shortfall_bins_v1.csv`

The first-pass national target-price diagnostic finds that eventual parents are
not more down-payment constrained than eventual childless renters in this
sample; they have higher income and are more often 20%-down-payment sufficient.

## Intergenerational Bequest Target Audits (July 15, 2026)

`audit_intergen_bequest_family_size_targets.R` repairs the old cross-sectional
late-life wealth construction. `build_intergen_bequest_retention_targets.R`
constructs survivor-conditional ten-year wealth changes by completed fertility
and lagged initial wealth, with a four-year sensitivity. Their self-contained
result packets are under:

- `output/intergen_bequest_family_size_audit/`
- `output/intergen_bequest_retention_targets/`

## First-Pass Wealth Findings (v1)

From `wealth_transition_regressions_v1.tex`:

- `asinh(networth2r_pre)` is positively associated with becoming an owner within 3 years after first birth (`p < 0.05`) among pre-birth renters.
- In the component model, `asinh(wlthsavetotr_pre)` shows a positive association (`p < 0.05`), while `asinh(wlthfundtotr_pre)` and `asinh(wlthodebtotr_pre)` are not significant in this first pass.
- For `moved_to_own_post3`, the net-worth coefficient is near zero in this specification.

From `wealth_tercile_outcomes_v1.csv`:

- Ownership by +3 years is higher in the top pre-birth wealth tercile than in the middle/lower terciles.
- Moving rates are lower in the top wealth tercile, consistent with some households crossing to ownership without relying as much on relocation.

These are descriptive/first-pass estimates and should be tightened with cleaner cohort restrictions and robustness checks before being treated as final evidence.

## No-Location Family-Space Packet (May 23, 2026)

`compile_no_location_family_space_packet.py` compiles existing PSID outputs and
a fresh split of the saved pre-birth-renter event sample without using
location-specific PSID. It writes:

- `output/no_location_family_space_packet_20260523/PSID_NO_LOCATION_FAMILY_SPACE_PACKET.md`
- `output/no_location_family_space_packet_20260523/psid_birth_rooms_event_summary.csv`
- `output/no_location_family_space_packet_20260523/psid_prebirth_renter_outcomes_by_income_tercile.csv`
- `output/no_location_family_space_packet_20260523/psid_prebirth_renter_outcomes_by_rooms_slack.csv`

The packet documents birth-related room responses, moved-for-size windows,
pre-birth renter heterogeneity by income and room slack, and the existing
non-monotone wealth-parenthood fact. It does not identify center-to-periphery
moves or metro family-size-premium effects.

## First-Pass Twins/Gender Findings

From `output/twins_gender_iv_v1/iv_first_stage_v1.tex`:

- `twin_firstbirth -> two_plus`: coefficient `0.316` (SE `0.050`), strong first stage.
- `same_sex_first2 -> three_plus`: coefficient `0.081` (SE `0.044`), weaker first stage.

From `output/secondborn_gender_v1/secondborn_first_stage_v1.tex`:

- `mixed_first2 -> three_plus`: coefficient `-0.085` (SE `0.044`), F-stat `3.7238`.
- Split by first-child sex, `second_girl` is not strongly predictive in either subsample.

Interpretation at this stage: twins design is much stronger than the same-sex/mixed-sex design in this sample; IV point estimates on ownership outcomes are imprecise and should be treated as exploratory.
