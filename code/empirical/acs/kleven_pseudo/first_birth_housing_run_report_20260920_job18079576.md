# First-birth housing run readout

This report is generated from the saved receipts and per-fit checkpoints for Torch job 18079576. It is a diagnostic matched pseudo-panel readout; the estimates are conditional on the constructed source-key matches and are not causal claims.

## Run and failure evidence

- Job 18079496: failed after 33 seconds at the driver source-key assertion: `true ACS source key is not SAMPLE:YEAR:SERIAL:PERNUM`. The candidate-support step had already written its compact outputs. The failure receipt is `failure_receipt_18079496.json`. The source packet contained observed whitespace-delimited tokens (`SAMPLE YEAR SERIAL PERNUM`), while relabeled CPS keys retained the `CPS:` form.
- Job 18079576: completed in 7:52 after the parser accepted the observed whitespace source-key form and retained the source-key guard. No new allocation is active.
- The corrected run passed dependency smoke, real matched-key smoke, source overlap validation, and the estimator interface before fitting.

## Saved-source and join evidence

- Verified source packet: 2,784,052 rows; SHA-256 `edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e`; 9,919,999,546 bytes. Verified unique key intersection is 2,190,987 and shared source years are 2005--2019.
- V5 panel rows: 3,779,568; true ACS rows: 2,004,566; original CPS rows: 1,736,538; relabeled CPS rows: 38,464; total CPS rows in the join audit: 1,775,002.
- True ACS rows matched to the verified source packet: 1,774,257; unmatched true ACS rows: 230,309.
- Repeated source-household clusters: 1,573,728; distinct source-household clusters: 200,529.
- Weights, event time, and labor outcomes are unchanged by the source join.

## Housing coding and support

- ROOMS: 1,751,814 valid; 68 unknown code 28; 22,375 missing code. The primary outcome is capped at 9.
- BEDROOMS: 1,751,882 valid; codes 1--22 are transformed by x-1 and capped at 5; 22,375 missing code.
- OWNERSHP: 1,751,882 valid; codes 1 and 2 are retained as 1/0; 22,375 missing code.
- Across outcome-specific source-matched support cells (three outcomes are reported separately): 5,322,771 rows, 5,255,578 outcome-valid rows, and 67,193 missing-outcome rows. State-by-gender detail is in `state_gender_support.csv`; short-window cell detail is in `short_window_support_cells.csv` and its compact gate is in `short_window_support_gate.csv` (36/36 groups pass positive outcome support before the future source-year filter).
- The compact support receipt does not retain sum(w^2), so weighted effective sample size is not fabricated. Its definition for a future retained-weight receipt is `(sum(w)^2)/sum(w^2)`.

## Saved-fit inference diagnostics

- Primary fit: level coefficients with source-household clustered covariance; heteroskedastic covariance is a sensitivity comparison. `contrast_se_comparison.csv` records both standard errors for post-minus-pre and +3-minus-(-1).
- `pretrend_wald_clustered.csv` reports joint clustered-V Wald tests for event times -5, -4, -3, and -1, excluding the -2 reference event, with estimable rank and df. These are pretrend pattern diagnostics, not a validity proof.
- The compact fit has one gender-specific regression across the six states; the curves retain the common fit-level nobs and source-household-cluster counts. The reference event is normalized to zero, so the saved curves do not identify raw housing level means at event -2.

## V6 comparison

- The earlier V6 ownership continuation (receipt: `/scratch/td2248/projects/kleven_acs_pilot_20260917/overnight_ne_benchmark/ne_housing_v6_18047612/ne_housing_receipt.json`) used the same V5 panel lineage, a true-ACS observed-ownership sample with `n_ownership=962538`, percentage-point ownership levels, and heteroskedastic standard errors. It excluded 38,464 relabeled CPS rows and 22,375 unknown true-ACS ownership codes.
- The new readout uses the verified extract27 source packet, retains the full source-key join audit, and uses source-household clustering as primary with heteroskedastic standard errors as sensitivity. The sample/source-year overlap and variance estimator differ from V6; this is a specification comparison, not a quantitative decomposition or evidence of a failed replication.

## Prepared next sensitivity

- `short_window_sensitivity_recipe.json` prepares, but does not run, the common implied-event-cohort window [-2,+3] over the verified 2005--2019 overlap, with -2 as reference, -1 as the pre-event, +3 minus -1 using the full covariance matrix, and no rematching.
- The cohort is a support label `doiy - numeric(t_es_lw)`, not a biological birth year. The compact gate requires positive outcome-valid source-matched support in every event cell; the future invocation must reapply that gate after restricting true ACS source YEAR to 2005--2019 and reports missing outcomes separately.

## Second-birth candidate availability diagnostic

- The same allocation's bounded candidate-support step found 21919 target rows, 7598 target rows with any donor, and 7598 with any eligible donor under the declared exact demographic cells.
- It found 8855 fixed full-pre-window target rows and 1414 with any donor; the event-time receipt is `second_birth_candidate_support_20260920/candidate_support_by_event_time.csv`. This is candidate availability only, with no matching assignment, coarsening, or causal interpretation.


## Files

- `housing_event_curves_labeled.png` (saved-CSV figure with visible labels rooms, bedrooms, and ownership (pp), separated by state).
- `run_first_birth_short_window.R` is the reviewed no-submit driver for fixed implied cohorts 2007--2016 and true ACS source years 2005--2019; `test_first_birth_short_window.R` checks the six-cell support gate and Kish ESS on a tiny fixture.
- `state_gender_support.csv`, `short_window_support_cells.csv`, `short_window_support_gate.csv`, `pretrend_wald_clustered.csv`, and `contrast_se_comparison.csv` are compact machine-readable receipts.
- This report and its generator are durable source artifacts; the large raw panel, full `fits.rds`, and cluster submission are not copied or recomputed locally.
