# PSID wealth lifecycle data audit — September 28

Status: source and saved-table audit complete. No confirmed numerical or unit error. Microdata tail/uncertainty audit remains pending; this is not certification of the entire dataset. No model code, calibration target, input data or graph changed. No microdata read or model solve was performed on the Mac or Torch.

## Confirmed

The dashboard builder `code/model/tools/build_e5f_lifecycle_data_dashboard.py:44–65` selects survey years 2005/2007, reference persons (`RELTOHEAD_=10`), ages18–85, positive finite IW and finite NETWORTHR. Working ages18–65 additionally require nonnegative finite EARNINDRRC. It does not replace missing observations with zero. The selected sample reproduces11,324 family-years and aggregate wealth/earnings6.926583791073. Independently recomputing the aggregate from the saved17 age-cell numerators and denominators yields6.926583791072991. Each plotted scaled mean agrees with its saved dollar mean/82,881.954 denominator to<5e-14. See `saved_table_check.json`.

This matches the accepted July accounting convention in `docs/model/intergen_wealth_target_beta_audit_20260723.md` and the sample filter in `code/data/psid_followup_mar2026/audit_aggregate_wealth_earnings_ratio.R:99–140`. July pooled2005–2019; the current initial-environment diagnostic deliberately pools only2005/2007. The current figure is not a replay of July's four broad age-bin ratios: those ratios divide each age-bin wealth by that bin's earnings; the graph divides every age-bin wealth mean by the common overall working-age earnings mean. Neither convention is inherently erroneous, but the numbers must not be interchanged.

Construction files under `/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/Construction_Files/Code/` establish:
- `01 Collect wealth variables.do:289–340`: NETWORTH is total family-unit net worth, from S717(2005) and S817(2007), not nonhousing wealth.
- `03 Generate income measures.do:109–132`: earnings combine RP and spouse; label explicitly says tax year.
- `07 Calculate inflation adjustments.do:111–151,165–204`: earnings use previous tax-year deflator and wealth survey-year deflator, both to a common real-dollar year. This is not a nominal/real mismatch.
- `01 Collect survey identifiers.do:1046–1074`: IW is the individual's longitudinal weight. Using the reference person's IW is the accepted project convention, not evidence that a family cross-sectional weight was used.

## Unresolved or interpretation limits

1. Extreme-value checks are missing. Original wealth construction preserves finite top/bottom codes, and inflation code preserves selected positive codes. The dashboard's finite-value test does not explicitly distinguish these. Such topcodes are censored amounts, not automatically missing-value contamination; whether any occur in this actual sample, and their influence, is unverified. Earnings similarly preserves9999999. Do not claim observed contamination from source possibilities alone.
2. Mean wealth can be sensitive to a few rich families. Oldest82–85 cell has245 pooled family-years; ages66–69 has357,70–73 has331,74–77 has332,78–81 has309. These are not independent-person counts. No weighted effective sample size, medians, top shares, leave-one-family influence, or uncertainty bands exist in the saved packet.
3. The plot is a cross-section of surviving reference-person households in2005/07, not a cohort trajectory. Household composition, reference-person selection, mortality selection and calendar conditions are not eliminated by dividing by common earnings. A drop across ages cannot be read as an individual retirement drawdown rate.
4. Survey-date wealth and prior-tax-year earnings differ in timing. Common denominator limits age-specific denominator artifacts, but this remains a measurement approximation.
5. The6.3GB source was not independently rehashed or read in this audit. Earlier provenance records file size/mtime and aggregate replay. The source was not found at the checked standard Torch paths; no bulk copy was attempted.

## Small next check

Locate an existing remote shelf or narrow extract, then perform one bounded selected-column data job on Torch, not the Mac. Retain ID, family interview identifier, year, age, RP relation, IW/family weights, wealth, gross earnings and relevant flags. Verify unique family-wave records, exclusions and actual topcodes. Produce each2005/2007 wave separately and pooled, median/mean/p90, effective N, top1% wealth share, and delete-one-family influence. Bootstrap by original family/reference-person cluster with the common earnings denominator recomputed per draw; uncertainty must reflect repeated observations and ratio estimation. Compare alternative weights only as labeled diagnostics, preserving the accepted target. This would tell us whether the jagged retirement data profile is genuine, composition-driven or noisy before attributing the model gap to bequests.
