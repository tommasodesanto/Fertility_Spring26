# SPEC — bootstrap standard errors + covariance for the 14 calibration moments (T2)

Author: lead (design + identification). Implementer: Codex. Verifier: lead.
Date: 2026-07-05 overnight loop. This is a NEW data-side harness that only
READS the existing builders/microdata; it does NOT modify any builder or any
model code. Output is a calibration-input artifact (the SE/covariance file the
D10 optimal-weighting program has been blocked on).

## Goal

Produce, for the 14-moment target vector of target set
`candidate_replacement_post_audit_v1`, a standard-error vector and a
14x14 (block-diagonal by data source) covariance matrix, by nonparametric
resampling of the SOURCE microdata — with a HARD reproduction gate on every
point estimate before any SE is reported.

Plus two supplementary tabulations (same sources, same session):
- (T2b) childlessness by income (and by education if available), for the D12
  gradient check (model says childlessness FALLS in income; US completed
  childlessness is flat-to-increasing).
- (T2c) young-childless-renter wealth-to-income BY ERA (PSID sub-periods vs
  SCF 2022), for the D7 era-sensitivity note (pooled PSID 0.179/0.10 vs
  SCF 2022 1.14/0.385).

## The 14 moments, published targets, and presumed source (VERIFY via the gate)

| key | target | source |
|---|---:|---|
| tfr | 1.918000 | fertility source — LOCATE + reproduce; if not reproducible, flag SE as "source-unconfirmed", do NOT fabricate |
| childless_rate | 0.188000 | same fertility source |
| own_rate | 0.575472 | ACS/MMS |
| own_family_gap | 0.167662 | ACS/MMS |
| housing_increment_0to1 | 0.664435 | PSID event study (Sun-Abraham) |
| old_parent_childless_nonhousing_wealth_to_income_gap_6575 | 1.007450 | PSID |
| prime30_55_childless_renter_mean_rooms | 3.805288 | ACS/MMS |
| prime30_55_childless_owner_share_rooms_ge6 | 0.596131 | ACS/MMS |
| old_nonhousing_wealth_to_income_median_6575 | 2.230461 | PSID |
| young_childless_renter_liquid_wealth_to_annual_gross_income_2535 | 0.179226 | PSID |
| prime30_55_childless_owner_minus_renter_mean_rooms | 2.418762 | ACS/MMS |
| old_age_own_rate | 0.764261 | PSID |
| own_rate_2534 | 0.341166 | ACS/MMS |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.367700 | ACS/MMS |

## Builders + microdata (from the grounding pass — read these for the exact filters)

- PSID moments: `code/data/psid_followup_mar2026/build_intergen_oldage_wealth_targets.R`
  and the young-wealth logic reproduced in `code/model/intergen_housing_fertility/calibration.py`
  (constants `PSID_ENTRY_WEALTH_RATIO_NODES_2535`, weighted mean 0.17922556,
  median 0.09996729). Microdata: `~/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta`.
  Key vars: NETWORTH2R, INCFAMR, RELCHIREP/RELCHINUM (children), AGEREP, year, IW (weight),
  and a family/interview id + person id for clustering.
- PSID event study: `code/empirical/roundup/empirical_roundup_first_birth_by_wealth_v1.do`
  (Stata) and `code/data/psid_followup_mar2026/sa_rooms_first_birth_variants_v1.do`.
  Reimplement the Sun-Abraham estimator in R with `fixest::feols` + `sunab()` if
  Stata is unavailable. Vars: ACTUALROOMS_, HOMEOWN, RELCHI1BYEAR, NETWORTH2R, IW.
  Target 0.664435 = the K=3 (0->1 child) rooms coefficient; the K=5 plateau
  coef_p5 = 0.843 is also produced (needed for T1b — report it here too).
- ACS/MMS moments: `code/data/mms_center_periphery/reconstruct_target_moments.R`,
  `build_intergen_one_market_housing_targets.R`, `audit_ownership_targets.R`.
  Microdata: `code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta`
  and processed `.../processed_data/yearly_rds_v7/fertility_microdata_*.rds`
  (2012-2023). Vars: age, owner/renter, rooms, nchild, eldch, perwt, ownershp.
- Childlessness-by-income: `code/data/mms_center_periphery/analyze_income_fertility_cross_section.R`
  (ACS) and `code/data/psid_followup_mar2026/build_income_wealth_fertility_master_v1.py`.
- Era: `code/data/mms_center_periphery/validate_acs_home_value_scf.R` (SCF 2022);
  PSID sub-period splits from the young-wealth logic above.

## HARD reproduction gate (do this FIRST, before any bootstrap)

For each moment, recompute the POINT estimate from the source using the builder's
exact sample filters and weights, and assert it matches the published target to
tolerance: |reproduced - target| <= max(1e-6, 1e-4*|target|). Log every
(key, target, reproduced, abs_diff, PASS/FAIL) to `output/reproduction_log.txt`.
Precedent: phase11 reproduced 0.1792255599/N=7163 to 1.7e-9 — aim for that on
the young-wealth moment. If a moment FAILS the gate (e.g., tfr/childless source
not found, or a filter mismatch), do NOT bootstrap it: mark its SE as NA with a
reason. Report all failures prominently — a wrong point estimate means the SE is
meaningless.

## Bootstrap design (LEAD-SPECIFIED — implement exactly)

Efficiency: load each source ONCE, apply filters, and CACHE the small filtered
analysis samples to `.rds` (per moment / per source). Bootstrap on the cached
samples, never re-read the 6-10 GB `.dta` per replicate. B = 1000 replicates
(fall back to B = 500 with a logged note if wall-clock exceeds ~2h per source).
Use a FIXED seed (e.g., 20260705) and save it.

Clustering + weights (do not deviate without logging why):
- PSID moments (young wealth, old-age wealth level/gap, old_age_own_rate):
  cluster bootstrap at the PSID FAMILY level (the panel/PSU unit — resample
  family ids with replacement, take all person-year rows of drawn families),
  then recompute the survey-weighted moment with IW. This respects within-family
  dependence across waves. Report the effective number of clusters.
- PSID event study (housing_increment_0to1): the PRIMARY SE is the analytic
  cluster-robust SE of the Sun-Abraham K=3 coefficient, clustered at the family
  id (fixest `cluster = ~famid`). Also report the K=5 plateau coefficient and its
  clustered SE. A family-cluster bootstrap of the SA regression is OPTIONAL
  (only if it finishes < 30 min); if run, report both.
- ACS moments: household-level weighted resample. ACS here is a repeated
  cross-section with person/household weights (perwt) and no replicate weights,
  so resample the household analysis-sample rows with replacement (B draws),
  each draw recomputing the perwt-weighted moment. Document this as
  "survey-weighted nonparametric bootstrap of the weighted estimator". If a
  cluster/strata var (e.g., PUMA or metro) is present in the builder, cluster at
  the metro level instead and log that choice.
- Fertility moments (tfr, childless): once the source is located, use its
  natural PSU (person for a cross-section; family for PSID). Log the choice.

Covariance:
- Reuse the SAME bootstrap draws across all moments from the SAME source so the
  within-source cross-moment covariance is captured (compute each replicate's
  full moment vector for that source, then Cov over replicates).
- Set cross-source covariance blocks to 0 (independent surveys) and SAY SO in
  the README. Assemble the 14x14 block-diagonal covariance.
- Sanity: diagonal of the covariance == se^2 for every moment; all diagonal
  entries > 0; correlations in [-1,1].

## Outputs (write to `code/data/moment_standard_errors/output/`)

1. `moment_se.csv` — columns: key, published_target, reproduced_point,
   repro_abs_diff, repro_pass, boot_mean, boot_se, n_analysis, n_clusters,
   source, psu, weight_var, B, method_note.
2. `moment_covariance.csv` — 14x14 labeled matrix (block-diagonal by source).
3. `moment_correlation.csv` — the implied correlation matrix.
4. `reproduction_log.txt` — the hard-gate log.
5. `childlessness_by_income.csv` — childless share by income quintile/decile
   (+ by education if the source has it), PSID and/or ACS, with the model's
   0.48/0.56/0.17/0.10/0.07 (z=0.6..1.4) pasted in a comment for comparison.
6. `young_wealth_by_era.csv` — young-childless-renter networth/income mean+median
   by PSID sub-period (e.g., 1984-1999 / 2000-2009 / 2010-2019) and SCF 2022
   (1.14/0.385), with N per cell.
7. `README.md` — the design (this spec's choices), what reproduced vs failed,
   the exact clustering/weight decisions, B, seed, and any fallbacks taken.

## Validation the implementer runs before declaring done

- `Rscript build_moment_bootstrap_se.R` runs end to end (or a `--source psid` /
  `--source acs` split if memory forces it), writing all 7 outputs.
- Every reproduction-gate line logged; at least the young-wealth 0.179 and the
  ACS room/own moments PASS.
- Covariance diagonal == se^2 (assert, tol 1e-12); SEs all finite and > 0 for
  passed moments.
- Print a short console summary: n moments passed gate, min/median/max SE,
  and the two supplementary tables' head.

## What to hand back to the lead

The console summary + the reproduction_log + the moment_se.csv contents, and any
moment that FAILED the gate with the reason. Do NOT report SEs for
gate-failed moments as if valid.
