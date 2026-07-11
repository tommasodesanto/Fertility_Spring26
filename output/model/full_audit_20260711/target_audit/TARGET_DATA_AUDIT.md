# Target-System Audit — DATA Side (15 moments, combined spec, July 10 overnight)

Auditor: data-side target specialist, 2026-07-11.
Scope: for each of the 15 targets in `candidate_replacement_post_audit_v1` +
`aggregate_mean_occupied_rooms_18_85` (fixed-spec extra), locate the in-repo
builder, verify the numerical value, sample, and definitions, and check
cross-target consistency (room units, head definition, weights, denominators,
negative wealth) plus the new `code/data/moment_standard_errors/` harness.

All file references are repo-relative from
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.

## 0. Headline

- 13 of 15 point values are byte-exact reproducible from in-repo builders +
  microdata. The July 5-6 bootstrap harness
  (`code/data/moment_standard_errors/`) reproduced 12 of the 14 post-audit
  targets to < 5e-7 (`output/reproduction_log.txt`), and I verified the 15th
  (new aggregate rooms target) matches its builder CSV byte-exactly.
- The 2 exceptions are `tfr = 1.918` and `childless_rate = 0.188`: no in-repo
  builder or raw artifact; provenance was closed 2026-07-08 at CITATION level
  only (CPS June 2024 Fertility Supplement, Census Historical Table H2, women
  40-44:
  `output/model/fable_size_mapping_audit_20260701/staging_drivers_20260706/TARGET_PROVENANCE_DRAFT.md`).
  The published cell (1,918 CEB per 1,000; 18.8% childless) matches the targets
  exactly and the two are wave-consistent (both June 2024), but the table file
  itself is not checked in.
- The NEW target `aggregate_mean_occupied_rooms_18_85 = 5.779970481941968` is
  built by an UNCOMMITTED working-tree edit to
  `code/data/mms_center_periphery/build_intergen_one_market_housing_targets.R`
  (diff vs HEAD adds `age <= 85` and the `post_moment(...)` block at lines
  179-183), output regenerated 2026-07-10 19:50, ~3 hours before the overnight
  launch. Its own output CSV carries
  `status = "candidate; re-audit source/sample/formula before SMM use"`.
- It is NOT a national ACS aggregate: the sample is the matched MMS metro
  sample (42 metros, middle collapsed to center), owner-or-renter household
  heads (PERNUM==1), ages 18-85, literal ROOMS, HHWT. This is internally
  consistent with the other room targets (same builder, same run), but any
  "US aggregate" description would be wrong.
- The bootstrap SE file (`output/moment_se.csv`) implies inverse-variance
  weights that differ from the overnight's ad hoc weights by up to ~500x in
  relative terms (details in section 4). The SE file predates the overnight
  and was not used.

## 1. Per-target ledger

Notation: builder = generating script; artifact = saved output that pins the
value; SE = family/metro-cluster bootstrap SE from
`code/data/moment_standard_errors/output/moment_se.csv` (B=1000, seed
20260705).

### 1.1 aggregate_mean_occupied_rooms_18_85 = 5.779970481941968 (NEW)

- Builder: `code/data/mms_center_periphery/build_intergen_one_market_housing_targets.R:179-183`
  ("`rows <- post_moment(rows, "aggregate_mean_occupied_rooms_18_85",
  weighted_mean_safe(df$rooms, df$hhwt), ...`"), UNCOMMITTED (working tree
  only; `git diff HEAD` shows the block and the added `age <= 85` filter at
  line 145).
- Artifact: `code/data/mms_center_periphery/output_intergen_one_market_targets/intergen_one_market_acs_housing_targets.csv`
  row 1, value 5.779970481941968, n=3,986,589, weight=419,883,313, generated
  2026-07-10 19:50:40 EDT. Byte-exact match to the target constant in
  `tmp/overnight_combined_20260710/run_overnight_combined.py:27`
  (`ROOMS_TARGET = 5.779970481941968`) and `run_wave2_lateral.py` (weight 6.0).
- Sample: ACS/IPUMS `extract27.dta` (9.9 GB, code/data/Spatial_aggregate_withmicrodata/raw_data/),
  2012-2023, `met2013 > 0`, `gq in (1,2)`, `pernum == 1`, `age in [18,85]`,
  `hhwt > 0`, `rooms > 0`, joined to MMS PUMA lookups and restricted to
  `mms_location in (center, periphery)` with middle→center
  (script lines 139-170). Tenure restricted to owner or renter.
- Formula: `sum(HHWT * ROOMS) / sum(HHWT)`; mean, literal ROOMS.
- Top-code: no explicit handling anywhere in the builder; I tabulated the
  cached ACS analysis sample
  (`code/data/moment_standard_errors/cache/acs_analysis_samples.rds`): ROOMS
  ranges 1-26 with weighted share at max 7.7e-6, so top-coding is immaterial.
- Model object (uncommitted diff in
  `code/model/intergen_housing_fertility/calibration.py:1090-1092`):
  `aggregate_housing_demand / max(total_mass, 1e-12)`. Production runs the
  markov-income path (`production_profile.py:134-135`
  `"use_income_types": True, "income_type_transition": "markov"`), where
  `aggregate_housing_demand` is a direct sum over the joint distribution, so
  the permanent-type aggregation caveat (solver.py:1415/1466-1476, weights vs
  raw total_mass) does NOT bind for the overnight. Model ages: J=17 four-year
  periods from 18 → 18-85, matching the new data filter.
- SE: NONE — the target postdates the SE harness (14 moments only).
- Consistency: same room units, head definition, weights, and geography as the
  prime30_55 room targets (same builder run; the regenerated CSV reproduces
  the older prime30_55 constants to all printed digits, confirming the sample
  is unchanged by the `age <= 85` edit).

### 1.2 tfr = 1.918 and childless_rate = 0.188

- Builder: NONE in repo. SE harness verdict (2026-07-06,
  `code/data/moment_standard_errors/output/reproduction_log.txt`): "status=FAIL
  | source=source-unconfirmed ... exact values appear only as calibration
  target constants/prose, not a source artifact."
- Provenance closure (2026-07-08):
  `output/model/fable_size_mapping_audit_20260701/staging_drivers_20260706/TARGET_PROVENANCE_DRAFT.md`
  — CPS June 2024 Fertility Supplement, Census Historical Table H2, women
  40-44, all marital classes: CEB per 1,000 = 1,918 (exact) and percent
  childless = 18.8 (exact). Wave-consistent (both 2024). Raw H2 table NOT
  checked in; citation-level provenance only.
- Model mapping: `calibration.py:1094` `"tfr": 2.0 * household_parity` — the
  household-doubling halving the provenance note demands is implemented.
  `childless_rate` = completed `parity_dist[0]` (post-fertility), the correct
  ever-childless counterpart.
- Caveats from the provenance note itself: CEB is per WOMAN unconditional on
  household formation, model unit is a household; 2024 H2 folds 7+ into the
  5-6 bin; wave-dependence is large for childlessness (0.150 in 2018 → 0.188
  in 2024).
- SE: NA (no microdata in repo).

### 1.3 own_rate = 0.57547241; own_rate_2534 = 0.34116609; old_age_own_rate = 0.76426097; own_family_gap = 0.16766167

- Builder: `code/data/mms_center_periphery/audit_ownership_targets.R`.
  Artifact: `output_ownership_audit/ownership_window_targets_all_sources.csv`,
  sample `household_heads_hhwt_due_housing`:
  - 30_55 overall = 0.575472413107657 (n=1,806,067) → own_rate ✓
  - 25_34 overall = 0.341166094183662 (n=538,998) → own_rate_2534 ✓
  - 65_75 overall = 0.764260969085936 (n=696,028) → old_age_own_rate ✓
  - 30_55 newparent_minus_nochildren = 0.167661670600099 → own_family_gap ✓
- old_age_own_rate is ACS/MMS, NOT PSID: the SE-harness log records "source
  correction: published target is reproduced by ACS/MMS DUE-housing ownership
  audit; PSID old-age ownership check is 0.863730454 and does not match".
- Sample (script lines 206-259): extract27, 2012-2023, `met2013>0`, gq 1/2,
  ages 20-84 base, `ownershp in (1,2)`, heads = `pernum==1 & relate==1`,
  HHWT, DUE-housing restriction `unitsstr %in% 3:10` + `rooms > 0`, MMS
  matched metros, middle→center.
- newparent = `nchild > 0 & eldch < 4` (eldest child under 4, line 69);
  no-child = `nchild == 0`. Model counterpart
  `own_gap_newparent_nonparent_3055` — data "no child in household" vs model
  nonparent state distinction is documented in `TARGET_MOMENT_OBJECTS`
  (calibration.py:225-230).
- Owner = `ownershp == 1` (owned or being bought, incl. mortgaged).
- SEs: own_rate 0.0204, own_rate_2534 0.0182, old_age_own_rate 0.0126,
  own_family_gap 0.0065 (42 metro clusters).

### 1.4 prime30_55 room targets (3.8052881; 0.59613112; 2.41876173)

- Builder: `build_intergen_one_market_housing_targets.R` (lines 190-213).
  Artifact CSV (regenerated 2026-07-10) matches:
  renter_mean_rooms = 3.805288099342435 (n=333,362);
  owner_share_rooms_ge_6 = 0.5961311181352554 (n=409,475);
  owner_minus_renter_mean_rooms = 2.418761731003141 (n=742,837). All ✓ to the
  8-digit target constants (calibration.py:97-100, 142).
- Sample: same as 1.1 but ages 30-55, childless = `nchild == 0`; owner/renter
  by OWNERSHP; literal ROOMS; HHWT; MMS metros middle→center.
- IMPORTANT sample inconsistency vs ownership targets: the room builder uses
  heads = `pernum == 1` ONLY (no `relate == 1`) and does NOT impose the
  DUE-housing `unitsstr %in% 3:10` structure restriction that all four
  ownership targets impose. So the room-target sample includes mobile homes /
  nonstandard structures and non-relate-1 first persons that the ownership
  sample excludes. Within the room targets themselves everything is
  consistent (same builder).
- SEs: 0.0738, 0.0175, 0.0714.

### 1.5 prime30_55_parent_3plus_minus_1to2_mean_rooms = 0.36769955881

- Original provenance gap (July 2 audit,
  `output/model/fable_size_mapping_audit_20260701/CODE_AUDIT_REPORT.md` A14)
  is now CLOSED by reproduction: the SE harness
  (`build_moment_bootstrap_se.R:484-510`) recomputes it from extract27 as
  `weighted_mean_safe(parent_3plus$rooms, hhwt) - weighted_mean_safe(parent_1to2$rooms, hhwt)`
  with `parent_3plus = prime_parent[nchild >= 3]`,
  `parent_1to2 = prime_parent[nchild in 1..2]`, prime_parent = ages 30-55
  heads with `nchild>0 & yngch<18`; reproduced 0.367699558813823
  (abs diff 4.4e-7, PASS). The dedicated target builder still does not emit
  this moment — the SE harness is currently its only in-repo generator.
- Model object: high-child current-parent minus one-child current-parent mean
  housing services, 30-55 (TARGET_MOMENT_OBJECTS, calibration.py:243-248);
  the coarse model high-child state → empirical 3+ bin mapping is a documented
  approximation, not a bug.
- SE: 0.0666.

### 1.6 housing_increment_0to1 = 0.66443467 (PSID)

- Artifact: `code/data/psid_followup_mar2026/output/sa_rooms_first_birth_one_variant_v1/rooms_f_c_y_all_summary.csv`:
  `coef_p3 = .66443467, se_p3 = .14851184` (also coef_p5=0.8431), sample
  248,566 obs / 23,041 ids; Sun-Abraham (eventstudyinteract) K=3 rooms
  coefficient, cluster(ID), absorb(year), covariates i.AGEREP i.EDUYEAR.
  Builder: `sa_rooms_first_birth_variants_v1.do` /
  `code/empirical/roundup/empirical_roundup_first_birth_by_wealth_v1.do`.
- Sample: PSID NATIONAL, all ages, ACTUALROOMS_ (PSID rooms; pre-event mean
  6.0095). Units are PSID rooms applied to a model calibrated to ACS ROOMS —
  roughly comparable but not identical instruments/samples.
- Mean-type: regression coefficient (controlled event-study), not a raw mean.
- SE: analytic clustered 0.1485.

### 1.7 old_nonhousing_wealth_to_income_median_6575 = 2.23046078 and old_parent_childless_nonhousing_wealth_to_income_gap_6575 = 1.00744952 (PSID)

- Builder: `code/data/psid_followup_mar2026/build_intergen_oldage_wealth_targets.R`.
  Artifact `output/intergen_oldage_wealth_targets/intergen_oldage_wealth_targets.csv`:
  `all_old_65_75_nonhousing_nw_to_income_median = 2.23046078407069` and
  `old_parent_minus_childless_nonhousing_nw_to_income_mean_gap = 1.0074495242682238`
  (both n=11,449). Both ✓.
- Sample: PSIDSHELF_MOBILITY (outside repo, `~/Desktop/Projects/Fertility/PSID/`),
  1984-2019, ages 65-75, IW weights, ALL INDIVIDUALS (no reference-person
  restriction — couples count twice with family-level wealth/income),
  children = completed RELCHIREP (RELCHIREP-based childless, answering the
  lead's question; ownership cross-check HOMEOWN exists but the target uses
  ACS for old_age_own_rate).
- Statistic types: the LEVEL is a weighted MEDIAN; the GAP is a difference of
  weighted MEANS of individual ratios (mean-of-ratios, not ratio-of-sums).
- Denominator: NETWORTH2R / INCFAMR with `income > 1000` filter — ANNUAL
  family income, real 2022$; negative net worth INCLUDED (no trimming; the
  gap's bootstrap SE is 0.682, i.e. 68% of the target — this moment is nearly
  uninformative and heavy-tail sensitive).
- Model side: `solver.py:4895-4966` builds cell ratios `bg/yj` with
  `yj = annual_gross_income_at_state(...)` (solver.py:148-156 divides period
  income by period_years, grosses up by 1/(1-tau) pre-retirement) and takes
  `weighted_median_from_cells` for the level and a difference of
  `weighted_mean_from_cells` for the gap. So the previously flagged
  period-vs-annual and ratio-of-sums issues (TARGET_MOMENT_OBJECTS
  "needs-fix"; also still in the working-tree
  `docs/model/intergen_target_object_audit.md`) appear FIXED in the code the
  overnight ran; the ledger/doc entries are STALE. Residual model-side caveat:
  `yj` is evaluated at z=1.0 only (income heterogeneity in the denominator is
  ignored, relevant for the partially working 65-66 period) — for task #6.
- SEs: median 0.0886; gap 0.6819 (family-cluster, 4,263 clusters).

### 1.8 young_childless_renter_liquid_wealth_to_annual_gross_income_2535 = 0.17922556 (PSID)

- Original builder: `code/data/psid_followup_mar2026/entry_wealth_targets_v1.do`
  → `output/entry_wealth_v1/entry_wealth_candidate_targets_focus_v1.csv` row
  `young_childless_rent_25_35, liq_nw_to_inc, .1792255598714484, .061549962, 7163`
  (chain verified line-by-line in
  `output/model/fable_size_mapping_audit_20260701/WEALTH_MOMENT_DIAGNOSIS.md`).
- Independent reproduction: SE harness (`build_moment_bootstrap_se.R:300-311`)
  from PSIDSHELF_MOBILITY, ages 25-35, childless (RELCHIREP-based n_children==0),
  renter (`own == 0`, HOMEOWN==2), `income > 1000`, IW-weighted mean of
  NETWORTH2R/INCFAMR = 0.179225561595449 — agrees with the Stata builder to
  1.7e-9. N bookkeeping differs (Stata reports 7163; harness counts 5453
  non-missing-ratio rows) but the estimate is identical; not substantive.
- Definitions: weighted MEAN (the weighted median is 0.099967; the CSV median
  0.061550 is unweighted — documented in calibration.py:261-266); ANNUAL gross
  family income denominator; negative net worth INCLUDED (entry-node grid in
  calibration.py:33-36 includes -2.519).
- Model side: `solver.py:4456` uses `annual_gross_income_at_state` — the old
  period-income mismatch is fixed for this target key.
- Note: the same 0.17922556 number does double duty as (i) an SMM target and
  (ii) the external entrant-wealth distribution mean
  (`external_entry_wealth_overrides`, calibration.py:48-56) — the target is
  partially hard-wired into the model's entry condition, which weakens its
  identifying content (flag for the identification adjudicator).
- SE: 0.0955 (53% of the target).

## 2. Answers to the specific consistency questions

(b) Room-unit / head / weight consistency: all five room-level targets
(aggregate + 3 prime30_55 + 3plus gap) come from ONE builder pass on the same
sample (heads PERNUM==1, HHWT, literal ROOMS>0, MMS 42-metro sample,
middle→center, 2012-2023) — internally consistent. The four OWNERSHIP targets
use a stricter sample (heads PERNUM==1 & RELATE==1, UNITSSTR 3:10, rooms>0):
head definition and structure universe differ between the room block and the
ownership block. PSID event-study rooms are national PSID ACTUALROOMS, a
different survey/universe from the ACS room levels.

(c) Mean vs median: all room and ownership targets are means/shares;
old_nonhousing level is a weighted median; old gap is a difference of weighted
means; young wealth is a weighted mean; increment_0to1 is a regression
coefficient. Model statistics use matching functional forms per solver.py
(median via weighted_median_from_cells; means via weighted_mean_from_cells).

(d) Annual vs 4-year denominators: all three PSID wealth-ratio targets use
ANNUAL INCFAMR; the model now converts period income to annual gross via
`annual_gross_income_at_state` for exactly the three target keys in the active
set. The stale period-denominator keys (`young_liquid_wealth_to_income`,
`old_nonhousing_wealth_to_income_6575`) are no longer in the active set.
Doc `docs/model/intergen_target_object_audit.md` (even its working-tree
version) still labels the old gap "Needs fix" — stale relative to code.

(e) Negative wealth: included everywhere (no trimming, no winsorizing);
income>1000 is the only guard on the ratios. Consequence: the old gap SE is
68% of its target value.

(f) Provenance holes as of July 2 → status now: parent_3plus gap CLOSED
(reproduced from extract27 by the SE harness); tfr/childless CLOSED at
citation level only (CPS June 2024 H2 cell; no raw artifact in repo; SE file
still lists them as source-unconfirmed because it predates the July 8 note).
The NEW aggregate rooms target introduces a fresh, milder provenance issue:
its builder edit is uncommitted and its own artifact says "candidate;
re-audit ... before SMM use".

## 3. code/data/moment_standard_errors/ (untracked, new)

Bootstrap SE + covariance harness for the 14 post-audit moments, spec'd
2026-07-05 (`SPEC_moment_bootstrap_se.md`), outputs written 2026-07-06.
Read-only over the existing builders; hard reproduction gate (12/14 PASS,
tfr/childless FAIL-no-source). Design: B=1000, seed 20260705, PSID
family-cluster bootstrap, ACS metro(met2013)-cluster bootstrap (42 clusters),
event-study analytic clustered SE; 14x14 block-diagonal covariance by source
(`moment_covariance.csv`, `moment_correlation.csv`); plus
childlessness-by-income and young-wealth-by-era tables.

## 4. SEs vs the overnight's ad hoc weights

Objective is sum_m w_m (model_m - target_m)^2 in raw units. Inverse-variance
(diagonal-optimal) weights implied by moment_se.csv vs the weights actually
used (CANDIDATE_REPLACEMENT_POST_AUDIT_V1_WEIGHTS + fixed-spec extras):

| moment | boot SE | 1/SE^2 | used w | used/(1/SE^2) |
|---|---:|---:|---:|---:|
| own_family_gap | 0.00650 | 23,650 | 45 | 0.0019 |
| old_age_own_rate | 0.01259 | 6,313 | 160 | 0.025 |
| prime_owner_share_ge6 | 0.01752 | 3,258 | 25 | 0.0077 |
| own_rate_2534 | 0.01820 | 3,018 | 80 | 0.027 |
| own_rate | 0.02037 | 2,411 | 100 | 0.041 |
| 3plus_minus_1to2 rooms | 0.06658 | 226 | 8 | 0.035 |
| owner-minus-renter rooms | 0.07137 | 196 | 12 | 0.061 |
| renter mean rooms | 0.07379 | 184 | 6 | 0.033 |
| old nh median | 0.08863 | 127 | 0.8 | 0.0063 |
| young wealth | 0.09553 | 110 | 12 | 0.11 |
| increment_0to1 | 0.14851 | 45 | 14 | 0.31 |
| old nh gap | 0.68191 | 2.15 | 2.0 | 0.93 |
| tfr / childless | NA | — | 20 / 20 | — |
| aggregate rooms | NA | — | 6 | — |

Relative weighting differs from inverse-variance by a factor of ~500 across
moments (own_family_gap vs old gap). In particular the noisiest moment (old
gap, SE = 68% of target) is the one closest to its optimal weight, while the
most precisely measured moments (family gap, old own rate, 2534 own rate) are
relatively downweighted by 1-3 orders of magnitude. The covariance file also
shows within-source correlations that a diagonal W ignores. Consequence: the
overnight point estimate is a valid (consistent) but efficiency-arbitrary
GMM/SMM point; no standard errors or J-statistics can be quoted from this
loss, and the loss ranking across candidates embeds an arbitrary weighting
choice made before the SE file existed. The SE file was available (July 6)
and unused (July 10 spec).

## 5. VERIFIED items (data side)

1. aggregate rooms target value = builder CSV = overnight driver constant
   (byte-exact, 5.779970481941968).
2. 12/14 post-audit targets reproduce from in-repo builders to <5e-7
   (reproduction_log.txt) — independently of my checks.
3. tfr halving convention implemented (`"tfr": 2.0 * household_parity`).
4. ACS ROOMS top-coding immaterial (range 1-26, weighted mass at max ~8e-6).
5. Model ages 18-85 match the new data filter (J=17 x 4y from 18).
6. Prime30_55 room constants unchanged by the July 10 builder edit
   (regenerated CSV matches the old constants to all printed digits).
7. Annual-income denominators implemented model-side for the three active
   PSID wealth targets (`annual_gross_income_at_state`).
8. old_age_own_rate source is ACS/MMS DUE-housing heads 65-75 (not PSID;
   PSID cross-check 0.8637 does not match).

## 6. Open questions / blocked

- Cluster-side check that the torch snapshot contained the same uncommitted
  builder/calibration diffs: BLOCKED (SSH unavailable). Lead reproduced the
  candidate loss locally, which indirectly confirms the model-side moment
  code; the data-side CSV is local-only either way.
- Which MMS data suffix (`MMS_DATA_SUFFIX`) was set for the July 10 builder
  run is not recorded in the output README; the reproduced prime30_55 values
  imply the default `data/` lookups (same as all prior targets), but the env
  is not logged.
- Whether the paper text describes aggregate_mean_occupied_rooms_18_85 as
  "US aggregate" (would be wrong) — paper-staleness auditor's scope.
