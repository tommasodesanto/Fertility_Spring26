# Empirical Roundup Prompt

## Objective

Audit what is already empirically established for the spatial-fertility theory block, what is already exportable from the repo, what exists on disk but is not yet clean enough to cite, and what is genuinely missing.

Do **not** start by inventing new data merges. The first task is an inventory tied to the actual files in this repo.

## Important Corrections

- The working PSID file in the current pipeline is `PSIDSHELF_MOBILITY.dta`.
- In the current PSID shelf / workflow, there is **no usable MSA/CBSA panel** for the proposed first-birth-to-cheaper-MSA design. Do **not** propose or implement an HUD FMR / ACS B25064 merge keyed by MSA unless you first verify that a metro identifier exists in the actual working data. The current code does not use one.
- The spatial center/periphery evidence in this project comes from the ACS/PUMA-based `mms_center_periphery` build, not from PSID geography.
- Many PSID event-study figures already exist in legacy `../../Outputs/Graphs` and `../../Outputs/Tables`. Treat those as existing results, not as missing results.

## Theory Sources To Audit Against

Use these as the canonical theory-side source files:

- `../latex/theoryexps/simple_spatial.tex`
- `../latex/theoryexps/claude_slides.tex`
- `dynamic_housing_policy_theory_note_v9.tex`
- `../latex/april_20_project_presentation.tex`

The clean prediction set from `simple_spatial.tex` is:

1. Sorting:
   - Parent share falls with local rent, conditional on income.
   - First birth triggers relocation toward lower-price areas.
2. Elasticity attenuation:
   - Parents have lower housing income elasticity.
   - Parents have lower housing price elasticity.
3. Wealth gradient of flight:
   - First-birth relocation / family flight is decreasing in pre-birth wealth.
   - There is a wealth threshold above which parents stay central.

The clean mechanism prediction from `claude_slides.tex` is:

4. Low-slack nonlinearity:
   - The first-child housing cost becomes sharply nonlinear at low pre-child slack on fixed branches.

## Repo Orientation

### Main slide deck

- `../latex/april_20_project_presentation.tex`

### ACS center/periphery evidence

- `../code/data/mms_center_periphery/output_middle_center/`

Most useful files:

- `mms_location_by_parent_summary.csv`
- `mms_origin_transition_summary_allparent.csv`
- `mms_move_regressions_allparent.csv`
- `mms_origin_regressions_allparent.csv`
- `mms_four_way_move_shares_allparent.csv`
- `mms_location_by_parent.png`
- `mms_grouped_4way_decomposition_allparent.png`
- `mms_within_cbsa_origin_destination_allparent.png`

### PSID follow-up evidence

- `../code/data/psid_followup_mar2026/`
- `../code/data/psid_followup_mar2026/output/`

Most useful files:

- `output/sa_rooms_first_birth_one_variant_v1/rooms_f_c_y_all_summary.csv`
- `output/sa_rooms_first_birth_grouped_v1/rooms_f_c_y_one_kid_by3_summary.csv`
- `output/sa_rooms_first_birth_grouped_v1/rooms_f_c_y_two_plus_by3_summary.csv`
- `output/rooms_first_birth_one_vs_two_horizon_v1/rooms_first_birth_one_vs_two_horizon_v1.csv`
- `output/sa_rooms_second_birth_with_onechild_controls_v1/rooms_s_c_y_all_summary.csv`
- `output/fertility_wealth_v1/within_wealth_gap_v9.csv`
- `output/fertility_wealth_v1/lpoly_parenthood_v12.png`
- `output/fertility_wealth_v1/parenthood_by_quintile_v12.csv`
- `output/wealth_v1/wealth_tercile_outcomes_v1.csv`
- `output/wealth_v1/wealth_transition_regressions_v1.tex`
- `output/moved_to_own_state_change_v1/moved_to_own_state_change_summary_v1.csv`
- `output/moved_for_size_iv_alltenure_v1/moved_for_size_design_summary_v1.csv`
- `output/moved_for_size_iv_alltenure_v1/moved_for_size_second_stage_ols_iv_v1.tex`
- `output/moved_for_size_iv_alltenure_v1/moved_for_size_fstats_v1.csv`

### Legacy PSID graphs still used by the slide deck

- `../../Outputs/Graphs/own_f_c_y_all.png`
- `../../Outputs/Graphs/rooms_f_c_y_all.png`
- `../../Outputs/Graphs/mv_s_f_c_y_all.png`
- `../../Outputs/Graphs/mv_n_f_c_y_all.png`
- `../../Outputs/Tables/own_f_c_y_all_estimates.dta`
- `../../Outputs/Tables/rooms_f_c_y_all_estimates.dta`
- `../../Outputs/Tables/mv_s_f_c_y_all_estimates.dta`

There is also a cleaned ownership replication figure here:

- `../code/data/psid_followup_mar2026/output/sa_replication/own_f_c_y_all_repl.png`

## What Is Already Established And Exportable

Treat the following as the baseline empirical roundup unless inspection of the files contradicts it.

### 1. ACS: parents are less central than non-parents

Status: **tested, slide-ready**

Files:

- `../code/data/mms_center_periphery/output_middle_center/mms_location_by_parent_summary.csv`
- `../code/data/mms_center_periphery/output_middle_center/mms_location_by_parent.png`

Current headline numbers:

- Non-parents center share: `0.494`
- New parents center share: `0.416`

This is the cleanest empirical counterpart to the basic sorting prediction.

### 2. ACS: movers with children sort more toward the periphery within metros

Status: **tested, slide-ready**

Files:

- `../code/data/mms_center_periphery/output_middle_center/mms_origin_transition_summary_allparent.csv`
- `../code/data/mms_center_periphery/output_middle_center/mms_move_regressions_allparent.csv`
- `../code/data/mms_center_periphery/output_middle_center/mms_origin_regressions_allparent.csv`
- `../code/data/mms_center_periphery/output_middle_center/mms_within_cbsa_origin_destination_allparent.png`

Current headline numbers:

- Among center-origin movers, center destination share:
  - non-parents `0.670`
  - parents `0.580`
- Among periphery-origin movers, center destination share:
  - non-parents `0.374`
  - parents `0.263`
- Regression:
  - parent indicator lowers center destination among within-CBSA movers by about `11.4 pp`
  - parent indicator lowers across-CBSA moving by about `5.2 pp`

This already gets most of the way to the “parents re-sort toward peripheral family space” claim.

### 3. PSID: first birth raises housing demand sharply

Status: **tested, slide-ready**

Files:

- `../code/data/psid_followup_mar2026/output/sa_rooms_first_birth_one_variant_v1/rooms_f_c_y_all_summary.csv`
- `../code/data/psid_followup_mar2026/output/rooms_first_birth_one_vs_two_horizon_v1/rooms_first_birth_one_vs_two_horizon_v1.csv`
- `../../Outputs/Graphs/rooms_f_c_y_all.png`

Current headline numbers:

- First-birth rooms event study:
  - `+0.664` rooms at `k=+3`
  - `+0.843` rooms at `k=+5`
- Simple one-vs-two horizon comparison:
  - one kid by +2: `+0.259`
  - two-plus by +2: `+1.368`

### 4. PSID: the first-birth housing step is much larger when fertility continues quickly

Status: **tested, exportable**

Files:

- `../code/data/psid_followup_mar2026/output/sa_rooms_first_birth_grouped_v1/rooms_f_c_y_one_kid_by3_summary.csv`
- `../code/data/psid_followup_mar2026/output/sa_rooms_first_birth_grouped_v1/rooms_f_c_y_two_plus_by3_summary.csv`

Current headline numbers:

- one-kid by +3: `+0.289`
- two-plus by +3: `+1.031`

This is useful for documenting that the family-space response scales with realized fertility path.

### 5. PSID: second birth raises housing demand again

Status: **tested, slide-ready**

Files:

- `../code/data/psid_followup_mar2026/output/sa_rooms_second_birth_with_onechild_controls_v1/rooms_s_c_y_all_summary.csv`
- `../code/data/psid_followup_mar2026/output/sa_rooms_second_birth_with_onechild_controls_v1/rooms_s_c_y_all.png`

Current headline numbers:

- second-birth rooms event study:
  - `+0.566` rooms at `k=+3`
  - `+0.576` rooms at `k=+5`

### 6. PSID: ownership transition is the central fertility margin among young renters

Status: **tested, slide-ready**

Files:

- `../code/data/psid_followup_mar2026/output/fertility_wealth_v1/within_wealth_gap_v9.csv`
- `../code/data/psid_followup_mar2026/output/fertility_wealth_v1/childless_by_transition_v7_bar.png`
- `../latex/april_20_project_presentation.tex` around the ownership-gap slides

Current headline interpretation already in the slide deck:

- among young pre-birth renters, those who transition to ownership by 35 are about `32 pp` less likely to remain childless by 45
- the gap persists at every wealth quintile
- this is already the strongest “constraint margin” fact in the repo

### 7. PSID: unconditional wealth-fertility relationship is non-monotone / hump-shaped

Status: **tested, slide-ready**

Files:

- `../code/data/psid_followup_mar2026/output/fertility_wealth_v1/lpoly_parenthood_v12.png`
- `../code/data/psid_followup_mar2026/output/fertility_wealth_v1/parenthood_by_quintile_v12.csv`

This is already in the slide deck and is important because it warns against over-claiming a monotone wealth mechanism from reduced-form plots.

### 8. PSID: moves for size are already in the empirical stack

Status: **tested, but split across legacy and follow-up folders**

Files:

- `../../Outputs/Graphs/mv_s_f_c_y_all.png`
- `../../Outputs/Tables/mv_s_f_c_y_all_estimates.dta`
- `../code/data/psid_followup_mar2026/output/moved_for_size_iv_alltenure_v1/`

Interpretation:

- The event-study fact “birth raises moves for space” is already part of the slide deck.
- The IV follow-up exists as robustness / mechanism support, not as the main headline:
  - twin first-stage F-stat is strong (`124.7`)
  - same-sex first-stage is modest (`5.8`)

Do not treat the IV as the core empirical fact. The core fact is the event-study pattern already in the deck.

## What Exists But Should Not Be Treated As Clean Evidence Yet

### 1. First-birth hazard by wealth files

Files:

- `../code/data/psid_followup_mar2026/output/fertility_wealth_v1/first_birth_hazard_by_tercile.csv`
- `../code/data/psid_followup_mar2026/output/fertility_wealth_v1/first_birth_hazard_by_quintile.csv`

These currently collapse to zeros and should be treated as **not exportable** until checked.

### 2. Wealth post-birth transition summary

Files:

- `../code/data/psid_followup_mar2026/output/wealth_v1/wealth_tercile_outcomes_v1.csv`
- `../code/data/psid_followup_mar2026/output/wealth_v1/wealth_transition_regressions_v1.tex`

These are promising and informative, but some raw tercile patterns look odd enough that they should be treated as **exploratory** until interpreted carefully.

## What Is Not Currently Established

### 1. First birth triggers moves to cheaper MSAs

Status: **not established; not feasible in current PSID workflow**

Reason:

- the current PSID shelf / scripts do not give a usable metro panel for this design
- current geography evidence is ACS center/periphery, not PSID MSA-to-MSA moves

If you want a spatial-mobility claim, use the ACS within-metro mover evidence already on disk.

### 2. Parent vs childless housing price elasticity

Status: **not established**

Reason:

- no live MSA rent panel in the PSID workflow
- no existing cleaned price merge in the repo

Do not claim this is “easy” with current data.

### 3. Parent vs childless housing income elasticity

Status: **not yet done cleanly, but potentially feasible**

Reason:

- `INCFAMR`, `rooms`, and parental status are in the PSID workflow
- this could be estimated as a rooms or ownership elasticity proxy
- but it is still new work, not an already-completed result

### 4. Wealth gradient of spatial flight

Status: **not established in the spatial sense**

Reason:

- current data do not support a clean pre-birth-wealth × metro-direction design
- what *is* already established is wealth × ownership / fertility, not wealth × center-to-periphery relocation at birth

### 5. Slack-dependent first-birth hazard / power-law squeeze

Status: **not established**

Reason:

- rooms and fertility timing are available
- but there is no finished slack-hazard design on disk
- there is definitely no finished IV version on disk

## Feasible Next Tests With Current Data

If you are asked to extend beyond the audit, prioritize only designs that match the actual data.

### Priority A. Wealth gradient in first-birth housing adjustment

Best feasible extension.

Use pre-birth wealth bins and interact them with first-birth event studies for:

- ownership
- rooms
- move-for-size

This does **not** identify “flight to cheaper metro,” but it does test whether low-wealth households make larger housing adjustments at first birth.

### Priority B. Slack gradient in first-birth timing, descriptive only

Feasible with current PSID variables:

- define pre-birth slack from rooms among childless households
- estimate a descriptive first-birth hazard by lagged rooms / slack bins
- renter-only version is probably the right first pass

Do **not** oversell this as a structural validation unless the design is very clean.

### Priority C. State-change composition around first birth

Only if you need an extra mobility fact.

Files already show state variables exist:

- `GEOSTATE`
- `STATEFIPS_`

You can study:

- whether birth-related moves are mostly within-state vs interstate
- whether move-to-own is mostly local-ish rather than interstate

But this is **not** a substitute for a cheaper-location or center/periphery spatial test.

### Priority D. Balance-sheet response around first birth

Potentially useful and data-consistent, using:

- `NETWORTHR`
- `NETWORTH2R`
- `HOMEEQUITYR`
- `HOMEMORTOTR`
- `WLTHSAVETOTR`

This could sharpen the “ownership as financing technology for family space” argument.

### Priority E. Income-elasticity attenuation using PSID rooms

Medium priority.

If attempted, be explicit that:

- this is an income elasticity of a proxy for housing services
- price elasticity remains untested in the current setup

## Deliverables

Produce a single audit document with one row per prediction and these columns:

1. theory prediction
2. current status
3. best existing file(s)
4. exact headline number if already available
5. slide-ready / exportable / exploratory / broken / missing
6. best next step

Then add a short narrative with four buckets:

- already in slides and fully supported
- on disk and exportable with minimal cleanup
- on disk but not yet trustworthy
- genuinely missing

## Working Rules

- Start with the existing files above. Do not duplicate work that is already on disk.
- If a legacy graph already exists, note it and only re-export it into a cleaner folder if needed.
- Do not claim a PSID spatial design that requires MSA identifiers.
- Do not treat the zero-valued hazard summaries as valid without debugging.
- If you extend the analysis, stop after the first useful extension and report back before branching out.

## Short Version

The corrected bottom line is:

- ACS already gives the spatial sorting facts.
- PSID already gives the first-birth housing-step, second-birth housing-step, ownership-constraint, and move-for-size facts.
- The true missing pieces are elasticity attenuation, slack-hazard nonlinearity, and any wealth gradient of **spatial** flight.
- With the current data, the best next test is not an MSA-rent merge. It is wealth-heterogeneity in first-birth housing adjustment, plus possibly a descriptive slack-hazard exercise.
