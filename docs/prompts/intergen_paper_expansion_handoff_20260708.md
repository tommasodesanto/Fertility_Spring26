# Intergenerational Housing-Fertility Paper Handoff

Date: 2026-07-08

## Objective

Expand `latex/intergenerational_housing_fertility_paper_draft.tex` into a serious paper draft. Keep the pure empirical sections out for now, but include the non-empirical core:

- the Part I sufficient-statistic theory,
- the quantitative environment,
- the solution and equilibrium object,
- calibration and identification,
- baseline solution graphs,
- baseline fit and known misses,
- mechanism diagnostics,
- policy counterfactuals.

This is a paper-facing writing and correctness pass. Do not smooth over mathematical gaps. If a claim is not supported by the live model or by the live calibration status, either fix the statement or flag it.

Do not use Oracle for this handoff unless the user explicitly asks for it. This prompt already incorporates a separate model review.

## Mandatory Startup

Before editing, read these in this order:

1. `memory/AGENT_MEMORY.md`
2. the latest `memory/daily/YYYY-MM-DD.md`
3. `CALIBRATION_STATUS.md`
4. `docs/style/econ_writing_style_guide.md`

Then read the source files listed below. If the status files disagree with older notes, treat `CALIBRATION_STATUS.md` as the live source and report the discrepancy.

## Files To Edit

Primary file:

- `latex/intergenerational_housing_fertility_paper_draft.tex`

Likely supporting files:

- `latex/lit_review_extra.bib`, only if citation keys need to be added or corrected.
- `latex/intergenerational_housing_fertility_paper_draft.pdf`, rebuilt from the `.tex`.

Do not edit model code for this task unless the paper exposes a genuine code/equation mismatch that must be checked. This pass is about the paper draft.

## Source Files To Read

Current draft:

- `latex/intergenerational_housing_fertility_paper_draft.tex`
- `latex/intergenerational_housing_fertility_paper_draft.pdf`

Part I theory sources:

- `latex/intergenerational_housing_fertility_part1_sufficient_stats.tex`
- `latex/intergenerational_housing_fertility_part1_sufficient_stats.pdf`
- `latex/intergen_housing_fertility_short_note.tex`
- `latex/intergen_housing_fertility_short_note.pdf`

Quantitative model and live status:

- `CALIBRATION_STATUS.md`
- `code/model/intergen_housing_fertility/README.md`
- the active files named by `CALIBRATION_STATUS.md` if equations need verification.

July deck and figures:

- `latex/July_26_slides.tex`
- `latex/July_26_slides.pdf`
- `latex/July_26_slides_eq_policies.png`
- `latex/July_26_slides_eq_consumption_by_tenure.png`
- `latex/July_26_slides_eq_fertility.png`
- `latex/July_26_slides_eq_housing.png`
- `latex/July_26_slides_eq_wealth.png`
- `latex/July_26_slides_eq_tenure.png`
- `latex/July_26_slides_eq_rungs.png`
- `latex/July_26_slides_eq_market.png`
- `latex/July_26_slides_age_profiles.png`
- `latex/July_26_slides_buy_flow_by_age_children.png`

Current numerical memos:

- `output/model/fable_size_mapping_audit_20260701/MORNING_NOTE_20260708.md`
- `output/model/fable_size_mapping_audit_20260701/MORNING_MEMO_20260707.md`
- `output/model/fable_size_mapping_audit_20260701/MECHANICS_MEMO_20260706.md`

Live target-fit CSV:

- `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26_fable_size_mapping_audit_20260701/output/model/diag_packet_nb120_best_20260707/target_fit.csv`

Policy CSV to treat cautiously:

- `output/model/intergen_policy_counterfactuals_current/wealthstress_diagnostic/policy_summary.csv`

The policy CSV appears to be an older or diagnostic wealth-stress run, not necessarily the headline July deck policy result. Do not use it as the headline policy table unless you reconcile it against `CALIBRATION_STATUS.md`, the July deck, and the current output provenance.

## Critical Correctness Issues To Fix First

1. The quantitative state must include tax basis when the capital-gains block is active.

   The current draft defines sale proceeds, adjustment wealth, and lock-in using a basis object, but the formal state vector omits it. Owner values, policies, and the invariant distribution are not well defined unless the state includes \(B^p\), or unless the gains block is explicitly switched off. Fix the environment and Bellman equations before adding policy results that rely on lock-in.

2. The quantitative entry and local-birth closure cannot be overstated.

   The live status says the published intergenerational runs are fixed-population. Entry and scale margins are driver-level or off-model accounting, not an embedded equilibrium entrant problem with an outside option and entrant law. If Part II keeps local births, entrants, or market scale, either write down that full quantitative closure or explicitly label it as off-model accounting / appendix material.

3. The active quantitative supply benchmark is unit-elastic, not fixed stock.

   `CALIBRATION_STATUS.md` on 2026-07-08 says the benchmark uses unit-elastic supply with `H0=4` and `xi=1`. The fixed-stock theory is a useful simplification or limiting case, not the active quantitative benchmark. Do not write the quantitative section as if fixed stock is the baseline.

4. Clarify physical housing units versus service units.

   The current draft says there is one aggregate housing-services market, but owners consume \(\chi_O H_k\) while market clearing uses physical \(H_k\). Decide and state the convention. The simplest repair is likely: the market clears in physical housing units, and \(\chi_O\) is a tenure-specific service shifter in preferences. Then all clearing and supply equations should be physical-unit equations.

5. Do not reinterpret the completed-fertility model as a sequential hazard model.

   The live model is one-shot completed fertility. Housing responses around first birth or second birth discipline housing demand and child-cost blocks; they are not separate first-birth or second-birth hazard equations. Use "completed fertility," "childlessness," and "number of children" precisely.

6. Be honest about the sharp fertility-choice limit.

   `kappa_fert` is a denominator and is currently at the sharp-choice lower bound. Poor households choosing zero births in the calibration is not a bug by itself. Do not describe \(\kappa_F\) as an ordinary well-identified taste-shock parameter without noting that it is bound-pinned and load-bearing.

## Separate Model Review Incorporated Here

A separate GPT-5.4 review identified the following major issues.

Critical:

- The formal quantitative state omits \(B^p\), even though sale proceeds and lock-in depend on it.
- Part II talks about entry, local births, and market scale without a quantitative entry problem or demographic closure.
- The market object mixes physical housing and housing services inconsistently.

Major:

- Part I needs a numbered equilibrium definition listing objects and conditions.
- Part I should restore the older note's switching sign conditions and corollaries, either in text or appendix.
- The compact model should explicitly say that \(\vartheta \log n_i\) gives an intensive-margin model among entrants, while childlessness is handled outside that particular block.
- The quantitative section after calibration reads like a scaffold. Replace "should report" prose with actual paper text, tables, and figures.
- The supply discussion should reflect the live unit-elastic benchmark.
- The calibration section must reflect the live numerical status: `Nb=120` search, `Nb=240` verification, and the main remaining misses.
- The parameter table should separate state/institutional objects from estimated parameters.
- The literature review is now structurally close, but spatial references should be positioned as allocation logic, since the current draft is not a full spatial equilibrium model.

Minor:

- Avoid using superscript \(O\) both for old-stage objects and owner objects.
- Tie "births generated by a cohort" to Part I or to off-model accounting, not to the embedded quantitative model unless the closure is added.
- Add prose around what is fixed externally, what is estimated, what closure is used, and what is deliberately off-model.

## Writing And Structure Instructions

Follow `docs/style/econ_writing_style_guide.md`.

Use an economics-paper structure:

1. Introduction
2. Literature review
3. Part I: compact theory / sufficient statistics
4. Part II: quantitative model
5. Solution and equilibrium computation
6. Calibration and identification
7. Baseline fit and mechanism diagnostics
8. Policy counterfactuals
9. Conclusion

Do not add a table of contents.

For model sections:

- Start with environment primitives.
- Then write the household problem.
- Then define equilibrium with numbered objects and conditions.
- Put one sentence of economics before and after every important display.
- Use explicit scoping sentences for local formulas, fixed-price results, off-model accounting, and policy exercises.

## Part I Tasks

Keep the Part I sufficient-statistic material because it is the closest existing theory draft to a real paper object. But clean the structure.

Required repairs:

1. Add a proper equilibrium definition.

   Define equilibrium objects first, then conditions. Include prices, young choices, old choices, entrant mass, distributions, clearing, and the closure \(N+\bar R=\bar M\) if that closure remains in the compact model.

2. State the scope of fertility in Part I.

   The \(\log n_i\) block is an intensive-margin model among entrants, not a childlessness model. Make this a sentence in the text, not an apologetic aside.

3. Restore the policy-relevant switching results.

   The current draft has the boundary-integral object, but not enough sign content. Bring back the older note's conditions and corollaries from `latex/intergenerational_housing_fertility_part1_sufficient_stats.tex`, roughly around lines 1183-1373:

   - when owner-access policies add fertility through tenure switching,
   - when rental-access policies are offset by households switching back into renting,
   - how the property-tax / lock-in capitalization logic maps into the sufficient statistic.

4. Put demography-literature details in a footnote.

   The literature section should have one broad fertility bucket. Demography references can be mentioned in a footnote along the lines of "The issue also appears in the demography literature..." rather than becoming a separate strand that reads like a duplicate.

## Part II Tasks

Rewrite the quantitative section as a complete environment, not as a placeholder.

Required objects:

1. Agents, timing, and lifecycle.

   State the ages / stages, death or aging transitions, initial distribution, and whether the benchmark is fixed-population.

2. State vector.

   Include assets, income state, age, fertility/children state, tenure, housing rung, and \(B^p\) when the tax basis / gains block is active. If the code uses a compressed basis convention, explain that convention exactly.

3. Preferences.

   Write utility over consumption, housing services, fertility / children, and any child costs. Interpret \(\chi\), child costs, and tenure service shifters.

4. Housing and tenure.

   Define renter options, owner rungs, the rental cap, owner service shifter, transaction/sale terms, mortgage/down-payment terms, and property-tax / capital-gains institutions.

5. Markets and supply.

   State that the live benchmark has unit-elastic supply. If using physical units, write supply and clearing in physical units and treat \(\chi\) as a service shifter.

6. Household problem.

   Present the Bellman problem compactly but completely. Include fertility choice, tenure choice, owner/renter continuation values, sale/adjustment wealth, and borrowing/down-payment constraints. Avoid overlong code-like notation, but do not omit state variables needed for correctness.

7. Equilibrium.

   Define an equilibrium with objects and conditions:

   - prices,
   - household value and policy functions,
   - invariant distribution,
   - transition law / KFE,
   - housing market clearing under unit-elastic supply,
   - fixed-population closure for the live benchmark,
   - any off-model entry or scale accounting separated from the embedded equilibrium.

## Solution And Verification Section

Add a short paper section that says how the model is solved.

Include:

- Bellman iteration / backward or stationary recursion as appropriate.
- KFE or invariant-distribution computation.
- Price search / market-clearing object.
- Grid protocol: search at `Nb=120`, verify at `Nb=240`.
- Diagnostic packet: policy functions, consumption by tenure, fertility policy, housing choices, wealth distribution, tenure by age, owner rungs, market clearing, age profiles, buy flow by age and children.
- Verification caveats from live status: the old `Nb=60` canonical point is not the live benchmark; fixed-theta comparisons across numerical regimes are not valid.

Use the July deck figures as paper figures where appropriate. The graph files are currently in `latex/` root, not in `latex/figures/`.

## Calibration And Identification Section

Replace the placeholder calibration discussion with actual content.

Lead with the live scalar loss and then show the full target-fit table. Do not report selected moments only.

Live target-fit table from:

`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26_fable_size_mapping_audit_20260701/output/model/diag_packet_nb120_best_20260707/target_fit.csv`

| moment | target | model | gap | weight | loss contribution |
|---|---:|---:|---:|---:|---:|
| old_age_own_rate | 0.764 | 0.893 | 0.128 | 160.0 | 2.635 |
| own_rate_2534 | 0.341 | 0.226 | -0.115 | 80.0 | 1.055 |
| prime30_55_childless_owner_minus_renter_mean_rooms | 2.419 | 2.151 | -0.268 | 12.0 | 0.862 |
| own_family_gap | 0.168 | 0.298 | 0.130 | 45.0 | 0.765 |
| young_childless_renter_liquid_wealth_to_annual_gross_income_2535 | 0.179 | 0.353 | 0.173 | 12.0 | 0.361 |
| own_rate | 0.575 | 0.531 | -0.045 | 100.0 | 0.200 |
| prime30_55_childless_renter_mean_rooms | 3.805 | 3.655 | -0.150 | 6.0 | 0.136 |
| childless_rate | 0.188 | 0.260 | 0.072 | 20.0 | 0.103 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.368 | 0.269 | -0.099 | 8.0 | 0.079 |
| old_parent_childless_nonhousing_wealth_to_income_gap_6575 | 1.007 | 0.842 | -0.166 | 2.0 | 0.055 |
| housing_increment_0to1 | 0.664 | 0.618 | -0.046 | 14.0 | 0.030 |
| tfr | 1.918 | 1.949 | 0.031 | 20.0 | 0.019 |
| old_nonhousing_wealth_to_income_median_6575 | 2.230 | 2.328 | 0.097 | 0.8 | 0.008 |
| prime30_55_childless_owner_share_rooms_ge6 | 0.596 | 0.593 | -0.003 | 25.0 | 0.000 |

The contributions sum to about 6.31.

Parameter table:

- Use `output/model/fable_size_mapping_audit_20260701/MORNING_MEMO_20260707.md`, section 1b.
- Report fixed/external parameters separately from estimated parameters.
- Include estimate, lower bound, upper bound, and near-bound flag.
- Current status: eight of thirteen estimated parameters are at or against bounds. This is not a throwaway caveat; it belongs in the calibration discussion.

Live parameter/bounds note:

- `beta_annual = 0.9400`, lower bound, declared restriction.
- `alpha_cons = 0.620`, interior.
- `c_bar_0 = 1.279`, upper bound.
- `c_bar_n = 0.116`, interior.
- `h_bar_0 = 1.000`, lower bound.
- `h_bar_jump = 2.120`, interior-high.
- `h_bar_n = 1.382`, interior.
- `psi_child = 0.349`, upper bound.
- `kappa_fert = 1.019`, lower bound.
- `tenure_choice_kappa = 0.000`, lower bound.
- `chi = 1.150`, upper bound.
- `theta0 = 0.0015`, lower bound.
- `theta_n = 0.972`, interior.

Also report:

- `phi = 0.80` fixed externally.
- `sigma = 2` fixed externally.
- `hR_max = 6.0` fixed externally from the DUE / p90 rule, with cap-profile robustness.

Cap profile:

| `hR_max` | best loss | status |
|---:|---:|---|
| 5.5 | 7.69 | converged |
| 6.0 | 6.31 | baseline external rule |
| 6.5 | 6.24 | still improving in the July 7 memo |
| 7.0 | 7.15 | wave A only |

Explain that the model mildly prefers 6.5 in the July 7 profile, but the external 6.0 rule remains defensible and within half a room of the data-preferred value.

## Baseline Fit And Mechanism Diagnostics

Add a paper-facing section that reports what the model gets right and where it misses.

Core honest takeaways:

- TFR and the first-child housing increment are no longer the main failure at resolved grids.
- The main misses are old-age ownership, ownership at ages 25-34, young childless renter wealth, and the owner-renter room gap.
- The model overpredicts childlessness relative to the target.
- The current calibration has several parameters at bounds, so some quantitative claims should be presented as disciplined but not final.

Mechanism facts from the July deck / status:

- Fertile renters at the cap have \(\zeta/q \approx 1.28\).
- The model's exposed fertility-response mass is about 0.734 with \(\zeta>0\).
- The model produces a birth-triggered deadline purchase margin: about 29 percent of the age-42 cohort buys, about 70 percent of buyers are parents, roughly half have a newborn, and buyers are close to the down-payment threshold.
- Large homes are disproportionately held by older/post-fertile households. For rooms `>= 6`, the data have about 72 percent held by ages 46+, while the model has about 87 percent; the data have about 55 percent held by households with no children at home, while the model has about 72 percent; fertile ages 30-45 are about 25 percent in the data and 13 percent in the model.

Use graphs to carry these facts where possible:

- `latex/July_26_slides_eq_policies.png`
- `latex/July_26_slides_eq_consumption_by_tenure.png`
- `latex/July_26_slides_eq_fertility.png`
- `latex/July_26_slides_eq_housing.png`
- `latex/July_26_slides_eq_wealth.png`
- `latex/July_26_slides_eq_tenure.png`
- `latex/July_26_slides_eq_rungs.png`
- `latex/July_26_slides_eq_market.png`
- `latex/July_26_slides_age_profiles.png`
- `latex/July_26_slides_buy_flow_by_age_children.png`

## Policy Counterfactuals

Add a policy section, but keep embedded-model GE policies separate from off-model entry/scale accounting.

Headline policy facts from the July deck / status:

- A 2 percent large-home tax / property-tax instrument lowers large-home prices by about 12.3 percent, but raises TFR only about 0.004 in the live headline.
- Birth-linked purchase grants have much larger fertility effects and limited price movement:
  - grant size 0.1: \(\Delta\)TFR about 0.01,
  - grant size 0.4: \(\Delta\)TFR about 0.13,
  - grant size 0.58: \(\Delta\)TFR about 0.26,
  - grant size 1.0: \(\Delta\)TFR about 0.49.
- The package of a large-home tax plus a birth-linked purchase grant gives about \(\Delta\)TFR \(=0.15\) and about an 11 percent decline in large-home prices.

Important scoping:

- Do not claim the old-owner lock-in wedge has a large direct TFR effect in the live quantitative model. The tax price effect is robust, but the TFR effect is small at the resolved grid.
- Do not present local-birth scale effects as an embedded equilibrium result unless the entry/scale closure is written into Part II.
- Report completed fertility, childlessness, ownership, prices, retention, and housing allocation if those objects are available.
- If exact policy tables come from `policy_summary.csv`, first reconcile that file against the July deck and status notes. It may be diagnostic rather than the live headline table.

## Literature Review Direction

The literature review should read like an economics paper:

> This paper contributes to several strands of literature. First, ...

Use three main buckets:

1. Fertility decline and the causes of fertility choices.

   Put housing-cost effects on fertility, completed fertility, and demography references in this bucket. The demography literature can be in a footnote if it would otherwise fragment the review.

2. Intergenerational allocation, lock-in, and household balance-sheet mechanisms.

   This bucket must cite the papers that are close to the paper's core logic, including Coven et al. on intergenerational allocation and Fonseca and Liu. Do not omit this literature.

3. Structural housing, tenure, spatial, and urban allocation models.

   Position spatial references carefully: the current draft borrows allocation logic and housing-market equilibrium discipline, but it is not yet a full spatial equilibrium model.

The review should not be a list of paper summaries. Each bucket needs to say what the existing literature does, what margin it leaves open, and how this paper contributes.

## Verification Before Final Handoff

After editing:

1. Compile the paper from `latex/`:

   ```bash
   latexmk -pdf -interaction=nonstopmode -halt-on-error intergenerational_housing_fertility_paper_draft.tex
   ```

2. Inspect the resulting PDF.

3. Check for:

   - undefined references,
   - missing citations,
   - overfull boxes above roughly 2pt,
   - figures that do not load,
   - notation collisions,
   - any sentence saying "should report," "will be disciplined," or similar scaffold language.

4. Run `git status -sb` and report only the files changed for this pass. Do not revert unrelated dirty files.

## Final Deliverables

Return:

- the updated `.tex`,
- the rebuilt PDF,
- a short change memo,
- any unresolved correctness concerns,
- the compile command and result.

Do not bury limitations. This paper is strongest when it is explicit about what is mathematically closed, what is calibrated, and what is only an accounting extension.
