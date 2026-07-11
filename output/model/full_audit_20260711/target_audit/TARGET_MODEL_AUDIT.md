# Target-System and Objective Audit (MODEL side) — July 10 Overnight Combined Spec

Auditor: target/objective specialist (task 6). Date: 2026-07-11.
Audited object: dirty tree at HEAD fe092ca (the dirty tree IS the audited object).
Active target system: `candidate_replacement_post_audit_v1` (14 moments,
`code/model/intergen_housing_fertility/calibration.py:170-179` targets,
`:494-503` weights) plus overnight extra target
`aggregate_mean_occupied_rooms_18_85 = 5.779970481941968` weight 6.0
(`tmp/overnight_combined_20260710/run_overnight_combined.py:27,157-158`).
15 moments, 14 free parameters (13 in `PRODUCTION_SEARCH_BOUNDS`,
`production_profile.py:39-53`, plus `H0` in [1,10],
`run_overnight_combined.py:156`).

## 0. Age-grid convention (verified)

`J=17`, `age_start=18`, `da=period_years=4`. Model age of index `j` is
`18 + 4j`: 18,22,26,30,34,38,42,46,50,54,58,62,66,70,74,78,82 (each index is a
4-year period, so index 16 covers 82–86).
`age_to_index` (`solver.py:4614-4616`):
```python
idx = int(round((float(age) - float(P.age_start)) / max(float(P.da), 1e-12)))
return int(np.clip(idx, 0, P.J - 1))
```
Window resolutions used by the active moments (period-start ages):

| label | indices | period-start ages | period coverage |
|---|---|---|---|
| 25–34 | 2..4 | 26,30,34 | 26–38 |
| 25–35 | 2..4 | 26,30,34 | 26–38 |
| 30–55 | 3..9 | 30,…,54 | 30–58 |
| 65–75 | 12..14 | 66,70,74 | 66–78 |
| 18–85 | 0..16 (all) | 18,…,82 | 18–86 |

No off-by-one bug relative to the rounding rule (no .5 rounding ties are hit:
round(1.75)=2, round(4.0)=4, round(4.25)=4, round(9.25)=9, round(11.75)=12,
round(14.25)=14). But every window is shifted right relative to its data label:
"25–34" contains no age-25 mass and includes households up to age 37+;
"65–75" starts at 66 and runs to 77+. This is a systematic window-drift
(≈ +1 to +3 years) shared by all age-windowed moments, not a coding error in
any single one. Index 12 (age 66) is exactly the first retired index
(`J_R = round((66-18)/4) = 12`), so the 65–75 window is 100% retirees.

Fertility window: births occur for `j+1 in [A_f_start, A_f_end] = [1,7]`
(`solver.py:3763`), i.e. decision ages 18–42 inclusive. Completed-fertility
statistics condition on `g[:, :, :, A_f_end:, ...]` = indices 7..16 = ages
46–82, strictly post-fertility (`solver.py:4662-4665`). No overlap. Verified.

Distribution used for statistics: `g_current` = post-decision realized
cross-section (`use_postdecision_current_distribution=True`,
`production_profile.py:166`; `solver.py:4014-4029`), except wealth moments
which use `asset_current` = beginning-of-period `b` conditioned on current
tenure (`solver.py:4030-4041`, `stats.wealth_moment_timing =
"beginning_of_period_b_conditioned_on_current_tenure"`, `solver.py:4094`).
Population mass is normalized to `N_target=1` (`parameters.py:128`,
`solver.py:4005-4008`); entry mass `1/J` per period and no mortality gives a
uniform age distribution.

"Childless" conditioning: `current_child_bin_dt(nn, cs, dep_last)`
(`solver.py:5739-5747`) returns bin 2 when `cs == 0 or cs > dep_last` — i.e.
never-parents AND empty-nest parents whose (stochastically maturing,
`use_stochastic_aging=True`, `parameters.py:37`) child has left. All
"childless" ROOM and YOUNG-WEALTH moments are childless-IN-HOUSEHOLD.
By contrast, the OLD-AGE wealth gap uses `nn > 0` = ever-parent vs
`nn == 0 and cs == 0` = never-parent (`solver.py:4920,4929`). Two different
"childless" definitions coexist by design; see M7/M10 caveats.

## 1. Moment-by-moment worksheet (15 rows)

Model counterpart line references: `calibration.py` `extract_moments`
(:1081-1190) maps solver attributes to moment names; solver statistics are in
`compute_statistics` (:4619-5101, applied to the income-collapsed `g_total`),
corrected/extended by `markov_renter_room_moments` (:4331-4395) and
`add_annual_gross_liquid_wealth_moments` (:4398-4463) on the full
income-resolved state inside `compute_markov_statistics` (:4466-4553).

| # | moment | target | weight | model counterpart (exact object) | age window (model) | conditioning | mean/median | denominator / units |
|---|---|---|---|---|---|---|---|---|
| M1 | tfr | 1.918 | 20 | `2.0 * sol.mean_completed_fertility` = `2*mean_parity`, `parity_dist` over `g[.., A_f_end:, nn, :]` (calibration.py:1093; solver.py:4661-4665) | ages 46–82 | all households post-fertility; `n_parity=3` so n∈{0,1,2}; each model child counts ×2 | mean (mass-weighted) | children ("completed-fertility-equivalent", ledger :207-212; NOT a period TFR) |
| M2 | childless_rate | 0.188 | 20 | `parity_dist[0]` same block (solver.py:4664; calibration.py:1100) | ages 46–82 | n==0 post-fertility (completed childlessness) | share | rate ∈[0,1] |
| M3 | own_rate | 0.57547241 | 100 | `own_rate_3055` = owner mass / total mass (solver.py:4829-4831; calibration.py:1094) | indices 3–9 (30–54) | all n, all cs; owner = tenure≥1 contemporaneous | mass share | rate |
| M4 | own_family_gap | 0.16766167 | 45 | `own_gap_newparent_nonparent_3055` = own(nn≥1 & cs∈{1}) − own(nn==0 & cs==0) (solver.py:4846-4866); `newparent_cs=[1]` = the single dependent-child stage (:4715) | 30–54 | "newparent" = parent w/ dependent child at home; "nonparent" = never-parent; empty-nesters in neither | share difference | rate diff |
| M5 | housing_increment_0to1 | 0.66443467 | 14 | `housing_increment_0to1_eventstudy_t3` = `(birth_es3_post_sum − birth_es3_control_post_sum)/birth_es3_mass` (solver.py:4070-4072); horizon 0 (`housing_event_horizon=0`, production_profile.py:22); birth cohort vs same-mass childless control, both advanced 0 periods and realized post-decision (solver.py:3779-3858) | birth ages 18–42 | new-birth cohorts, all tenures; realized rooms (renter hR or owner rung) | mass-weighted mean diff | rooms |
| M6 | old_nonhousing_wealth_to_income_median_6575 | 2.23046078 | 0.8 | `weighted_median_from_cells(old_nonhousing_ratio_vals, wts)` of `bg/yj`, `yj = annual_gross_income_at_state(P,i,jj,1.0)` = pension/4 (solver.py:4894-4916, 4965-4967; utils.py:237-240) | indices 12–14 (66–74), all retired | all tenures, parents+childless; b includes negatives (mortgage debt) | **median, uninterpolated grid value** | liquid b / ANNUAL pension income |
| M7 | old_parent_childless_nonhousing_wealth_to_income_gap_6575 | 1.00744952 | 2.0 | `parent_nonhousing_ratio − childless_nonhousing_ratio`, each a weighted MEAN of `bg/yj` (solver.py:4938-4961) | 66–74 | parent = nn>0 (any cs); childless = nn==0 & cs==0 (never-parent) | mean-of-ratios (both sides) | annual pension income |
| M8 | prime30_55_childless_renter_mean_rooms | 3.80528810 | 6 | `renter_rooms_3055_childless / renter_mass_3055_childless`, renter hR on income-collapsed policy (exact for means) (solver.py:5040-5050, 5080-5082) | 30–54 | ten==0, `current_child_bin_dt==2` (childless-in-household incl. empty-nest), hr>0 finite | mass-weighted mean | rooms (renter cap 6) |
| M9 | prime30_55_childless_owner_share_rooms_ge6 | 0.59613112 | 25 | `owner_rooms_ge6_3055_childless / owner_mass_3055_childless`; owner rung `H_own[ten-1] >= 6-1e-8` (solver.py:5060-5070, 5092-5094) | 30–54 | owners, childless-in-household | mass share (exact — rungs discrete per tenure index) | rate; H_own=[2,4,6,8,10] so ≥6 = top 3 rungs |
| M10 | young_childless_renter_liquid_wealth_to_annual_gross_income_2535 | 0.17922556 | 12 | mean of `bg/y`, `y = annual_gross_income_at_state` = period income/4/(1−tau_pay), z-resolved; dist = `asset_current` (solver.py:4419-4459, esp. :4427,4439-4441,4456-4458) | indices 2–4 (26–34) | ten==0 after current-choice realization; child_bin==2 | weighted mean of ratios | liquid b / annual GROSS labor income (grossed up by 1/(1−0.179)) |
| M11 | prime30_55_childless_owner_minus_renter_mean_rooms | 2.41876173 | 12 | `prime30_55_childless_owner_mean_rooms − prime30_55_childless_renter_mean_rooms` (solver.py:5086-5088) | 30–54 | childless-in-household, contemporaneous tenure | mean diff | rooms |
| M12 | old_age_own_rate | 0.76426097 | 160 | `old_age_own_rate_6575` = owner mass/total mass (solver.py:4868-4870; calibration.py:1098) | 66–74 | all n/cs | mass share | rate |
| M13 | own_rate_2534 | 0.34116609 | 80 | `own_rate_2534` (solver.py:4832-4834) | 26–34 | all n/cs | mass share | rate |
| M14 | prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.36769955881 | 8 | mean realized rooms bin4 (nn≥2 dependent) − bin3 (nn==1 dependent), renters+owners pooled (solver.py:5051-5079, 5098-5101) | 30–54 | current-parent (dependent child at home) only | mean diff | rooms; model nn=2 (top parity ≈ "4 kids" under the ×2 convention) maps to data 3+; nn=1 maps to data 1–2 |
| M15 | aggregate_mean_occupied_rooms_18_85 | 5.779970481941968 | 6 | `sol.aggregate_housing_demand / total_mass`; demand = Σ g_current·hR (renters) + Σ g_current·H_own (owners), normalizer N_target=1 (calibration.py:1090-1092; solver.py:5477-5494, 4065) | all indices 0–16 (18–82 starts, covering 18–86 vs data label 18–85) | everyone | mass-weighted mean | rooms per household |

Income-object summary: model `P.income[i,j] = (1−tau_pay)·w·profile[j]` for
workers, `pension` for retirees (`parameters.py:414-420`).
`annual_gross_income_at_state` (`solver.py:148-156`) = period income / 4,
grossed up by 1/(1−tau_pay) for workers only. So M10's denominator is annual
GROSS labor earnings (z- and age-varying), M6/M7's denominator is annual
pension (identical across all retirees since pension is flat and z-free after
`J_R`; z=1.0 hard-coded at `solver.py:4896` is innocuous because the whole
65–75 window is retired). Wealth = liquid `b` only for "nonhousing"/"liquid";
housing enters only via home equity `(1−psi)·p·H_own` in "total" variants
(`solver.py:4898`). Negative `b` (mortgages, b_min=−12) is included in
nonhousing wealth.

## 2. Objective recomputation (from stored moment vectors, no model call)

`diagnostic_loss` (`calibration.py:1244-1261`) is exactly
`sum_i w_i·(m_i − t_i)²` in raw units — no normalization, no scaling by
target, plus a +100 penalty iff `market_residual` is non-finite or > 5e-3, and
`inf` if any targeted moment is non-finite. `run_local_panel_case` stores this
as `rank_loss` (`local_panel.py:857`); candidate selection additionally
requires `strict_converged` (residual ≤ tol_eq=1e-4 AND solver strict flag,
`local_panel.py:862-867, 1151-1155`), so the stored rank_loss of any selected
record is penalty-free.

Script: `output/model/full_audit_20260711/scripts/recompute_loss.py` (run
2026-07-11, output at
`output/model/full_audit_20260711/target_audit_recompute_output.txt`).
Result: ALL CHECKS PASSED.

- current_bound_best: recomputed sum = 14.780020699972585, exactly equal
  (diff 0.0) to stored `loss` and `contribution_sum`; every stored
  `gap == model − target` and `loss_contribution == w·gap²` exactly;
  `diagnostic_loss(stored moments, residual 2.906e-05)` returns the identical
  float; +100 branch verified to fire at residual 6e-3.
- housing_relaxed_best: recomputed sum = 13.90346531862862 = stored `loss`
  (diff 0.0). Stored `contribution_sum` = 13.903465318628617 differs from
  stored `loss` by 1 ulp (summation-order noise only).
- Stored fit targets and weights match `get_target_set("candidate_replacement_
  post_audit_v1")` + the rooms extra target exactly, for all 15 rows.

## 3. Weight vector and loss shares

Weight vector (raw-unit squared gaps): own-rate block {old_age 160, own_rate
100, 2534 80, family gap 45}, room-tail 25, fertility 20/20, event-study 14,
{roomgap, young wealth} 12, {3plus rooms} 8, {renter rooms, agg rooms} 6,
old wealth {gap 2.0, median 0.8}.

Loss shares (from recomputation):

| moment | current_bound share | housing_relaxed share |
|---|---|---|
| old_age_own_rate | 36.5% | 45.5% |
| prime30_55_childless_owner_minus_renter_mean_rooms | 18.0% | 4.6% |
| young_childless_renter_liquid_wealth_2535 | 10.8% | 5.4% |
| prime30_55_childless_renter_mean_rooms | 9.2% | 2.4% |
| old_nonhousing_wealth_to_income_median_6575 | 7.8% | 19.7% |
| own_rate_2534 | 4.4% | 0.02% |
| own_rate | 4.1% | 10.5% |
| owner_share_rooms_ge6 | 3.6% | 2.9% |
| aggregate_mean_occupied_rooms_18_85 | 2.1% | 1.6% |
| own_family_gap | 1.2% | 0.6% |
| parent_3plus_minus_1to2 | 0.8% | 0.0% |
| old wealth gap 6575 | 0.8% | 0.8% |
| tfr | 0.6% | 4.6% |
| housing_increment_0to1 | 0.05% | 0.9% |
| childless_rate | 0.04% | 0.2% |

Scale-dependence assessment: the mix is ad hoc (hand-set ledger constants that
evolved across target-set versions; no variance/SE-based weighting is
documented anywhere in `calibration.py`). Under the active vector, a
1-percentage-point ownership error costs 0.016 loss (w=160) while a 0.1-room
error costs 0.12 (w=12) and a 0.1 wealth-ratio error costs 0.008–0.12
(w=0.8–12). Consequences visible in both candidates: `old_age_own_rate` alone
is 36–46% of the loss (an 18–20pp overshoot at w=160 dominates everything),
while the old wealth MEDIAN — missed by a factor of 2–6 — costs only 1.2–2.7
raw loss at w=0.8, so the optimizer rationally sacrifices old-age wealth
entirely. Both candidates share essentially the same dominant-miss structure;
the loss difference (14.78 vs 13.90) is small relative to the economic
difference in theta_n (0.77 vs 1.11) and the housing block, which is itself a
weak-identification symptom. Not recommending dropping targets (per rules);
but any interpretation of "fit improved by X" is weight-vector-specific.

## 4. Identification map (14 parameters × 15 moments)

Prior evidence: `docs/model/intergen_sensitivity_jacobian_audit_20260618.md`
(J=17 core point, PREDATES the combined spec: no H0 search, no Rouwenhorst
income, no bequest normalization, no aggregate-rooms moment — directions
indicative only).

| parameter | plausibly informative moments | status at current_bound_best |
|---|---|---|
| beta_annual | M6, M10 (wealth ratios), M3/M12/M13 lifecycle ownership | interior 0.9577; Jacobian: collinear with kappa_fert (cos 0.892) |
| alpha_cons | M8, M11, M15 (rooms), M3 | interior 0.673 |
| c_bar_0 | M10 (young saving), M1/M2 via resources | **AT upper bound 1.28** — not interior-identified |
| c_bar_n | M1, M2, M14 | interior 0.456 |
| h_bar_0 | M8, M15 | **AT lower bound 1.0** |
| h_bar_jump | M5, M14, M4 | interior 1.476 |
| h_bar_n | M14, M5 | interior 0.984 |
| psi_child | M1, M2 | interior 0.269; cos(psi_child,kappa_fert)=−0.984, cos(psi_child,theta_n)=0.975 at roomcost point |
| kappa_fert | M1 vs M2 split | interior 1.891 |
| tenure_choice_kappa | M3/M13 smoothing | **AT lower bound 0** (deterministic tenure; per user rule this bound is acceptable/expected) |
| chi | M3, M12, M9, M11 | **AT upper bound 1.15** (housing_relaxed, with bound moved to 1.6, takes 1.20 → the production bound binds) |
| theta0 | M6, M7 (old wealth) | interior 0.1318 — D11's theta0≈0 pathology is gone |
| theta_n | M7 ONLY (parent−childless old wealth gap, weight 2.0, <1% of loss) | interior 0.768 but 0.768 vs 1.109 across the two near-equal-loss candidates → **still weakly identified**; June 18 audit: theta0/theta_n columns "tiny" and old_age_own_rate loads on chi/alpha/beta, not bequests. D11 is at best partially resolved. |
| H0 | M15 (aggregate rooms) and price level → M3/M12/M9 | **9.99967 ≈ upper bound 10** (near_bound flag True in the report CSV) while M15 still undershoots (5.55 < 5.78) — the supply-scale parameter is bound-constrained, its identifying moment unreached from below |

Count check: 15 moments ≥ 14 parameters, so nominally just-identified plus
one. BUT at current_bound_best, 4 parameters sit at bounds (c_bar_0, h_bar_0,
chi, H0≈; tenure_choice_kappa=0 is an accepted corner) — these are
bound-constrained, not moment-identified, and the interior problem is
effectively ~9-10 parameters against 15 moments with at least one direction
(theta_n) still disciplined by a single 2%-of-loss moment.

Near-redundancies: (i) M3/M12/M13 are three points of one ownership lifecycle
profile; useful jointly but strongly co-moving through chi/beta.
(ii) M8 + M11 jointly pin owner mean rooms, so M9 adds only tail information
given the 5-rung grid; (iii) M15 is close to a mass-weighted combination of
the tenure-conditional room means and ownership rates over age — its marginal
information is precisely the supply/H0 scale, which is why H0-at-bound
matters.

## 5. Findings ledger (this task's scope)

1. MAJOR — hard-targeted grid-discrete median. M6 is an uninterpolated
   weighted median over ratios `bg/pension` (`utils.py:228-231` returns
   `values[pos]`, a raw b-grid node value scaled by a constant), so the moment
   moves in discrete steps of the b-grid (Nb=120 search) and is
   piecewise-constant in parameters. This violates the standing project rule
   "never hard-target median moments (grid-discrete)" and the June 22
   surrogate finding that discrete-median moments break smooth-search
   behavior. Weight is only 0.8, but at housing_relaxed it is 19.7% of the
   loss. Decisive test: perturb one parameter (e.g. beta) by ±0.5% at the
   candidate and plot M6 — it will be a step function.
2. MAJOR — H0 at its search bound. The one parameter added for the combined
   spec ends at 9.9997 of [1,10] while its identifying moment (M15) still
   misses from below; the extra parameter bought no interior optimum, and the
   14-parameter/15-moment count is not the effective identification story
   (4 parameters at bounds).
3. MAJOR (inherited, confirmed still present) — theta_n leans on a single
   low-weight moment (M7, 2% of loss); two near-equal-loss candidates differ
   by 45% in theta_n. D11 "theta_n unidentified" is not resolved by theta0
   being interior.
4. MINOR — systematic age-window drift: every windowed moment excludes the
   lower labeled edge and extends past the upper (25–34→ model 26–37;
   65–75→66–77; 30–55→30–57; 18–85→18–86). Consistent across moments, so
   comparative statements are fine, but level targets are matched against
   slightly older model windows than the data labels.
5. MINOR — "childless" means childless-in-household (includes stochastic
   empty-nesters) for room/young-wealth moments, but never-parent for the
   old-age wealth gap; PSID young "childless renters 25–35" are plausibly
   never-parents, so M10's model cell includes early empty-nest parents the
   data cell excludes (small mass at 26–34 given maturity hazard 1/4.5 per
   period, but nonzero).
6. MINOR — stale measurement ledger: `TARGET_MOMENT_OBJECTS` entries for M6
   and M7 (`calibration.py:267-278`) still say "divided by model period
   income … needs-fix", but the code divides by ANNUAL income
   (`solver.py:4896`); the ledger lags the implementation (the fix is in
   HEAD, not a dirty-tree change).
7. VERIFIED — objective arithmetic, target values, weights, gaps,
   contributions, and both stored losses reproduce exactly from calibration.py
   constants; the extra rooms target entered via `additional_targets/weights`
   into `run_local_polish`'s rank system (`local_panel.py:558-561`) with key
   equality enforced.
8. NOTE — M14's mapping (model top parity nn=2 ↦ data 3+ children; nn=1 ↦
   data 1–2) is coherent under the model's each-child-counts-double TFR
   convention (M1 = 2×mean parity) and is documented in the ledger
   (`calibration.py:243-248`), but it is a convention, not a measurement.
