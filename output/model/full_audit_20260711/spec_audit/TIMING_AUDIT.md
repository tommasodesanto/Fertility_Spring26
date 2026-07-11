# Timing / Distribution Audit (Task D) — July 9 repair re-audit from first principles

Auditor: timing/distribution specialist, 2026-07-11.
Object: DIRTY working tree at HEAD fe092ca (the tree the July 10 overnight snapshot was cut from).
Active evaluation chain: `intergen_housing_fertility.local_panel.run_local_panel_case` ->
`run_model_cp_dt` -> `solve_markov_income_equilibrium` -> `solve_bellman_full_markov_income`
+ `forward_distribution_markov_income` -> `compute_markov_statistics` -> `extract_moments`.

Decisive test script: `output/model/full_audit_20260711/scripts/timing_identity_test.py`
(run 2026-07-11, tiny model J=17, Nb=28, 2 income states, 2 owner rungs; all checks PASS —
output reproduced at the end of this file's summary section).

---

## 1. Bellman: flow utility at age j uses the NEWLY chosen tenure/housing — VERIFIED

`solve_bellman_core` (solver.py:2333) computes, per age j, conditional values `Vd[:, ten, ...]`
whose flow utility comes from THAT tenure's housing this period (renter: chosen `hR` with rent
`ri`, solver.py:2491-2546; owner rung ten: services `hsv = P.H_own[ten-1]` this period,
solver.py:2548-2597 "`ht_o = owner_service_premium * max(hsv - owner_h_bar_scale*SD.hb_flat, ...)`").
The tenure choice at j (solver.py:2686-2752) maximizes over `Vd` evaluated at transaction wealth
(buy: `bab = b_grid - hc` with infeasibility `inf_m = (b_grid < dpn) | (bab < bmn)`,
solver.py:2713-2731; sell: `ba = np.maximum(b_grid + sp, 0.0)`, solver.py:2705). So the period-j
choice delivers period-j flow utility. Order of nesting at age j: savings/housing inside tenure
inside location logit inside fertility logit (solver.py:2780-2791) — identical to the KFE order
(fertility split -> location -> tenure -> savings, solver.py:3762-3966).

Fertility timing: at fertile ages the chosen family size applies at j itself —
`Vfa[:, :, :, nn] = VI[:, :, :, nn, 1]` (solver.py:2784), i.e. the birth state (n, cs=1) with
its child costs `c_bar`, housing floor `h_bar`, and `psi_child` binds in the SAME period, and the
tenure choice inside that branch uses the child-state down payment `dp_choice[id_, tn, nn, cs]`
(solver.py:2717). Matches KFE: the birth split mutates `g` at age j in place
(solver.py:3776-3777 "`g[:, :, :, j, zz, 0, 0] = 0.0`" / "`g[:, :, :, j, zz, :, :] += gpf`").

## 2. Post-decision distribution: same object for moments AND clearing — VERIFIED

- Flag read at solver.py:4014: `use_postdecision_current = bool(getattr(P, "use_postdecision_current_distribution", True))`.
- `g_current = realize_current_cross_section(g, ...)` (solver.py:4015-4029): applies location
  (lmm map, home-equity cash-out, solver.py:2967-2985), tenure transactions (tmx map built by
  `build_forward_tenure_transition_maps`, solver.py:2875-2941, mirroring the Bellman accounting
  incl. buy cost, sale proceeds `np.maximum(bg + sale_proceeds, 0.0)`, borrowing floor
  `-financed_share*purchase_cost`, and the birth-grant/waiver precedence at solver.py:2916-2920),
  at EVERY age including j=J-1, WITHOUT aging/saving.
- Moments: `compute_markov_statistics(g_current, ..., asset_g=asset_current)` (solver.py:4049-4058).
- Market clearing: the equilibrium loop's demand comes from the SAME construction — fast path
  `compute_markov_eq_stats(g_current, ...)` (solver.py:4047; called with fast_stats=True inside
  `solve_markov_income_equilibrium`, solver.py:1019-1025 uses `sol_it.housing_demand`); the scalar
  Brent refine uses the same (solver.py:1097-1104); the final full solve at best_p (solver.py:1057)
  recomputes both moments and `best_max_abs_rel_excess` on `g_current`
  (pack_solution_markov_income, solver.py:5474-5501, demand recomputed from the returned `g`
  which IS `g_current`). Decisive test T5: recomputed demand from returned g == sol.housing_demand
  to 1.8e-15. The new target `aggregate_mean_occupied_rooms_18_85 = aggregate_housing_demand /
  total_mass` (calibration.py:1090-1092) is the same post-decision object divided by
  `stats.total_mass = sum(g_current)` (solver.py:4065).
- Renter demand integrand `g_current[:,0,...] * hR_pol[:,0,...]` is coordinate-consistent: the tmx
  sell branch places just-sold renters at `max(b+sale, 0)` (solver.py:2927), exactly the wealth at
  which the Bellman evaluated the renter branch (solver.py:2705), so `hR_pol` read at the realized
  coordinate is (up to linear interpolation) the housing actually chosen.

Mass conservation: in-solver assert is AGGREGATE only (solver.py:4042-4045, atol=1e-10).
Decisive test T1 checked it age-by-age and by (age, n, cs): max gap 8.7e-17. Conserved.

Wealth moments deliberately use a THIRD object `asset_current =
assign_current_cross_section_to_beginning_assets(g, ...)` (solver.py:4030-4041; function at
3105-3177): beginning-of-period b conditioned on the current period's location/tenure choice
(decision evaluated at transaction wealth, mass recorded at original b) — declared at
solver.py:4094 `stats.wealth_moment_timing = "beginning_of_period_b_conditioned_on_current_tenure"`.
This is a survey-timing convention, internally consistent; used by
`add_annual_gross_liquid_wealth_moments`, `liquid_wealth_4555`, old-age wealth ratios,
`owner_neg_liquid_share_*` (solver.py:4727, 4771, 4901, 4987-4993).

## 3. Births — VERIFIED (with one diagnostic-only caveat, finding T-2)

- n increments at j (in-place split, solver.py:3763-3777); child costs and the child housing floor
  bind at j through the (nn, cs=1) column of `SD.cb_flat/hb_flat` in the same-period Bellman flow.
- Birth-period purchases enter the horizon-0 statistic: `event_horizon = 0`
  (production_profile.py:22, solver.py:3728); `advance_cohort_horizon_markov_income` with
  horizon 0 is the identity (solver.py:4146-4150, loop `range(1, 0+1)` empty); the observed birth
  cohort is `realize_current_choices_markov_income(birth_cohort, j+0, ...)` under the flag
  (solver.py:3806-3821) — so a renter who gives birth AND buys at j is measured as an owner at
  H_own in `post_h`. The moment is a DiD against a same-period childless counterfactual cohort
  (solver.py:3829-3857): `housing_increment_0to1_eventstudy_t3 = (birth_es3_post_sum -
  birth_es3_control_post_sum) / birth_es3_mass` (solver.py:4070-4072), mapped to the calibration
  moment at calibration.py:1134.
- Decisive test T4: in the tiny model the newborn-state owner share at the birth age is strictly
  larger under the post-decision distribution than under inherited tenure (0.0059 vs 0.0002 at
  j=6), i.e. birth-triggered same-period purchases are counted.

## 4. Age-specific ownership moments — windows checked; July 2 finding STILL PRESENT

All tenure/room moments in `compute_statistics` integrate `g` = `g_current` (post-choice tenure):
`own_rate_2534` (solver.py:4832-4834), `own_rate_3055`, `own_rate_3544`,
`old_age_own_rate_6575` (solver.py:4868-4870), newparent/nonparent gaps (4846-4866), all
prime30_55 room moments (5035-5101), renter cap shares (5119-5177), `markov_renter_room_moments`
(4331+, on the full income-resolved g_current). Decisive test T2: reported own-by-age equals the
manual application of the age-j tenure policy to the age-j beginning mass to 1.5e-15 at every age,
and the legacy (flag False) series is EXACTLY the one-period lag — which is precisely the
pre-repair defect, now demonstrably gone.

Windows (age_to_index, solver.py:4614-4616, round((age-18)/4) clipped):
- `own_rate_2534`: j in [2,4] -> model ages 26,30,34 (4-year bins [26,30),[30,34),[34,38)).
- `old_age_own_rate_6575`: j in [12,14] -> model ages 66,70,74 (bins [66,70),[70,74),[74,78)).
Under the "model index = 4-year bin starting at 18+4j" reading, the 65-75 window still DROPS age
65 (lives in bin j=11, [62,66)) and still INCLUDES 75-77 (bin j=14 spans 74-78) — the July 2
finding is NOT fixed. Under the "model index = point age" reading, {66,70,74} is the correct
nearest-point set for [65,75] and {26,30,34} for [25,34]. Which reading matches the ACS/HRS
measurement is a data-side question the lead must adjudicate; the code has not changed either way.
These two moments carry weights 80 (`own_rate_2534`) and 160 (`old_age_own_rate`)
(calibration.py:489-492 CANDIDATE_REPLACEMENT_ROOMGAP_14MOMENT_V1_WEIGHTS inherited by
POST_AUDIT_V1), and both feed `old_minus_young_owner_rate_6575_2534` (calibration.py:1103).

## 5. Market clearing — VERIFIED (see section 2; T5 decisive)

Demand aggregates the identical `g_current` the moments use, in every branch of the equilibrium
loop (damped iteration, bracketing, Brent refine, final full solve). Supply side
`P.H0 * (user_cost / r_bar) ** xi_supply` (solver.py:5522-5523, 1025) is out of scope here
(task C).

## 6. Remaining old-timing paths — enumerated

The pre-repair code is fully retained behind the flag: `g_current = ... if use_postdecision_current
else g` (solver.py:4015-4029 markov; 3602-3629 non-markov mirror). BUT the default is TRUE in
`parameters.py:101` ("`P.use_postdecision_current_distribution = True`"; also
`P.propagate_birth_entry_grant = True` at :100), so every tool that does not pass the switch runs
REPAIRED timing. The task brief's hypothesis (default False -> tools silently on old timing) is
FALSE on this tree. Explicit flag-False consumers, all deliberate comparison arms:
`production_profile.FROZEN_SOURCE_REPRO_SWITCHES` (production_profile.py:75-80) used by
`tools/run_intergen_code_repair_comparison.py`, `tools/compare_intergen_combined_specification.py`,
`tools/build_intergen_mechanics_packet.py`, `tools/compare_intergen_income_processes.py`,
`tools/compare_intergen_bequest_normalization.py` via `comparison_arm_switches`. Residual risk:
any tool that REPLAYS a pre-repair record's saved overrides verbatim (e.g.
`run_policy_counterfactuals_from_record.py`) would faithfully reproduce old timing if the stored
record embeds the switch — that is reproduction semantics, not a defect, but the lead should keep
pre-July-9 records out of headline policy numbers.

The July 10 overnight itself: `tmp/overnight_combined_20260710/run_overnight_combined.py:135-166`
calls `lp.run_local_polish(..., profile_name=PRODUCTION_PROFILE_NAME, ...)`;
`run_local_polish` merges `production_profile_overrides()` (local_panel.py:565-579) which pins
`"use_postdecision_current_distribution": True`, `"propagate_birth_entry_grant": True`,
`"housing_event_horizon": 0` (production_profile.py:159-166); final merge order in
`run_local_panel_case` is `{base_overrides, **extra_overrides, **income, **theta}`
(local_panel.py:848-853), and no income/theta key collides with the timing switches. The overnight
ran repaired timing.

## 7. tests/test_code_repair_20260709.py — what it proves and what it does not

Proves (unit level, synthetic fixtures, I=1, 2 tenures):
- Transaction-map math: zero-grant map identical to historical map; grant applied only renter->buy;
  waiver precedence over grant matches the Bellman branch (lines 349-426).
- `realize_current_choices` records a same-period birth+purchase in the same current period, and
  `compute_eq_stats` on the realized vs inherited cross-section yields 6.0 vs 2.0 demand
  (lines 429-462) — the flag's intent, on a hand-built 1-point cohort.
- Terminal age gets current choices applied (464-480); beginning-asset conditioning keeps b while
  switching tenure label (482-502); one-period advance applies transaction, saving, child stage,
  and income transition exactly once with mass 1.0 conserved (504-572).
- Profile pinning: production overrides/theta/switch sets exact; DE panel passes the two switches
  (278-346).

Does NOT prove:
- Nothing end-to-end: no test solves the model and checks that the RETURNED distribution, the
  moments, and the market-clearing demand are all the post-decision object (my T2/T5 close this).
- Mass conservation only on 1-point cohorts, never age-by-age on a solved model (T1 closes this).
- No test of the horizon-0 DiD statistic construction (birth vs control cohort), none of the
  age-window indices (T6), none of the Bellman-side flow-utility timing, and no cross-check that
  the non-markov and markov KFEs define the same statistics (they do not — finding T-3).

## Findings

### T-1 (MINOR, pre-existing, high-weight targets): 2534/6575 age windows unchanged since July 2
solver.py:4706-4712 with age_to_index rounding gives own_rate_2534 = model ages {26,30,34} and
old_age_own_rate_6575 = {66,70,74}. Under the bin reading the old window drops 65 and includes
75-77 — the July 2 finding is still present, on moments weighted 80 and 160. Decisive test:
timing_identity_test.py T6 output. Consequence: level bias in two of the three highest-weighted
ownership moments if the data windows are bin-aligned; calibration compensates through theta.

### T-2 (MINOR, diagnostic-only): fert_by_age / mean_age_first_birth use POST-split childless mass
The KFE mutates g at fertile ages before stats (solver.py:3776-3777), so
`stats.fert_by_age` (solver.py:4691-4704) and `mean_age_first_birth` (solver.py:4796-4809,
"`pb = 1 - fp[nz, ten, i, j, 0]`" applied to `gs = g[:, ten, i, j, 0, 0]`) weight age j by
(pre-split childless mass)*fp0*(1-fp0) instead of the true birth flow (pre-split mass)*(1-fp0).
Not in the active 15-moment target set (mean_age_first_birth was dropped in all *_no_timing and
replacement sets, calibration.py:397-399), but it IS in extract_moments (calibration.py:1133) and
in the reported diagnostic `full_old_nonlocation_loss` (local_panel.py:858, weight 12). Pre-dates
the July 9 repair. Decisive test: place unit childless mass at one fertile age, compare
fert_by_age against `E[n' choice]` computed from pre-split mass.

### T-3 (MINOR, path inconsistency): housing_increment_0to1_eventstudy_t3 has two definitions
Markov KFE (ACTIVE path): `(birth_es3_post_sum - birth_es3_control_post_sum)/birth_es3_mass`
(solver.py:4070-4072), a same-period DiD against a childless counterfactual cohort. Non-markov KFE:
`post_sum/mass - pre_sum/mass` (solver.py:3659-3661) — post minus PRE-birth housing, no control
cohort exists in that function at all (no `birth_es3_control_post_sum` accumulator before
solver.py:3683). The comment at solver.py:3726-3727 ("at horizon 0 the control is the pre-birth
state") is not exact under repaired timing: the control cohort gets CURRENT-period realized choices
(solver.py:3836-3851), the pre-birth mean uses inherited holdings. Any tool running the
income-type/non-markov path reports a different object under the same moment name. Active
calibration unaffected (markov path).

### T-4 (SUSPECTED, reproduction semantics): replaying stored pre-repair records replays old timing
Records that embed `use_postdecision_current_distribution: False` in saved overrides (frozen-source
repro arm, any pre-July-9 result JSON) will faithfully re-run legacy timing through tools that
replay stored overrides. Not a code defect; a bookkeeping hazard for headline numbers.

## Verified items (checked and correct)

1. Bellman flow utility at j uses newly chosen tenure/housing and same-period child state
   (solver.py:2491-2546, 2548-2597, 2686-2752, 2780-2791).
2. KFE order (fert -> location -> tenure -> savings -> aging) matches Bellman nesting
   (solver.py:3762-4003 vs 2468-2793).
3. `g_current` (post-decision) is the single object behind moments, clearing, packed `sol.g`,
   `total_mass`, and the new rooms target (solver.py:4014-4058, 1019-1025, 1097-1104, 5474-5501).
   Decisive T5: 1.8e-15.
4. Mass conservation age-by-age and by (age,n,cs): 8.7e-17 max gap (decisive T1; in-solver assert
   at 4042-4045 is aggregate-only, now independently verified at cell level).
5. Contemporaneous-tenure identity: reported own-by-age == manual choice-application, max dev
   1.5e-15; legacy arm is exactly the one-period lag (decisive T2/T3) — the pre-repair defect is
   gone and was exactly what the July 9 report described.
6. Birth-triggered same-period purchases counted in the horizon-0 statistic and in the stationary
   moments (decisive T4; solver.py:3806-3824, 4070-4072).
7. Transaction map == Bellman accounting incl. grant/waiver precedence, borrowing floor, sale
   nonnegativity (solver.py:2875-2941 vs 2698-2743; unit tests 349-426).
8. Flag defaults TRUE in parameters.py:100-101; overnight drivers inherited repaired switches via
   production_profile_overrides (production_profile.py:159-166, local_panel.py:565-579, 848-853).
9. `advance_cohort_horizon*` with horizon 0 is the identity (solver.py:4146-4150).
10. Wealth moments' beginning-b-conditional-on-current-choice convention implemented as declared
    (solver.py:3105-3177, 4094).

## Open questions for the lead

- Which age convention (bin [18+4j, 18+4j+4) vs point age 18+4j) do the ACS 25-34 and 65-75 data
  moments assume? That decides whether T-1 is a real misalignment or the correct nearest-point map.
- Cluster-side confirmation that the torch snapshot's parameters.py/production_profile.py match
  this tree (SSH blocked; local tmp/intergen_combined_specification_20260710 copy does match on
  the flag lines).
