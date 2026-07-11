# Diagnostics packet — July 10 combined-spec candidates (audit 2026-07-11)

Scope: numerical diagnostics for `current_bound_best` (GE re-solve, Nb=120, loss
reproduced exactly at 14.780020699972585, p_eq=0.6897624769313723) and a light
PE pass on the relaxed/housing candidate at its stored p_eq=0.7564643263592232.
All artifacts produced by the isolated drivers in
`output/model/full_audit_20260711/scripts/` (`audit_diag_packet.py`,
`audit_diag_probe.py`, `audit_diag_extra.py`, `audit_diag_locate.py`). No
production files touched. The plots are the package's STANDARD set from
`intergen_housing_fertility/diagnostics.py::write_diagnostics` (no new graphs).

## Figures — current_bound (`diagnostics/current_bound/`)

- ownership_by_age.png
- fertility_by_age.png
- housing_market.png
- market_clearing_by_market.png
- market_clearing_residuals.png
- tenure_services.png
- owner_rungs.png
- housing_prices.png
- income_state_outcomes.png
- ownership_by_age_income_state.png
- liquid_wealth_by_age_income_state.png
- housing_by_age_income_state.png
- fertility_policy_by_age_income_state.png
- wealth_dist_childless_renter_age30.png
- wealth_dist_childless_renter_age42.png
- policy_childless_renter_age30.png
- policy_childless_renter_age42.png
- summary.json (standard scalar summary)
- numeric_checks.json, probe_checks.json, extra_checks.json,
  locate_invalid_j0.json (audit numeric packs)

## Figures — relaxed light pass (`diagnostics/relaxed_light/`)

Same standard plot set (PE at stored p_eq, so market plots reflect the stored
price, residual 6.6e-05), plus numeric_checks_light.json.

## Numeric findings (current_bound, Nb=120 GE)

1. Value monotonicity: 168 pairwise decreases in wealth among "valid"
   (V>-1e9) nodes; concentrated at age indices 11-14, owner-parent states.
   Population mass at violating nodes 2.5e-05 (share 2.5e-05). The apparent
   violations are driven by sentinel contamination: valid-node V ranges down to
   -9.997e8 (just above the -1e9 validity cutoff); 7,612 valid nodes have
   V < -1e6 carrying 0.57% of population mass. The Bellman propagates the
   -1e10 infeasibility sentinel through expectations instead of restricting
   the choice set (kernels.py NEG_INF=-1e10; solver logsumexp blocks).

2. Populated infeasible mass (MAJOR/FATAL, see audit findings): 1.72% of total
   population mass sits on nodes with V=-1e10, almost all renters (0.01712 of
   0.01719), decaying with age (j=0: 0.01199 = 20.4% of the entering cohort;
   j=1: 0.0032; j=2: 0.00095). Cause: PSID entry-wealth ratio distribution
   (calibration.py external_entry_wealth_overrides) injects entrants over 27
   wealth nodes in [-2.907, +3.512]; at the candidate's fixed floor
   c_bar_0=1.28, low-z / negative-wealth entrant nodes are infeasible even as
   renters. Households on those nodes still transit through the KFE with
   default policies (c in [1.32, 3.08], b' in [-1.10, 2.74]).

3. Uniform-random fertility at infeasible nodes (FATAL-candidate): at fully
   infeasible childless nodes fert_probs = (1/3, 1/3, 1/3) exactly (verified:
   0.333333485667725 per option; logsumexp over three -1e10 values). The KFE
   applies these, so 2/3 of the 20.4% infeasible entrants receive a "birth"
   at age 18-22 (1/3 receive two children). Artifact births at j=0 alone:
   0.003995*1 + 0.003995*2 = 0.011986 births vs total cohort births
   ~0.0583 (TFR/2 * 1/17), i.e. ~21% of all model births are mechanical.
   TFR net of the j=0 artifact would be roughly 1.57 vs the reported 1.983
   (target 1.918). Childlessness, fertility-income gradient, age at first
   birth are contaminated in the same direction.

4. p_birth = 0 cells: 8.9% of childless fertile-window cells, carrying 40.3%
   of childless fertile-window mass (0.0980 of 0.2434), all on FEASIBLE nodes
   (invalid-node contribution zero). Mechanism: the birth option is
   infeasible (child floors) so the Gumbel/logit assigns exactly zero. Mean
   wealth of the p_birth=0 childless mass is 0.0002 vs 0.1232 for all
   childless: the extensive fertility margin is dominated by hard feasibility,
   not smoothed choice.

5. Consumption/saving: no negative consumption on populated states (min
   c=1.32); no bunching at the wealth-grid edges (zero mass at b_min=-12 and
   b_max=30; populated b' spans [-5.52, 8.20]); 47.2% of population mass has
   negative liquid wealth.

6. Tenure/housing: ownership rises from 9.3% at 18 to 94.8% at 66+, largest
   jump 38->42 (+19.8pp); aggregate own rate 67.8%. Owner rung shares
   (H=2,4,6,8,10): 1.0%, 26.7%, 38.2%, 28.9%, 5.2% — interior, no top-rung
   pile-up. Renter-cap (hR_max=6) bunching: 13.2% of renter mass.

7. Market clearing: demand 5.5511 vs supply 5.5509, max rel. residual
   2.906e-05 < tol_eq 1e-4 (strict).

8. Mass accounting: total mass 1.0 (exactly), per-age mass exactly uniform at
   1/17 = 0.0588235 (max abs deviation 5.3e-16): no leakage across ages. The
   fert_probs rows that sum to 0 (58.8% of cells) all lie at ages j>=7 where
   the KFE never applies fertility (operative window is j=0..6 via the
   j+1 in [A_f_start, A_f_end] condition), so they do not destroy mass; note
   diagnostics.py plots the fertility window as j=0..7 (off-by-one vs the
   KFE), so the age-46 point in fertility_policy_by_age_income_state.png is
   an artifact of never-computed probabilities.

## Relaxed candidate (light pass, PE at stored p_eq)

Same pathologies at similar magnitude: populated-infeasible mass 1.74%;
p_birth=0 childless mass share 43.2% (operative window); 153 monotonicity
decreases; renter-cap bunching 10.7%; mass exactly conserved; residual
6.60e-05; own rate 71.9%; min populated consumption 1.32.
