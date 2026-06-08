# One-Market Calibration Structure Memo

Date: 2026-06-08

## Purpose

This memo fixes the calibration-order problem in the current one-market
intergenerational housing/fertility scaffold. The diagnostic random screen is
not a formal calibration. Before another long run, the model needs a parameter
inventory, a target inventory, and an explicit parameter-to-moment map.

Startup note: `CALIBRATION_STATUS.md` still describes the old workhorse as the
live model. For the current pivot, the controlling implementation records are
`code/model/intergen_housing_fertility/IMPLEMENTATION_STATUS.md` and
`docs/model/intergen_housing_fertility_audit_20260608.md`.

## Current Model Object

The current implementation is a runnable scaffold, not a production calibrated
model. The economically active one-market ingredients are:

- One aggregate housing-services market, with `I=1`.
- Owner housing is a discrete rung choice, `H_own`.
- Renter housing is continuous up to `hR_max`.
- The code stores the period interest rate as `P.q`; owner user cost is
  `(P.q + P.delta + P.tau_H) * P_asset`.
- The one-market price clears
  \[
  H^D(P) = H_0 \left(\frac{(r+\delta+\tau^p)P}{\bar q}\right)^\eta .
  \]
- New owner purchases face a down-payment constraint
  \[
  b \ge (1-\phi)Ph,
  \]
  with `phi` interpreted as the financed share.
- New owner purchases also face a simplified payment-to-income screen
  \[
  (r\phi+\tau^p)Ph \le \psi^{PTI} y_{a,z}.
  \]
- Fertility is one-shot completed parity for childless fertile households, not
  a sequential parity hazard.
- The dependent-child state is geometric with expected duration 18 years.
- The default period is 4 years, with fertile ages represented by period-start
  ages 26, 30, 34, 38, and 42.

## Parameter Inventory

The current random-search driver varies 19 knobs:

| Parameter | Current status | Calibration status |
|---|---|---|
| `alpha_cons` | Structural preference | Should be externally fixed or tied to a housing-expenditure/share target before calibration |
| `phi` | Finance constraint, financed share | Prefer externally fixed from LTV/down-payment evidence unless a credit-access moment is added |
| `b_entry_fixed` | Entrant liquid wealth | Candidate internal parameter |
| `pti_limit` | Finance constraint | Prefer externally fixed from underwriting/payment evidence unless a PTI/denial/share-constrained moment is added |
| `c_bar_n` | Child goods cost slope | Candidate internal parameter |
| `chi` | Owner housing-service premium | Candidate internal parameter |
| `kappa_fert` | Fertility logit scale | Candidate internal parameter, but not enough by itself to identify timing |
| `psi` | Sale/transaction haircut in retained owner equity | Meaningful but should be externally fixed unless targeting turnover/down-sizing |
| `psi_child` | Child utility shifter | Candidate internal parameter |
| `theta0` | Bequest intensity | Candidate internal parameter |
| `theta_n` | Bequest child shifter | Candidate internal parameter |
| `theta1` | Bequest shift | Externally fixed normalization, not a target-calibrated parameter in this scaffold |
| `h_bar_jump` | First-child housing need jump | Candidate internal parameter |
| `h_bar_n` | Additional-child housing need slope | Candidate internal parameter |
| `hR_max` | Renter menu cap | Menu/market normalization unless disciplined by renter size distribution |
| `owner_h_bar_scale` | Owner effective subsistence scale | Menu/utility normalization unless separately disciplined |
| `owner_size_cost` | Extra large-owner-unit cost | Menu/supply device unless disciplined by owner size distribution |
| `owner_size_cost_ref` | Reference size for owner size cost | Menu/supply device, not structural without a size-distribution target |
| `H0` | Aggregate housing supply shifter | Market normalization, not a household preference parameter |

This is already the core failure: 19 varied knobs are being scored against 10
hard moments, and several varied knobs are normalizations or menu devices.

Additional economically meaningful parameters currently fixed:

| Parameter | Role | Recommended treatment |
|---|---|---|
| `beta` | Saving, wealth accumulation, tenure timing | Either externally fixed or internally calibrated to young liquid wealth and ownership age profile |
| `sigma` | CRRA curvature | Externally fixed |
| `c_bar_0` | Baseline childless consumption need | Fixed normalization unless a baseline consumption/poverty object is added |
| `h_bar_0` | Baseline childless housing need | Fixed or calibrated only if childless housing-size targets are hard moments |
| `q`, `delta`, `tau_H` | Interest, depreciation, property tax | Externally fixed for benchmark; policy experiments can move `tau_H` |
| `income_age_breaks`, `income_age_values` | Lifecycle income profile | Externally estimated |
| `z_grid`, `z_weights`, `income_shock_persistence`, `Pi_z` | Markov income process | Externally estimated before production |
| `tau_pay`, `pension_mode`, `pension` | Retirement income system | Externally fixed/accounting closure |
| `tenure_choice_kappa` | Type-I-EV smoothing over tenure/rungs | Numerical smoothing; keep fixed and report sensitivity |

Lifecycle and state-space parameters are model-design choices, not calibration
parameters: `period_years`, `J`, `age_start`, `J_R`, `A_f_start`, `A_f_end`,
`A_m`, `stage_durations`, `n_parity`, `use_stochastic_aging`,
`n_child_stages`, and `n_child_states`.

Market/menu parameters are normalizations unless matched to distributional
housing targets: `H_own`, `n_house`, `h_own_min`, `h_own_max`, `hR_max`,
`H0`, `r_bar`, `eta_supply`, `xi_supply`, `w_hat`, `N_0`, `N_target`,
`entry_by_loc`, and `entry_shares`.

Inherited one-market-degenerate parameters should not be calibrated in the new
one-market setup: `kappa_loc`, `eps_loc`, `E_loc`, `mu_stay`, `mu_move`,
`mu_move_parent`, `mu_age_decay`, `location_choice_form`, `kappa_entry`,
`outside_value`, `outside_entry_flow`, `outside_entry_shares`,
`local_birth_entry_weight`, `renewal_retention`,
`renewal_calibrate_outside_flow`, and related population-closure fields. They
remain in code because the scaffold starts from the old workhorse architecture.

Disabled policy toggles are not production parameters unless the corresponding
policy experiment is implemented and documented: `parent_dp_waiver`,
`parent_dp_waiver_phi`, `birth_dp_grant`, `birth_entry_grant`,
`birth_entry_grant_amount`, and their location/rung filters.

Numerical parameters should never enter the economic calibration vector:
`Nb`, `b_min`, `b_max`, `b_grid_power`, `max_iter_eq`, `tol_eq`,
`lambda_eq`, scalar price-refinement controls, adaptive damping controls,
price bounds, entry floors, interpolation flags, and kernel flags.

## Current Hard Moment List

The active one-market diagnostic target set is `old_nonlocation` in
`code/model/intergen_housing_fertility/calibration.py`.

| Moment | Target | Current model extractor |
|---|---:|---|
| `tfr` | 1.70 | `2 * mean_completed_fertility` |
| `childless_rate` | 0.15 | completed parity mass at zero |
| `mean_age_first_birth` | 26.0 | probability-weighted first-birth age over period-start ages |
| `own_rate` | 0.57547241 | `own_rate_3055`, not aggregate ownership |
| `own_family_gap` | 0.16766167 | new-parent minus non-parent ownership, ages 30-55 |
| `housing_increment_0to1` | 0.66443467 | birth-cohort housing response at horizon 3 periods |
| `housing_increment_1to2` | 0.56581378 | additional-child proxy at horizon 3 periods |
| `young_liquid_wealth_to_income` | 0.60 | childless renters, roughly ages 25-35 |
| `old_age_own_rate` | 0.76426097 | ownership ages 65-75 |
| `old_age_parent_childless_gap` | 0.070 | parent minus childless ownership ages 65-75 |

Important target inconsistency: `latex/main_note.tex` contains older benchmark
values for some ownership and old-age targets. The current one-market code uses
the values above. Production calibration should first decide which empirical
target file is canonical for the new project.

## Old Workhorse Moments That Disappear

The one-market model cannot identify moments that require multiple markets or
market-to-market movement. These old workhorse targets should be dropped or
replaced:

- Center population share / inversion share.
- Center-over-periphery rent ratio / unit-rent inversion ratio.
- Fertility gradient across markets.
- Ownership gradient across markets.
- Center share of non-parents and new parents.
- Market-switching or migration rate.
- Any counterfactual claim that depends on differential local supply or sorting.

With `I=1`, these objects are either undefined, mechanical, or zero by
construction. Keeping them in the objective would be false identification.

## Replacement Non-Location Moments

Available in the current code and usable after empirical validation:

- Ownership by age: `own_rate_2534`, `own_rate_3544`, `own_rate_3055`,
  `old_age_own_rate_6575`.
- Ownership by family status: `own_gap_newparent_nonparent_3055`,
  `old_age_parent_childless_gap_6575`.
- Housing-size distribution diagnostics: childless renter median rooms,
  childless owner median rooms, renter cap shares, owner rung shares, and owner
  room-bin shares for ages 25-45.
- Family housing responses: `housing_increment_0to1_eventstudy_t3`,
  `housing_increment_1to2_proxy_t3`,
  `housing_increment_0to1_onechild_eventstudy_t3`, and
  `housing_increment_0to2plus_eventstudy_t3`.
- Financial-state diagnostics: `young_liquid_wealth_to_income`,
  young owner negative-liquid share, and liquid wealth by age/income state.
- Income-state gradients: ownership, fertility, and housing demand by Markov
  income state.

Moments that should be built before production:

- A model-equivalent first-birth timing object: 4-year binned first-birth
  hazards or first-birth age from the same household-unit definition used in
  the model.
- A direct credit-constraint target if `phi` or `pti_limit` is internally
  calibrated: down-payment shortfall, LTV, denial/PTI bind rate, or first
  purchase by age conditional on income/wealth.
- Old-retention and underuse moments if an old-retention wedge is added:
  old-owner large-home share, downsizing/turnover rate, or old-vs-young family
  occupancy of family-sized units.
- Housing affordability ratios by age/tenure if the property-tax capitalization
  mechanism becomes quantitative rather than illustrative.

## Parameter Count

The current diagnostic screen fails the counting rule:

\[
19 \text{ varied knobs} > 10 \text{ hard moments}.
\]

After removing menu normalizations and externally fixed objects, a disciplined
production calibration can be at or below 10 internal parameters. If the
first-birth timing target remains hard, one missing timing parameter must be
added; otherwise the model is trying to fit timing through housing and finance
parameters that are supposed to identify other mechanisms.

## Identification Map

A disciplined one-market map should look like this:

| Parameter | Primary moment | Secondary moments/checks |
|---|---|---|
| `b_entry_fixed` | `young_liquid_wealth_to_income` | young ownership and owner-entry thresholds |
| `psi_child` | `tfr` | childlessness |
| `c_bar_n` | `childless_rate` and completed fertility mix | TFR and young consumption feasibility |
| `kappa_fert` | parity dispersion / childlessness | timing only weakly, not as the sole timing margin |
| new timing shifter, not implemented | first-birth timing object | age profile of first births |
| `h_bar_jump` | `housing_increment_0to1` | family ownership gap |
| `h_bar_n` | `housing_increment_1to2` | intensive completed fertility |
| `chi` | `own_rate_3055` | family ownership gap and owner/renter housing services |
| `theta0` | `old_age_own_rate` | lifecycle ownership slope |
| `theta_n` | `old_age_parent_childless_gap` | parent-childless wealth and old housing retention |

Finance parameters should be handled deliberately:

- If `phi` and `pti_limit` are externally fixed, the current 10 moments can
  discipline the 10-parameter vector above once the timing shifter exists.
- If either `phi` or `pti_limit` is internally calibrated, add at least one
  direct credit moment per added parameter and remove or externally fix a
  different parameter. Otherwise the financial constraint is not separately
  identified from entry wealth, discounting, or the ownership premium.

## First-Birth-Age Target

The target `mean_age_first_birth = 26.0` is not yet verified as comparable to
the model object.

The model statistic is
\[
\frac{\sum_{a \in \{26,30,34,38,42\}} a \cdot m_a
      \Pr(\text{birth at }a \mid \text{childless at }a)}
     {\sum_{a \in \{26,30,34,38,42\}} m_a
      \Pr(\text{birth at }a \mid \text{childless at }a)} ,
\]
using period-start ages for a unitary household agent.

The local source label in `latex/main_note.tex` is NVSS. That is not yet a
documented local construction showing the same unit, period convention, and
conditioning set. NVSS is plausibly a female period statistic over observed
first births; the model is a stationary household-agent hazard-weighted object.
Those are not automatically the same target.

There is also a mechanical support issue. With the current fertility window,
26 is the earliest model birth age. A target of exactly 26 puts the target at
the lower bound of the model's timing support. If period midpoints were used,
the earliest implied age would be 28, making the gap worse.

Production options:

1. Verify and reconstruct the target as a 4-year binned household first-birth
   hazard from PSID/ACS/NSFG/NVSS with the same unit convention.
2. Replace the scalar mean with a small age-profile target, for example
   first-birth shares or hazards in bins 22-25, 26-29, 30-33, 34-37, and 38-42,
   if the model fertility window is expanded to include the first bin.
3. If `26.0` is truly the desired hard target, add an explicit age-specific
   fertility shifter or child-cost schedule. Do not force the target through
   housing preferences, down-payment constraints, or bequests.

## Smallest Disciplined Production Vector

Under the current 10 hard moments, the smallest coherent production vector is:

\[
\Theta =
\{b_0,\psi_{\mathrm{child}},\bar c_n,\kappa_n,
\lambda_a,\bar h_{\mathrm{jump}},\bar h_n,\chi,\theta_0,\theta_n\},
\]

where \(\lambda_a\) is a new age-specific timing shifter or compact schedule
that is not implemented yet.

Externally fix for this first production calibration:

- `beta`, `alpha_cons`, `sigma`;
- `q`, `delta`, `tau_H`;
- `phi`, `pti_limit`;
- income process and lifecycle income profile;
- `H_own`, `hR_max`, `H0`, `r_bar`, `eta_supply`;
- `theta1`;
- `tenure_choice_kappa`;
- all numerical solver controls.

If the project wants to calibrate `phi` or `pti_limit` internally, the target
system must add direct credit moments first. In that case the production vector
should not grow silently; it should replace a weaker internal parameter or add
the corresponding moment.

## Next Required Steps

1. Decide the canonical empirical target file for the new one-market project.
2. Verify or replace the first-birth timing target before treating it as hard.
3. Add the timing shifter if the timing target remains hard.
4. Freeze normalizations and externally fixed parameters.
5. Write a production calibration driver with the final internal vector,
   explicit bounds, exact moment map, and rejection rules.
6. Only then run a new calibration search.
