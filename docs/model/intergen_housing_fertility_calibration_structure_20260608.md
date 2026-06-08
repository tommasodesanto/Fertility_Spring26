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

The code inventory is larger than the economic parameter vector. The active
`setup_parameters()` function defines 145 default `P` fields. Most are derived
objects, state-space choices, market normalizations, inherited one-market
degenerate fields, or numerical controls. The diagnostic random-search driver
currently varies 19 keys:

\[
\{\alpha,\phi,b_0,\psi^{PTI},\bar c_n,\chi,\kappa_n,\psi,
\psi_{\mathrm{child}},\theta_0,\theta_n,\theta_1,\bar h_{\mathrm{jump}},
\bar h_n,h_R^{\max},s_h,\omega_h,h_h^{ref},H_0\}.
\]

That 19-key list is not the production parameter vector. The production vector
must be built from the economic blocks below.

### Intertemporal Saving

| Parameter | Current code name | Production treatment |
|---|---|---|
| Discount factor | `beta` | Internally calibrated if the gross return `R_gross=1+q` is fixed externally. Fixing both \(R\) and \(\beta\) while targeting wealth removes the saving lever. |
| Gross return / period interest | `R_gross`, `q` | Externally fixed or policy-specified. If `q` is internally calibrated, then `beta` must not be treated as separately free without a second saving/return moment. |
| Curvature | `sigma` | Can be externally fixed if the project is not targeting risk/precautionary-saving gradients. If income-risk moments become central, it needs a clear moment or sensitivity band. |
| Entrant liquid wealth | `b_entry_fixed` | Internally calibrated or replaced by an externally estimated entry-wealth distribution. With a single entry wealth scalar, it is naturally tied to young liquid wealth and early owner entry. |

### Stone-Geary Consumption-Housing Preferences

| Parameter | Current code name | Production treatment |
|---|---|---|
| Cobb-Douglas housing share away from subsistence | `alpha_cons` | Internally calibrated if housing expenditure shares are targets. Fixing it is not coherent when Stone-Geary subsistence terms prevent shares from matching by construction. |
| Baseline consumption commitment | `c_bar_0` | Normalization only if a separate consumption floor is imposed and documented. Otherwise internally calibrated jointly with expenditure-share/consumption feasibility moments. |
| Baseline housing need | `h_bar_0` | Internally calibrated if childless housing services or housing expenditure shares are targets. |
| Child goods cost slope | `c_bar_n` | Internally calibrated to fertility composition and childlessness, jointly with child utility. |
| First-child housing jump | `h_bar_jump` | Internally calibrated to the \(0\to1\) housing response. |
| Additional-child housing slope | `h_bar_n` | Internally calibrated to the \(1\to2\) housing response and intensive fertility. |
| Owner service premium | `chi` | Internally calibrated to ownership and owner/renter service wedges. |
| Owner subsistence scale | `owner_h_bar_scale` | Not a primitive structural parameter unless the owner/renter service mapping is explicitly modeled. Treat as a fixed menu/measurement conversion or discipline it with owner-size distribution moments. |

The Stone-Geary block requires at least one direct housing expenditure/share
moment. Without it, \(\alpha\), \(\bar c_0\), and \(\bar h_0\) are weakly
separated from each other.

### Fertility Preferences And Timing

| Parameter | Current code name | Production treatment |
|---|---|---|
| Child utility shifter | `psi_child` | Internally calibrated to completed fertility/TFR. |
| Fertility logit scale | `kappa_fert` | Internally calibrated to parity dispersion and childlessness. It should not be the only first-birth-timing parameter. |
| Age timing shifter | not implemented | Required if a first-birth timing target remains hard. This can be a compact age shifter \(\lambda_a\) or an age-specific child-cost schedule. |
| Fertility window and child duration | `A_f_start`, `A_f_end`, `A_m`, `stage_durations` | Model-design/biological timing inputs. Do not tune them to fit the calibration objective without changing the model statement. |

### Housing Finance And Transaction Frictions

| Parameter | Current code name | Production treatment |
|---|---|---|
| Financed share | `phi` | First-stage external from LTV/down-payment evidence, or internally calibrated only with direct credit moments. It is not identified by fertility/ownership moments alone. |
| Payment-to-income limit | `pti_limit` | Same: first-stage external or internally calibrated with PTI/denial/bind-rate moments. |
| Transaction/sale haircut | `psi` | First-stage external from transaction costs or internally calibrated with turnover/downsizing/tenure-duration moments. |
| Property tax | `tau_H` | Externally fixed for the benchmark and moved in policy experiments. |
| Depreciation | `delta` | Externally fixed. |
| Payment formula | `use_pti_constraint`, current simplified formula | Model specification. If Coven-style amortized payments are implemented, the parameter map changes. |

### Bequests And Old Retention

| Parameter | Current code name | Production treatment |
|---|---|---|
| Bequest intensity | `theta0` | Internally calibrated to old-age ownership/wealth. |
| Child shifter in bequest utility | `theta_n` | Internally calibrated to old parent-childless ownership/wealth gaps. |
| Bequest utility shift | `theta1` | Normalization/curvature support. Fix unless separately identified. |
| Old retention wedge | not cleanly implemented | Required if incumbent retention is a central mechanism. Calibrate with old-owner large-home retention, downsizing, or turnover moments. |

### Income And Retirement

| Parameter | Current code name | Production treatment |
|---|---|---|
| Lifecycle income profile | `income_age_breaks`, `income_age_values` | First-stage external estimate. |
| Persistent income states | `z_grid`, `z_weights`, `income_shock_persistence`, `Pi_z` | First-stage external estimate from income moments. They should not be chosen by the fertility/housing SMM unless those income moments are included. |
| Payroll tax and pension closure | `tau_pay`, `pension_mode`, `pension` | External/accounting closure. |

### Housing Market, Menu, And Normalizations

| Parameter | Current code name | Production treatment |
|---|---|---|
| Aggregate supply level | `H0` | Market normalization or calibrated market-clearing shifter tied to a price/rent/quantity level. It is not a household preference parameter. |
| Rent/user-cost anchor | `r_bar` | Normalization paired with `H0` and price units. |
| Supply elasticity | `eta_supply`, `xi_supply` | External supply estimate or calibrated only with a supply elasticity/quantity response target. |
| Owner rung menu | `H_own`, `n_house` | Measurement/menu choice based on room-unit support. Do not search it as a preference parameter. |
| Renter cap | `hR_max` | Measurement/menu choice or calibrated only with renter-size distribution moments. |
| Large owner unit cost | `owner_size_cost`, `owner_size_cost_ref`, `owner_size_cost_power` | Menu/supply device. Include only if owner-size distribution moments are in the target system. |

### Not Production Economic Parameters

Inherited one-market-degenerate fields should not enter the new production
calibration: `kappa_loc`, `eps_loc`, `E_loc`, `mu_stay`, `mu_move`,
`mu_move_parent`, `mu_age_decay`, `location_choice_form`, `kappa_entry`,
`outside_value`, `outside_entry_flow`, `outside_entry_shares`,
`local_birth_entry_weight`, `renewal_retention`,
`renewal_calibrate_outside_flow`, and related population-closure fields.

Policy toggles should enter only when the policy is part of the model:
`parent_dp_waiver`, `parent_dp_waiver_phi`, `birth_dp_grant`,
`birth_entry_grant`, `birth_entry_grant_amount`, and their filters.

Pure numerical controls are not economic calibration parameters: grids,
iteration limits, tolerances, damping controls, price bounds, interpolation
flags, and kernel flags.

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

The current diagnostic screen fails before any moment-counting exercise:

\[
19 \text{ varied knobs} \ne \text{production economic parameter vector}.
\]

The production count depends on first-stage decisions:

- If \(R\), income process, finance limits, supply elasticity, housing menu,
  and bequest shift are externally/first-stage fixed, the implemented internal
  vector is roughly
  \[
  \{\beta,\alpha,b_0,\bar c_0,\bar c_n,\bar h_0,
  \bar h_{\mathrm{jump}},\bar h_n,\psi_{\mathrm{child}},\kappa_n,
  \chi,\theta_0,\theta_n\},
  \]
  which is 13 parameters.
- If one Stone-Geary baseline term is imposed as a normalization rather than
  calibrated, this falls to 12.
- If first-birth timing remains a hard target, one timing parameter must be
  added, taking the count to 13-14.
- If `phi`, `pti_limit`, or `psi` are internally calibrated, add one parameter
  and one direct credit/turnover moment for each. They should not be identified
  only from fertility and ownership.
- If `H0` or `r_bar` is used to match a price/rent level, treat it as a
  separate market-normalization equation, not as a household preference
  parameter.

## Identification Map

A disciplined map starts from the parameters:

| Parameter | Primary moment | Secondary moments/checks |
|---|---|---|
| `beta` | liquid wealth / saving profile | tenure timing; cannot be fixed jointly with \(R\) if wealth is targeted |
| `alpha_cons` | housing expenditure share | housing services by age/tenure |
| `b_entry_fixed` | `young_liquid_wealth_to_income` | young ownership and owner-entry thresholds |
| `c_bar_0` | baseline consumption/housing expenditure feasibility | consumption floor diagnostics |
| `h_bar_0` | childless housing services | childless renter/owner size distribution |
| `c_bar_n` | `childless_rate` and completed fertility mix | TFR and young consumption feasibility |
| `psi_child` | `tfr` | childlessness |
| `kappa_fert` | parity dispersion / childlessness | timing only weakly, not as the sole timing margin |
| new timing shifter, not implemented | first-birth timing object | age profile of first births |
| `h_bar_jump` | `housing_increment_0to1` | family ownership gap |
| `h_bar_n` | `housing_increment_1to2` | intensive completed fertility |
| `chi` | `own_rate_3055` | family ownership gap and owner/renter housing services |
| `theta0` | `old_age_own_rate` | lifecycle ownership slope |
| `theta_n` | `old_age_parent_childless_gap` | parent-childless wealth and old housing retention |
| `phi`, if internal | LTV/down-payment/bind-rate moment | ownership and fertility responses |
| `pti_limit`, if internal | PTI/bind-rate or denial moment | owner entry by income |
| `psi`, if internal | turnover/downsizing/tenure-duration moment | old retention and owner mobility |

Finance parameters should be handled deliberately:

- If `phi` and `pti_limit` are first-stage/external, they do not enter the
  internal production vector.
- If either `phi` or `pti_limit` is internally calibrated, add at least one
  direct credit moment per added parameter. Otherwise the financial constraint
  is not separately identified from entry wealth, discounting, or the ownership
  premium.

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

## Production Internal Vector

Given the current simplified one-market model, the production baseline should
start from the parameter vector

\[
\Theta_{\mathrm{base}} =
\{\beta,\alpha,b_0,\bar c_0,\bar c_n,\bar h_0,
\bar h_{\mathrm{jump}},\bar h_n,\psi_{\mathrm{child}},\kappa_n,
\chi,\theta_0,\theta_n,\lambda_a\},
\]

where \(\lambda_a\) is a required timing object if first-birth timing is a hard
target. This is a 14-parameter vector if \(\bar c_0\) is calibrated. If
\(\bar c_0\) is imposed as a normalization, the vector is 13 parameters.

First-stage or externally set objects:

- \(R\)/`q`, `delta`, and benchmark `tau_H`;
- `sigma`, unless risk/precautionary-saving gradients are targeted;
- income profile and Markov income process;
- `theta1`;
- supply elasticity and housing menu, unless size/supply moments are added;
- numerical solver controls.

Conditional internal additions:

- Add `phi` only with direct down-payment/LTV/bind-rate moments.
- Add `pti_limit` only with direct PTI/bind-rate/denial moments.
- Add `psi` only with turnover/downsizing/tenure-duration moments.
- Add an old-retention wedge only after implementing it cleanly and adding an
  old-home retention or downsizing moment.

This is the production logic: write down the internal parameter vector first,
then add at least as many moments as the vector requires, with the moments tied
to the corresponding mechanisms.

## Candidate No-Timing Trial Ledger

After the June 8 calibration discussion, the code defines a documented trial
target set named `candidate_no_timing_v0`. This is not a finalized empirical
target system. It is a disciplined trial ledger for learning how the simplified
one-market model behaves when the badly mapped first-birth-age target is removed
from the objective.

The internal parameter vector for this trial is
\[
\Theta_{\mathrm{no\ timing}} =
\{\beta,\alpha,b_0,\bar c_0,\bar c_n,\bar h_0,
\bar h_{\mathrm{jump}},\bar h_n,\psi_{\mathrm{child}},\kappa_n,
\chi,\theta_0,\theta_n\},
\]
which has 13 parameters. The 13 target moments are:

| Moment | Target | Status / source note |
|---|---:|---|
| `tfr` | 1.700 | existing old-target carryover |
| `childless_rate` | 0.150 | existing old-target carryover |
| `own_rate` | 0.57547241 | prime-age ownership, existing carryover |
| `own_family_gap` | 0.16766167 | new-parent minus non-parent ownership, existing carryover |
| `housing_increment_0to1` | 0.66443467 | first-child housing response, existing carryover |
| `housing_increment_1to2` | 0.56581378 | additional-child housing response proxy, existing carryover |
| `young_liquid_wealth_to_income` | 0.600 | selected young-balance-sheet target |
| `old_age_own_rate` | 0.76426097 | existing current one-market code target |
| `old_age_parent_childless_gap` | 0.070 | existing PSID completed-fertility old-age gap |
| `liquid_wealth_to_income` | 1.200 | archived SCF 45--55 liquid wealth/income target; candidate reactivation |
| `housing_user_cost_share` | 0.240 | archived alpha-share target; candidate reactivation |
| `prime_childless_renter_median_rooms` | 4.000 | archived room-size target; candidate reactivation |
| `prime_childless_owner_median_rooms` | 6.000 | archived room-size target; candidate reactivation |

The extra four moments are deliberately labeled as candidate reactivations,
because their empirical construction has not been re-audited for the new
one-market project. They are suitable for a documented diagnostic run, not for
claiming final production calibration.

Parity composition is kept diagnostic in this trial. Under the project
reporting convention
\[
TFR = 2E[n],
\]
the targets \(TFR=1.70\) and \(\Pr(n=0)=0.15\) imply \(E[n]=0.85\) and
\(\Pr(n>0)=0.85\). Since the smallest positive model family size is \(n=1\),
those two targets already imply zero mass on \(n\ge2\) if both are hit exactly.
Therefore a positive `parity_share_2plus` hard target is mechanically
inconsistent unless the fertility measurement convention or fertility state is
changed. The code still reports `parity_share_0`, `parity_share_1`, and
`parity_share_2plus` ex post.

## Next Required Steps

1. Decide the canonical empirical target file for the new one-market project.
2. Verify or replace the first-birth timing target before treating it as hard.
3. Add the timing shifter if the timing target remains hard.
4. Freeze normalizations and externally fixed parameters.
5. Write a production calibration driver with the final internal vector,
   explicit bounds, exact moment map, and rejection rules.
6. Only then run a new calibration search.
