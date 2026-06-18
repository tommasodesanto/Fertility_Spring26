# Current One-Market Intergenerational Calibration Identification Ledger

Date: 2026-06-18

## Scope

This note revisits the identification map for the current one-market
intergenerational housing-fertility model in
`code/model/intergen_housing_fertility/`. The comparison class is the older
spatial center-periphery model. The current model intentionally drops location
and should be assessed on its own one-market terms.

The live June 17--18 frontier runs used `J=16`, `Nb=60`, five Markov income
states, `n_house=5`, and `max_iter_eq=3`. They are diagnostic searches, not a
final SMM calibration.

## Internal Parameter Vector

The global-DE diagnostic varies 13 internal parameters:

\[
\Theta =
\{\beta,\alpha,b_0,\bar c_0,\bar c_n,\bar h_0,\bar h_{\mathrm{jump}},
\bar h_n,\psi_{\mathrm{child}},\kappa_n,\chi,\theta_0,\theta_n\}.
\]

In code these are:

| Parameter | Code name | Main role |
|---|---|---|
| Discount factor | `beta` | saving and liquid wealth profile |
| Housing share | `alpha_cons` | housing demand and housing expenditure share |
| Entry liquid wealth | `b_entry_fixed` | young balance sheet and ownership access |
| Baseline goods floor | `c_bar_0` | nonhousing subsistence / residual goods feasibility |
| Child goods cost | `c_bar_n` | child cost, childlessness, completed fertility |
| Baseline housing floor | `h_bar_0` | childless housing services |
| First-child housing jump | `h_bar_jump` | \(0 \to 1\) housing response |
| Extra-child housing slope | `h_bar_n` | \(1 \to 2\) housing response |
| Child utility shifter | `psi_child` | completed fertility level |
| Fertility logit scale | `kappa_fert` | fertility dispersion / childlessness |
| Owner/renter service wedge | `chi` | ownership and owner-renter housing wedge |
| Bequest strength | `theta0` | old-age wealth/ownership retention |
| Child shifter in bequests | `theta_n` | parent-childless old-age gap |

This vector can be internally calibrated only if the target system contains at
least 13 informative moments. A raw count of 13 is not enough; the local
moment-gradient matrix must have rank and the moments must be empirically
comparable to model objects.

## Fixed Or First-Stage Objects

The current global-DE runs hold these outside the internal SMM vector:

| Object | Code name | Treatment |
|---|---|---|
| Interest/user-cost primitives | `q`, `delta`, `tau_H` | external/accounting inputs unless policy moves them |
| Finance limits | `phi`, `pti_limit` | external unless direct LTV/PTI/bind-rate moments are added |
| Transaction/sale friction | `psi` | external unless turnover or downsizing moments are added |
| Income process | `income_age_profile`, `z_grid`, `Pi_z` | first-stage income process |
| Housing menu | `H_own`, `hR_max` | measurement/menu support unless room-distribution moments identify them |
| Supply normalizations | `H0`, `r_bar`, `eta_supply` | market normalization/supply inputs |
| Bequest shift | `theta1` | fixed normalization/support term |
| Numerical controls | grids and solver controls | not economic parameters |

If any of these become free parameters, the target system needs additional
moments for them. In particular, `phi` and `pti_limit` cannot be identified
only from ownership and fertility; they need direct credit-access evidence.

## Candidate 13-Moment System

The live no-timing trial target set, `candidate_no_timing_v0`, has 13 target
moments:

| Moment | Target | Intended block |
|---|---:|---|
| `tfr` | 1.700 | child utility and child costs |
| `childless_rate` | 0.150 | fertility dispersion/extensive margin |
| `own_rate` | 0.57547241 | owner wedge, young access, savings |
| `own_family_gap` | 0.16766167 | family housing need and owner access |
| `housing_increment_0to1` | 0.66443467 | first-child housing need |
| `housing_increment_1to2` | 0.56581378 | additional-child housing need |
| `young_liquid_wealth_to_income` | 0.600 | entry wealth and savings |
| `old_age_own_rate` | 0.76426097 | bequests / old retention |
| `old_age_parent_childless_gap` | 0.070 | child shifter in bequests |
| `liquid_wealth_to_income` | 1.200 | lifecycle saving |
| `housing_user_cost_share` | 0.240 | housing share / price of housing services |
| `prime_childless_renter_median_rooms` | 4.000 | baseline housing floor |
| `prime_childless_owner_median_rooms` | 6.000 | owner service wedge and room menu |

This is exactly identified by count, but the June 17--18 frontier shows that
three moments need measurement or model-object audits before they can be used
together as hard targets: aggregate housing user-cost share, old-age ownership,
and owner median rooms.

## What The Frontier Says

The diagnostic frontier is informative because each target subset fails in a
different direction:

| Target set | Best fit pattern |
|---|---|
| Core feasibility | Fits broad fertility/ownership/room block, but old-age ownership is too high and housing costs are too high. |
| Cost test | Adding the housing cost target still leaves user-cost share around 0.36, above the 0.24 target. |
| Old-age test | Adding old-age targets improves old-age moments, but pushes housing cost share near 0.50. |
| Room-cost test | Gets housing user-cost share near target only when owner median rooms collapse to 4 instead of 6. |

The clean interpretation is not "drop moments." It is:

1. Check whether each target is the correct empirical counterpart.
2. If the code and target are correct, identify the model mechanism that makes
   the joint target unreachable.
3. Replace any invalid hard target with another moment that identifies the same
   parameter block, or fix the affected parameter externally.

## Block-By-Block Identification Assessment

| Block | Current identifying moments | Main weakness | Identification-preserving fix |
|---|---|---|---|
| `beta` | `liquid_wealth_to_income`, `young_liquid_wealth_to_income` | confounded with `b_entry_fixed` and borrowing rules | use an age profile of liquid wealth or fix \(\beta\) externally and keep wealth moments for entry wealth |
| `alpha_cons` | `housing_user_cost_share`, room levels | aggregate user-cost share may not match the model object | replace with tenure-specific rent/user-cost-to-income shares or a clearly measured housing-expenditure share |
| `b_entry_fixed` | young liquid wealth, young ownership access | also moves ownership and fertility through constraints | use young liquid wealth and first-purchase/young-owner moments jointly |
| `c_bar_0` | currently indirect through affordability and saving | no clean direct goods-consumption moment | either fix as a normalization or add a nonhousing consumption/residual-goods share moment |
| `c_bar_n` | `tfr`, `childless_rate` | confounded with `psi_child` and `kappa_fert` | add a valid parity-composition or child-cost moment after resolving the TFR/household-unit convention |
| `h_bar_0` | childless renter/owner rooms | owner median is a coarse rung statistic | use childless renter room distribution and owner-renter room gap/bin shares |
| `h_bar_jump` | `housing_increment_0to1` | relatively clean if event-study object is comparable | keep, but report event-time and household-unit construction |
| `h_bar_n` | `housing_increment_1to2` | proxy is housing demand, not a second-birth hazard | keep as housing-demand moment, not fertility-hazard moment |
| `psi_child` | `tfr` | confounded with child costs | combine with childlessness and a valid parity distribution target |
| `kappa_fert` | `childless_rate` | weak without parity dispersion | add parity shares if measurement convention is coherent |
| `chi` | `own_rate`, `own_family_gap`, owner-renter rooms | confounded with finance and entry wealth | add owner-renter room gap or ownership by income/age; keep finance fixed or directly targeted |
| `theta0` | `old_age_own_rate` | old-age ownership inherits earlier ownership and may not isolate bequest strength | use old liquid wealth, old owner downsizing/retention, or ownership age-profile moments |
| `theta_n` | `old_age_parent_childless_gap` | ownership gap alone may miss wealth-bequest channel | add parent-childless old wealth/liquid wealth or housing-retention gap |

## Recommended Calibration Discipline

There are two coherent ways forward.

### Option A: Keep 13 Internal Parameters

Keep the 13-parameter vector only if we add or validate 13 informative moments.
The weakest current hard moments should be audited or replaced as follows:

| Current weak target | Why weak | Replacement that preserves identification |
|---|---|---|
| `housing_user_cost_share` | aggregate object may not match measured expenditure share | tenure-specific housing cost share, rent-to-income for renters, or owner user-cost-to-income with imputed owner costs |
| `prime_childless_owner_median_rooms` | median rung is discrete and fights cost share | owner room-bin shares or owner-renter room gap conditional on childless prime-age households |
| `old_age_own_rate` | may mostly inherit prime-age ownership | old-age ownership decline, old liquid wealth, old downsizing/retention, or old owner room distribution |
| `liquid_wealth_to_income` | must define net liquid wealth relative to model borrowing state | age-binned liquid-wealth profile with the same asset definition |
| no direct moment for `c_bar_0` | baseline goods floor is otherwise indirect | nonhousing consumption share or fix `c_bar_0` externally |

### Option B: Fix One Or More Parameters Externally

If we do not trust a direct moment for a parameter, reduce the internal vector.
For example, fix `c_bar_0` as a normalization and run a 12-parameter SMM with
12 hard moments. This is not "dropping a bad target"; it is replacing an
unidentified internal parameter with an external restriction.

## Immediate Next Check

Before another cluster calibration, compute a local sensitivity/Jacobian table
at one or two economically plausible candidates:

\[
J_{mp}=\frac{\partial m_m(\Theta)}{\partial \theta_p}.
\]

The table should report signs, normalized elasticities, rank, and near-collinear
parameter pairs. This is the fastest way to separate:

1. a search failure,
2. a measurement mismatch, and
3. a genuine model-impossibility result.

No new SMM objective should be called a production calibration until this
identification audit is complete.
