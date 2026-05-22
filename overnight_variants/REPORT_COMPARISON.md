# Overnight Variant Comparison

Date: 2026-05-22

## Implemented Scope

Test 1 copied the benchmark Python discrete-time model and added a finite
diagnostic augmentation for \(z\) and \(\mu\): a three-state persistent earnings
grid and a five-state mortgage-account grid with good/bad credit, balance bins,
default-hazard diagnostics, and bad-credit accounting. The copied benchmark
solver ran to an accepted equilibrium, but \(z\) and \(\mu\) are not yet inside
the Bellman recursion. Treat this as a state-accounting smoke test, not a full
HANK/mortgage implementation.

Test 2 copied the benchmark model and patched the copied solver so owner
budgets, maintenance, down-payment constraints, and sale proceeds can use
type-specific \(p_{iq(k)}\). Renters use a two-type rent approximation keyed to
the state's housing need. The smoke test pinned ordinary prices and updated
only middle-unit prices from the developer inverse supply curve. It is a partial
type-clearing test, not a full \(2\times2\) general-equilibrium fixed point.

## Moment Comparison

| Moment | Target | Current Benchmark | Income/Mortgage | Developer Middle |
|---|---:|---:|---:|---:|
| `tfr` | 1.700 | 1.898 | 1.890 | 1.838 |
| `childless_rate` | 0.150 | 0.145 | 0.148 | 0.159 |
| `mean_age_first_birth` | 26.000 | 33.535 | 33.553 | 33.346 |
| `tfr_gradient` | 0.133 | 0.119 | 0.119 | 0.100 |
| `own_rate` | 0.627 | 0.643 | 0.647 | 0.000 |
| `own_gradient` | 0.170 | 0.139 | 0.145 | 0.000 |
| `own_family_gap` | 0.110 | 0.114 | 0.106 | 0.000 |
| `prime_childless_renter_median_rooms` | 4.000 | 6.365 | 6.359 | 6.552 |
| `prime_childless_owner_median_rooms` | 6.000 | 6.800 | 6.800 | 6.800 |
| `housing_increment_0to1` | 0.664 | 0.441 | 0.442 | 0.740 |
| `housing_increment_1to2` | 0.566 | 0.192 | 0.176 | 0.314 |
| `young_liquid_wealth_to_income` | 0.600 | 0.527 | 0.511 | 0.600 |
| `center_share_nonparents` | 0.494 | 0.405 | 0.402 | 0.421 |
| `center_share_newparents` | 0.416 | 0.382 | 0.378 | 0.394 |
| `migration_rate` | 0.032 | 0.035 | 0.035 | 0.038 |
| `old_age_own_rate` | 0.863 | 0.947 | 0.836 | 0.000 |
| `old_age_parent_childless_gap` | 0.070 | 0.062 | -0.018 | 0.000 |
| `inv_pop_share_C` | 0.450 | 0.441 | 0.439 | 0.450 |
| `inv_rent_ratio_C_over_P` | 1.140 | 1.182 | 1.185 | 1.182 |

Notes: current benchmark values are from `CALIBRATION_STATUS.md`. Variant
values are smoke-test outputs and are not globally recalibrated.

## Dimensions And Loops

| Object | Current Benchmark | Income/Mortgage | Developer Middle |
|---|---:|---:|---:|
| household state | \((b,d,i,a,n,s)\) | diagnostic \((b,d,i,a,n,s,z,\mu)\) | \((b,d,i,a,n,s)\) |
| \(b\) states | 80 | 80 | 40 in smoke |
| tenure states | 7 | 7 | 7 |
| locations | 2 | 2 | 2 |
| ages | 60 | 60 | 60 |
| fertility states | 4 | 4 | 4 |
| child-age states | 7 | 7 | 7 |
| \(z\) states | 0 | 3 | 0 |
| \(\mu\) states | 0 | 5 | 0 |
| price dimensions | \(p_i\), 2 scalars | \(p_i\), 2 scalars | \(p_{iq}\), 4 prices with ordinary pinned |
| fixed-point loop | scalar GE prices plus entry closure | copied scalar GE only | partial PE outer loop for middle prices |

## Verdicts

| Branch | Verdict | Reason |
|---|---|---|
| Income risk + mortgage account | yellow | Copied benchmark equilibrium solved and diagnostic \(z,\mu\) distributions are non-degenerate, but the new states are not yet structural. Default hazards are essentially inactive. |
| Developer missing-middle supply | red | Partial middle-price loop stabilized, but ordinary prices were pinned, type excess remained large, the room screen did not improve, and ownership collapsed. |

## Second Pass

After the first pass failed, I added second-pass drivers rather than overwriting
the original evidence.

Test 1 V2 (`run_income_mortgage_risk_v2_scenarios.py`) solves 15
partial-equilibrium scenario economies. In those solves \(z\) scales earnings,
bad credit changes \(\phi_g\), and active mortgage accounts add the compact
owner-payment wedge \(\rho_gm\). This is still not a joint stationary
distribution over \((z,\mu)\), but the new objects now enter choices.

Test 2 V2 (`run_developer_missing_middle_v2.py`) fixes the main mechanical
error in V1. It uses three types, \(S/M/L\), with \(M\) bounded to the 5--6.5
room interval, and it calibrates type-specific supply curves to pass through
baseline scalar-price type demand before moving \(M\) and \(L\) prices.

| Moment | Target | Current Benchmark | Income/Mortgage V2 | Developer V2 |
|---|---:|---:|---:|---:|
| `tfr` | 1.700 | 1.898 | 1.532 | 1.956 |
| `childless_rate` | 0.150 | 0.145 | 0.388 | 0.128 |
| `mean_age_first_birth` | 26.000 | 33.535 | 34.333 | 33.353 |
| `tfr_gradient` | 0.133 | 0.119 | 0.144 | 0.119 |
| `own_rate` | 0.627 | 0.643 | 0.531 | 0.606 |
| `own_gradient` | 0.170 | 0.139 | 0.101 | 0.025 |
| `own_family_gap` | 0.110 | 0.114 | 0.177 | 0.062 |
| `prime_childless_renter_median_rooms` | 4.000 | 6.365 | 5.972 | 6.232 |
| `prime_childless_owner_median_rooms` | 6.000 | 6.800 | 6.870 | 8.200 |
| `housing_increment_0to1` | 0.664 | 0.441 | 0.505 | 0.382 |
| `housing_increment_1to2` | 0.566 | 0.192 | 0.228 | 0.259 |
| `young_liquid_wealth_to_income` | 0.600 | 0.527 | 0.443 | 0.560 |
| `center_share_nonparents` | 0.494 | 0.405 | 0.340 | 0.421 |
| `center_share_newparents` | 0.416 | 0.382 | 0.293 | 0.399 |
| `migration_rate` | 0.032 | 0.035 | 0.031 | 0.034 |
| `old_age_own_rate` | 0.863 | 0.947 | 0.845 | 0.778 |
| `old_age_parent_childless_gap` | 0.070 | 0.062 | -0.058 | -0.050 |
| `inv_pop_share_C` | 0.450 | 0.441 | 0.373 | 0.442 |
| `inv_rent_ratio_C_over_P` | 1.140 | 1.182 | 1.182 | 1.182 |

Second-pass verdicts:

| Branch | V2 Verdict | Reason |
|---|---|---|
| Income risk + mortgage account | yellow | \(z\) and \(\mu\) now affect solved choices in scenario economies; spreads are large enough to matter, but aggregate childlessness, geography, and old-age ownership gaps deteriorate. |
| Developer missing-middle supply | yellow | Type clearing is now numerically stable and ownership no longer collapses; still not green because owner rooms worsen and \(H01\) falls. |

## Recommendation

Do another diagnostic before coding either branch into the live model. The
second pass makes Branch 2 the better next diagnostic target, not because it is
ready, but because the corrected type map gives a stable price system without
destroying ownership. The next Branch 2 diagnostic should vary the type
boundaries and supply calibration against ACS room-bin targets before any live
solver merge.

Do not code Branch 1 into the live model from the scenario mixture. The margins
have bite, but the scenario mixture is not a coherent stationary model and it
breaks too many core moments. Branch 1 should proceed only as an explicit
state-space expansion over \(z\) and \(\mu\).

## Failure Modes

- Test 1 V1 failure mode: \(z\) and \(\mu\) were bookkeeping diagnostics, not
  state variables in household optimization. The default/bad-credit margin was
  essentially inactive, with maximum diagnostic default hazard around
  \(4.15\times10^{-4}\).
- Test 1 V2 remaining failure mode: \(z\), credit status, and mortgage-payment
  wedges now move choices, but only through separate scenario solves. The
  output is not a coherent joint distribution over \((z,\mu)\), and it breaks
  childlessness, geography, and the old-age parent-childless ownership gap.
- Test 2 V1 failure mode: the two-type map treated too much of the owner ladder
  as middle housing, so middle prices became high enough to wipe out ownership.
- Test 2 V2 remaining failure mode: the corrected \(S/M/L\) type clearing is
  numerically stable, but ordinary prices remain pinned, owner rooms move away
  from the target, and \(H01\) falls relative to the current benchmark.

## Files Outside `overnight_variants/`

None touched by this task. The starting git status already had modified live
code, docs, and LaTeX files; those were read but not edited.

## Rerun Commands

Test 1:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_smoke.py --quiet --max-iter-eq 35
```

Test 2:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_developer_missing_middle
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_developer_missing_middle_smoke.py --quiet --iterations 4 --nb 40
```

Test 1 V2:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_v2_scenarios.py --quiet --nb 40
```

Test 2 V2:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_developer_missing_middle
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_developer_missing_middle_v2.py --quiet --iterations 5 --nb 50
```
