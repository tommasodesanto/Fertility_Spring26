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

## Recommendation

Do another diagnostic before coding either branch into the live model. If one
branch must be prioritized after this, Branch 2 is still the more relevant
economic direction because it targets the live room-distribution failure, but
the prototype shows the full type-price GE loop and renter pricing need to be
designed more carefully before production work. Do not code Branch 1 first
unless the next prototype actually expands the Bellman and distribution arrays
over \(z\) and \(\mu\).

## Failure Modes

- Test 1 failure mode: \(z\) and \(\mu\) are bookkeeping diagnostics, not state
  variables in household optimization. The default/bad-credit margin is mostly
  inactive, with maximum diagnostic default hazard around \(4.15\times10^{-4}\).
- Test 1 economic warning: the smoke solve moves old-age ownership and the
  old-age parent-childless gap noticeably relative to the status table even
  before structural risk is active, so this branch needs strict target checks.
- Test 2 failure mode: middle prices became high enough to wipe out ownership
  in the partial PE allocation.
- Test 2 market-clearing warning: ordinary demand was nearly zero while ordinary
  supply stayed pinned, so the reported type wedge is not an equilibrium wedge.
- Test 2 mechanism warning: the H2/starter-family margin is not alive in the
  smoke output; new-parent feasible-but-rent share is essentially one.

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
