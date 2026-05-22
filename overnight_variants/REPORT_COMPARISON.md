# Overnight Variant Comparison

Date: 2026-05-22

## Implemented Scope

Test 1 is now the standard HANK-style income-risk branch, not a two-new-state
mortgage branch. `run_income_mortgage_risk_v5_hank_z_outside_closure.py` puts
one finite idiosyncratic earnings state \(z\) directly into the copied Bellman
arrays, policy arrays, fertility/location probabilities, tenure choices, and
forward distribution, then solves the copied price and entry fixed point under
the paper-facing outside-option scale closure. The current serious grid is
Rouwenhorst \(N_z=7\), \(\rho_z=0.95\), unconditional
\(\sigma_z=0.35\). The mortgage-account object \(\mu\) is not structural in
this decision run. The 2026-05-22 engineering pass added a compiled
HANK-\(z\) forward-distribution kernel for the GE loop; the final
full-statistics pass is retained for reporting.

Test 2 is now a full type-price developer GE prototype. `run_developer_missing_middle_v3_ge.py`
updates all six \(p_{iq}\) and \(r_{iq}\) objects for \(i\in\{C,P\}\) and
\(q\in\{S,M,L\}\), plus entry shares. The copied solver prices owner rungs by
\(p_{iq(k)}\) and, after the V3 fix, prices renters by realized room interval
\(r_{iq(h^R)}h^R\). Fixed costs \(F_{iq}\) are reported through entry
thresholds but are not used for discrete shutdown.

## GE Moment Comparison

Current benchmark values are from `CALIBRATION_STATUS.md`. Variant values are
coarse smoke outputs, not full recalibrations.

| Moment | Target | Current Benchmark | HANK-z V4 GE | HANK-z V5 Outside | Developer V3 GE |
|---|---:|---:|---:|---:|---:|
| `tfr` | 1.700 | 1.898 | 1.858 | 1.562 | 1.899 |
| `childless_rate` | 0.150 | 0.145 | 0.181 | 0.351 | 0.148 |
| `mean_age_first_birth` | 26.000 | 33.535 | 34.318 | 35.279 | 33.698 |
| `tfr_gradient` | 0.133 | 0.119 | 0.087 | -0.194 | 0.110 |
| `own_rate` | 0.627 | 0.643 | 0.639 | 0.608 | 0.578 |
| `own_gradient` | 0.170 | 0.139 | -0.066 | -0.127 | 0.018 |
| `own_family_gap` | 0.110 | 0.114 | 0.154 | 0.354 | 0.070 |
| `prime_childless_renter_median_rooms` | 4.000 | 6.365 | 6.180 | 6.001 | 6.500 |
| `prime_childless_owner_median_rooms` | 6.000 | 6.800 | 6.800 | 6.800 | 8.200 |
| `housing_increment_0to1` | 0.664 | 0.441 | 0.415 | 0.589 | 0.418 |
| `housing_increment_1to2` | 0.566 | 0.192 | 0.402 | 1.494 | 0.252 |
| `young_liquid_wealth_to_income` | 0.600 | 0.527 | 0.935 | 1.873 | 0.815 |
| `center_share_nonparents` | 0.494 | 0.405 | 0.375 | 0.233 | 0.412 |
| `center_share_newparents` | 0.416 | 0.382 | 0.372 | 0.371 | 0.394 |
| `migration_rate` | 0.032 | 0.035 | 0.035 | 0.033 | 0.035 |
| `old_age_own_rate` | 0.863 | 0.947 | 0.772 | 0.753 | 0.767 |
| `old_age_parent_childless_gap` | 0.070 | 0.062 | 0.004 | 0.263 | -0.044 |
| `inv_pop_share_C` | 0.450 | 0.441 | 0.433 | 0.385 | 0.447 |
| `inv_rent_ratio_C_over_P` | 1.140 | 1.182 | 1.186 | 1.156 | 1.178 |

## Dimensions And Loops

| Object | Current Benchmark | HANK-z V5 Outside | Developer V3 GE |
|---|---:|---:|---:|
| household state | \((b,d,i,a,n,s)\) | \((b,d,i,a,n,s,z)\) | \((b,d,i,a,n,s)\) |
| \(b\) states | 80 | 30 | 30 |
| tenure states | 7 | 7 | 7 |
| locations | 2 | 2 | 2 |
| ages | 60 | 60 | 60 |
| fertility states | 4 | 4 | 4 |
| child-age states | 7 | 7 | 7 |
| \(z\) states | 0 | 7 | 0 |
| \(\mu\) states | 0 | 0 | 0 |
| price dimensions | \(p_i\), 2 prices | \(p_i\), 2 prices | \(p_{iq}\), 6 prices |
| rent dimensions | scalar user cost by \(i\) | scalar user cost by \(i\) | \(r_{iq}\), 6 rents |
| closure / fixed-point loop | renewal-valve prices plus entry | outside-option scale closure, scalar prices plus entry | all type prices plus entry |
| accepted GE run | live benchmark | yes, strict tolerance | yes, coarse type-price tolerance |
| cost category | benchmark reference | project-scale | expensive |

## Verdicts

| Branch | Verdict | Reason |
|---|---|---|
| Branch 1: HANK-z income risk | yellow | Full GE clears with a real Rouwenhorst \(z\) state and the outside-option closure, but the un-recalibrated moments are not acceptable: gradients flip sign, fertility is too low/late, and young liquid wealth overshoots sharply. |
| Branch 2: developer missing-middle supply | yellow | Full type-price GE clears after the realized-rent fix and middle demand is non-degenerate, but the room screen does not improve, \(H01\) falls, and the ownership gradient is nearly flat. |

## Recommendation

Recommendation: **code Branch 1 first, but only as a controlled live branch with
immediate recalibration diagnostics**.

The HANK-z branch is the cleaner next live branch because it is the canonical
one-state income-risk extension the user asked for, it uses a standard
Rouwenhorst process, and it now clears in the copied GE loop under the
paper-facing outside-option closure. When moving it out of
`overnight_variants/`, carry over exactly one new state \(z\) first. Do not add
\(\mu\) in the same live branch; the current evidence says the next problem is
recalibrating the HANK-z economy and disciplining transitions, not expanding the
state space again.

Follow-up borrowing-wedge diagnostics reinforce that recommendation. With the
existing \(\phi\)-based down-payment and borrowing-floor wedge, high-\(z\)
renters buy much more often than low-\(z\) renters, 0.091 versus 0.015, and are
more often feasible for the starter down payment, 0.590 versus 0.336. Low-\(z\)
owners are more likely to sit near the borrowing floor, 0.089 versus 0.029.
This means the current model already has a meaningful mortgage-like
collateral/liquidity channel once \(z\) is structural.

Do not merge the developer branch into live code yet. The corrected full-GE
prototype is useful and no longer fails mechanically, but it does not solve the
room target problem: prime childless renter median rooms are 6.500 against a
target of 4.000, prime childless owner median rooms are 8.200 against a target
of 6.000, and middle owner mass is essentially zero. The next Branch 2 step
should be a narrower diagnostic on room-bin/type mapping and the owner rung
ladder before another solver merge.

## Failure Modes

- Branch 1 fixed closure issue: the first outside-option attempt reused
  \(\kappa_\ell\) for the entry/outside margin. That made \(q^E\) jump to
  zero or one because entry values are lifetime-utility objects. The accepted
  V5 run uses a separate \(\kappa_E=10^6\), keeps \(q^E=0.90008\), and has
  final scale \(S=1.00023\).
- Branch 1 fixed income-grid issue: the Rouwenhorst transition was correct,
  but the stationary-distribution helper initially clipped a signed eigenvector
  and returned uniform weights. It now uses power iteration and reports the
  correct binomial stationary weights.
- Branch 1 remaining failure mode: the \(z\) state is economically active and
  the closure works, but the ownership gradient is \(-0.127\) versus target
  \(0.170\), the fertility gradient is \(-0.194\) versus target \(0.133\),
  and young liquid wealth/income is 1.873 versus target 0.600. This is a
  recalibration and transition-discipline problem.
- Branch 1 implementation boundary: \(\mu\) is deliberately absent from the GE
  decision run. The earlier \(\mu\) work remains diagnostic and should not be
  interpreted as a solved mortgage/default account.
- Branch 2 fixed first-pass failure: the earlier renter block priced by
  \(\bar h_n\), which mismatched the demand diagnostic by realized \(h^R\).
  V3 fixes this by solving renter choices piecewise over \(S/M/L\) intervals.
- Branch 2 remaining failure: full type-price GE clears, but almost all owner
  mass is still in large units, renter room mass is almost entirely in 7--8
  rooms, \(H01\) is 0.418 versus target 0.664, and the ownership gradient is
  only 0.018 versus target 0.170.
- Branch 2 fixed costs: \(F_{iq}\) is documented through entry-threshold
  diagnostics only; it is not active as a discrete type-shutdown condition.

## Files Outside `overnight_variants/`

None touched by this task. The repo had pre-existing modified live code, docs,
and LaTeX files outside `overnight_variants/`; those were not edited or staged.

## Rerun Commands

Branch 1 full GE HANK-z:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_v5_hank_z_outside_closure.py --quiet --nb 30 --nz 7 --rho-z 0.95 --sigma-z 0.35 --kappa-entry 1000000 --baseline-max-iter-eq 35 --max-iter-eq 60 --tol-eq 5e-4
```

Branch 2 full GE developer supply:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_developer_missing_middle
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_developer_missing_middle_v3_ge.py --quiet --nb 30 --iterations 12 --price-damp 0.25 --entry-damp 0.25
```

Branch 1 borrowing-wedge diagnostic:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_hank_z_borrowing_wedge_diagnostics.py --quiet --nb 30 --nz 3 --max-iter-eq 35
```
