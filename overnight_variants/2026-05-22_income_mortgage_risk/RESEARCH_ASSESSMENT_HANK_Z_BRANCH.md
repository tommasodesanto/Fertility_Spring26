# Research Assessment: HANK-z Income-Risk Branch

Date: 2026-05-23

Scope: isolated prototype under
`overnight_variants/2026-05-22_income_mortgage_risk/`. This note does not
evaluate or modify live code under `code/model/dt_cp_model/`.

## Bottom Line

The HANK-\(z\) branch is mechanically encouraging but economically incomplete.
It now solves a real structural earnings-risk state in GE, uses a standard
Rouwenhorst process, and implements the paper-style outside-option benchmark
normalization correctly:
\[
S E_0(p)=q^E(p)\left[M+S B_0(p)\right],
\qquad
S=1,\qquad
M=E_0(p)/q^E(p)-B_0(p).
\]
That is a useful numerical achievement.

It is not yet evidence of a calibrated benchmark. The best formal global-search
candidate has objective loss `30.5893`, no penalty, and strict GE convergence,
but it obtains the low loss by sacrificing the tenure and spatial price
structure:

| Moment | Target | Best formal candidate |
|---|---:|---:|
| `own_rate` | 0.627 | 0.213 |
| `old_age_own_rate` | 0.863 | 0.297 |
| `inv_rent_ratio_C_over_P` | 1.140 | 0.854 |
| `mean_age_first_birth` | 26.000 | 31.606 |
| `prime_childless_owner_median_rooms` | 6.000 | 9.519 |

This is not an acceptable economics result. A model in which the calibrated
center becomes cheaper than the periphery and ownership collapses cannot be the
paper benchmark, even if the scalar objective falls.

The right next step is targeted cross-sectional diagnostics, followed by a
guardrailed/reweighted calibration only if those diagnostics show useful
mechanisms. Do not add a structural mortgage-account state \(\mu\) yet.

## What Should Count As Success

Success for this branch should be defined relative to the current
renewal-valve benchmark, not relative to the scalar loss alone. The branch is
worth keeping only if structural earnings risk adds disciplined
cross-sectional content that the renewal-valve benchmark cannot deliver.

Minimum success criteria:

1. **Numerical validity.** Strict GE convergence at the validation grid
   (`Nz=7`), finite outside-option scale, \(S=1\) in the benchmark,
   \(q^E=0.9\), positive residual outside-born flow \(M\), and stable
   economics moving from `Nz=5` to `Nz=7`.

2. **Aggregate guardrails.** The search cannot be allowed to trade away core
   objects. Before comparing losses, candidate benchmarks should satisfy hard
   admissibility floors such as:
   \[
   \text{own}_{30-55}>0.5,\quad
   \text{own}_{65-75}>0.7,\quad
   \nabla \text{own}>0,\quad
   p_C/p_P>1,
   \]
   plus finite-scale closure and strict GE acceptance. These are not final
   targets; they are screens against economically invalid optima.

3. **Cross-sectional mechanism.** Conditional on age and wealth, the \(z\)
   state should generate meaningful heterogeneity in fertility, timing,
   tenure, housing demand, location, and migration. The mechanism should not
   be only a relabeling of wealth.

4. **Counterfactual heterogeneity.** Housing-cost counterfactuals should have
   different effects by \(z\), liquid wealth, tenure, and parent status. If all
   aggregate effects are unchanged relative to the renewal-valve benchmark, the
   state-space expansion is not earning its cost.

5. **Validation against data.** The cross-sectional predictions need empirical
   counterparts, even if first used as untargeted diagnostics. The branch
   should face facts from ACS/PSID/SCF-style data on fertility, ownership,
   rooms, wealth, and location by income or earnings rank.

## Cross-Sectional Facts The Branch Should Match Or Generate

The branch should be judged on a set of facts that map directly to the added
state \(z\). The empirical object should be persistent or current earnings
rank, not education unless education is explicitly treated as a proxy for
permanent income.

Recommended diagnostic targets:

| Domain | Cross-sectional fact to discipline |
|---|---|
| Fertility | TFR, childlessness, and age at first birth by income/earnings rank and wealth rank. |
| Fertility timing | Birth hazards over ages 22-45 by \(z\) and wealth. A successful model should explain why low-resource households delay or accelerate relative to high-resource households. |
| Tenure | Ownership by age, \(z\), wealth, parent status, and location. |
| Housing demand | Rooms by \(z\), tenure, and parity; housing response to first and second births by \(z\) or income rank. |
| Spatial sorting | Center share by \(z\), parity, and tenure; parent versus nonparent sorting within each \(z\). |
| Migration | Center-periphery moves by \(z\), birth status, tenure, and wealth. |
| Finance wedge | Down-payment feasibility, renter-to-owner transition shares, owners near the borrowing floor, and liquid-wealth buffers by \(z\). |
| Counterfactuals | \(\Delta\)TFR, \(\Delta\)childlessness, \(\Delta\)ownership, \(\Delta\)rooms, and \(\Delta\)center share after housing-cost shocks by \(z\) and wealth. |

These are not all calibration targets. The immediate task is to measure them
in the model and, where possible, construct matching empirical bins. If the
model gives the wrong signs before calibration, more search is unlikely to
solve the research question.

## What HANK-z Adds Relative To Renewal Valve

The live benchmark closure in `CALIBRATION_STATUS.md` is
`renewal_valve_calibrated`. Its population scale channel is accounting-based:
the renewal valve changes city scale through \(M\), \(B_0(p)\), and a retention
object, but it does not create persistent idiosyncratic earnings histories.

The HANK-\(z\) branch adds three distinct objects:

1. **Persistent earnings heterogeneity.** Households differ in current and
   expected future income through a Markov \(z\) state. This can affect the
   timing and quantity of children, not just current budget feasibility.

2. **Precautionary saving and liquidity heterogeneity.** With persistent risk,
   liquid wealth is partly an insurance buffer. That matters for fertility and
   tenure because the housing block has down-payment and borrowing-floor
   constraints.

3. **Income-risk sorting and incidence.** Housing-price and rent changes can
   have different effects on low-\(z\), middle-\(z\), and high-\(z\)
   households. This creates a natural HANK-style incidence decomposition for
   fertility and spatial sorting.

The branch already shows that the existing collateral wedge is active once
\(z\) is structural. In the borrowing-wedge diagnostic, high-\(z\) renters buy
at a share of `0.091` versus `0.015` for low-\(z\) renters; starter down-payment
feasibility is `0.590` versus `0.336`; low-\(z\) owners are more often near the
borrowing floor (`0.089` versus `0.029`). The `Nz=5` figure packet also shows
strong monotone sorting in the uncalibrated prototype: ownership rises from
about `0.218` at the lowest \(z\) to `0.720` at the highest \(z\), and center
share rises from about `0.224` to `0.523`.

That is the useful part of the branch. It is not proof that the branch is
calibrated.

## Current Calibration Evidence

The evidence is mixed:

- **Encouraging mechanically.** The V5 `Nz=7` validation run accepted strict GE
  with final GE error `0.000222`, \(S=1\), \(q^E=0.9\), positive \(M\), and
  the benchmark-normalized outside-option closure. The `Nz=5` exploratory grid
  solves fast enough for local search.

- **Encouraging locally.** The directional probe
  `alpha_high_fertility_mid_finance` moved loss from `311.30` to `166.72` and
  improved 16 of 19 target gaps. The model is not mechanically stuck.

- **Discouraging as a calibration.** The global search's best formal candidate
  is not economically admissible. Ownership collapses to `0.213`, old-age
  ownership to `0.297`, and the center/periphery rent ratio falls to `0.854`.
  Among the formal top 12 candidates, ownership is generally around
  `0.11`-`0.21` and old-age ownership around `0.17`-`0.30`.

- **Discouraging under ownership screens.** Filtering to penalty-free accepted
  candidates with `own_rate > 0.4` and `old_age_own_rate > 0.5`, the best
  objective rises to `78.96`. The best filtered candidate has ownership levels
  but a wrong-signed ownership gradient (`-0.151`), late first birth
  (`33.11`), high childlessness (`0.216`), and high wealth (`0.993`). No
  candidate in the run passed a loose combined story screen requiring positive
  fertility and ownership gradients, acceptable ownership levels, first-birth
  age below 33, nonextreme childlessness, and nonpathological geography.

- **Objective design is too permissive.** The current objective is a relative
  SMM loss. This treats a five-to-eight-year first-birth age miss as a modest
  relative deviation, while allowing large raw errors in age and rooms to sit
  beside severe ownership failures. This is useful for exploratory search but
  not adequate as the final acceptance rule.

My read: the current evidence is **mechanically encouraging, economically
inconclusive, and calibration-discouraging**. It justifies diagnostics; it does
not justify promoting the branch.

## Recommended Next Step

Do **targeted cross-sectional diagnostics first**. This is option (b), followed
by option (a) only if the diagnostics are promising.

Priority order:

1. **Build a candidate comparison packet** for:
   - live renewal-valve benchmark,
   - V5 `Nz=5` baseline,
   - directional best eval `2`,
   - formal best eval `956`,
   - filtered ownership-floor candidates `554`, `612`, `351`, `208`, and
     `803`.

2. **For each candidate, report cross-sectional objects by \(z\), age,
   wealth, tenure, parent status, and location.** The key question is whether
   \(z\) produces interpretable gradients conditional on wealth and age.

3. **Only then run a constrained/reweighted calibration.** Add hard screens
   for ownership, old-age ownership, \(p_C/p_P>1\), sign of ownership gradient,
   strict GE convergence, finite outside-option accounting, and perhaps a
   maximum acceptable first-birth age. Use stronger penalties for age in years
   and room gaps rather than relying only on relative deviations.

4. **Do not add structural \(\mu\) yet.** A mortgage-account state should be
   added only if the \(z\)-diagnostics show that the income-risk channel is
   promising but the missing account history - amortization, locked-in debt,
   delinquency/default, or refinancing - is the clear reason for failure.
   Adding \(\mu\) now would expand the state space before the one-state branch
   has passed the research test.

5. **Do not change the closure now.** The benchmark-normalized outside-option
   closure is the right paper prototype. The problem is not the closure
   equation; it is whether this state extension can produce admissible
   cross-sectional economics under a disciplined objective.

6. **Do not demote yet, but set a kill criterion.** If the diagnostic packet
   shows that \(z\) only improves aggregate loss by destroying ownership or
   spatial sorting, demote the branch to an appendix/diagnostic result:
   "income risk is mechanically implementable and active, but does not improve
   the benchmark mechanism under admissible tenure/geography constraints."

## Figures And Tables That Would Convince Us

The branch is worth keeping if the following packet looks coherent and improves
on the renewal-valve benchmark.

1. **Aggregate scorecard with guardrails.** A table comparing targets, live
   renewal-valve benchmark, V5 baseline, formal best, ownership-floor bests,
   and any constrained-calibration best. Include target values next to model
   values.

2. **Income-risk mechanism table.** By \(z\) quintile and wealth tercile:
   mass, income, liquid wealth/income, TFR, childlessness, mean age at first
   birth, ownership, renter-to-owner transition, rooms, center share, migration,
   down-payment feasibility, and owner mass near borrowing floor.

3. **Lifecycle profiles by \(z\).** Age profiles for fertility hazards,
   ownership, rooms, center share, and liquid wealth. This should reveal
   whether low-\(z\) households delay fertility because of liquidity/risk or
   whether the model just pushes all births too late.

4. **Tenure and housing ladder heatmaps.** Ownership and owner-rung mass by
   age, \(z\), parent status, and location. The branch should not concentrate
   ownership only in implausible late-life states or oversized rungs.

5. **Spatial sorting heatmaps.** Center share by \(z\), parity, tenure, and age.
   In a Rosen-Roback housing-cost mechanism, high-income households should have
   greater ability to absorb center prices, while parents should face stronger
   space-cost pressure to sort outward.

6. **Birth-response event-study heterogeneity.** Housing response to first and
   second births by \(z\) or income/wealth rank. The model needs to show whether
   the space-cost channel differs across liquidity states.

7. **Counterfactual heterogeneity table.** For a center housing-cost shock or
   broad price/rent shock, report \(\Delta\)TFR, \(\Delta\)childlessness,
   \(\Delta\)age at first birth, \(\Delta\)ownership, \(\Delta\)rooms,
   \(\Delta\)center share, and welfare by \(z\)/wealth group. This is the
   clearest test of whether HANK-\(z\) contributes to the paper rather than
   only to fit.

8. **Grid validation table.** Rerun the accepted constrained candidate at
   `Nz=7` and compare all core moments and key cross-sectional slopes to the
   `Nz=5` result. A branch that changes qualitative economics across grids is
   not ready.

## Decision Rule

Continue this branch if targeted diagnostics show:

- positive and interpretable \(z\)-gradients in fertility/tenure/housing,
- admissible ownership and old-age ownership under constrained calibration,
- center prices above periphery prices,
- meaningful counterfactual heterogeneity relative to the renewal-valve
  benchmark,
- and stable results at `Nz=7`.

Demote it if:

- the best admissible candidates still require wrong-signed fertility or
  ownership gradients,
- the model matches fertility only by collapsing ownership or geography,
- \(z\)-heterogeneity disappears after conditioning on wealth and age,
- or counterfactual effects are not materially different from the simpler
  benchmark.

The current evidence does not warrant another blind global search. It warrants
a cross-sectional mechanism audit.
