# Intergen Housing-Block Audit Plan

Date: 2026-06-26

This note is a plan for a housing-block audit before further fertility
calibration. The current symptoms are not primarily fertility symptoms. They
are coming from wealth grids, tenure/rung discreteness, purchase constraints,
owner retention, and housing/accounting conventions.

## Coven Reference Check

The relevant paper is Coven, Golder, Gupta, and Ndiaye, "Property Taxes and
Housing Allocation Under Financial Constraints," June 19, 2025. The paper is
useful for mechanism discipline, but it does not provide a visual benchmark for
individual consumption or housing policy functions.

What the paper does say:

- The quantitative model includes tenure, housing quantity, savings, and
  consumption choices, with down-payment and payment-to-income constraints.
- Housing markets clear in general equilibrium, with house prices and rents
  tied by a user-cost condition.
- The solution method is backward induction for value functions, forward
  induction for the distribution, and iteration until taxes and rebates are
  stationary.
- The paper validates the model mainly with aggregate lifecycle objects:
  homeownership/renting by age and location, housing quantity by age, total
  wealth by age, migration, and counterfactual changes.

What the paper does not appear to do:

- It does not plot individual consumption policy functions over wealth.
- It does not plot housing policy functions over liquid wealth or total wealth.
- It does not discuss sawtooth policy functions, grid-induced kinks, or
  interpolation artifacts.

Implication: Coven et al. are not evidence that our jagged policy functions are
acceptable. Their mechanism supports the broad down-payment/capitalization
logic, but our current graphs need their own numerical and economic audit.

## Current Housing Symptoms

The solved June 26 packet shows several housing-block issues that should be
treated as first-order until explained:

1. The liquid-wealth grid is poorly occupied. Most stationary mass sits in a
   narrow region near zero, while the grid extends far into states with little
   or no mass.
2. There is a very large stationary mass at or near the zero liquid-wealth
   node. This can be economic, but it is also exactly where purchase and
   borrowing constraints create branch-switching cliffs.
3. The mass-weighted average consumption function by liquid wealth remains
   visibly wiggly, even after integrating over age, income, tenure, parity, and
   child states.
4. Replotting against total wealth using the live convention
   \(W=b+(1-\psi)pH\) for owners does not remove the jaggedness. Some of the
   total-wealth sawtooth is mechanical because owner equity shifts interleave
   several discrete owner grids, but the object is still too hard to read.
5. The owner housing menu is coarse: \(H^{O}\in\{2,4,6,8,10\}\), while renter
   housing is continuous up to \(h_R^{\max}=6\). This mixes one continuous
   choice with a coarse discrete product ladder.
6. The dormant \(H=2\) owner rung has near-zero residual service and is floored
   rather than cleanly masked. Even if mass is tiny, this is a bad numerical
   precedent.
7. Old ownership remains near absorbing, while young ownership is still too
   low. This is a housing lifecycle failure, not a fertility failure.
8. Total-wealth conventions have been inconsistent across tools. Live
   statistics use home equity \(b+(1-\psi)pH\), while at least one older helper
   used accounting net worth \(b+pH\).

## Audit Goal

The audit should answer one question:

> Does the model have a coherent, numerically stable housing block before
> fertility incentives are allowed to matter?

The answer should be based on plots, invariants, and fixed-price experiments,
not on a new broad SMM calibration.

## Phase 1: Minimal Reproducible Visual Packet

Build one command that regenerates the same housing-block packet from any saved
solution cache or source record.

Required plots:

1. \(E[c\mid b]\): mass-weighted average consumption by liquid wealth, one line.
2. \(E[c\mid W]\): mass-weighted average consumption by total wealth, one line,
   with \(W=b+(1-\psi)pH\) for owners.
3. Stationary mass by liquid wealth.
4. Stationary mass by total wealth.
5. Ownership probability by liquid wealth and by total wealth.
6. Average next liquid wealth \(E[b'\mid b]\) and \(45^\circ\) line.
7. Housing services by liquid wealth and total wealth, realized after tenure.
8. Owner-rung shares by age.
9. Renter cap share by age.
10. Old-owner same-rung retention, downsizing, upsizing, and sell-to-rent rates.

Deliverable:

- One script under `code/model/tools/`.
- One output folder under the diagnostic packet.
- One markdown readout with the exact wealth definitions.

## Phase 2: Grid Audit

This phase determines whether the observed wiggles are grid-amplified.

Diagnostics:

1. Re-solve the best point at fixed price with `Nb` in `{60, 120, 240}`.
2. Run the same packet under the default grid and under dense-core grids
   centered on the occupied support.
3. Compare \(E[c\mid b]\), \(E[c\mid W]\), \(E[b'\mid b]\), tenure entry
   probabilities, and owner-rung choices.
4. Report grid occupancy: mass by wealth decile, mass at zero, mass below zero,
   and mass outside each proposed core region.
5. Compare `linear` and `monotone_cubic` interpolation at fixed price, but do
   not treat the interpolation swap as a final fix unless moments are stable.

Decision rule:

- If wiggles shrink materially with finer or denser grids, redesign the wealth
  grid before calibration.
- If wiggles survive, treat them as discrete housing-choice economics and move
  to the product-ladder audit.

## Phase 3: Housing Product-Ladder Audit

This phase tests whether the owner menu is creating artificial thresholds.

Diagnostics:

1. Sweep owner menus at fixed price:
   \[
   \{2,4,6,8,10\},\quad
   \{3,4,5,6,8,10\},\quad
   \{4,5,6,7,8,10\},\quad
   \{4,5,6,7,8,9,10\}.
   \]
2. Compare the same visual packet plus owner-rung shares.
3. Remove or explicitly mask the \(H=2\) rung and check moment changes.
4. Add a temporary continuous-owner-housing diagnostic at fixed price, if
   feasible, to see whether the discrete owner menu is the source of the
   consumption sawtooth.

Decision rule:

- If the coarse menu is the source, move to a better owner product menu or a
  nested discrete-continuous owner problem.
- If the menu is not the source, focus on tenure accounting and constraints.

## Phase 4: Tenure And Purchase Accounting Audit

This phase checks whether owner entry and exit are being represented in the
right state variable.

Diagnostics:

1. For each high-mass state near \(b=0\), tabulate the value gap between
   renting and each owner rung.
2. Compute the exact cash need for purchase:
   \[
   (1-\phi)pH.
   \]
3. Compute the branch liquid resources after each tenure transition.
4. Plot the branch value functions around purchase thresholds.
5. Check whether consumption dips occur at the same points where next-period
   tenure or owner rung changes.
6. Compare total-wealth definitions:
   - live net-equity wealth \(b+(1-\psi)pH\);
   - accounting gross-housing wealth \(b+pH\);
   - liquid wealth \(b\).

Decision rule:

- If the model is effectively using liquid debt as the main state while housing
  equity is too discrete, consider adding a cleaner housing-equity or mortgage
  accounting state.
- If the accounting is right but thresholds dominate, the model may need a
  smoother owner product market or a starter-owner product.

## Phase 5: Old Ownership And Downsizing Audit

This phase isolates the old-owner absorbing-margin problem.

Diagnostics:

1. Decompose old owner transitions into same rung, downsize, upsize, and sell.
2. Repeat after setting sale costs to zero at fixed price.
3. Repeat after setting owner service premium `chi` to zero at fixed price.
4. Repeat after changing terminal bequests from gross housing wealth to net
   equity wealth.
5. Add a temporary old-age maintenance or downsizing cost/benefit diagnostic
   and check whether old ownership moves without breaking young ownership.

Potential fixes:

- Explicit old-age moving/downsizing shock.
- Maintenance or property-tax burden that rises with owner size.
- Net-equity bequest convention.
- Old-age liquidity or health expense shock.
- Separate old-owner retention wedge from young-owner access wedge.

## Phase 6: Renter Cap And Rental Supply Audit

This phase tests whether the renter cap is doing too much work.

Diagnostics:

1. Sweep \(h_R^{\max}\in\{5,6,7,8,10,\infty\}\) at fixed price.
2. Compare renter mean rooms, renter share rooms \(\ge 6\), owner-renter room
   gap, and housing increments.
3. Replace the hard renter cap with a rental size price wedge in a diagnostic
   solve, if feasible.

Potential fixes:

- Keep the cap only if it is a measurement/product-support object.
- Otherwise use a rental-size price wedge or segmented rental supply curve.

## Phase 7: Fertility-Off Housing Test

This phase prevents fertility from hiding housing bugs.

Diagnostics:

1. Freeze fertility choices or force parity to zero.
2. Run the housing packet under the same theta and fixed price.
3. Check whether consumption, housing, ownership, and grid occupancy still show
   the same kinks.

Decision rule:

- If the kinks remain, they are housing/block/grid issues.
- If the kinks disappear, fertility/family-space thresholds are interacting
  with housing in a way that needs separate study.

## Phase 8: Acceptance Criteria

Do not return to broad fertility calibration until these conditions hold:

1. The wealth grid has at least 95 percent of mass in a region with dense
   support, and no major policy jumps are driven by empty support.
2. The definition of total wealth is explicit and consistent across statistics
   and plots.
3. \(E[c\mid b]\) and \(E[c\mid W]\) are explainable from tenure/rung
   thresholds or become smoother under grid/product changes.
4. The owner ladder has no utility-dead rung used only because of numerical
   floors.
5. Old owner retention is decomposed and has at least one tested mechanism that
   can move it independently of young owner access.
6. Young owner entry is decomposed into feasibility and value-ranking margins.
7. Renter size support is an economic object, not an accidental cap used to hit
   room moments.
8. The housing packet can be regenerated by one command.

## Immediate Next Step

Start with a fixed-price grid and product-ladder audit on the June 26 best
packet. This is the fastest way to answer whether the weird graphs are mostly:

1. wealth-grid coarseness;
2. discrete owner-product coarseness;
3. tenure-transition accounting;
4. old-owner retention;
5. or a real but ugly economic threshold-saving object.

## Current Issue Ledger: 2026-06-27

This ledger separates three objects that should not be mixed:

1. **Numerical pathology:** a feature created by grid support, interpolation,
   masks, or KFE injection.
2. **Economic miss:** a model mechanism that is coherent but empirically wrong.
3. **Calibration fit:** the weighted SMM loss under the current target system.

For the next audit cycle, judge each issue by its own diagnostic object. Do not
use SMM loss as the primary criterion unless the issue is explicitly a
calibration-fit issue.

| Priority | Issue | Pathology Or Miss | Current Evidence | Correct Validation Object | Candidate Fix |
|---|---|---|---|---|---|
| `1` | Entrant liquid-wealth atom | Mostly closure/numerical support, with an economic point-mass assumption underneath. | The old point-entry run had a large atom near \(b_{entry}=0.1465\). The live calibration base now uses an external PSID income-ratio distribution with five weighted quintile-bin means \([-2.519,-0.079,0.102,0.353,3.040]\), mapped to model wealth by entrant income and scattered linearly over grid nodes. | Model-vs-data liquid-wealth distributions by age, tenure, and parent status; atom-origin heatmap by age/income/tenure/parity; not the SMM loss. | Treat point-entry as a legacy switch. Validate the empirical-entry \(G_0(b\mid z)\) against PSID/SCF/SIPP liquid-wealth distributions before recalibration. |
| `1b` | Young wealth target denominator | Measurement/statistic mismatch, fixed for the active target key. | The old model statistic `young_liquid_wealth_to_income` divides liquid wealth by 4-year period after-tax income, while the PSID candidate targets are annual-income ratios. The active 14-moment target now uses `young_childless_renter_liquid_wealth_to_annual_gross_income_2535 = 0.17922556`; historical old-key losses are stale. | Annual-gross model moments compared to PSID wealth/income ratios by age, tenure, and child status; distribution plots, especially medians and near-zero mass. | Keep the active key; validate the whole model/data young-wealth distribution, not only the mean. |
| `2` | Total-wealth jaggedness | Partly mechanical from discrete owner rungs; possibly amplified by coarse owner menu/grid. | Total wealth \(W=b+(1-\psi)pH\) remains sawtoothed even after entry spreading because each owner rung shifts the liquid grid by a different equity amount. | Total-wealth densities conditional on owner rung \(H\), then smoothed aggregate; rung transition matrices by age/parity. | Do not smooth the plot alone. Test denser owner ladders or adjacent-rung lotteries only after conditional-rung diagnostics show artificial bunching. |
| `3` | Poor grid occupancy and threshold sensitivity | Numerical support problem. | Much of the liquid grid is unused, while mass and purchase thresholds sit in the near-zero/dense region. Prior `Nb=60 -> 240` diagnostics moved young ownership and other moments materially. | Fixed-price nested-grid convergence for densities, affected mass, \(E[c\mid b]\), \(E[b'\mid b]\), and owner-entry probabilities. | Redesign the grid around occupied support and pin economic thresholds: \(0\), entry support, \((1-\phi)pH\), relaxed thresholds, and \(-\phi pH\). |
| `4` | Young ownership too low and old ownership too high | Economic lifecycle miss, not primarily the entry atom. | Entry spreading raises young ownership, but old ownership remains near absorbing. Earlier grid diagnostics showed old ownership barely moves with `Nb`. | Young owner feasibility/value-gap tables; old-owner transition decomposition into same rung, downsize, upsize, sell-to-rent. | Separate young access from old retention: starter-owner mapping, rental premium/supply, old-age maintenance/downsizing/health/liquidity mechanism, or net-equity bequest convention. |
| `5` | Tenure/product segmentation is hand-built | Economic measurement/model-object miss. | Current hard renter cap \(h_R^{max}=6\) and owner ladder \(\{2,4,6,8,10\}\) are convenient but not yet disciplined by the empirical joint distribution. | ACS/AHS tables for \(F_R(h)\), \(F_O(h)\), \(\Pr(O\mid h)\), rent gradients by size, owner value gradients by size, and old-owner retained-size distributions. | Replace pure hand segmentation with either empirically chosen owner rungs plus overlapping renter support, or a soft rental-size premium/availability wedge. |
| `6` | Parent-credit counterfactuals have tiny fertility effects | Mechanism/affected-mass issue, not a calibration-loss issue. | Model affected-mass tables suggested raw threshold mass can exist while birth-weighted and owner-entry-weighted mass is tiny. | Sufficient-stat table: raw, birth-weighted, space-constrained, owner-entry-weighted, and \(\Psi\)-weighted affected mass; PSID liquid assets plus ACS/AHS housing thresholds. | Do not infer policy failure from \(\Delta\)TFR alone. First measure whether the empirical susceptible mass exists and whether the model places parent-margin households there. |
| `7` | Housing event responses are weak or unstable | Economic transition miss, possibly product/tenure mapping. | The one-child-to-two-plus housing increment remains small in current runs; entry spreading helped it slightly but did not solve it. | Birth-cohort transition tables: pre-birth tenure/rung, post-birth tenure/rung, housing services change, and wealth state, weighted by birth probability. | Fix only after tenure/product mapping is disciplined. Candidate mechanisms include better family-sized rental scarcity, starter-owner/family-owner ladder, or child-dependent moving/adjustment costs. |
| `8` | Infeasible cliffs and masks | Numerical hygiene issue. | Diagnostics found small positive mass near \(-10^{10}\) continuation interpolation, shrinking with finer grids. Owner residual service floors remain a bad precedent. | Occupied-branch mask audit and Bellman residuals weighted by stationary mass. | Replace penalty values/floors with explicit feasibility masks where feasible; verify moment invariance before using as a production change. |

### Recommended Order

1. Validate the wealth denominator and entry-wealth distribution against data:
   model and data liquid wealth by age and tenure, including near-zero and
   negative mass, with annual-gross income denominators reported explicitly.
2. Run the atom-origin heatmap in the model: \(b\) bins by age, faceted by
   tenure/parity/income, comparing the legacy point-entry switch to the
   external PSID entry distribution.
3. Diagnose total wealth by owner rung before changing owner products.
4. Run fixed-price nested-grid diagnostics with threshold-pinned grids.
5. Only then choose the economic redesign: empirical entry distribution,
   tenure/product mapping, and old-owner retention mechanism.

The external five-node entry distribution should therefore be recorded as
**progress on Issue 1 only**. It does not solve Issues 2--8, and it should not
be judged primarily by its current calibration loss.
