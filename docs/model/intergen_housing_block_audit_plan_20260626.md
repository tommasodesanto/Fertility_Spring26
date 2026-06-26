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

