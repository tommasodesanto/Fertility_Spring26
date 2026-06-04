# Implementation Status: Intergenerational Housing Fertility

Updated: 2026-06-04

## Mandate

Implement a quantitative OLG model of intergenerational housing mismatch and
fertility without geography. The code lives in this separate folder so the
current calibrated center/periphery implementation remains intact.

The first-pass mechanism is:

\[
\text{old retention / lock-in} \rightarrow \text{scarce housing services}
\rightarrow \text{young finance constraint} \rightarrow \text{lower fertility}.
\]

The target is now a Coven-style scaffold with fertility added. The model should
stay small until calibration actually requires extra state variables.

## No Hidden Simplifications Rule

Every simplification must be listed in this file before or in the same commit
as the code that uses it.

Allowed labels:

- `INTENDED`: part of the current target model.
- `SIMPLIFICATION`: deliberately simpler than a richer target model.
- `NOT IMPLEMENTED`: target object not yet coded.
- `DIAGNOSTIC ONLY`: computed for inspection, not part of equilibrium or
  calibration.

Do not leave a shortcut in code with only a local comment. Put it here as well.

## Current First-Pass Target

The first pass is runnable model code, not a calibrated quantitative result.

## Coven Alignment

This is the implementation map used for the current code. The source model is
Coven, Golder, Gupta, and Ndiaye, "Property Taxes and Housing Allocation Under
Financial Constraints," Section 3.

- `INTENDED`: copy the tenure menu structure conceptually. Households choose
  between renting, buying/adjusting, or keeping owned housing. The current code
  represents this with a single tenure/quantity index.
- `INTENDED`: copy the rental-price logic. Rents are not a separate market
  price. They are pinned by the user-cost relation
  \[
  R=(r+\delta+\tau^p)P=q
  \]
  after setting expected capital gains to zero in the stationary first pass.
- `INTENDED`: copy aggregate housing clearing. Total renter plus owner housing
  services clear against
  \[
  H^S(P)=cP^\eta.
  \]
  There is no separate landlord state and no rental-stock-by-size clearing
  condition in the current target.
- `INTENDED`: copy the owned-housing minimum-size idea. The owner grid starts
  at \(h^O_{\min}>h^R_{\min}\). Renters may choose smaller housing quantities.
- `INTENDED`: copy the homeownership utility term. Coven includes an ownership
  benefit \(\Xi^O\). The no-geography code has one scalar
  `owner_utility_bonus`, normalized relative to renter tenure.
- `INTENDED`: copy the financial-friction channels. New owner choices face
  down-payment and payment-to-income screens.
- `SIMPLIFICATION`: renters and buyers choose from short discrete housing
  quantity grids. Coven's renter and adjuster housing choice is continuous.
- `SIMPLIFICATION`: the code has no purchase-price balance-sheet accounting,
  mortgage debt state, amortization, or capital-gains tax. It uses a
  collateral-constrained user-cost shortcut: \(q h\) enters the flow budget and
  \((1-\phi)P h^O\le b\) enters the ownership feasibility screen.
- `SIMPLIFICATION`: the payment-to-income screen is \(q h^O\le \psi y\), not
  Coven's amortized mortgage-payment plus property-tax condition.
- `NOT IMPLEMENTED`: two-region location choice, migration costs, regional tax
  systems, bequests, estate taxation, capital-gains taxation, and transition
  dynamics.

### State Variables

- `INTENDED`: age \(a\).
- `INTENDED`: idiosyncratic income/productivity \(z\).
- `INTENDED`: liquid assets \(b\).
- `INTENDED`: tenure/quantity \(o\). Indices \(0,\ldots,N_R-1\) are renter
  quantities and indices \(N_R,\ldots,N_R+K-1\) are owner quantities.
- `INTENDED`: completed children \(n\).
- `INTENDED`: renter and owner housing choices use separate quantity grids,
  \(h^R\in\mathcal H^R\) and \(h^O\in\mathcal H^O\), with the owner grid
  starting at a larger minimum size.
- `SIMPLIFICATION`: the grids are small and discrete. Coven's model has a
  continuous housing quantity choice; this first pass approximates it with a
  short quantity ladder.
- `SIMPLIFICATION`: fertility is a one-shot completed-fertility choice at one
  fertile age. There is no child-age vector or sequential birth hazard yet.
- `INTENDED`: first pass uses extreme-value smoothing over discrete tenure and
  fertility alternatives so aggregate demand is continuous enough for price
  clearing.

### Income

- `INTENDED`: lifecycle income with heterogeneity,
  \[
  y(a,z)=W e_a z.
  \]
- `INTENDED`: income enters both the household budget and the flow
  affordability constraint.
- `SIMPLIFICATION`: smoke mode uses a small finite \(z\)-grid.

### Housing and Finance

- `INTENDED`: one aggregate housing-services market with flow user cost \(q\)
  and asset price
  \[
  P=\frac{q}{r+\delta+\tau^p}.
  \]
- `INTENDED`: renters and owners pay the same per-unit user cost \(q\) for
  housing services. This is the no-arbitrage rent/user-cost condition in the
  Coven-style block.
- `INTENDED`: owners receive a utility term \(\Xi^O\), implemented as
  `owner_utility_bonus`. This is not a hidden tuning device: Coven has this
  object, and it is needed because owners face extra constraints and adjustment
  costs.
- `INTENDED`: owner purchase is constrained by a down-payment condition,
  \[
  (1-\phi)Ph^O \le b,
  \]
  where \(\phi\) is the financed share.
- `INTENDED`: owner purchase also has a flow affordability screen,
  \[
  qh^O \le \psi y(a,z).
  \]
- `SIMPLIFICATION`: the household budget uses the flow user cost \(qh\), while
  the down-payment requirement enters only as a liquidity constraint. The code
  does not subtract the down payment from liquid wealth and then carry housing
  equity as a separate state. This is a collateral-constrained user-cost
  shortcut.
- `SIMPLIFICATION`: there is no mortgage balance, coupon, amortization,
  refinancing, or mortgage-duration state in the first pass.
- `SIMPLIFICATION`: incumbent owners who keep the same tenure do not requalify
  for the down-payment and flow affordability screens each period.
- `INTENDED`: rental quantities are part of aggregate housing demand. There is
  no separate rental-stock state because the supply side provides aggregate
  housing services, as in the simple Coven-style supply block.

### Old Retention and Lock-In

- `INTENDED`: old owners may retain housing because of a moving or downsizing
  wedge.
- `SIMPLIFICATION`: first pass uses a reduced-form age-dependent sale or
  downsizing cost.
- `SIMPLIFICATION`: the first-pass old-retention wedge is a private
  resource-equivalent adjustment cost in the household budget. It has no fiscal
  or lender counterpart yet.
- `NOT IMPLEMENTED`: present-value mortgage lock-in from coupon gaps.
- `NOT IMPLEMENTED`: lender-side accounting for coupon lock-in.

### Fertility

- `INTENDED`: fertility responds to household housing quantity \(h(o)\) and to
  the finance constraints that restrict access to larger owner quantities.
- `SIMPLIFICATION`: fertility is discrete in code even though the compact theory
  uses continuous \(n\) for clean first-order conditions.
- `DIAGNOSTIC ONLY`: parity progression or second-birth hazard moments are not
  equilibrium targets unless the sequential fertility state is implemented.

### Bequests and Estate Tax

- `INTENDED`: eventual model includes estates, inheritance transmission, estate
  taxation, revenue rebates, and bequest-tax counterfactuals.
- `NOT IMPLEMENTED`: inheritance kernel \(\Gamma\), estate-tax revenue,
  transfer schedule \(T^b\), and gross bequest principal adding-up in code.

### Equilibrium

- `INTENDED`: solve one scalar housing-services market-clearing condition,
  \[
  H^R(q)+H^O(q)=H^S(P(q)).
  \]
- `INTENDED`: first pass uses a static upward-sloping housing supply,
  \[
  H^S(P)=cP^\eta.
  \]
  This follows Coven's simple competitive supply convention. It is not a
  dynamic construction or transition block.
- `SIMPLIFICATION`: the default supply shifter \(c=3.40\) is a smoke-test
  normalization. It is not a calibrated housing stock or supply estimate.
- `INTENDED`: renter demand, owner demand, and supply are measured in service
  units per normalized adult in the lifecycle cross-section.
- `SIMPLIFICATION`: the lifecycle distribution is normalized by entrant mass
  and does not yet feed fertility choices back into cohort size or entry.
- `INTENDED`: price iteration reports the excess-demand metric even when it
  converges.
- `NOT IMPLEMENTED`: transition dynamics.

### Calibration

- `NOT IMPLEMENTED`: calibration, SMM objective, and counterfactual tables.
- `DIAGNOSTIC ONLY`: first pass reports moments for inspection only. They are
  not calibrated estimates.

## Current Coding Plan

1. Keep the state as \((a,z,b,n,o)\).
2. Use renter and owner quantity grids with one aggregate housing-services
   price.
3. Use property-tax capitalization in \(P=q/(r+\delta+\tau^p)\).
4. Use down-payment and flow affordability constraints, with no amortized
   mortgage object.
5. Use a reduced-form old-retention wedge.
6. Solve fixed-price smoke tests before any price iteration.
7. Solve one scalar housing-services market-clearing condition.
8. Produce diagnostics for ownership by age, fertility by age, aggregate
   housing clearing, tenure services, quantity demand, and prices.

## Open Decisions

- Whether to move from one-shot completed fertility to sequential hazards.
  Current recommendation: keep one-shot for the first runnable version.
- Whether to add bequests before calibration. Current recommendation: defer
  until the Coven-style baseline solves robustly.
