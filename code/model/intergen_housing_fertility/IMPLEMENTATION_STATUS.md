# Implementation Status: Intergenerational Housing Fertility

Updated: 2026-06-04

## Mandate

Implement a quantitative OLG model of intergenerational housing mismatch and
fertility without geography. The code lives in this separate folder so the
current calibrated center/periphery implementation remains intact.

The first-pass mechanism is:

\[
\text{old retention / lock-in} \rightarrow \text{scarce family homes}
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

### State Variables

- `INTENDED`: age \(a\).
- `INTENDED`: idiosyncratic income/productivity \(z\).
- `INTENDED`: liquid assets \(b\).
- `INTENDED`: tenure \(o\in\{R,O\}\).
- `INTENDED`: completed children \(n\).
- `SIMPLIFICATION`: there is one scarce family-home asset \(h_O\), not a
  housing-size ladder.
- `SIMPLIFICATION`: renter housing \(h_R\) is an outside rental service. It is
  not part of the scarce family-home stock.
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

- `INTENDED`: one national family-home asset with flow user cost \(q\) and
  asset price
  \[
  P=\frac{q}{r+\delta+\tau^p}.
  \]
- `INTENDED`: owner purchase is constrained by a down-payment condition,
  \[
  (1-\phi)Ph_O \le b,
  \]
  where \(\phi\) is the financed share.
- `INTENDED`: owner purchase also has a flow affordability screen,
  \[
  qh_O \le \psi y(a,z).
  \]
- `SIMPLIFICATION`: there is no mortgage balance, coupon, amortization,
  refinancing, or mortgage-duration state in the first pass.
- `SIMPLIFICATION`: incumbent owners who keep the same tenure do not requalify
  for the down-payment and flow affordability screens each period.
- `SIMPLIFICATION`: renters consume the outside service \(h_R\) at exogenous
  user cost \(R^R\). Since renters do not use the scarce family-home stock,
  there is no rental-market clearing condition in the first pass.
- `NOT IMPLEMENTED`: if future code lets renters occupy the scarce family-home
  stock, then the model must add an explicit rental-stock clearing block. That
  block is deliberately absent from the current Coven-style scaffold.

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

- `INTENDED`: fertility responds to the difference between \(h_R\) and \(h_O\)
  and to the finance constraints that restrict access to \(h_O\).
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

- `INTENDED`: solve one scalar family-home market-clearing condition,
  \[
  H^O(q)=H^S(q).
  \]
- `SIMPLIFICATION`: first pass uses a static upward-sloping family-home supply,
  \[
  H^S(q)=\bar H(q/q_0)^\eta.
  \]
  This is not a construction or transition block.
- `SIMPLIFICATION`: the default stock normalization \(\bar H=1.34\) is
  chosen to avoid placing the scalar market exactly on a coarse-grid tenure
  threshold. It is not a calibrated housing stock.
- `INTENDED`: owner demand and supply are measured in service units per
  normalized adult in the lifecycle cross-section.
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
2. Keep one scarce family-home asset \(h_O\) and one renter outside option
   \(h_R\).
3. Use property-tax capitalization in \(P=q/(r+\delta+\tau^p)\).
4. Use down-payment and flow affordability constraints, with no amortized
   mortgage object.
5. Use a reduced-form old-retention wedge.
6. Solve fixed-price smoke tests before any price iteration.
7. Solve one scalar family-home market-clearing condition.
8. Produce diagnostics for ownership by age, fertility by age, scarce-stock
   clearing, tenure services, and prices.

## Open Decisions

- Whether to move from one-shot completed fertility to sequential hazards.
  Current recommendation: keep one-shot for the first runnable version.
- Whether to make \(h_R\) endogenous. Current recommendation: keep it exogenous
  until the single-home ownership channel is stable.
- Whether to add bequests before calibration. Current recommendation: defer
  until the Coven-style baseline solves robustly.
