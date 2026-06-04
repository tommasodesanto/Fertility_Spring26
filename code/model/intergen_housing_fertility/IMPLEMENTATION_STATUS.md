# Implementation Status: Intergenerational Housing Fertility

Updated: 2026-06-04

## Mandate

Implement a quantitative OLG model of intergenerational housing mismatch and
fertility without geography. The code should start from the existing Python
model conventions where useful, but it must live in this separate folder so the
current calibrated center/periphery implementation remains intact.

The model mechanism is:

\[
\text{old retention / lock-in} \rightarrow \text{scarce family housing}
\rightarrow \text{young finance constraint} \rightarrow \text{lower fertility}.
\]

The main imported mechanisms are:

- Coven et al.: property-tax capitalization and financial constraints. A higher
  recurring tax lowers the asset price \(P_k\) for a given user cost and can
  relax down-payment constraints for high-income, low-liquid-wealth young
  households.
- Fonseca, Liu, and Mabille: mortgage lock-in as a moving or downsizing wedge
  for incumbent owners.

## No Hidden Simplifications Rule

Every simplification must be listed in this file before or in the same commit
as the code that uses it.

Allowed labels:

- `INTENDED`: part of the target model.
- `SIMPLIFICATION`: deliberately simpler than the target model.
- `NOT IMPLEMENTED`: target object not yet coded.
- `DIAGNOSTIC ONLY`: computed for inspection, not part of equilibrium or
  calibration.

Do not leave a shortcut in code with only a local comment. Put it here as well.

## Current First-Pass Target

The first pass is a runnable model, not a calibrated model. It should solve and
produce diagnostics for a stripped-down version of the full blueprint.

### State Variables

- `INTENDED`: age \(a\).
- `INTENDED`: idiosyncratic income/productivity \(z\).
- `INTENDED`: liquid assets \(b\).
- `INTENDED`: tenure \(o\in\{R,O\}\).
- `INTENDED`: occupied housing size \(k\in\{S,M,L\}\) or a small size ladder.
- `INTENDED`: fertility/children state.
- `SIMPLIFICATION`: first pass may retain the existing one-shot completed
  fertility architecture from `dt_cp_model` rather than implementing sequential
  parity hazards immediately.
- `INTENDED`: first pass uses extreme-value smoothing over discrete
  tenure/fertility alternatives. This keeps aggregate owner demand continuous
  enough for market clearing and matches the existing model's use of smoothed
  discrete choices.
- `NOT IMPLEMENTED`: mortgage coupon, remaining maturity, and debt duration
  states \((\iota,m^d)\).
- `NOT IMPLEMENTED`: inheritance receipt state \(B\) and estate-transmission
  kernel.

### Income

- `INTENDED`: lifecycle income with heterogeneity,
  \[
  y(a,z)=W e_a z.
  \]
- `INTENDED`: payment-to-income constraint uses current income:
  \[
  \mathcal P(d',R^m,m_0)\le \psi y(a,z).
  \]
- `SIMPLIFICATION`: first pass can use a small finite \(z\)-grid or even a
  deterministic \(z=1\) smoke mode, but the production solver should be written
  so income heterogeneity is active.

### Housing and Tenure

- `INTENDED`: no location choice and no location-specific prices.
- `INTENDED`: national housing sizes \(k\), with size-specific owner prices
  \(P_k\).
- `INTENDED`: renters and owners choose from the same size ladder in the
  current code, with renter tenure indices followed by owner tenure indices.
- `INTENDED`: owner purchase constrained by
  \[
  d'\le \phi P_k h_k,
  \qquad
  \mathcal P(d',R^m,m_0)\le \psi y(a,z).
  \]
- `INTENDED`: property-tax capitalization affects owner asset prices through
  the user-cost relation.
- `SIMPLIFICATION`: the first pass uses the compact draft's
  collateral-constrained user-cost representation. Owner choices pay a flow
  user cost, while asset prices enter the down-payment and payment-to-income
  constraints. The code does not yet track mortgage balances, home equity, or
  amortization.
- `SIMPLIFICATION`: LTV and payment-to-income constraints apply when a renter
  buys or an owner changes size. Incumbent owners who keep the same size are
  not forced to requalify each period.
- `INTENDED`: first pass includes a competitive landlord residual by size:
  landlords hold the non-owner-occupied stock and rent it to tenants.
- `SIMPLIFICATION`: landlord zero profit is imposed mechanically by tying
  rental user costs to owner user costs,
  \[
  R_k=(1+\upsilon^R)q_k+\iota^R,
  \]
  rather than solving a separate landlord portfolio problem.
- `SIMPLIFICATION`: the default rental markup \(\upsilon^R\) is a first-pass
  ownership-margin device. It is not calibrated and should be replaced by a
  rent/landlord-cost target or a richer tenure block before quantitative
  claims.
- `SIMPLIFICATION`: landlords do not have balance sheets, taxes, default risk,
  or financing constraints in the first pass.
- `SIMPLIFICATION`: first pass uses a static upward-sloping owner supply curve
  by size,
  \[
  H_k^S(q_k)=\bar H_k(q_k/q_{k0})^{\eta_k}.
  \]
  This is not a transition construction block and does not yet track durable
  investment \(I_{k,t}\).

### Old Retention and Lock-In

- `INTENDED`: old owners may retain housing because of a moving/downsizing
  wedge.
- `SIMPLIFICATION`: first pass uses a reduced-form age-dependent sale or
  downsizing cost.
- `SIMPLIFICATION`: the first-pass old-retention wedge is a private
  resource-equivalent adjustment cost in the household budget. It has no fiscal
  or lender counterpart yet.
- `NOT IMPLEMENTED`: present-value mortgage lock-in from coupon gaps.
- `NOT IMPLEMENTED`: explicit lender-side accounting for coupon lock-in.

### Fertility

- `INTENDED`: fertility responds to housing access and the cost of child-related
  housing needs.
- `SIMPLIFICATION`: first pass may use the existing discrete completed-fertility
  choice architecture to preserve solver continuity.
- `SIMPLIFICATION`: first pass lets fertility be chosen once at a fixed fertile
  age. Children then remain as a completed-fertility state; there is no
  child-age vector or period-by-period birth hazard.
- `NOT IMPLEMENTED`: sequential birth hazards and child-age vector.
- `DIAGNOSTIC ONLY`: parity progression or second-birth hazard moments are not
  equilibrium targets unless the sequential fertility state is implemented.

### Bequests and Estate Tax

- `INTENDED`: eventual model includes estates, inheritance transmission, estate
  taxation, revenue rebates, and bequest-tax counterfactuals.
- `SIMPLIFICATION`: first pass omits bequest-tax policy.
- `NOT IMPLEMENTED`: inheritance kernel \(\Gamma\), estate-tax revenue,
  transfer schedule \(T^b\), and gross bequest principal adding-up in code.

### Equilibrium

- `INTENDED`: solve for national housing prices or user costs that clear total
  physical occupancy by size,
  \[
  H^O_k(q)+H^R_k(q)=H_k^S(q).
  \]
- `SIMPLIFICATION`: current smoke/default GE clears aggregate housing services
  with one common user-cost shifter across size rungs. Size-specific excess
  demand is diagnostic unless the `by-size` clearing mode is requested.
  This avoids over-interpreting a coarse three-rung menu as three fully
  separate physical submarkets before the richer supply block is implemented.
- `INTENDED`: first-pass owner supply and demand are measured in service units
  per normalized adult in the lifecycle cross-section; owner demand is
  normalized by total lifecycle mass before comparing to supply.
- `INTENDED`: first-pass owner, renter, landlord, and total stock quantities
  are measured in service units per normalized adult in the lifecycle
  cross-section.
- `SIMPLIFICATION`: rental prices are tied mechanically to owner user costs in
  the first pass.
- `SIMPLIFICATION`: the lifecycle distribution is normalized by entrant mass
  and does not yet feed fertility choices back into cohort size or entry.
- `INTENDED`: price iteration reports the best excess-demand metric even when
  it converges, so failures are visible rather than hidden.
- `INTENDED`: first-pass GE uses dependency-free price search routines. The
  default is aggregate bisection over a common owner user-cost shifter;
  by-size coordinate search remains available for diagnostics.
- `INTENDED`: the by-size coordinate search minimizes the sum of squared
  relative excess demands and reports the maximum relative excess demand
  separately. This is a numerical device for the discontinuous discrete
  housing-size menu, not a calibrated error structure.
- `DIAGNOSTIC ONLY`: in the current smoke model, short by-size searches improve
  but do not reliably clear all size markets. The default aggregate-clearing
  solution is the maintained smoke equilibrium until the supply and tenure
  blocks are refined.
- `DIAGNOSTIC ONLY`: with renters and owners on the same size ladder, smoke
  solves are slower than the earlier owner-only scaffold. Treat long by-size
  price searches as diagnostics, not routine smoke tests.
- `INTENDED`: full size-specific rental/owner stock decomposition
  \(H_k^O+H_k^R=H_k\) is now represented in the static equilibrium accounting.
- `NOT IMPLEMENTED`: transition dynamics.

### Calibration

- `NOT IMPLEMENTED`: calibration, SMM objective, and counterfactual tables.
- `DIAGNOSTIC ONLY`: first pass should report moments next to placeholders for
  future targets, but these are not calibrated estimates.

## Initial Coding Plan

1. Create a minimal parameter module with no location objects.
2. Create grids for age, assets, income states, fertility states, tenure, and
   housing sizes.
3. Port the existing Bellman/distribution pattern only where it is structurally
   compatible.
4. Add Coven-style owner finance constraints: LTV and payment-to-income.
5. Add a reduced-form old-retention wedge.
6. Solve a fixed-price smoke test before any general-equilibrium price update.
7. Add owner housing-market clearing by size.
8. Produce diagnostics: ownership by age/income, housing size by age/children,
   fertility by tenure/housing access, old retention by age, and excess demand
   by size.

## Open Decisions

- Whether first-pass fertility should keep one-shot completed fertility or move
  immediately to sequential hazards. Current recommendation: keep one-shot for
  the first runnable version and mark it as a simplification.
- Whether renters choose only a fixed rental option or a capped continuous
  rental amount. Current recommendation: choose the smaller code change that
  preserves the current model's renter logic.
- Whether first-pass owner price clearing should use fixed supply by size or a
  simple upward-sloping supply curve. Current implementation: simple
  upward-sloping supply by size.
