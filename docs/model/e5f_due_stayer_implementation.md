# Native DUE-style existing-owner credit: conditional implementation

Author conditionally adopted this rule September 27, subject to a matched
numerical check. This isolated implementation is not yet the calibrated baseline.

Let b be beginning liquid wealth (negative is debt), L=-phi*p*H the current
collateral floor, and b_next end-period wealth. Existing owners keeping the
same house face b_next >= min(b,L). New buyers retain the existing origination
screen and b_next >= L. The budget remains R*b + income - costs = c + b_next:
the inherited balance is b, not R*b, so interest cannot be capitalized into
additional debt above the collateral ceiling. DUE reference equations 2.2–2.3:
`tmp/september_slides_review/greaney_reference.txt`, lines 440–460.

At any age with positive death probability (including certain terminal death),
a separate floor b_next >= -(1-selling_cost)*p*H enforces the author's
nonnegative-estate requirement. This is an explicit solvency extension: a large
price fall can require repayment for estate solvency even though LTV restoration
is not required. No default, insurance, amortization or payment-to-income rule
is added. All prices here are the current decision-date prices, preserving the
existing estate valuation timing.

The flag `native_due_stayer_credit` defaults off. It requires reviewed native
purchase-income accounting, one market, postdecision current distribution and
existing supported native modes. Natural-credit and legacy stayer switches are
incompatible. The legacy origination-only switch is NOT DUE: it prevents all
cash-out for indebted stayers even when below the current collateral ceiling.

## Payload and routing

- Saving-stage buyer policies remain `bp_pol`, `c_pol`.
- Stayer policies are `P._bp_pol_stay`, `P._c_pol_stay`, packed publicly as
  `solution.bp_pol_stay`, `solution.c_pol_stay`.
- Fast retained solutions carry both private fields; upgrade restores both
  before KFE/reporting. Existing dictionary/pickle serialization carries these
  fields; external serializers and dated policy adapters require explicit wiring.
- `realize_stayer_cross_section(g, loc_probs, tenure_choice, tenure_probs)`
  returns seven-axis mass after location/tenure choice and before saving for
  same-location owners choosing the same house. Renters are excluded.
- `solution.g_stay_distribution` is that mass. Other current mass is
  `solution.g - solution.g_stay_distribution`. Aggregate bequest flow uses the
  corresponding saving policy for each branch.
- Existing stationary/one-period forward stayer scatter is reused. A dated
  evaluator must bind the date's policies, not stale P fields left by backward
  induction. Companion accounting/transition changes are a separate work scope.

## Validation completed / pending

Eight tiny synthetic tests cover the floor, interest timing, death-solvency
extension, default-off behavior, owner-kernel branches, exact budgets, stayer
origin mass, aggregate estate mixtures, packed payload retention and fail-closed
rejection of the unsupported nonkernel saving fallback. All eight passed JIT on/off. Eight existing native-household fixtures also pass.
No full household or equilibrium solve was launched by this implementation.

Before promotion: lead reviews mathematical diff, default-off authenticated
native replay, matched same-price full targets/parameters/policy arrays, and a
short falling-price matched path with origin-specific saving/consumption/estate
checks. Full-grid equality with flag on is not expected: formerly infeasible
inherited-debt nodes become feasible, and grid interpolation can populate nodes
below the continuous collateral boundary even at stationary prices. Quantify
those effects without changing grids or silently treating them as economics.

Companion reporting prerequisite: the current transition runtime now requires an
estate-audit capability marker for DUE and rejects the frozen legacy audit for
this mode. A new authenticated audit/source contract is required before a
matched case can run; never relax the old pin or reuse an unmodified frozen
calibration evaluator with the new flag. The accounting companion covers dated
policy packaging, origin-specific budgets and estates, and explicit forward
stayer policies, with its separate test receipt.

## September27 integration status

Same-price comparison and permanent10% lower-price dated household/operator check
pass; all31 original parameters fixed. Main now contains the reviewed optional
core and origin-aware accounting/date routing.45 focused main tests pass.
Existing defaults and frozen calibration contracts remain unchanged; a new
production loader/contract must opt into the approved DUE rule and audits.
This is not a cleared DUE equilibrium or completed frictionless transition.
Full dated evidence: output/model/daytime_calibration_20260927/due_price_fall.
Existing-owner supplemental policies inspected; see the dated packet. They
remain distinct from standard buyer-conditional plot templates.
