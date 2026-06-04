# Implementation Status: Intergenerational Housing Fertility

Updated: 2026-06-04

## Rule

Do not leave simplifications implicit. If code uses a shortcut, record it here
with one of these labels:

- `INTENDED`: part of the current target model.
- `SIMPLIFICATION`: deliberately simpler than a richer target.
- `NOT IMPLEMENTED`: target object not coded yet.
- `DIAGNOSTIC ONLY`: produced for inspection, not a calibration target.

## Current Pass

The current code is a runnable first-pass scaffold, not a calibrated
quantitative model.

- `INTENDED`: one aggregate housing-services market. The old multi-market
  choice dimension from the workhorse code is kept with `I=1`, so it is
  mechanically present but economically degenerate.
- `INTENDED`: Coven-style user-cost closure in stationarity:
  \[
  q=(r+\delta+\tau^p)P.
  \]
- `INTENDED`: aggregate housing services clear against an upward-sloping supply
  curve:
  \[
  H^D(q)=H^S(P)=H_0\left(\frac{q}{\bar q}\right)^\eta.
  \]
- `INTENDED`: renter housing is continuous up to `hR_max`; owner housing is a
  discrete rung choice on `H_own`.
- `INTENDED`: owner adjustment includes transaction/sale wedges and the
  workhorse down-payment constraint. In this code, `phi` is the financed share,
  so the down-payment threshold is \((1-\phi)Ph\).
- `INTENDED`: new owner choices face a payment-to-income screen. The current
  implementation blocks new purchases/adjustments when
  \[
  (q_{\text{int}}\phi+\tau^p)Ph > \psi^{PTI} y_a.
  \]
  Incumbent owners who keep the same owner rung do not requalify each period.
- `INTENDED`: lifecycle income profile \(y_a\) enters the budget and the
  payment-to-income screen.
- `INTENDED`: 4-year decision periods. Defaults are `age_start=22`, `J=16`,
  `J_R=11`, so agents enter at 22, retire around 66, and die after the last
  age-82 decision period.
- `INTENDED`: one dependent-child stage, `stage_durations=[1.0]`. Since a model
  period is 4 years, this means one 4-year child-at-home stage.
- `INTENDED`: period-unit flows. Discounting, interest, depreciation, property
  tax, income, child goods costs, and minimum consumption are scaled to the
  4-year decision period in the default parameter file.

## Simplifications

- `SIMPLIFICATION`: no idiosyncratic income/productivity state \(z\) yet.
  Current income heterogeneity is only lifecycle income by age. Adding a
  Coven-style \(z\) process is a real state-space extension and has not been
  done in this pass.
- `SIMPLIFICATION`: fertility remains the workhorse one-shot completed-family
  choice for childless fertile households. It is not a sequential parity hazard.
- `SIMPLIFICATION`: the model uses a collateral-constrained user-cost shortcut:
  \(qh\) enters the flow budget and \((1-\phi)Ph\le b\) enters the new-owner
  feasibility screen. The down payment is not subtracted as a separate asset
  purchase and housing equity is not a separate continuous state.
- `SIMPLIFICATION`: no mortgage coupon, amortization, maturity, refinancing, or
  present-value mortgage lock-in state.
- `SIMPLIFICATION`: the payment-to-income screen uses a simple interest-plus-tax
  payment. It is not yet Coven's full amortized mortgage payment formula.
- `SIMPLIFICATION`: the housing supply shifter `H0` and rent anchor `r_bar` are
  smoke-test normalizations, not calibrated targets.
- `SIMPLIFICATION`: old retention currently comes only through the inherited
  workhorse owner state, transaction/sale wedge, and bequest utility. A clean
  policy-created old-retention wedge is not yet implemented.

## Not Implemented

- `NOT IMPLEMENTED`: calibration, SMM objective, counterfactual tables, and
  parameter search.
- `NOT IMPLEMENTED`: idiosyncratic income process \(z_t\) with transition
  probabilities.
- `NOT IMPLEMENTED`: estate-tax counterfactuals, inheritance kernels, bequest
  principal adding-up, and estate-revenue rebates.
- `NOT IMPLEMENTED`: mortgage-rate lock-in from coupon gaps.
- `NOT IMPLEMENTED`: transition dynamics.
- `NOT IMPLEMENTED`: landlord balance sheets by house size. The current model
  clears aggregate housing services; rents are pinned by user cost.

## Verification Targets

Every nontrivial change should pass:

```bash
python -m compileall -q code/model/intergen_housing_fertility
cd code/model && .venv/bin/python -m intergen_housing_fertility.cli smoke --quiet
```

For model inspection, also run:

```bash
cd code/model && .venv/bin/python -m intergen_housing_fertility.cli diagnostics --fixed-prices --quiet
```
