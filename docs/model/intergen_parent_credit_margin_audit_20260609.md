# Parent Credit Treated-Margin Audit

Date: 2026-06-09

## Scope

This note documents a diagnostic audit of the weak fertility response to the
parent-targeted credit-relief proof of concept.

The audit uses the final global-DE toy best:

- run directory:
  `results_intergen_housing_fertility_intergen_candidate_no_timing_v0_globalde_3g_20260609`
- task `27`, case `236`, label `de_g008_i011`
- baseline loss under `candidate_no_timing_v0`: `11.503191936648555`

The reproducible driver is:

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
code/model/.venv/bin/python code/model/tools/audit_intergen_parent_credit_margin.py
```

Primary output:

```text
output/model/intergen_parent_credit_margin_audit_20260609/
```

## Cases

| Case | Description |
|---|---|
| `baseline` | Final toy best |
| `parent_ltv95_birth` | Raises financed share to 95 percent for parent owner choices in child state \(c=1\) |
| `parent_ltv100_birth` | Raises financed share to 100 percent for parent owner choices in child state \(c=1\) |
| `parent_ltv100_all_child` | Raises financed share to 100 percent for all dependent-child states |

With the current one-stage child architecture, `parent_ltv100_all_child` is
mechanically identical to `parent_ltv100_birth`: there is only one
dependent-child state.

## Equilibrium Results

| Case | TFR | Childless | Prime own | Own 25--34 | Own 35--44 | Family own gap |
|---|---:|---:|---:|---:|---:|---:|
| `baseline` | `1.660` | `0.238` | `0.222` | `0.004` | `0.147` | `0.083` |
| `parent_ltv95_birth` | `1.677` | `0.236` | `0.227` | `0.004` | `0.151` | `0.091` |
| `parent_ltv100_birth` | `1.677` | `0.236` | `0.229` | `0.004` | `0.152` | `0.094` |
| `parent_ltv100_all_child` | `1.677` | `0.236` | `0.229` | `0.004` | `0.152` | `0.094` |

Going from LTV 95 to LTV 100 barely changes fertility. Extending the policy to
all dependent-child states also does nothing because the model currently has
one dependent-child state.

## Baseline-State Treated Margin

The audit evaluates policy functions on the baseline mass of fertile
childless renters:

\[
g(b,\text{renter},i,j,z,n=0,c=0), \qquad j\in\{26,30,34,38,42\}.
\]

This object is not aggregate TFR. It asks whether the households currently on
the active fertility margin change fertility or buy family housing when the
credit policy is applied.

| Case | TFR-equivalent on states | Birth probability | Birth and owner entry | Owner entry conditional on birth |
|---|---:|---:|---:|---:|
| `baseline` | `0.352` | `0.169` | `0.0158` | `0.0939` |
| `parent_ltv95_birth` | `0.353` | `0.168` | `0.0163` | `0.0966` |
| `parent_ltv100_birth` | `0.353` | `0.168` | `0.0165` | `0.0978` |
| `parent_ltv100_all_child` | `0.353` | `0.168` | `0.0165` | `0.0978` |

The fertility probability on baseline fertile childless renter states is
essentially flat. The policy raises the probability of buying conditional on a
birth only from `9.39%` to `9.78%`.

## Feasibility Exposure

For one child, the current calibrated housing need implies the first
family-feasible owner rung is `3.6` rooms on the temporary owner ladder
\((2,3.6,5.2,6.8,8.4,10)\).

On the baseline fertile-childless-renter mass:

| Case | DP feasible | PTI feasible | Both feasible | DP fail |
|---|---:|---:|---:|---:|
| `baseline` | `0.818` | `1.000` | `0.818` | `0.182` |
| `parent_ltv95_birth` | `0.884` | `1.000` | `0.884` | `0.116` |
| `parent_ltv100_birth` | `0.985` | `1.000` | `0.985` | `0.015` |
| `parent_ltv100_all_child` | `0.985` | `1.000` | `0.985` | `0.015` |

This rules out a simple PTI-blocking story for the first-child family rung:
the first-child PTI screen is slack for this baseline-state population. LTV
relief expands down-payment feasibility, but many households were already
feasible and the policy functions barely move.

For the two-child family rung, PTI becomes more relevant:

| Case | DP feasible | PTI feasible | Both feasible |
|---|---:|---:|---:|
| `baseline` | `0.818` | `0.932` | `0.807` |
| `parent_ltv95_birth` | `0.884` | `0.847` | `0.832` |
| `parent_ltv100_birth` | `0.985` | `0.716` | `0.714` |

But the two-plus probability on the baseline fertile-childless-renter margin is
only about `0.0074` in the baseline, so this intensive-margin PTI channel does
not move aggregate fertility much.

## Interpretation

The weak fertility effect is not primarily caused by the policy leaving a 5
percent down payment, and it is not caused by a first-child PTI block.

The current toy baseline has the wrong treated margin:

1. Early ownership is nearly zero: ownership ages 25--34 is `0.004`.
2. Among fertile childless renters, the baseline already has high feasibility
   for the first family owner rung, but owner entry conditional on birth is
   only `9.4%`.
3. The policy raises feasibility more than it raises choice probabilities.
4. Fertility probabilities on baseline renter states are almost unchanged, and
   young-age birth probabilities even fall slightly while late-age probabilities
   rise slightly.

So this is not evidence that credit constraints are irrelevant. It is evidence
that, at this toy parameter vector and temporary owner ladder, the active
fertility margin is not a clean down-payment-constrained young-parent
homeownership margin.

The next model improvement should focus on the baseline: lifecycle ownership,
the owner-entry policy shape, and the owner ladder. A parent-credit
counterfactual is only interpretable once the calibrated baseline puts
substantial young households near the housing-purchase constraint.
