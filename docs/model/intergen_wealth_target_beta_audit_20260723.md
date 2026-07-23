# Wealth-target and conditional-beta audit

Date: 2026-07-23

Scope: the timing-repaired **current one-shot model**, not the E-series. The
conditional profile deliberately leaves the current fourteen-moment system
unchanged, including the tenure-Brier row and free `tenure_choice_kappa`, so
that it isolates the beta question. It is not a promoted calibration.

## 1. What the published target establishes

[De Nardi and Yang (2014)](https://users.nber.org/~denardim/research/De-Nardi_Yang_EER.pdf)
report a wealth/after-tax-earnings target of 6.90 and match it with annual
beta 0.96 (Table 2). They attribute the target to Hendricks (2007).

That citation does not, by itself, give this project a matched target.
[Yang (2009)](https://crr.bc.edu/wp-content/uploads/2009/01/wp_2009-6-508.pdf),
in a closely related model and citing the same Hendricks paper, instead
reports 4.9. The same paper describes Hendricks's earnings as labor income of
the household head and spouse net of both income-tax payments and Social
Security contributions. The difference between 4.9 and 6.9 may reflect a
changed tax or wealth normalization, but neither paper provides enough detail
to reconcile the two numbers mechanically.

## 2. Why the current code is not a matched denominator

The timing-repaired model statistic is internally coherent:

\[
  \frac{\sum_{\text{all living households}}(b_t+pH_t)}
       {\sum_{\text{working households}} y_t^{\mathrm{after\ payroll}}}.
\]

Mortgage debt is carried in \(b_t\), so \(b_t+pH_t\) is net worth rather than
gross housing wealth. Pensions and lump-sum transfers are excluded from the
denominator.

The denominator, however, deducts only `tau_pay = 0.179`, the Greaney
payroll/pension wedge used in the model. It does not reproduce Hendricks's
income-tax-plus-Social-Security construction. The metadata label
`borrowed-target-matched-denominator` is therefore too strong. The defensible
status is **borrowed target with an unresolved tax normalization**.

## 3. Transparent PSID diagnostic

`code/data/psid_followup_mar2026/audit_aggregate_wealth_earnings_ratio.R`
constructs a ratio that avoids the unresolved tax convention. It uses one
reference-person row per family, ages 18--82, `NETWORTHR` in the numerator,
RP/spouse combined `EARNINDRRC` for households whose reference person is at
most 62 in the denominator, and the longitudinal weight `IW`.

| PSID vintage | Gross wealth / gross earnings | Equivalent after only the model's 17.9% payroll wedge |
|---|---:|---:|
| 1984--2003 | 5.209733 | 6.345594 |
| 2005--2019 | 7.024080 | 8.555518 |

These are diagnostics, not adopted targets. They do not reproduce Hendricks's
TAXSIM calculation, correct PSID top-wealth undercoverage, or settle the
appropriate calibration vintage. They establish two points:

1. The normalization matters materially.
2. There is no evidence-based case for lowering the target merely to obtain a
   conventional beta; the recent-vintage gross ratio is higher, not lower.

## 4. Conditional annual-beta profile

Six bounded Torch searches reoptimized the other thirteen parameters at three
fixed annual betas. Every chain completed 320 evaluations, passed the strict
equilibrium check, and produced two exact strict repeats. The lower loss from
the two chains in each cell is:

| Annual beta | Strict loss | Wealth / current after-payroll earnings | Living-old p90/p50 | TFR |
|---:|---:|---:|---:|---:|
| 0.980 | 0.398933 | 3.871255 | 2.040787 | 2.031539 |
| 0.990 | 0.287679 | 5.016297 | 1.925948 | 2.066726 |
| 0.995 | 0.244511 | 5.665811 | 1.916880 | 2.064915 |

For comparison, the unrestricted timing-repaired result has loss 0.237479
and annual beta 0.999364. Thus fixing beta at 0.995 costs only 0.00703 in the
current objective, but beta 0.99 or 0.98 produces progressively larger wealth
shortfalls that nuisance-parameter reoptimization does not undo.

As a score-only sensitivity, replacing the current row by the 1984--2003
gross target 5.209733 lowers the beta-0.995 loss from 0.244511 to 0.223994 at
the same parameter vector. Using the 2005--2019 gross target 7.024080 raises
it to 0.326599. These are not reoptimized calibrations; they show why the
vintage must be chosen before interpreting beta.

Complete per-cell target fits, parameter tables, and strict-chain checks are
in
`output/model/intergen_new_moment_beta_reasonable_probe_20260723/report/`.

## 5. Decision

The beta result is real **conditional on the current 6.90 row**, but that row
is not yet a defensible matched target. The next clean specification should:

1. replace the ambiguous after-tax ratio with aggregate net worth divided by
   aggregate gross labor earnings;
2. choose and document the empirical vintage;
3. measure the data and model objects identically; and
4. rerun the beta profile before changing the saving mechanism.

If beta remains implausibly high after that repair, the next structural
diagnosis is the model's thin precautionary- and intergenerational-saving
channels. De Nardi and Yang combine much larger earnings risk, inherited
ability, accidental and intended inheritances, a consumption floor, Social
Security, and pensions; their beta 0.96 is not generated by the wealth target
alone.
