# Fable review prompt: intergenerational bequest calibration and old-age balance sheet

You are reviewing the live intergenerational extension of the Fertility Spring 2026 model. Act as an independent senior quantitative-macro/household-finance referee and model diagnostician. Do not merely summarize the previous work: verify its accounting, challenge its conclusions, and recommend the smallest economically credible next step.

Repository root:

`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

## Scope and constraints

- This is a review and decision memo, not an authorization to edit model code, launch cluster jobs, or redesign the entire calibration.
- Use the local files linked below. No broad repository scan is needed.
- Do not conclude that bequest motives are absent merely because the current model cannot match the targets.
- The user is concerned that three bequest/estate moments now absorb too much calibration attention. Treat the number of hard targets as a substantive choice. Diagnostics do not automatically need to become targets.
- Preserve identification discipline. If you recommend dropping or demoting a target, state which parameter it currently disciplines and what moment or external restriction replaces it. If there is no replacement, say that the proposed system is underidentified.
- Distinguish a measurement/accounting problem, a numerical problem, and a missing economic mechanism.
- If you propose a model change, write the relevant equation or state transition precisely and name the empirical variation that would identify every new parameter.
- Recommend at most one new bounded experiment, and only if the existing evidence cannot settle the next decision. Give its acceptance and rejection criteria in advance.

## Mandatory orientation

Read these first, in order:

1. [`memory/AGENT_MEMORY.md`](../../memory/AGENT_MEMORY.md)
2. The latest note in [`memory/daily/`](../../memory/daily/)
3. [`CALIBRATION_STATUS.md`](../../CALIBRATION_STATUS.md)
4. [`code/model/README.md`](../../code/model/README.md)

The active implementation is the Python one-market model under:

[`code/model/intergen_housing_fertility/`](../../code/model/intergen_housing_fertility/)

Do not use archived MATLAB or obsolete two-location results as the live benchmark.

## The decision we need from you

The narrow question is:

> Does the poor fit of the old-age estate distribution mainly reflect the current bequest-preference specification, or does it reveal a more basic failure of the model's late-life balance sheet—especially its absence of positive nonhousing wealth among most old households? Given that diagnosis, what is the minimal credible calibration recipe?

We need a concrete recommendation about:

1. Which of the three estate moments should remain hard calibration targets.
2. Which should instead be diagnostics.
3. Whether \(\theta_0\), \(\theta_1\), and \(\theta_n\) should all remain internally estimated, or whether any should be externally restricted.
4. Whether the next task should be a measurement fix, a debt/portfolio mechanism repair, a bequest-specification change, a target redesign, or no further run until the model is repaired.

## Current target system

The live SMM system has 15 targets for 14 estimated parameters. Eleven established non-bequest targets are retained, and three old targets were replaced by three PSID estate targets:

| Estate target | Data target | Bootstrap SE | Intended parameter margin |
|---|---:|---:|---|
| Median total estate / annual income, ages 76–84 | 6.5013 | 0.2320 | overall bequest strength \(\theta_0\) |
| p90/p50 of total estate / annual income, ages 76–84 | 3.4481 | 0.1325 | curvature/threshold \(\theta_1\) |
| Median estate/income gap, 2+ children minus 1 child, ages 65–75 | 0.1011 | 0.5630 | child scaling \(\theta_n\) |

The empirical estate object is total net worth: nonhousing net worth plus home equity. The intended model counterpart is

\[
W=b+pH,
\]

with annual gross income as the denominator and no house-sale wedge in the measurement equation.

The external restriction is `tenure_choice_kappa = 0`. The remaining 14 parameters are estimated internally.

Please verify this system from the active target and calibration code rather than accepting the table on faith:

- [`code/model/intergen_housing_fertility/calibration.py`](../../code/model/intergen_housing_fertility/calibration.py)
- [`code/model/intergen_housing_fertility/parameters.py`](../../code/model/intergen_housing_fertility/parameters.py)
- [`code/model/intergen_housing_fertility/solver.py`](../../code/model/intergen_housing_fertility/solver.py)

## Current bequest specification

The intended specification is a warm-glow bequest term with linear child scaling, approximately

\[
B(W,n)=\theta_0\max\{1+\theta_n n,0\}
\frac{(\theta_1+W)^{1-\sigma}-\theta_1^{1-\sigma}}
{1-\sigma}.
\]

Inspect the exact implementation in `bequest_utility_vec` and the terminal Bellman logic. In particular, verify:

- whether \(W\) is gross housing plus liquid assets, net estate, or another object at each point in the code;
- where negative estates are clamped to zero;
- whether taxes, liquidation, housing, and debt are handled consistently between terminal utility and measured moments;
- whether the model moment uses the same timing convention as the empirical estate measure.

## Evidence already produced

### 1. Full internal recalibration

Main report:

[`output/model/intergen_internal_bequest_recalibration_20260715/report/`](../../output/model/intergen_internal_bequest_recalibration_20260715/report/)

The strict repeated winner has loss 166.65 and equilibrium residual \(1.22\times10^{-6}\). Its bequest parameters are approximately:

\[
\theta_0=0.310,
\qquad \theta_1=0.536,
\qquad \theta_n=0.710.
\]

All estimated parameters are interior, but the estate fit is poor:

| Moment | Target | Model | Approx. loss contribution |
|---|---:|---:|---:|
| Median estate/income | 6.501 | 6.267 | 1.02 |
| p90/p50 estate/income | 3.448 | 1.881 | 139.93 |
| 2+ minus 1 child family gap | 0.101 | 1.770 | 8.79 |

Inspect the complete 15-moment table and all parameter restrictions in the report. Do not judge the calibration from these three rows alone.

### 2. Local Jacobian

[`output/model/intergen_internal_bequest_recalibration_20260715/identification/`](../../output/model/intergen_internal_bequest_recalibration_20260715/identification/)

The full \(15\times14\) SMM-weighted Jacobian has numerical rank 9/14 at relative tolerance \(10^{-2}\), 12/14 at \(10^{-3}\), and condition number about 5,499. The target-scaled version is weaker. Its weakest direction mixes \(\theta_n\), `h_bar_0`, \(\theta_1\), and `h_bar_n`.

Assess whether this is acceptable local weakness, a sign that the three estate parameters are not separately identified by these moments, or evidence that the estate targets load on the wrong mechanisms.

### 3. Wide \(\theta_1\) reachability frontier

[`output/model/intergen_bequest_reachability_20260715/report/`](../../output/model/intergen_bequest_reachability_20260715/report/)

This diagnostic fixes the 11 non-bequest parameters, profiles \(\theta_1\) over 12 cells from 0.02 to 16, and reoptimizes \(\theta_0\) and \(\theta_n\) against the median and family-gap targets.

- No cell matches both the median and family gap within one bootstrap SE.
- Around \(\theta_1=1.2\), the median remains near target, but p90/p50 is only about 1.74 and the family gap is about 1.45.
- At \(\theta_1=8\) or 16, the distribution can become sufficiently dispersed, but the median collapses to roughly 2.52 or 2.19.

Interpret exactly what this rules out. It may rule out fixing the problem by moving \(\theta_1\) alone within the current equilibrium, but it does not by itself rule out bequest motives or alternative balance-sheet mechanisms.

### 4. Retirement-income heterogeneity diagnostic

[`output/model/intergen_retirement_income_dispersion_20260716/README.md`](../../output/model/intergen_retirement_income_dispersion_20260716/README.md)

[`output/model/intergen_retirement_income_dispersion_20260716/report/`](../../output/model/intergen_retirement_income_dispersion_20260716/report/)

This diagnostic adds persistent retirement-income heterogeneity,

\[
y_j(z)=\bar y_j[1+s_R(z-1)],
\]

profiles \(s_R\in\{0,0.25,0.5,0.75,1,1.5,2\}\), fixes the 11 core parameters, and reoptimizes the three bequest parameters.

- The best strict cell is \(s_R=0\).
- No cell matches all three estate targets within one bootstrap SE.
- The high-dispersion cell \(s_R=2\) raises p90/p50 only to about 2.37 while lowering the median to about 5.51; it also puts \(\theta_0\) at zero and is not a strict equilibrium solution.

This experiment is diagnostic-only. Decide whether it cleanly rejects income dispersion as the first-order missing mechanism or whether its design is too restricted to do so.

## New balance-sheet decomposition—the central evidence

Empirical construction and outputs:

- [`code/data/psid_followup_mar2026/diagnose_intergen_bequest_distribution.R`](../../code/data/psid_followup_mar2026/diagnose_intergen_bequest_distribution.R)
- [`code/data/psid_followup_mar2026/output/intergen_bequest_distribution_diagnostic/`](../../code/data/psid_followup_mar2026/output/intergen_bequest_distribution_diagnostic/)

Matched model construction and outputs:

- [`code/model/tools/diagnose_intergen_bequest_distribution.py`](../../code/model/tools/diagnose_intergen_bequest_distribution.py)
- [`output/model/intergen_bequest_distribution_diagnostic_20260716/README.md`](../../output/model/intergen_bequest_distribution_diagnostic_20260716/README.md)
- [`output/model/intergen_bequest_distribution_diagnostic_20260716/`](../../output/model/intergen_bequest_distribution_diagnostic_20260716/)

The matched model diagnostic uses a fresh strict solve with residual \(1.19\times10^{-6}\) and reproduces the three calibrated estate moments to numerical precision.

All values below are normalized by each dataset's median annual income at ages 76–84:

| Distribution, ages 76–84 | PSID p50 | PSID p75 | PSID p90 | Model p50 | Model p75 | Model p90 |
|---|---:|---:|---:|---:|---:|---:|
| Annual income | 1.000 | — | 2.957 | 1.000 | — | 1.000 |
| Total estate | 6.867 | 17.101 | 36.335 | 6.267 | 8.866 | 11.789 |
| Nonhousing net worth | 2.368 | 10.336 | 26.272 | 0.000 | 0.000 | 0.000 |
| Home equity | 3.192 | 6.396 | 11.609 | 6.267 | 8.866 | 11.789 |
| Estate/income | 6.501 | 12.685 | 22.417 | 6.267 | 8.866 | 11.789 |

Additional facts:

- Estate-income correlation is 0.393 in the PSID and 0 in this model diagnostic.
- The top estate decile's housing share is about 0.186 in the PSID and 0.999 in the model.
- The raw total-estate p90/p50 is 5.291 in the PSID and 1.881 in the model. Therefore income in the denominator compresses the empirical estate/income ratio; denominator noise is not creating the missing upper tail.
- In the current model diagnostic, the normalized liquid/debt state \(b\) has p50/p75/p90 approximately \(-6.44,-5.14,-3.84\), while gross housing \(pH\) has p50/p75/p90 approximately 11.97, 15.96, and 19.95.
- The diagnostic maps owner home equity to \(pH+\min\{b,0\}\), owner nonhousing wealth to \(\max\{b,0\}\), and renter \(b\) to nonhousing wealth. Verify whether that mapping is dictated by the model or merely a diagnostic convention.

This decomposition suggests that the model is matching the estate median almost entirely with housing while generating virtually no positive nonhousing wealth through the p90. That is potentially more fundamental than a poor choice of \(\theta_1\).

### Family-size decomposition, ages 65–75

| Group | PSID estate/income p50 | PSID p75 | PSID p90 | Model p50 | Model p75 | Model p90 |
|---|---:|---:|---:|---:|---:|---:|
| 1 child | 4.804 | 11.607 | 24.936 | 11.968 | 16.336 | 20.975 |
| 2+ children | 4.905 | 9.989 | 18.675 | 13.738 | 18.285 | 22.832 |

In the PSID, the 2+-minus-1-child median gap is about \(-0.886\) among married households and \(+0.884\) among nonmarried households, producing the small aggregate target of \(+0.101\). The model has no marital-status state.

Assess whether the aggregate family-gap target is a credible hard target for \(\theta_n\), given this cancellation and the model's lack of marriage. Do not interpret the small aggregate gap as evidence that parents do not alter saving motives with family size; child costs, selection, marital composition, and bequest preferences all enter the reduced-form gap.

## Required code/accounting audit

Trace the exact old-age balance-sheet logic in the active code. Focus on:

1. **Meaning of \(b\).** Is it liquid wealth, net financial wealth, mortgage balance folded into a single net asset, or a mixture? Is the diagnostic split into mortgage debt and nonhousing wealth valid?
2. **Collateral and debt limits.** How are borrowing caps tied to housing? Can households carry large negative \(b\) late in life, refinance implicitly, or avoid amortization?
3. **Timing.** Is the reported distribution beginning-of-period wealth, post-choice wealth, or an internally inconsistent mixture? Are empirical age bins compared with the correct model ages?
4. **Terminal incentives.** What happens to debt, housing, taxes, and liquidation at death? Does the bequest function reward gross housing while inadequately penalizing debt?
5. **Policy and distribution numerics.** Could interpolation, the asset grid, KFE transitions, boundary mass, or market-clearing prices mechanically create the negative-\(b\) distribution?
6. **Income process.** Why is old-age model income degenerate in the matched packet? Is that an intended restriction, a lost persistent state, or a diagnostic aggregation error?
7. **Housing ladder.** Do house-size discreteness, tenure persistence, transaction costs, and the existing loan-to-value rule push essentially all retirement wealth into housing?

Useful implementation anchors include:

- `income_at_state`, `annual_gross_income_at_state`, `add_annual_gross_estate_moments`
- `bequest_utility_vec` and terminal Bellman code
- owner feasibility floors and collateral constraints
- `forward_distribution_markov`
- the moment construction in [`code/model/tools/run_intergen_bequest_exit_chain.py`](../../code/model/tools/run_intergen_bequest_exit_chain.py)

If the empirical/model component comparison is not apples-to-apples, identify the exact mismatch and explain how it changes the diagnosis.

## Questions your recommendation must answer

1. Is the failure to match p90/p50 primarily a taste-curvature failure, an asset-composition failure, or both? What evidence separates them?
2. Why is \(b<0\) through the model p90 at ages 76–84? Is that economically intended and empirically plausible?
3. Does the terminal bequest object treat housing and debt coherently?
4. Is the p90/p50 target still useful as a hard target for \(\theta_1\), or should it remain a diagnostic until the balance sheet is repaired?
5. Is the aggregate children gap useful as a hard target for \(\theta_n\), or is it too composition-driven without marital status? If demoted, how should \(\theta_n\) be identified internally—or should it be externally restricted?
6. Is the median estate/income target sufficient to retain as the single bequest target for overall strength, conditional on an external restriction for the other taste parameters?
7. Which existing non-bequest moments, if any, already contain identifying information about retirement saving, leverage, or portfolio composition?
8. What is the minimal defensible target system that does not devote a disproportionate fraction of calibration effort to this block?

## Deliverable

Write the review to:

`docs/model/intergen_bequest_balance_sheet_fable_review_20260716.md`

Use this structure:

1. **Executive verdict** — no more than six bullets.
2. **Verification table** — each important prior claim marked verified, qualified, or rejected, with file/code evidence.
3. **Economic and numerical diagnosis** — identify the causal chain producing the old-age balance sheet.
4. **Minimal calibration recipe** — a table listing every current estate moment as keep hard target / diagnostic / drop, the parameter it disciplines, and any replacement restriction or moment.
5. **Parameter recommendation** — explicit treatment of \(\theta_0\), \(\theta_1\), and \(\theta_n\).
6. **Next action** — either no new run, or one bounded experiment with exact design and ex ante pass/fail criteria.
7. **Open risks** — what remains genuinely unresolved.

Lead with your own judgment. The useful answer is not “the model fits poorly”; it is a defensible explanation of why, and a small recipe for what to do next.
