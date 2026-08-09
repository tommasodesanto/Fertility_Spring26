# PSID Fertility-IV Full-Panel Diagnostic

## Bottom line

Keeping the full event-window panel does not solve the precision problem. The
rooms sample contains 2,334 mothers but only 38 same-birth-year first-birth
proxies, versus 34 in the five-year endpoint design. The moved-for-size panel
contains 35 proxy-positive mothers, versus 34 before. The binding restriction
is therefore observing the mother both before and after the birth, not the
earlier requirement for an interview specifically at event time +4 or +5.

The one suggestive result is the twins rooms response. The weighted panel IV
estimate is 0.706 rooms (SE 0.424, p=0.096), close to both the first-birth
event-study magnitude and the five-year endpoint twins estimate. Without PSID
weights it is 0.671 rooms (SE 0.302, p=0.027). This is not a robust 5-percent
result: the weighted estimate misses that threshold, and the weighted dynamic
reduced form has a marginal joint pre-period test (p=0.082). It should be
reported as suggestive triangulation, not selected as a new headline estimate.

No other weighted average panel IV is statistically significant. The same-sex
rooms response remains negative, rather than supporting the positive
first-birth rooms response.

## Design and estimand

The twins design is indexed at the first birth and the same-sex design at the
second birth. Each design retains event times -4 through +5 and requires at
least one usable pre-event and one usable post-event observation for the
particular outcome. The same-sex sample excludes same-year multiple-birth
proxies at the first and second births.

The average panel specification is

\[
Y_{it}=\alpha_i+\lambda_t+\gamma_k+
\beta\left(Z_i\times Post_{it}\right)+X_i'\Pi Post_{it}+u_{it},
\]

where $\alpha_i$ is a mother fixed effect, $\lambda_t$ is a calendar-year fixed
effect, and $\gamma_k$ is a common event-time profile. Baseline age, age
squared, education, and tenure are allowed to shift the post-event response;
the rooms regression also uses baseline rooms, and the same-sex design uses
first-child sex. Move outcomes control for the length of the interview gap.
Standard errors are clustered by mother. Weighted specifications use the last
available pre-event PSID individual weight.

For 2SLS, current additional-child status is instrumented by
$Z_i\times Post_{it}$. Thus the panel IV is an average effect of currently
having the additional child on the housing state or interview-interval outcome
over event times 0 through 5. It is not identical to the five-year cumulative
endpoint estimand.

## Weighted average results

| Design | Outcome | Mothers | $Z=1$ mothers | First stage | $F$ | Reduced form (SE) | 2SLS (SE) | $p$ for 2SLS |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Twins proxy | Rooms | 2,334 | 38 | 0.691 | 751.4 | 0.488 (0.289) | 0.706 (0.424) | 0.096 |
| Twins proxy | Moved for size | 2,211 | 35 | 0.691 | 722.1 | -0.023 (0.072) | -0.033 (0.104) | 0.751 |
| Twins proxy | Ownership, baseline renters | 1,192 | 17 | 0.692 | 310.1 | 0.165 (0.180) | 0.238 (0.251) | 0.344 |
| Twins proxy | Moved to ownership, baseline renters | 776 | 8 | 0.724 | 242.9 | -0.016 (0.049) | -0.022 (0.067) | 0.747 |
| Same sex | Rooms | 2,173 | 1,081 | 0.039 | 6.6 | -0.121 (0.076) | -3.114 (2.353) | 0.186 |
| Same sex | Moved for size | 2,027 | 1,001 | 0.039 | 6.1 | -0.002 (0.014) | -0.048 (0.355) | 0.892 |
| Same sex | Ownership, baseline renters | 1,179 | 568 | 0.054 | 5.9 | 0.011 (0.036) | 0.210 (0.661) | 0.751 |
| Same sex | Moved to ownership, baseline renters | 820 | 402 | 0.057 | 4.6 | 0.011 (0.016) | 0.189 (0.298) | 0.526 |

The very large twins first-stage statistics arise because a same-birth-year
first-birth proxy almost mechanically switches additional-child status on at
the birth. They do not overcome the fact that only 8--38 instrument-positive
mothers identify the housing reduced forms. Conversely, the same-sex panel
first stage is weaker than the five-year endpoint first stage because current
third-child status remains zero during early post-second-birth observations.

## Dynamic reduced forms

The weighted twins rooms reduced form rises to 0.859 rooms at +3 (SE 0.497,
p=0.084) and 1.129 at +5 (SE 0.608, p=0.064), relative to event time -2. The
joint pre-period test has p=0.082, so the dynamic path is directionally useful
but not clean enough to claim a significant causal effect.

The same-sex rooms path reaches -0.318 rooms at +5 (SE 0.145, p=0.028). This
opposite-signed result may reflect sampling noise or a direct sex-composition
effect on room sharing; it does not validate a positive additional-child rooms
effect. Joint post-period rejections for some twins move outcomes are driven by
coefficients with alternating signs and should not be described as supportive
evidence. The twins ownership path also fails its joint pre-period test.

## Files and reproduction

- `iv_housing_panel_estimates.csv`: weighted and unweighted first stages,
  reduced forms, fixed-effect OLS, and 2SLS estimates.
- `iv_housing_panel_samples.csv`: outcome-specific observation and mother
  counts.
- `iv_housing_panel_eventstudy.csv`: weighted dynamic reduced forms.
- `iv_housing_panel_tests.csv`: joint pre- and post-period tests.
- `eventstudy_*.png`: diagnostic plots generated locally by the driver.

Run:

```bash
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -b do \
  /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/iv_housing_panel_reaudit_20260809.do
```

This packet is a robustness diagnostic. It does not replace the five-year
mother-level audit, create a new calibration target, or resolve the active
first-birth rooms-code hold.
