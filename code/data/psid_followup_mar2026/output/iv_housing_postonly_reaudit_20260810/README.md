# PSID Fertility-IV Post-Only Diagnostic

## Bottom line

Dropping the pre-birth-observation requirement increases the twins sample, but
does not produce a robustly significant rooms effect. The rooms regression now
contains 3,521 mothers and 52 same-birth-year first-birth proxies, compared with
2,334 mothers and 38 proxies in the full-panel design and 1,974 mothers and 34
proxies in the baseline-adjusted five-year endpoint design.

The post-only twins rooms IV remains positive and economically similar to the
main first-birth response: 0.776 rooms (SE 0.828, p=0.349) with exact maternal-
age and race fixed effects. Its large standard error reflects the fact that 52
proxy-positive mothers are still a small effective sample.

Post-birth ownership produces the most favorable result. Twins increase the
average ownership rate over event times 0 through 5 by 0.330 in the flexible-
age weighted IV (SE 0.179, p=0.065). With the more conventional quadratic age
control, the estimate is 0.392 (SE 0.178, p=0.028). Because significance is
sensitive to the maternal-age specification, this is suggestive robustness,
not a new causal headline.

## Design

The design deliberately removes the pre-birth housing requirement. It keeps
one observation per mother, requires an interview around event time +4 or +5,
and constructs outcomes over event times 0 through 5:

- rooms and ownership are the mother's average observed housing state;
- moved-for-size and moved-to-ownership indicate whether the transition is
  ever observed during the window.

The twins design is indexed at the first birth and instruments a second child
within five years with the same-birth-year sibling proxy. The same-sex design
is indexed at the second birth, excludes same-year multiple-birth proxies, and
instruments a third child within five years with same sex of the first two.

The preferred robustness absorbs exact maternal-age-at-event indicators,
event-birth-year indicators, and race indicators. The same-sex design also
controls for first-child sex. A quadratic maternal-age specification is
reported because it matches the usual empirical implementation. Regressions
use the PSID individual weight observed around +4/+5 and heteroskedasticity-
robust standard errors; unweighted estimates are included as sensitivity.

## Weighted 2SLS results: exact age and race controls

| Design | Outcome | Mothers | $Z=1$ | First-stage $F$ | 2SLS (SE) | $p$ |
|---|---|---:|---:|---:|---:|---:|
| Twins proxy | Rooms mean, 0:5 | 3,521 | 52 | 156.4 | 0.776 (0.828) | 0.349 |
| Twins proxy | Moved for size, 0:5 | 3,282 | 46 | 145.3 | 0.124 (0.220) | 0.571 |
| Twins proxy | Ownership mean, 0:5 | 3,475 | 52 | 154.4 | 0.330 (0.179) | 0.065 |
| Twins proxy | Moved to ownership, 0:5 | 2,971 | 43 | 124.4 | -0.094 (0.089) | 0.287 |
| Same sex | Rooms mean, 0:5 | 2,684 | 1,333 | 4.4 | -1.895 (1.823) | 0.299 |
| Same sex | Moved for size, 0:5 | 2,494 | 1,233 | 3.9 | 0.171 (0.468) | 0.714 |
| Same sex | Ownership mean, 0:5 | 2,658 | 1,318 | 4.4 | 0.081 (0.372) | 0.827 |
| Same sex | Moved to ownership, 0:5 | 2,312 | 1,136 | 4.9 | 0.165 (0.224) | 0.460 |

The corresponding weighted twins ownership reduced form is 0.117 (SE 0.063,
p=0.065). Under quadratic age controls it is 0.143 (SE 0.065, p=0.028), giving
2SLS 0.392 (SE 0.178). The unweighted exact-age twins ownership estimate is
0.255 (SE 0.144, p=0.076).

The same-sex first stages remain weak. Its unweighted rooms reduced form is
negative and significant, but the weighted reduced form and 2SLS are not; this
does not support a positive additional-child rooms response and retains the
direct room-sharing concern.

## Balance and interpretation

Conditional only on event-birth year, twins-proxy mothers are about 2.3 years
older (weighted p=0.033). This is expected because multiple births rise with
maternal age, but it makes flexible age adjustment essential. First-child sex
is balanced. Post-only housing levels are also more vulnerable than within-
mother changes to unobserved socioeconomic differences and fertility treatment.

Accordingly, this exercise is useful corroboration of positive twins effects
on rooms and ownership, but it should not replace the baseline-adjusted IV or
the first-birth event study.

## Reproduction

Run:

```bash
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -b do \
  /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/iv_housing_postonly_reaudit_20260810.do
```

Compact outputs are `iv_housing_postonly_estimates.csv`,
`iv_housing_postonly_samples.csv`, and `iv_housing_postonly_balance.csv` in
this directory. This is a diagnostic, not a calibration target.
