# Earnings and pensions: transition accounting check

September 12, 2026. Scope: the six-date 2019–2039 presentation patch, not a
certified infinite-horizon solution or a completed 2007–2023 historical fit.
No new model solve or calibration change was performed.

## Economic identity

The May production interpretation is consistent with the maintained unit wage:
competitive firms have $Y_t=AL_t$, $w=A=1$. Let $e_a$ denote effective labor
supplied over the four-year model period, and $z_t$ the combined permanent-group
and persistent earnings state. Then $L_t$ aggregates $e_a z_t$ across working
households and the gross wage bill $W_t=wL_t$ equals output under this
interpretation. In code, $e_a$ is four times the annual age profile.

With payroll tax $\tau^{SS}=0.179$ and actual retired-household mass $H_t^R$,
flat benefits satisfy $\varpi_t H_t^R=\tau^{SS}W_t$. Consequently aggregate
net labor income plus pensions equals the gross wage bill:

$$(1-\tau^{SS})W_t+\varpi_t H_t^R=W_t.$$

This identity remains valid as the working and retired populations change;
the benefit adjusts, rather than remaining at its initial value.
Here $H_t^R$ is only an audit notation for retirees, not rented housing.

## Code and units

The quantitative task checked the actual frozen checkpoint: `w_hat=[1.]`,
`period_years=4`, `retirement_income_z_scale=0`, permanent income groups enabled,
15 combined income states, and payroll tax 0.179. The initial pension is overridden
by the dated path. Initial checkpoint SHA-256:
`120ffc45c0fb8756f4182f999c96b7c0236adf315cb938190ec31cd2068c87c2`.

- Frozen `tools/e5f_social_security.py:29–49` binds workers' after-tax annual
  earnings times four, and an already-period pension for retirees. Pension
  binding is idempotent: no second multiplication by four.
- `intergen_eqscale_seq_optimized/solver.py:226–233` multiplies worker income
  by the current income state. Retirement income is common because its
  income-state scale is zero. Property-tax transfers are added separately.
- `e5f_social_security.py:78–124` uses the actual current household distribution,
  including its age and income composition. It does not use resident persons,
  stationary reference weights, or fixed initial worker/retiree ratios.
- `run_e5f_perfect_foresight_transition.py:541,625` and the person-demography
  forward path apply the same dated fiscal values before computing household
  policies. The backward solution therefore incorporates the pension path
  that the forward simulation actually pays. The surprise-path observer also
  checks the dated pension and tax explicitly.
- The joint path root checks housing and pension residuals and their replay;
  accepting a generic solver termination is insufficient.

Local independent fiscal checks:
`/opt/anaconda3/bin/python -B -m unittest test_e5f_social_security -q`
in `tmp/e5f_matched_pf/code/model/tools`: **23 tests passed**. These are pure
accounting tests, including hand-calculated income, actual retiree exposure,
period scaling, population rescaling, transfers excluded from payroll, and
zero-exposure safeguards. They do not solve the model.

## Saved transition receipts

All values below are in the saved normalization and four-year flow units.
The lead independently recomputed tax × payroll base, benefit × retiree mass,
and net labor income plus pensions minus the gross wage bill at every date.
The maximum absolute pension-budget residual is $1.589\times10^{-10}$;
the maximum relative residual is $2.386\times10^{-10}$.

| Year | Gross wage bill | Retiree mass | Pension per retiree | Payroll receipts | Pension payments | Receipts minus payments |
|---|---:|---:|---:|---:|---:|---:|
| 2019 | 2.963205891 | 0.259198526 | 2.046361383 | 0.530413855 | 0.530413855 | 9.297e-12 |
| 2023 | 3.546519182 | 0.272458481 | 2.329995128 | 0.634826933 | 0.634826934 | -1.127e-10 |
| 2027 | 3.569002079 | 0.306245684 | 2.086074693 | 0.638851372 | 0.638851372 | 6.729e-11 |
| 2031 | 3.598614112 | 0.329915927 | 1.952472956 | 0.644151926 | 0.644151926 | 1.058e-10 |
| 2035 | 3.663189303 | 0.336839647 | 1.946655899 | 0.655710885 | 0.655710885 | 3.240e-11 |
| 2039 | 3.719692940 | 0.339373335 | 1.961925018 | 0.665825036 | 0.665825036 | 1.588e-10 |

Source: `source/expected_transition.csv`, with the root receipt in
`source/root_receipt.json`. `source/verification.json` records an exact replay
(maximum aggregate difference zero) and the same source identity. In 2023,
resident persons are 2.9659779384 and household heads are 1.1292476712; those
counts must not be interchanged. Worker and retiree head masses sum to the latter.

## Limits of the conclusion

The baseline property-tax regime is unrebated. For example, 2023 property-tax
receipts are 0.2152850541 and equal-transfer payments are zero. The separate
`government_budget_residual` therefore is not a Social Security failure, but
neither is it evidence of a balanced total government budget. No government
spending allocation should be invented to close that field.

The labor-production formula interprets fixed wages; it does not establish a
separately solved closed goods or asset market. No complete aggregate resource
constraint, external-asset position, or allocation of the unrebated property-tax
receipts was verified in this check. These remain outside the certification.

The finite path clears housing and balances pensions. Its terminal-distance and
horizon checks fail, and the final historical fertility window remains unmatched.
Thus the check validates the implemented dated earnings/PAYGO accounting, not
long-run convergence or a fully closed general equilibrium.

Slide corrections: define effective labor in model-period units and replace
“both fiscal budgets hold” with the verified pay-as-you-go pension condition.
