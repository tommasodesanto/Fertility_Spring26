# Verified initial candidate

Loss 179.2984252281, versus rebated seed 182.6491468669 (1.83% lower). Two exact numerical reproductions passed. The search stopped after repeated numerical root evaluation limits; this is a verified improvement, not a converged optimizer or a certified final calibration.

Twelve scored restrictions and one separate completed-fertility normalization; nine free structural parameters. This table retains the original measurement and weighting caveats in the linked full CSV.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.1 | 8.24545e-07 | — | — |
| Childless women, ages 40–44 | 0.198279 | 0.196258 | -0.00202037 | 35532.3 | 0.145039 |
| Exactly one child among mothers, ages 40–44 | 0.213655 | 0.216439 | 0.00278364 | 26952.8 | 0.208848 |
| Period mean first-birth age | 25.9763 | 26.0952 | 0.118982 | 139.828 | 1.9795 |
| First births at age 30+ | 0.249278 | 0.232155 | -0.0171228 | 13866.1 | 4.06539 |
| Wealth / annual gross labor earnings | 6.14586 | 4.89343 | -1.25243 | 7.5951 | 11.9136 |
| Annual bequests / aggregate wealth | 0.0088 | 0.00858241 | -0.000217591 | 5.16529e+06 | 0.244554 |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 4.4998 | 0.983865 | 10.6164 | 10.2765 |
| Mean occupied rooms, capped at 9 | 5.5611 | 6.42465 | 0.863555 | 128.021 | 95.4685 |
| Ownership, heads 30–55 | 0.648334 | 0.539896 | -0.108438 | 2339.36 | 27.508 |
| First-birth room response, −1 to +3 | 0.720246 | 0.980806 | 0.26056 | 137.565 | 9.3395 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.172043 | -0.175024 | 280.528 | 8.59354 |
| Recent-parent ownership gap | 0.162896 | 0.144103 | -0.0187928 | 27055.8 | 9.55534 |

| Parameter | Estimate | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| beta_annual | 0.99 | 0.94 | 0.99 | True | diagnostic candidate; not a certified estimate |
| kappa_fert | 0.33773423411167025 | 0.02 | 50.0 | True | diagnostic candidate; not a certified estimate |
| kappa_fert_continuation | 0.39785645171756406 | 0.02 | 50.0 | True | diagnostic candidate; not a certified estimate |
| chi | 1.0496534423047694 | 0.1 | 5.0 | False | diagnostic candidate; not a certified estimate |
| H0 | 8.11210048056786 | 0.2 | 80.0 | False | diagnostic candidate; not a certified estimate |
| theta0 | 0.08105103333987912 | 0.0 | 8.0 | False | diagnostic candidate; not a certified estimate |
| theta1 | 0.08520356663830632 | 0.02 | 16.0 | True | diagnostic candidate; not a certified estimate |
| first_birth_fixed_cost | 0.26576477618628114 | 0.0 | 8.0 | False | diagnostic candidate; not a certified estimate |
| h_P | 2.3 | 0.1 | 2.3 | True | diagnostic candidate; not a certified estimate |
| hbar_child_rooms | 0.0 | — | — | False | zero restriction |
| psi_child | 0.1489153142205283 | — | — | False | normalized to 2.1 |
| payroll_tax | 0.179 | — | — | False | externally fixed |
| pension_period | 2.0463613896121218 | — | — | False | budget derived |
| housing_supply_elasticity | 0.63 | — | — | False | externally fixed |
| tenure_choice_kappa | 0.005 | — | — | False | externally fixed |
| alpha_cons | 0.733 | — | — | False | externally fixed |
| sigma | 2.0 | — | — | False | externally fixed |

The displayed beta cap is the enforced search restriction. The raw scorer table records an obsolete 0.9995 upper bound and is preserved separately. Other near-bound flags are the frozen scorer convention, which uses the full search range; they do not imply a constraint binds exactly.

Housing, household, source/target and accounting gates passed. Pension relative gap: 3.83e-10; property-tax rebate relative gap: 3.84e-7. Historical fertility and policy paths remain separate work.
