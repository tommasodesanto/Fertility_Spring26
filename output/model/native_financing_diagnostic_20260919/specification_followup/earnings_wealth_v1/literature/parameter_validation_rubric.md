# Parameter validation rubric for the earnings-and-wealth calibration

**Purpose.** This is a pre-calibration plausibility audit for the frozen scorer's
17-row `external_inputs/parameters.csv` table and the live nine-parameter search
bounds. It defines what can be compared with external evidence and what must be
judged in the model's own normalized units. It does not assess the current fit,
select parameters, or change any bound.

**Lead definition correction, September 22:** the table below was checked against
`e5f_parenthood_utility.parenthood_utility_metadata`, the frozen solver's housing
supply rule, fertility logit, and `bequest_utility_vec`. Earlier labels incorrectly
called H0 a utility scale, theta1 a child-dependence parameter, and h_P a per-child
requirement. These are documentation corrections only; no parameter, equation,
target or bound changed.

## Evidence classes

Use one status for every row in the eventual full table:

| Status | Meaning | Required evidence |
|---|---|---|
| Directly comparable | Same economic object, time unit, income concept, and normalization. | Source definition, sample/estimate, and conversion to the model period documented. |
| Comparable after conversion | External object is portable only after an explicit, reproducible conversion. | Annual versus four-year units, gross/net treatment, and aggregation formula recorded. |
| Model-specific | Utility coefficient, taste shifter, or normalized technology object without a portable literature level. | Bounds, sign, monotonicity, and implied behavior audited internally. |
| Fixed or derived | Externally fixed contract value or algebraic consequence of another object. | Source contract or derivation; no independent plausibility claim. |

The reviewed income literature supports external risk comparisons only after
matching period and income concept: Sommer (2016) supplies annual persistent and
iid wage-risk inputs; De Nardi (2004) estimates five-year PSID cells directly;
Bick (2016) constructs three-year gross-income cells; De Nardi and Yang (2016)
uses a five-year Markov process; Kolasa (2024) is a four-year fertility model
but its Warsaw working-paper earnings entries are explicitly annual. These are
measurement precedents, not portable values for the project's normalized
household process. See [multiyear review](multiyear_review.md).

## Seventeen-row audit map

| Row | Object and class | Audit rule |
|---|---|---|
| `beta_annual` | Annual discount factor; **directly comparable after period conversion** | Record \(\beta_4=\beta_a^4\). At the live bounds [0.94, 0.99], the implied four-year factors are 0.78075 and 0.96060. Check the four-year factor and the implied patience against saving, ownership, housing, and old-wealth profiles. Do not compare annual beta directly with a four-year Bellman coefficient. |
| `kappa_fert` | First-birth extreme-value taste-shock scale; **model-specific utility units** | Audit positivity and the smoothness of first-birth responses by age, current income and wealth. Literature can motivate the fertility mechanism, but cannot port this normalized coefficient. |
| `kappa_fert_continuation` | Subsequent-birth extreme-value taste-shock scale; **model-specific utility units** | Report its relation to `kappa_fert` and the subsequent-birth hazard; distinct margin-specific scales do not impose an ordering restriction in this sequential specification. Treat near-bound values as an identification warning, not as evidence for a new bound. |
| `chi` | Multiplier on owner housing services relative to rental services; **model-specific** | Check the tenure response and housing services conditional on income, wealth, age, and children. It is not a portable housing-price elasticity. |
| `H0` | Housing supply intercept; **model-specific quantity normalization** | The supply rule is \(H^s=H_0(r/\bar r)^\xi\). Check its units jointly with rents, prices, rooms demand and market clearing; it is not a utility coefficient. |
| `theta0` | Bequest utility intercept/scale; **model-specific** | Inspect bequest probability and bequest flow by age, wealth, and number of children. Compare the resulting aggregate bequest flow to the target only after verifying the target's timing and denominator. |
| `theta1` | Shift inside bequest utility in wealth units; **model-specific** | The bequest utility uses \(\theta_1+b\) (or per-child estate under equal division). Audit marginal bequest incentives and low-estate behavior. Child dependence is governed separately by the bequest specification and any child-scale parameter. |
| `first_birth_fixed_cost` | First-birth utility cost; **model-specific one-time utility cost** | Inspect first-birth timing, childlessness, and first-birth age response. Do not compare its level to an observed dollar cost without the complete utility and period normalization. |
| `h_P` | Housing-services floor whenever dependent children are present; **model-specific physical model unit** | The live upper bound is 2.3 and the maintained rule is \(\bar h(m)=h_P\mathbf{1}\{m>0\}\), with zero additional per-child slope. Inspect this first-child kink. It should be audited jointly with room caps and the first-birth housing response. |
| `hbar_child_rooms` | Child-space restriction; **fixed restriction** | Record whether zero means inactive. Verify the constraint in policy functions and the rooms-by-child plots; do not interpret zero as an estimated preference. |
| `psi_child` | Fertility normalization to completed fertility 2.1; **derived/fixed normalization** | Recompute the normalization from the stated fertility target and preserve the target contract. It is not an independent preference estimate. |
| `payroll_tax` | Payroll tax; **directly comparable external fiscal object** | Compare rate, tax base, incidence, and period treatment with the authoritative source contract. Keep separate from income-risk or utility normalization. |
| `pension_period` | Pension replacement/period object; **derived or externally fixed contract value** | Recompute from the pension rule and annual-to-four-year timing. Audit retirement income and wealth accumulation; do not compare the scalar without the benefit formula. |
| `housing_supply_elasticity` | Housing supply elasticity; **directly comparable only with same supply object** | Require the same geographic market, supply definition, and price/rent measure as the source. Inspect price, rent, quantity, and market residual plots together. |
| `tenure_choice_kappa` | Tenure extreme-value taste-shock scale; **model-specific utility units** | Check ownership probability, renter/owner housing services, and sensitivity around the financing threshold. It is not a portable ownership rate or mortgage coefficient. |
| `alpha_cons` | Consumption exponent in the consumption-housing composite; **model-specific normalization** | Verify utility curvature and expenditure shares over age, wealth, tenure, and children. Do not compare its level with a demand-system share unless the utility aggregator is identical. |
| `sigma` | Relative-risk-aversion/utility curvature; **model-specific but economically interpretable** | Check consumption smoothing, saving, wealth concentration, and bequest behavior. Literature can provide a broad risk-aversion context, but the normalized level is not directly portable across aggregators. |

The nine live search coordinates are the first nine active rows through `h_P`.
The remaining rows are restrictions, normalizations, externally fixed objects,
or derived quantities. A full audit must retain all 17 rows even when only nine
are searched.

## Boundary and implied-quantity checks

For each active coordinate, record the estimate, lower and upper bounds, distance
to each bound in both raw and normalized units, and whether the optimum is within
the project's declared near-bound tolerance. A boundary flag is diagnostic only;
it does not justify widening a bound. Repeat the flag under a one-sided local
perturbation when the point is at a bound.

Always report the economic quantity implied by the parameter, not only the raw
coefficient:

- Convert annual discounting and interest to four-year objects. With the
  inherited annual benchmark, \(R=1.02^4=1.08243216\) and
  \(\beta=0.99^4=0.96059601\). A candidate at the lower live beta bound implies
  \(0.94^4=0.78074896\). These quantities govern saving and purchase timing.
- For `h_P`, tabulate the minimum housing requirement for one through the maximum
  modeled number of resident children. Inspect the discrete kink where a family
  moves from feasible to infeasible housing, including renter and owner states.
- For `theta0` and `theta1`, report bequest utility and the implied bequest
  resource threshold at zero, one, and several children across old-age wealth.
  The annual bequest-flow target is a flow-to-aggregate-wealth ratio, so it must
  not be treated as a direct utility target.
- For fertility costs and tastes, report the change in first-birth probability,
  first-birth age, childlessness, and completed fertility along income and wealth
  slices. Use children currently at home and children ever born consistently.
- For housing and tenure tastes, report implied rooms, ownership, mortgage
  eligibility, and renter/owner consumption at the financing boundary.

## Standard 17-plot inspection

The rubric is satisfied only when the standard 17 diagnostic plots are inspected
for lifecycle shape, income sorting, and boundary behavior. At minimum, record
whether each plot shows monotone or economically explained non-monotone patterns:

1. fertility by age;
2. fertility by age and income state;
3. housing by age and income state;
4. housing market quantities;
5. housing prices;
6. income-state outcomes;
7. liquid wealth by age and income state;
8. market-clearing quantities by market;
9. market-clearing residuals;
10. owner policy rungs;
11. ownership by age;
12. ownership by age and income state;
13. childless-renter policy at age 30;
14. childless-renter policy at age 42;
15. tenure services;
16. childless-renter wealth distribution at age 30;
17. childless-renter wealth distribution at age 42.

Inspect poorest and richest income states, youngest and oldest working ages,
renter and owner states, childless and parent states, and the highest child-count
state wherever the plot supports them. Flag mass at wealth-grid boundaries,
discontinuous room or ownership jumps without a contractual kink, negative or
non-monotone values that lack an economic explanation, and market residuals that
are large relative to the declared solver tolerance.

## Decision record template

For every row, the eventual calibration table should append: `status`, source or
derivation, period and unit, bound distances, implied economic quantity, plots
checked, and unresolved concern. A literature value may support a prior plausibility
range, but promotion into a production specification requires matching the income
concept, timing, normalization, target contract, and identification count. This
rubric deliberately leaves all parameter choices and bounds unchanged.

**References:** Sommer (2016), Section 4.1/Table 2;
De Nardi (2004), Appendix A.3/Table A.1;
De Nardi and Yang (2016), pp. 132–133;
Bick (2016), Appendix C.1/Table C.1;
Kolasa (2024), JEDC and Warsaw working-paper Table 1;
Doepke and Kindermann (2019), three-year fertility timing;
Storesletten, Telmer, and Yaron (2004), persistent/transitory income-risk
decomposition.
