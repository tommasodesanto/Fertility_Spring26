# Permanent lowest-productivity credit diagnostic

This is a fixed-price household comparison at the October 3 working post-interest,
soft-financing chain-13 parameter vector, without Estate A. Both arms hold the
reference house price at `0.7760569760205563`; they are neither general
equilibria nor recalibrations. The only difference **between the two arms** is
the uniform financed share `phi`: `0.8` versus `1.0`.

Relative to the named chain-13 reference, the experimental economic changes are:

* Every entrant receives the original lowest productivity value
  `z=0.10346813919312171`, and every transition leads to that value. The
  usual gross age-earnings profile is retained and multiplied by this `z`.
  The nine numerical income nodes remain distinct, but only the lowest is
  occupied. No productivity-mean renormalization is applied.
* The original **marginal** distribution of entrant financial wealth over the
  120 wealth nodes is preserved exactly. Original wealth-productivity
  correlation cannot survive a single occupied productivity level. The entrant
  conditional wealth column for the low node is therefore set to the original
  weighted marginal; other columns are unreachable.
* The payroll tax remains `0.08028070961950022`. The native fixed-payroll
  PAYGO rule recomputes the balanced four-year pension as
  `0.09496140757220382` in both arms. Retirement income is not multiplied by
  productivity under the retained `retirement_income_z_scale=0`.
* The two `phi` values are the intended financing comparison. All other
  economic primitives, the birth menu, entry timing, preferences, targets,
  transfer rules, and housing price are held at the chain-13 inputs. No target
  loss is interpreted for this counterfactual population.

The [driver](run_fixed_price.py) uses `production.equilibrium.solve_at_price`.
It runs one core, with a 600-second and 8 GiB supervisor per arm and no
fallback. The source/input hashes and exact income and entry-wealth contracts
are in each arm's `input_contract.json`. The [exporter](export_children_inputs.py)
reads saved native distributions and prepares input for the established
children-by-age plotter without another solve. Rerun command, from the project
root, using the authenticated environment:

```sh
output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python output/model/experiments/low_productivity_credit_v1/run_fixed_price.py --phi 0.8 --preflight
```

Each complete arm already has a directory; the driver deliberately refuses
to overwrite it. Preserve those receipts if further diagnosis is needed.

## Verification and reporting status

| Arm | Solver and fiscal checks | Full aggregate reporter | Saved artifacts |
|---|---|---|---|
| `phi_08` | Unit population; zero mass outside low `z`; scaled PAYGO residual `1.90e-14`; supervisor source hash unchanged | Passed | Complete native arrays, plain arrays, 17 standard figures, summary |
| `phi_10` | Unit population; zero mass outside low `z`; scaled PAYGO residual `2.09e-14`; supervisor source hash unchanged | **Failed** at one buyer-policy value entry with `2.681e-15` population mass | Complete native arrays, plain arrays, 17 standard figures, failed supervisor receipt |

The flagged `phi_10` cell is at age 82, financial wealth
`-2.6279069767441863`, owner tenure index 2, lowest productivity, and a
one-child-ever-born state. Its saved value is `-9.540e9`, below the aggregate
reporter's infeasible threshold. The saved solution has not been altered, and
no aggregate-reporter gate was relaxed. Consequently, the `phi_10` full
wealth/consumption aggregate packet is not validated. The children-by-age
comparison uses only the unchanged current and beginning distributions and
passes its own mass checks; see [`children_by_age/`](children_by_age/).

The 17 standard figures include policy lines for all nine numerical income
nodes, although eight nodes are unoccupied. Their housing-market residual
views evaluate the held price and should not be read as market-clearing tests.
