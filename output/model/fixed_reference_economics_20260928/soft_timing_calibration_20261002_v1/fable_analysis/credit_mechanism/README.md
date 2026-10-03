# Supplemental native-grid housing and tenure mechanism

This packet reads the pinned revised-interest-timing, soft-purchase-constraint
chain 13 selected repeat. It does not solve or change the model. The exact
`solution_arrays.npz` path and SHA-256 are in `native_grid_data.json` and are
checked against `../revision1/saved_array_source_receipt.json` by the script.

`housing_tenure_fertility_native.png` compares conditional policy schedules for
ages 22–25, inherited renters, location 0, and discrete income states 4 and 6.
At identical beginning net financial wealth \(b\), it plots the first-birth
attempt probability when childless, rental rooms if the renter branch is chosen,
ownership probability, and expected physical rooms across tenure choices.
Family schedules use \((n,m)=(0,0),(1,1),(2,2)\), where \(n\) is children ever
born and \(m\) children currently at home. They are not a forced-birth
counterfactual. The executed saved housing floor is zero at \(m=0\) and
\(2.593759507\) rooms for every feasible \(m\geq1\); it does not rise again
for the second child. Renters can choose up to six rooms. Owner products have
2, 4, 6, 8, or 10 rooms.

`childless_mass_native.png` shows the native-node pre-fertility childless mass
used to identify the occupied wealth range. Saved
`g_beginning_distribution` is post-fertility and pre-tenure. For childless
renters at a given node and income state, pre-fertility mass equals saved
post-fertility mass divided by \(1-\pi a\), where \(a\) is the saved first-birth
attempt probability and \(\pi=1-0.02\exp(0.134(22-18))\) is fecundity.
Mass is normalized within each income state. There are 28 positive-mass nodes
per state, from \(b=0\) to \(15.474\). Figures show nodes through \(b=5.5\),
covering more than 99.9% of this conditional mass in each state; the full
native-node values are in `native_grid_data.json`. Dotted vertical lines mark
the 99.9% cutoffs: \(b=2.535\) in state 4 and \(5.465\) in state 6.

Income states 4 and 6 are discrete model states, not income quantiles. Their
saved multipliers are approximately 0.470 and 1.287. With the age-22–25
four-year profile value 2.650830657 in `code/model/production/inputs.py`,
the corresponding pre-tax period earnings are approximately 1.245 and 3.413
model units. The asset coordinate is beginning **net financial wealth**, not
house value or total wealth.

At \(b=0\), income state 4 has a first-birth attempt probability of 0.000061,
with 2.376 childless rental rooms versus 3.966 rooms in the one-child schedule;
ownership probabilities are 0.003 and approximately zero. In income state 6,
the corresponding values are 0.697, 4.781 versus 6.000 rooms, and 0.389 versus
0.281 ownership probability. These are same-coordinate conditional policies,
not observed changes for a household following a birth. `selected_native_nodes.csv`
holds exact values at \(b=0\) and \(b=1\).

The purchase and borrowing mechanism is **not diagnosed** by these plots.
Owner saving policies use a transaction-adjusted asset coordinate, and the
revised-interest-timing purchase map must be applied to reconstruct an ending
debt slack at fixed beginning \(b\). The housing and tenure panels use direct
saved policy slices at their correct inherited-renter coordinate.

Regenerate:

```sh
MPLBACKEND=Agg code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/plot_credit_mechanism.py
```
