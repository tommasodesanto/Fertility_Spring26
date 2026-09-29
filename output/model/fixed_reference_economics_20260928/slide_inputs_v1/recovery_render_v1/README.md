# Three-price recovery figure and table

**Historical render:** actual job 18815186 passed. The lead's visual check found
overlapping footer lines, corrected by the formatting-only [v2 successor](../recovery_render_v2/README.md)
on Torch 18817838. Use its final figure; v1 numerical inputs and table values
are unchanged and this packet remains preserved.

**Synthetic Torch verification passed:** job **18814969**, including rendered
test-only output and rejection of wrong hashes, missing factors, wrong regimes,
comparison/elasticity disagreement, and a false five-price completion claim.
Evidence is in `torch_validation/`. No synthetic number is an economic result.

The lead submitted actual zero-solve rendering job **18815186**, dependent on
successful recovery job **18815133**; invalid dependencies cancel the render.
It uses `render_completed.sh` (one CPU, 4 GiB, five minutes) and requires a
passed, authenticated three-price completion receipt. Do not submit a duplicate.
Frozen remote renderer:
`/scratch/td2248/projects/fixed_reference_economic_figures_20260929/recovery_render_v1/source_v1/build_three_price_inputs.py`,
SHA-256 `7e0d3a806dca7c0f27fce677921de5fef105a3aec3f11efa5e6beb9b5d06d766`.
Inputs will come from
`/scratch/td2248/projects/fixed_reference_elasticity_recovery_20260929/solve_v1/`;
actual outputs will be in
`/scratch/td2248/projects/fixed_reference_economic_figures_20260929/recovery_render_v1/results/actual_v1/`.
Log and error files are `render.log` and `render.err` under that renderer root.
Actual output remains pending until that job succeeds and the lead visually
checks the chart.

This renderer consumes only the completed three-price comparison; it does not run or import the model. It requires explicit paths for `comparison.csv`, `elasticities.csv`, `completed.json`, and a new output directory:

```sh
python build_three_price_inputs.py \
  --comparison /path/to/comparison.csv \
  --elasticities /path/to/elasticities.csv \
  --completed /path/to/completed.json \
  --output /path/to/new_output
```

The completion receipt must say `status: passed`, `complete_three_price: true`, `complete_five_price: false`, carry the exact frozen reference label, list only factors `[0.99, 1.0, 1.01]` and step size `[0.01]`, and pin both CSV SHA-256 hashes. The comparison must cover each of the three prices for both regimes and both scopes, with no extra price factors. Regimes are `reference` and `credit`; scopes are `impact` and `cohort`. The elasticity table reads the central log elasticity at the 1% step for total, first, and second births on impact and completed fertility by cohort.

Outputs are one two-panel PNG/PDF figure, a four-row CSV/TeX elasticity table, and a manifest with source paths and hashes. The graph plots percent changes from each regime's own baseline over ±1%; it does not extrapolate to ±2%. Notes state that house price and mapped rent move together, all preferences including \(\psi_{\mathrm{child}}\) stay fixed, and impact uses the same inherited households. The cohort distribution has unit household mass; completed fertility is free to change and the cohort calculation is not a transition. The exact frozen reference label appears on every output. Existing 14-moment, 31-parameter fits and standard diagnostic plots remain linked through their run packet and are not copied into a report here.

Run the synthetic contract tests only on Torch. They create clearly branded `TEST ONLY` figures under a separate test output folder and check a complete three-price render plus rejection of missing-price, wrong-regime, wrong-hash, and five-price-claim inputs:

```sh
python test_three_price_inputs.py /path/to/new/synthetic_test_only
```
