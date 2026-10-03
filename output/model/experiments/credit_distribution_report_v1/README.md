# Saved household distributions under two fixed-price credit comparisons

The [report manifest](manifest.json) indexes twelve paired figures, four concise comparison tables, full underlying CSVs, source hashes, validation checks, and the 17 original diagnostic plots for each of four saved cases. The companion PDF is assembled from this manifest. Rebuild figures and CSVs without a model solve using:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python \
  output/model/experiments/credit_distribution_report_v1/build_figures.py
```

Both comparisons hold the house price at `0.7760569760205563` and use the October 3 post-interest chain-13 parameters without Estate A. Within each population, the only primitive changed between arms is the uniform financed share \(\phi=0.8\) versus \(\phi=1\). The latter removes the nominal down payment but retains the owner collateral floor and other purchase screens. These are fixed-price household outcomes, not market-clearing general equilibria, transitions, or recalibrations; no counterfactual target loss is reported. The existing housing-market residual figures evaluate the held price.

The permanent-low-productivity population additionally places every entrant and future income state at the original lowest productivity value, preserves the original *marginal* entrant financial-wealth distribution, and rebalances the pension under the retained payroll tax. Its levels should not be interpreted as an isolated credit effect relative to the heterogeneous-productivity population. [The low-productivity source README](../low_productivity_credit_v1/README.md) records the complete economic change list.

## Distribution timing and units

`g` is the realized post-birth, post-tenure cross-section used for children ever born, children currently at home, ownership, owner room sizes, and realized rental housing services. `g_beginning_distribution` is post-fertility but pre-tenure. It weights beginning net financial assets \(b\) and the explicitly constructed pre-transaction gross net worth \(b+p h\), where inherited renters have \(h=0\) and inherited owners have their fixed room size. The latter excludes selling costs and is an accounting display, not the empirical wealth-target observer. All stock units retain the original mean annual gross-earnings normalization, including the low-productivity cases; consumption is a four-year flow in those same units. Ages are starts of four-year model cells: 18 means ages 18–21. The last ever-born count state means three or more children.

The policy aggregation follows `code/model/tools/model_policy_tools.py`: buyers use mass `g - g_stay_distribution`, owner stayers use `g_stay_distribution`, and their respective consumption and next-asset policies remain separate. The policy CDF plots use exact occupied mass sorted by value, sampled at probability grid points; complete histogram and sampled CDF CSVs are retained. The whole-grid financial-wealth and net-worth CSVs retain tail mass even where a figure displays the central 1–99% range. Next financial assets are not national-account saving.

## Validity boundary

The native aggregate-policy reporter passes for both \(\phi=0.8\) arms. It fails its unmodified initialized-infeasible-value gate for both \(\phi=1\) arms. In the heterogeneous population, 234 flagged buyer cells carry total mass `5.25755947452896e-38`; in the low-productivity population, one flagged buyer cell carries `2.680992916328218e-15`. The latter is the already documented terminal-age owner cell. No mass was filtered and no gate was weakened. Accordingly, the consumption and next-financial-asset figures show only the validated \(\phi=0.8\) arms. The \(\phi=1\) policy entries are missing in the concise quantile table, rather than coded as zero.

All four `g` and beginning distributions have unit mass to floating-point tolerance, nonnegative finite cells, matching age mass, and complete within-age child-count shares. Renter housing is used only where occupied renter mass has finite value above the native infeasibility threshold, strictly positive room demand at most six, and at least the active `2.593759507364224` room parent floor when children are at home. All four arms pass these separate checks; owned rooms follow realized tenure and need no renter-policy inference. Source hashes, exact mass residuals, and the policy-gate locations are in `manifest.json`.

In the permanent-low population, households with one child ever born have pooled mass below `1.1e-12` in either arm, while two-or-more-child states have zero saved mass. The conditional ownership point for the one-child group is therefore based on negligible mass; the age-by-count CSV records its denominator. Empty two-or-more groups have no ownership rate.

The four concise `tables/comparison_*.csv` files are formatted for the PDF. The detailed case files include every age-cell count share and CDF, ownership denominator, occupied productivity-state rate, owned and rented room distribution, beginning-asset grid and age quantiles, gross-net-worth grid/tenure/age distributions, and validated policy histograms and CDFs. `manifest.json` lists each path and all 68 original standard diagnostics.
