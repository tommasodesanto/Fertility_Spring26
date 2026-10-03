# Children ever born when everyone has permanently low productivity

The [four-arm age comparison](reference_vs_low_productivity_by_age.png) shows
the share with no children and the share with two or more children ever born.
The [low-productivity count shares](figures/child_count_shares_by_age.png) and
[CDFs](figures/child_count_cdf_by_age.png) show all four count states and both
financed shares. Full values are in [comparison_by_age.csv](comparison_by_age.csv),
the [low-productivity age table](figures/children_by_age.csv), and its
[difference table](figures/differences_by_age.csv).

In the saved low-productivity household distributions, essentially all
households have no children ever born at every age. The largest share with any
children is \(1.73\times10^{-12}\); the two-or-more share is zero at the saved
precision. The \(\phi=0.8\) and \(\phi=1.0\) lines therefore overlap at the
scale of the figures. By contrast, the heterogeneous-productivity reference
has 18.888% childless and 67.044% with two or more at ages 42–45 under
\(\phi=0.8\). This is a stark outcome of this particular fixed-price low-income
experiment, not an estimated population forecast.

The experiment moves every household to the original lowest productivity
node, \(z=0.1034681392\), at entry and keeps it there by replacing the income
transition. The marginal entrant wealth distribution is held fixed, while its
conditional allocation across income types changes. The payroll tax remains
0.0802807096 and the pension is rebalanced from the reference 0.917784 to
0.0949614 per period. Price 0.7760569760 and the other loaded preferences are
held fixed. These income and fiscal changes move together; the comparison does
not isolate a single channel or solve a new general equilibrium. The two
low-productivity arms differ only in the financed share \(\phi\).

All age shares use the raw saved post-birth `g` distribution, including every
household and without deleting small cells. The age points are starts of
four-year model cells for the reproductive household member, rather than
single-year interview ages. [Source checks](figures/checks_and_provenance.json)
verify finite and nonnegative distributions, unit total mass, within-age share
sums, monotone CDFs, and current/beginning age-count margin agreement. The
[four-arm comparison receipt](comparison_provenance.json) pins both input CSVs.
The \(\phi=1\) aggregate policy reporter separately stopped on an infeasible
buyer-value cell with only \(2.68\times10^{-15}\) mass; the native distribution
and fiscal checks passed. That reporter issue was not filtered or waived to
make these age plots.

Regenerate the low-productivity age plots with the reusable
[saved-array script](../../../fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/credit_relaxation/children_by_age/plot_children_by_age.py),
passing the two files `phi_08/solution_arrays.npz` and
`phi_10/solution_arrays.npz` under this directory, and `--output figures`.
Then run [compare_with_reference.py](compare_with_reference.py) to regenerate
the four-arm figure and CSV. Neither plotting step runs the model.
