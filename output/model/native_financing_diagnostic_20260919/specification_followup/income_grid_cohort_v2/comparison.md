# Income-grid cohort collection

Jobs `18079046` (smoke) and `18079047` (production) completed successfully. The three readout arms are defined by joint income-state counts: 5x3 has 15 states, 9x3 has 27, and 15x3 has 45. The saved aggregate comparison is in [comparison.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/income_grid_cohort_v2/comparison.csv). The three rows report lifetime explicit births, first-birth probability, and conditional first-birth age directly from each saved `receipt.json`:

| arm | joint income states | lifetime explicit births | first-birth probability | conditional first-birth age |
|---|---:|---:|---:|---:|
| 5x3 | 15 | 1.857720411166992 | 0.792264752756171 | 24.70009336511201 |
| 9x3 | 27 | 1.814527049505432 | 0.782737380377278 | 25.39634456590778 |
| 15x3 | 45 | 1.806882123522255 | 0.781607775405935 | 25.57706921192440 |

All three saved entry wealth marginals have 120 positions and unit mass. Pairwise L1 distances, computed from the saved `entry_distribution_summary.json` files, are: 5x3 versus 9x3 `0.354668498907649`, 5x3 versus 15x3 `0.390308140567585`, and 9x3 versus 15x3 `0.395733510646080`. A wealth mean is unavailable because the saved summaries contain no wealth-node values; no checkpoint or new household solve was loaded.

The 17-row `cohort_by_age.csv` profiles remain supplemental artifacts for each arm. They are actual entry-cohort profiles, not stationary ACS-fit moments. The comparison is a grid specification readout; it does not identify a pure policy effect, convergence result, GE effect, or calibration change.

The remote PNG hash manifest covers 51 graphs, with 17 per arm. The 17 production 15x3 standard plots remain copied under `results/plots/production_15x3_standard_diagnostics/`.

See [collection_receipt.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/income_grid_cohort_v2/collection_receipt.json) and [plot_manifest.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/income_grid_cohort_v2/results/plot_manifest.json).
