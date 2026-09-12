# Fertility age-support check

The age-support sensitivity reduces the observed decline across the matched four-year windows from **16.67% to 15.11%**. This matters modestly; it does not eliminate the historical decline or settle the target convention. No calibration target, weight, model equation or production input was changed. No equilibrium was solved; the full continuation remains paused.

| Birth window | Published TFR | Ages 18–44 | Approximate ages 18–45 |
|---|---:|---:|---:|
| 2008–2011 | 1.9749 | 1.9118 | 1.9124 |
| 2012–2015 | 1.8610 | 1.8191 | 1.8199 |
| 2016–2019 | 1.7554 | 1.7264 | 1.7273 |
| 2020–2023 | 1.6458 | 1.6229 | 1.6239 |

Annual 2007-to-2023 declines are 23.54% for published TFR and 22.04% for ages 18–44. These differ from the four-year-window comparisons above.

## Construction and limitations

The partial index is (2 × rate18–19 + 5 × sum of rates20–24 through40–44)/1000, in births per woman. Four-year observations average the annual indices equally, preserving the existing empirical block convention. The model decision in 2007 maps to births2008–2011; decision2019 maps to births2020–2023.

The model covers ages18–45. The cached tables do not isolate age45. The approximate column adds one year at the published oldest-group rate; this is an explicit approximation, not an exact age45 rate. NCHS computes its oldest rate using births to women45+ divided by women45–49. The full older-group contribution is saved for scale, not used as a mathematical bound on the age45 contribution.

The published15–19 rate need not equal a duration-weighted average of the15–17 and18–19 rates because female exposures differ. Annual CSV separates the young-age contribution, the teen regrouping difference, and the older contribution; their sum reproduces the published-minus-partial difference exactly. No female exposure counts were inferred from rounded rates.

The existing initial first-birth timing contract collapses boundary ages into model endpoint cells and preserves all first births. It is a different moment and was not changed. Treating one model household as one representative potential mother remains an approximation, not a newly established arithmetic error. Moving to a support-restricted production target would require a consistent initial-level and transition observation convention; this audit alone does not authorize that change.

## Saved model check

All three short-path selected mappings exactly match their saved final mappings. Recomputed age-specific rates from top-code-adjusted birth flows divided by model age mass reproduce the saved period fertility indices. No extra factor four is needed: annualization cancels the four-year age-cell width. No births occur outside the modeled fertility cells.

The three saved short paths imply declines of 8.21%, 13.24%, and 23.59% across these windows. These are shock diagnostics under the earlier unrestricted calibration, not a fitted shock path or results for the new beta cap. Their complete prior calibration fit and parameter tables remain in ../return_home_20260911/READOUT.md. `blocks.csv` reports all four observations for all three paths and gaps against both empirical definitions. Long-horizon acceptance is not established by these short paths.

## Reproduce and sources

Run `python output/model/e5f_matched_pf_20260909a/current_candidate_transition/check_fertility_age_support.py` from the project root. This reads existing cached source text and saved model outputs only. See verification.json for source hashes and checked identities.

- NCHS, Births: Final Data for2015, Table4 (PDF page20): https://www.cdc.gov/nchs/data/nvsr/nvsr66/NVSR66_01.pdf —2007–2009 inputs.
- NCHS, Births: Final Data for2023, Table2 (PDF page14; oldest-age footnote page15): https://www.cdc.gov/nchs/data/nvsr/nvsr74/nvsr74-1.pdf —2010–2023 inputs.
