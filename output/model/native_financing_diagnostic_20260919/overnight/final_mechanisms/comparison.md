# Final mechanisms collection

All three production families completed the same dose design: 2 controls plus 48 unique grid arms, with 17 primary graph records per arm. The refit smoke completed all 3 cases.

| Family | Price | Population source | Checkpoint SHA-256 | Baseline exact |
|---|---:|---|---|---|
| original | 0.702608328515 | stationary | `3322a61994fb3654d67f4b1d6cf2d0f7cacbb3668d06a417e192ee363c174993` | yes |
| stationary_new_income | 0.678030033681 | saved_evaluation | `79f9dd5e56351bb6ced3a2b31d10657b1bfa8acd59c1c68436bfb671c99bb675` | yes |
| refit_new_income | 0.546936972019 | saved_evaluation | `b3491eedcee6250cf94833067646d3e6463496cbf5a64bdcabe7b13bc7e89eb2` | yes |

The CSV reports, for each family, levels and changes in birth flow, first-birth flow, explicit lifetime cohort births, ownership (including percentage-point change), rooms, and first-birth mean age. It includes baseline comparisons, mortgage increments at rental caps 6, 8, and 10, and credit-dose comparisons at cap 6.

The flow and lifetime cohort measures need to be read separately. In both new-income families, credit-dose arms raise snapshot birth flow relative to baseline while lowering explicit lifetime cohort births and first-birth probability; the original family has positive lifetime-cohort changes (about 3.02% at lambda 5). Lambda 1 and lambda 5 are numerically identical in both new-income grids for the reported outcomes, an observed plateau without a causal interpretation.

Population audit: the within-family pre-choice mass is fixed across all arms. Original uses the stationary source; both new-income families use the saved evaluated source. The stored contract-level raw/evaluated differences are floating-point scale and remain below the (10^{-12}) absolute checks.

Graph manifest: 2550 remote PNG records with relative paths, sizes, and SHA-256 hashes; only the refit baseline 17 PNGs were downloaded.

Lead accounting check: the refit baseline and credit doses 0.25, 1 and 5 have exactly identical saved entry distributions. Summing the 17 age records exactly reproduces each reported lifetime birth total. Credit raises births at age 18 but lowers later-age births; lifetime first births also decline. The pattern therefore is not caused by changing the initial cohort or by summation error, and it is not merely a relabeling of birth timing. Its underlying economic or numerical cause remains to be diagnosed. [Saved accounting check](cohort_accounting_check.json).
