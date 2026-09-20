# Exhaustive rental-wedge v1 collection

Jobs `18080796` (smoke) and `18080797` (production) both completed with exit `0:0`; runtimes were 4:50 and 12:04. Six cases produced one saved household solve each, 17 standard graphs each, and complete cohort receipts. Source verification passed before and after both jobs; the local source-manifest hash matches the launch pin, and source trees were unchanged. No policy arrays were downloaded.

| case | phi | slope | cap | snapshot birth flow | snapshot first-birth flow | lifetime cohort births | lifetime first-birth probability | cohort first-birth age | solve seconds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| cap6zero | 0.8 | 0 | 6.0 | 0.115590387830670 | 0.050649833856012 | 1.872410763655674 | 0.820460038842370 | 24.09524572744819 | 31.663 |
| cap10zero | 0.8 | 0 | 10.0 | 0.116957157814593 | 0.051618977888043 | 1.883434886344046 | 0.822684429925007 | 23.99790311889151 | 23.561 |
| cap10s005 | 0.8 | 0.05 | 10.0 | 0.115595104549409 | 0.050653538368573 | 1.872415500649061 | 0.820460776282149 | 24.09477867765188 | 171.869 |
| cap10s02 | 0.8 | 0.2 | 10.0 | 0.115590387804873 | 0.050649833810685 | 1.872410767298959 | 0.820460040145778 | 24.09524574318061 | 152.354 |
| cap10s1 | 0.8 | 1 | 10.0 | 0.115590387830670 | 0.050649833856012 | 1.872410766751898 | 0.820460040088875 | 24.09524573599539 | 164.964 |
| cap10s02_phi1 | 1.0 | 0.2 | 10.0 | 0.116388093829978 | 0.051179789196443 | 1.876793059575202 | 0.820121941452310 | 24.02538525867143 | 155.354 |

All six saving audits passed the (10^{-7}) maximum-value-gain gate; budget excess masses were zero, population raw-to-evaluation L1 gaps were zero, and all six standard-plot counts were 17. The 102-entry remote plot hash manifest is `results/plot_hash_manifest.json`; the 17 central `cap10s02` plots are copied under `results/cap10s02/standard_diagnostics/`.

The matched φ=1 case is retained as a direct comparison to `cap10s02` (φ=0.8). The saved-array supplement records fixed-price value-choice-set inequalities and descriptive renter mass above six rooms; these are artifact checks and descriptive quantities only. No GE, frictionless benchmark, or calibration/adoption claim is made.

The canonical reproducible saved-array readout is `lead_saved_array_review.json`, documented in `lead_saved_array_review.md`; the earlier hand-entered supplement remains in `results/saved_array_supplement.json` with its pre-axis-fix audit copy.
