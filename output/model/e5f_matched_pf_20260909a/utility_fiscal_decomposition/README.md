# Fixed-preference utility and pension decomposition

All four cases in array **17358470** completed with exit 0:0. This is a stationary diagnostic at a fixed preference intercept, not calibrated SMM, a historical perfect-foresight result, or a policy benchmark. No early empirical targets, weights or loss are applied.

The source remains immutable **c6dd3508**. Every case uses the same nine structural coordinates and child preference intercept **0.2900515293650047**, with one stationary equilibrium solve, no fertility renormalization, a 600-second internal limit and a 12-minute Slurm limit. Supply elasticity is 0.63 in all four cases, using the same one-time seed supply rebasing. The utility change only replaces the old jump-plus-slope floor with its first-child sum and a zero slope. The pension comparison uses the inherited pension versus the pension implied by actual stationary payroll revenue and retiree exposure.

| Case | Utility | Pension rule | Asset price | Legacy `tfr` | Childlessness | Loop seconds |
|---|---|---|---:|---:|---:|---:|
| old_old | old_jump_plus_slope | old_unbalanced_diagnostic | 0.625524371 | 2.100004856 | 0.156808632 | 74.76 |
| old_balanced | old_jump_plus_slope | actual_budget_balanced | 0.634497242 | 2.101396359 | 0.156535128 | 104.99 |
| new_old | parenthood_only | old_unbalanced_diagnostic | 0.624406264 | 2.113468776 | 0.154399732 | 103.38 |
| new_balanced | parenthood_only | actual_budget_balanced | 0.633328602 | 2.114974717 | 0.154101695 | 75.25 |

The table retains the solver's legacy stationary observer names. In particular, `tfr` is not an annual female fertility-rate path or the separately certified early observer. The full set of 97 finite scalar diagnostic comparisons is in [raw_decomposition.csv](raw_decomposition.csv) and [raw_decomposition.json](raw_decomposition.json); all original moments remain in each collected summary. Nonscalar/unavailable metrics are listed explicitly in the JSON.

The period pension parameter rises from 1.7184 to 2.0463613896 under actual-budget balance. The old-pension cases intentionally leave about 16.0266% of payroll revenue unspent. Their fiscal invalidity is retained, not waived; these two cases cannot be selected as valid fiscal benchmarks. Balanced cases have relative fiscal residuals -1.41e-9 and 1.82e-10. All four pass the unchanged housing, household-budget, population/operator, probability and occupied-value screens. Budget-excess mass and occupied negative value steps are zero.

The comparison holds the nine coordinates exactly equal across cases, including the mapped first-child requirement 0.7168500088318734. Full values, existing bounds, status and near-bound flags are in [parameters_all_cases.csv](parameters_all_cases.csv); separate original parameter tables remain under `collected/<case>/repetition_01/parameters.csv`. Differences separate utility effects at each pension rule, pension effects under each utility, and their interaction. No parameter search or new normalization was performed.

## Reproducibility and receipts

- [contract.json](contract.json): new contract; SHA256 `65c7b48a7d26a8c0ca67632acf499aeb53705d2b7436469833e61b6df8ece6a9`.
- [submit_array.sh](submit_array.sh) and [submission_receipt.json](submission_receipt.json): four one-CPU/16-GB jobs, account `torch_pr_570_general`, partition `cpu_short`; no duplicate was queued at submission.
- [case_summary.csv](case_summary.csv): per-case job ID, runtime, fiscal and market residual, remote checkpoint path and SHA256.
- [sacct_final.txt](sacct_final.txt): Slurm completion, CPU, memory and runtime evidence.
- [collect_remote.py](collect_remote.py), [archive_receipt.json](archive_receipt.json), and [collected/collection_manifest.json](collected/collection_manifest.json): all four completion receipts and original checkpoint hashes verified remotely, then 132 light artifacts checked byte-for-byte after collection.
- [build_decomposition.py](build_decomposition.py): deterministic receipt/parameter validation and raw tabulation. Run with `/usr/bin/python3 build_decomposition.py` from this folder or using its absolute path; it performs no model solve.
- [original_graphs_receipt.json](original_graphs_receipt.json): all 68 original standard diagnostic graph hashes; the PNGs are under `collected/<case>/repetition_01/standard_diagnostics/`. All 68 decode successfully; [graph_decode_verification.json](graph_decode_verification.json) records dimensions. This collector did not perform a visual economic review.

The 625 source pins and inherited seed checksum were verified before submission and again by every case driver. The four large checkpoints remain on Torch; they were not downloaded. Remote output root: `/scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508/output/utility_fiscal_decomposition_20260911a/`. Previous source, smoke results, thresholds and tests were untouched. No failures, retries or additional jobs occurred. The lead will combine these raw comparisons with the separately certified early observers.
