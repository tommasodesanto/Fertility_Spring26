# Experimental initial age-profile pilot

Recovered from saved native outputs without another equilibrium solve. Five primary cases verified, one numerical candidate failed, and the selected point passed two exact numerical repetitions. The original reporting failure remains preserved remotely; it compared checkpoint provenance hashes rather than only numerical results.

The augmented objective changes from306.587264 to306.192916, only0.13%. The original component improves182.649147 to180.187063 while the six extra age rows worsen123.938117 to126.005852. This is not evidence of an improved initial age-profile fit, and no extra targets or pilot candidate were promoted to the main history.

- `selected_original_target_fit.csv`: all13 original rows, including12 scored and the separate2.1 normalization.
- `selected_extra_target_fit.csv`: all6 experimental age-profile rows.
- `selected_augmented_target_fit.csv`: the18 scored original and experimental rows.
- `selected_candidate_parameters.csv`: all17 estimates/restrictions, with the actually enforced beta cap0.99. The original raw scorer metadata remains separately saved.
- `all_cases.json`: truthful seven-case ledger, retaining the failed proposal.
- `all_*`: complete original/extra fit and parameter/profile comparisons for every valid case.
- `exact_score_discrepancies.json`: numerical equality and separately recorded checkpoint, receipt and solve-time differences.
- `source_and_checkpoint_verification.json`: saved source/input/checkpoint verification.

Extra weights use5% of each target as a synthetic scale, not an empirical standard error; overlapping covariance is ignored. This pilot is experimental and not production eligible.
