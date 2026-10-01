# Matched credit diagnostic from saved q0 data

The saved-policy native forward decomposition is also complete, without a
household solve. [FORWARD_RESULTS.md](FORWARD_RESULTS.md) reports the exact
three-corner accounting by age, inherited net financial assets and origin
children states; its identity and add-up receipts are retained alongside it.

The separately reviewed v2 analytical reduction is complete. See
[RESULTS.md](RESULTS.md), `saved_reduction_v2.json` and
`completion_receipt_v2.json`. Its first-birth attempt response is positive at
every fecund age on common recoverable support; eligible excluded mass is
reported explicitly. No model acceptance gate or engine was changed.

Lead-reviewed remote execution was attempted once and failed the strict
action-probability sum gate before producing a reduction. No numerical
comparison is available. `failure_receipt.json` and `remote_stderr.log`
preserve the failure; no retry or tolerance change was made.
`reduce_saved.py` reads the two completed q0 solution archives and the exact
baseline pre-fertility checkpoint on Torch. It imports NumPy and standard
library modules only. It neither solves the model nor aggregates housing or
borrowing policies. Only its compact JSON output is collected locally.

The reference is frozen chain 11 case 0026, loss 146.4246168285, rather than the
newer calibration candidate. The comparison holds each baseline pre-fertility
childless state and its probability weight fixed across the reference-credit
and native solvency-credit regimes. Both use the same logit scale
0.17614503485474298. For interior action probabilities, the recovered value gap
is \(\kappa\log(p_{\rm try}/p_{\rm wait})\). Action levels are the saved
inclusive value plus \(\kappa\log p_a\). These are wait/try action values,
including their continuation consequences; they do not separate current
utility from continuation utility or identify a raw birth-success value.

Each age row reports matched-weight probability and value changes, positive
and negative change fractions, and excluded mass. Zero fecundability,
non-fertile ages, dead values and zero/one action probabilities are excluded
from value recovery without clipping. The report is conditional on common
recoverable support, and its weighted exclusions are part of the result.

`check_reduce.py` passed synthetic action/gap inversion, matched weighting,
sign fractions, zero-fecundability-age exclusions, dead/corner exclusions and
strict rejection of invalid probabilities. `toy_test_receipt.json` records
the zero-model-call test. A reviewed remote execution will use one process,
one numerical thread and a 180-second external limit. No Slurm job or model
call is required.
