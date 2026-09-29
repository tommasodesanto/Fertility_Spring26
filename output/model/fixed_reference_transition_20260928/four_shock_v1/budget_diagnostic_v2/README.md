# Pension-correction diagnostic: reporting fix

Retry job **18820811** follows the author's request to address the failed
reporter. The v1 job stopped after 6m44s while saving its first six-date mapping;
it never tested the pension correction. Both attempts remain separately retained.

The source fix supports NumPy arrays/scalars in JSON and converts terminal arrays
to numeric lists before repeat comparison. Verification **18820759** passed 66
tests, including the real mapping/reporting/repeat loop with a fake native solver.
A separate saved-data regression reproduced the old writer crash and verified an
exact round trip with the fix; native queue arrays were restored from the saved
terminal receipt. That verification performed zero model solves.

The numerical diagnostic is unchanged: hold the reference preference and prices
fixed, propose pension × recorded revenue/outlays, and evaluate fresh household
policies and distributions. At most three six-date smoke mappings precede two
104-date mappings. Full work starts only after smoke succeeds. Original economics,
queues, endpoints, grid, cache and all gates are unchanged; this does not fit shocks.

Budgets: smoke 45 minutes; full stage four hours, three hours per mapping; Slurm
five hours. One numerical thread, four allocated CPUs, 96-GiB memory and 64-GiB
cache. See configs/ and launch_receipt.json for exact pins and budgets. The new
reporting_regression.json must match the executed driver before native work.

Evidence: verification_tests.log, reporting_regression.json, run.sh and
verify_reporting.sh; later smoke/ and full/ receipts live on Torch. The monitor
checks every thirty minutes, read-only, and stops on job termination. No automatic
repair/restart/estimation. Success still requires longer-horizon and changed-
preference validation before full estimation can resume.
