# Fixed-price pension-correction diagnostic

Torch job **18818674** tests one numerical correction after both launch_v3 fits
failed before fitting any shock. It is not a new estimation launch. All 64 pure
tests and both actual-config preflights passed, but the job FAILED after 6m44s
while writing the first smoke receipt: a NumPy array was not JSON serializable.
The pension correction and full stage were never reached. The reporter fix and
author-authorized retry are retained in the adjacent budget_diagnostic_v2 packet.

At the unchanged September 28 block0506 preference and fixed house-price path,
set each dated pension to its old value times recorded payroll revenue divided
by pension outlays. This balances the recorded old-distribution accounts; a fresh
native backward/forward solution must establish whether it also clears the
updated economy. The frozen economics, original distribution and both entry
queues, terminal endpoint, wealth grid and acceptance tolerances are preserved.

The only change to existing numerical source propagates TimeoutError through
policy-cache serialization. Cache keys and successful results are unchanged.
No short-horizon Jacobian is used in this diagnostic.

The gated job runs the pure tests and validates both actual configurations before
native work. It then runs a six-date smoke with a 0.001% pension-path perturbation,
one accounting correction and a fresh repeat. Only smoke success permits the
104-date correction and repeat. The full stage reuses the authenticated initial
mapping from launch_v3 as accounting input, not as a converged equilibrium.

Limits: three smoke mappings, two full mappings, at most 452 backward/forward
policy calls before caching; 45 minutes for smoke, four hours for the full stage,
three hours per full mapping, five hours Slurm. Prior full changed-path mappings
took 89–106 minutes; total expected runtime is approximately three to four hours,
subject to queue time and numerical outcome. Four CPUs are allocated for the
96-GiB memory request; numerical threads remain one and cache capacity 64 GiB.

The native gates remain housing 2e-4, pension 1e-6, terminal 1e-3, replay 1e-10.
The driver records latest/best errors, fertility changes and process heartbeats.
It stops at an uncertified trial. Even success still requires a longer-horizon
check and a changed-preference candidate before restarting estimation.

Small evidence: launch_receipt.json, configs/, tests.log, preflight.json, then
smoke/complete.json and full/complete.json or failure.json on Torch. The retained
run.sh invokes the source snapshot under source/code/ in the Torch packet, with
the original frozen project mounted read-only. The monitor checks every thirty
minutes, reports meaningful failures or full success, and stops at termination.
There are no automatic repairs, restarts, or estimation submissions.
