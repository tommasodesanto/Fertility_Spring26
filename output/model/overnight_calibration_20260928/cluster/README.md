# September 28 overnight cluster preparation

Preparation only. No numerical job is authorized by these files alone.
New stage: `/scratch/td2248/projects/fertility_night_calibration_20260928_v1`.
`stage_night_remote.sh` physically copies the authenticated evening dependencies
on Torch, excluding accrued search trees, then adds the selected block0347 case,
its two fresh repeats, and selected export. The old stage is preserved.

New wrappers are `code/cluster/submit_e5f_night_calibration.sh` and
`code/cluster/submit_e5f_night_gated.sh`. Both passed `bash -n` locally without
model imports. They request 24 CPUs/192GB, single library threads, explicit JIT,
and the exact original project-path Apptainer bind. Override Slurm wall time
with the ceiling of remaining minutes to epoch1790596800 (08:00 EDT September28).
Search cutoff is1790593200 (07:00 EDT); waiting for approval never extends it.
The six-smoke gate must pass full source/target/parameter/plot review before the
lead writes an approval file. The launcher never writes its own approval.

The agent's ordinary SSH attempt was denied by its restricted execution
sandbox: control socket `Operation not permitted`, followed by unavailable host
resolution, exit255. No remote staging or tests were performed by this agent,
and no credential or connectivity workaround was attempted. The coordinating
lead can execute the saved remote script from its authorized SSH environment.

Pending: controller/builder/test freeze, their narrow remote refresh, manifest
pins, focused tests in a real Slurm allocation, and zero-solve `prepare`.
The new controller confirms `EXPECTED_E5F_NIGHT_SHA256`.
No unverified placeholder hash belongs in a submission.

Gated submission argument order remains:
`STAGE CONTRACT CONTRACT_SHA CONTROLLER CONTROLLER_SHA SEARCH_CUTOFF RUN_ROOT INNER_WRAPPER_SHA`.
Inner wrapper argument order remains:
`STAGE CONTRACT CONTRACT_SHA CONTROLLER CONTROLLER_SHA SEARCH_CUTOFF OUTPUT MODE APPROVAL APPROVAL_SHA`.

All three lanes start from the same authenticated selected block0347 parameters;
The exact evening normalization guess and step are retained; the child-benefit
level remains normalized for every proposal.
Scientific runtime, bounds, targets, weights and gates remain unchanged.

## Resumed preparation after access restoration

The ordinary hostname probe succeeded, but this agent's subsequent streamed SSH
staging command was again denied by its sandbox. The lead therefore owns the
actual remote copy. Do not run another copy over an existing stage.

After the controller agent confirms final freeze, the lead may run
`bash output/model/overnight_calibration_20260928/cluster/refresh_and_test.sh`.
This transfers only three new night tools, two wrappers and the test script;
records their hashes on Torch; and submits a five-minute, one-CPU, 8GB Slurm
preflight. The preflight checks the transfer manifest, runs the synthetic night
controller tests (no household solves), builds the immutable contract with
start epoch1790563080, and invokes the controller's zero-solve prepare stage.
Both scripts pass shell syntax checks. Results remain pending until actual
Slurm receipts and complete logs are collected.
