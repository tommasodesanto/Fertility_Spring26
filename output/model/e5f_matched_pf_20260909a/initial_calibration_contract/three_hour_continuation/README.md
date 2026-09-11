# Autonomous continuation after housing-gate candidate rejections

Author explicitly requests the job continue after a single candidate failure.
Previous17375937 finished24 candidates in8m16s:22 passed,2 failed the unchanged
initial housing-equilibrium gate. The previous controller stopped all adaptive
stages; no refinement or final selected repetitions ran. Its results are preserved
under `../three_hour_refinement/collected_17375937/` and remotely in that batch.

This separate continuation starts from the saved best candidate md_joint_10.
It first reproduces that candidate twice through the unchanged scored loop, then
runs up to2 rounds of18 derivative probes plus12 jointly adjusted proposals,
followed by two selected repetitions. Maximum62 candidate jobs/64repetitions,
512GE at the existing8GE cap per repetition; expected45–90minutes excluding
queue.24CPUs/192GiB, one numerical thread per child;3h Slurm hard cap.
Search has7800seconds and final verification/reporting reserves3000seconds.

Only the exact source-verified error `Initial housing equilibrium failed its
unchanged strict gate`, raised in stationary_equilibrium with a completed
preflight, becomes an inadmissible candidate. It receives no score and does not
stop other proposals. Missing preflight, source changes, unexpected exceptions,
observer/target errors and accounting failures remain fatal. A derivative may
use a saved one-sided difference if just one side fails; both sides failing
stops that derivative stage without inventing an observation. More than half
of a batch failing equilibrium also stops adaptive search for review.

Eight tests verify the parallel loop, transformed bounds, full residual vector,
failure preservation, exception propagation and exact failure classification.
Original27wrapper/scorer startup tests remain. Every candidate verifies all641
source files and the complete fixed early target/weight fingerprint. No model,
numerical tolerance, parameter bound, target, weight or measurement change.

The entire controller runs on Torch. No powered-on laptop or AI heartbeat is
required. All completed-case/best-so-far summaries,30second heartbeat, model
checkpoints, full fit and parameter tables and standard17graphs are retained.
The selected result remains provisional unless its exact repetitions pass;
this is an initial calibration search, not a historical/policy certification.

Remote:
`/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/three_hour_continuation_20260911/`.
Immutable configuration is in plan.json; job identifier in submission.json.
