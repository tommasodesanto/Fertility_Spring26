# Matched estate-A calibration deployment

The October 4 20-chain, one-birth continuation is in
[`binary_continuation_v2/`](binary_continuation_v2/README.md). Its pinned
October 4 stage and job receipts are under
`output/model/experiments/birth_count_choice/estate_a_binary_continuation_20261004_v2/deployment/`.
The historical scripts below belong to the failed 19127370 array.

**Current status, October 4:** all ten tasks in array **19127370** failed at
the search time boundary before final native selected-point verification.
The saved per-case and best-so-far checkpoints remain provisional. See the
[failure diagnosis](../../../output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/failure_diagnosis_20261004.md).
The dependent continuation controllers cannot release production. No repair
or restart was made.

The isolated search uses the net estate `bp+(1-psi)*P*h_prime` in utility and deceased-estate accounting, the experimental PSID wealth target 4.45838713455674, and the original numerical weight. All other target values and weights remain unchanged. The SCF recipient and wealth-scope mismatch remains provisional. No production default or parameter file changes.

`prepare_plan.py` writes the chain-13 anchor, the preserved unverified new-target chain-6 checkpoint, and three deterministic bounded chain-13 perturbations shared by binary cap 1 and count cap 3. Start 1 is seed-only: paused array 19111687, chain 6, case 0060_nm, saved loss 22.141841386410267; it has no final native verification and is not adopted. Its target identity and source SHA are pinned. Tasks 0–4 are binary chains 0–4; tasks 5–9 are count-3 chains 0–4. Each task has one core, 24 GiB, six hours, at most 500 calls, and 1,800 seconds reserved for a fresh native check. Every proposal calls full native GE with the unchanged 32-lifecycle cap. The complete contract, timing estimate and conservative storage budget are in `output/model/experiments/birth_count_choice/estate_a_calibration_v1/start_plan.json`.

`build_stage.py` reuses the authenticated parent archive and packages the isolated engine, recovered exact-source dependencies, input snapshot, anchor receipt and entrypoints. `stage_torch.sh` stages and initializes both arms with zero solves. Rebuild only after the engine and workflow are frozen; once staged, use a fresh attempt location for any changed source.

After lead review, submit the two smoke tasks with `sbatch --array=0,5%2 --time=01:30:00 --export=ALL,ESTATE_RUN_MODE=smoke launch_torch.sh`. Smoke evaluates two optimizer calls per arm, then a fresh child native verification and an exact native repeat with 14 targets, 31 parameters and 17 plots. `verify_smoke_gate.py` fails closed unless both pass under the current inventory. Invoke `submit_torch.sh` only following explicit lead production release.

`collect_torch.py --stage REMOTE --mode smoke|production --out NEW_JSON` produces complete tables and the best verified result by arm while rejecting mixed fingerprints. Copy this compact JSON, not the retained solution-array tree, for ordinary readouts.

The cancelled arrays 19111687 and 19112020 and their checkpoints remain preserved. No automatic resume, retry, fallback or scientific adoption occurs.
