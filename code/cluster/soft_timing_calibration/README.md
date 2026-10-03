# Matched soft-purchase timing calibration: Torch deployment

This package runs four independent chains under the original interest clock and
four under the experimental post-interest transaction clock. Both use the same
soft purchase rule, ten free coordinates, target and weight contracts, financed
share 0.8, entrant rule, earnings, floors and utility. It is an experimental
search; no fitted point becomes a paper baseline without a fresh native check.

The attempt-2 immutable source stage is built by `build_stage.py`. `stage_torch.sh` creates
the isolated `/scratch/td2248/projects/soft_timing_calibration_20261002_v2`
folder, uploads the source archive and runs a zero-solve preflight under the
existing authenticated Apptainer runtime. The preflight hashes every staged
source file, the original 234-file pin list, all four paired timing files, the
selected soft checkpoint, the target and weight contracts, then constructs the
native evaluator once for each arm without solving. It relies on the existing
frozen repository and prior authenticated base/floor overlays. In particular,
the frozen remote copies of `e5f_exact_policy_cache.py` and its test have the
required historical hashes; the current checkout's newer copies must not be
overlaid.

The author approved the first Torch upload. That archive and its inventory,
status, and failed smoke receipts are preserved locally under
`output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/deployment/attempt1/`.
The remote attempt-1 root `/scratch/td2248/projects/soft_timing_calibration_20261002_v1`
is immutable. Its smoke array `19085105` completed two exploratory calls per
arm, then both tasks failed when the full native postcheck was constructed in
the same interpreter. The repair moves postcheck into a fresh child interpreter.
The attempt-2 stage is prepared locally for lead review (SHA-256
`f7a8fec4ff370fd3690c0d0068ca595b75a17dd8aaac3bd47f6009ef73ecd68b`);
it has not been uploaded or submitted. The 323 source hashes differ from
attempt 1 only for `calibrate.py` and `driver_plan.json`. The target, weight,
selected-checkpoint and original 234-file source pins are unchanged. The
synthetic gate checks include the fresh child's search-receipt hash and
fast/full target and parameter comparison. This is a smoke repair, not a
production restart.

After the lead reviews the repaired driver and stage, run
`bash code/cluster/soft_timing_calibration/stage_torch.sh` once. If its host,
container, and two zero-solve evaluator checks pass, submit only the two smoke
tasks with:

```bash
ssh torch 'cd /scratch/td2248/projects/soft_timing_calibration_20261002_v2 && sbatch --parsable --array=0,4 --time=01:30:00 --export=ALL,SOFT_RUN_MODE=smoke launch_torch.sh'
```

Smoke task 0 is original chain 0 and task 4 is alternative chain 0. Each runs
the exact two-call optimizer loop and fresh native selected-point verification
under a 90-minute actual-start deadline, including the 30-minute native reserve.
Inspect both `results/smoke_<arm>_chain_0/run/completed.json`, launch terminal
receipts, checkpoints, standard 17 plots, and the target/parameter tables before
production. `submit_torch.sh` itself rejects a missing or failed smoke, mixed
fingerprints, incomplete checkpoints or reports, or an inventory changed since
smoke. The local `test_verify_smoke_gate.py` exercises those decisions with
clearly synthetic fixtures; no model values are used. Production requires
separate lead review and is submitted only by explicitly invoking
`submit_torch.sh` on Torch. It requests tasks 0–7, at most
eight concurrent, with one core, 24 GiB and six hours per chain. The driver
caps each chain at 250 objective calls and reserves 1,800 seconds for native
verification. There is no automatic repair, fallback or restart.

`collect_torch.py --mode production --root <local-output> --fetch` copies and
validates all eight chain folders, rejecting mixed target, weight or selection
fingerprints. It requires a passing native postcheck, 14 target rows, 31
parameter rows and 17 standard plots for each selected result. A result whose
optimizer lacks convergence remains provisional despite a passing numerical
postcheck.
