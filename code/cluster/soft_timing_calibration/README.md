# Matched soft-purchase timing calibration: Torch deployment

This package runs four independent chains under the original interest clock and
four under the experimental post-interest transaction clock. Both use the same
soft purchase rule, ten free coordinates, target and weight contracts, financed
share 0.8, entrant rule, earnings, floors and utility. It is an experimental
search; no fitted point becomes a paper baseline without a fresh native check.

The immutable source stage is built by `build_stage.py`. `stage_torch.sh` creates
the isolated `/scratch/td2248/projects/soft_timing_calibration_20261002_v1`
folder, uploads the source archive and runs a zero-solve preflight under the
existing authenticated Apptainer runtime. The preflight hashes every staged
source file, the original 234-file pin list, all four paired timing files, the
selected soft checkpoint, the target and weight contracts, then constructs the
native evaluator once for each arm without solving. It relies on the existing
frozen repository and prior authenticated base/floor overlays. In particular,
the frozen remote copies of `e5f_exact_policy_cache.py` and its test have the
required historical hashes; the current checkout's newer copies must not be
overlaid.

The first staged archive, SHA-256
`a10ec1dff0b22c26535d30450882736bf71ae9aba7b6020a6baf41b83694ae8f`,
passed the host 325-file inventory check but failed before Apptainer because
`launch_torch.sh` transposed two characters in the expected hash of the
unchanged `arrays.npz` input. Its remote archive is preserved as
`stage_attempt1.tar.gz`. The corrected local archive is recorded in
`output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/deployment/stage_receipt.json`.
At the last update, automatic approval review rejected uploading this corrected
archive, so **remote preflight and both native smokes remain pending**. There
are no Slurm jobs from this deployment and `submit_torch.sh` has not run.
The corrected package excludes two incidental `.DS_Store` files and hashes 323
source and contract files.
Once that upload is authorized, `bash code/cluster/soft_timing_calibration/restage_torch.sh`
is the exact guarded recovery command. It requires the preserved first archive,
no submission receipt and no result directory before replacing the isolated
stage, then runs the two-arm zero-solve preflight.

After the corrected archive is approved and remote preflight passes, a lead may
submit only the two smoke tasks with:

```bash
ssh torch 'cd /scratch/td2248/projects/soft_timing_calibration_20261002_v1 && sbatch --parsable --array=0,4 --time=01:30:00 --export=ALL,SOFT_RUN_MODE=smoke launch_torch.sh'
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
