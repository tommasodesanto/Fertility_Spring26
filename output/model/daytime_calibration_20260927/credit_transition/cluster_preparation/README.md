# Native transition cluster preparation

September 27: ordinary `code/cluster/torch.sh status` authentication passed;
Torch queue was empty. No numerical job was submitted and no remote source
was overwritten. Existing daytime remote contract has a different path-relocated
SHA; it cannot replace the native runtime's locally authenticated contract.

## Byte-preserving route

Apptainer is available on Torch. A read-only probe successfully mounted a scratch
directory at the exact original project path
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26` inside
`/share/apps/images/ubuntu-24.04.4.sif`. This avoids rewriting any contracts,
approvals, source literals or checkpoint bytes. With `module load
anaconda3/2025.06`, container execution of
`/share/apps/anaconda3/2025.06/bin/python` successfully imported NumPy 2.1.3,
SciPy 1.15.3, Numba 0.61.0 and Matplotlib 3.10.0. These Linux package versions
are a platform change; numerical cross-host checks are still required.

`code/cluster/prepare_e5f_native_transition_cluster.py` accepts an explicit
lead-pinned plan and fresh output directory. It writes exact file SHA/size
inventory; `--archive` additionally writes a tar preserving project-relative
paths. No simulation or submission occurs. Example (replace PLAN and NEW):

```
code/model/.venv/bin/python code/cluster/prepare_e5f_native_transition_cluster.py --plan PLAN --output NEW --archive
```

The existing native smoke plan inventory probe collected 752,215,451 bytes,
well below the 73 GiB portable tree; source/contract bytes are unchanged.
`inventory_probe_v1/manifest.json` is preparation evidence, not a root launch.
The first discovery attempt followed arbitrary historical JSON path strings and
exceeded the explicit 2GiB cap, failing before output/staging. The corrected
collector follows authenticated path/SHA records, source-pin keys and receipt
links, not arbitrary historical output references. This failure caused no
numerical run and remains documented here.

## Next bounded gate

Generate a fresh package from the final root plan (after local exact-loop
verification), transfer to a new scratch directory, verify every inventory hash,
then bind its project tree at the original absolute root. Run `setup` only first,
with one thread and a hard preparation timeout; fail on missing dependency or
pin mismatch. Only after that should a bounded cluster exact-loop smoke be
submitted with immutable source and explicit time/case caps. No full-transition
portability or numerical equivalence is certified by this packaging work.

Baseline source receipts include Mac `.venv` pure-Python files because they are
inside project `code/`. The package retains any pinned files as inert evidence;
Linux runtime imports use the external Anaconda environment. Do not execute a
Mac virtual environment on Torch or assert Linux library byte equivalence.

## Root-plan staging attempt

`root_smoke_package_v1/` packages the lead's root smoke plan SHA
`6991f512768e06f5f9bb87c33cde7c456b1540b88047dd22e515827e9cfc0126`.
Archive transfer is in progress to the newly created scratch directory
`/scratch/td2248/projects/fertility_native_transition_20260927_v1`.
Do not run this plan numerically: its local deadline is historical.

After transfer completes (manifest is copied last), run this *preparation only*
command using ordinary SSH:

```
ssh torch 'module load anaconda3/2025.06; base=/scratch/td2248/projects/fertility_native_transition_20260927_v1; python "$base/verify_staged_inputs.py" "$base" && bash "$base/run_in_environment.sh" "$base/native_setup_only.py" "$base/setup_only_v1"'
```

`verify_staged_inputs.py` refuses an existing extracted project and checks all
file hashes after extraction. `native_setup_only.py` sets a 180-second hard
cap, imports authenticated native runtime and loads the baseline checkpoint,
but performs zero numerical solves. The preparation result is written remotely
under `setup_only_v1/`; any failure remains in `failure.txt`. The command is
not yet executed at this note's writing. Wrapper uses external Linux Python,
single numerical threads and a separate Numba cache. A future numerical job
requires a NEW deadline/plan and independent source/hash review.

Staging completed: all **2,689** file hashes verified on Torch; native setup authenticated successfully with zero numerical solves. Collected receipts: `setup_only_v1/complete.json`, `environment.json`, `verification.json`. No auth/path blocker remains for this frozen dependency closure. Lead separately found a strict initial-distribution projection in the local root; this frozen root is **not approved** for a numerical transition. New corrected source/approval/plan must be separately pinned and staged before numerical execution.
