The launcher reuses the calibrated frozen-root, grid-resolution base overlays, continuation source overlays, and authenticated input bundle. It adds only an inventoried transition source and verified calibration handoff overlay. It never writes the calibration deployment.

`stage.py` snapshots the transition source, handoff, verified selected packet, source manifest, and any prepared plan and annual/blocks inputs into `stage/`. Upload that directory to `/scratch/td2248/projects/transition_readiness_v1/current_floor/` with `rsync -az`. The inventory authenticates every uploaded source before container entry. Rebuild the snapshot after any source or plan change; do not mutate an inventory used by a running job.

Zero-model-call target-architecture preflight (run on the Torch login node):

```sh
bash /scratch/td2248/projects/transition_readiness_v1/current_floor/floor_launch.sh --mode preflight --seconds 300 --label native_import_v1
```

After copying its `preflight.json` here as `native_import_preflight.json`, `make_plan.py` constructs the lead-specified diagnostic smoke plan. It takes the identity from the actual zero-call constructor, actual selected preference from the parameter table, and the retained complete annual target contract. It does not submit a job.

The lead must review and submit the final inventoried smoke plan. Pass `--time=01:31:00` and an explicit log output to `sbatch`, then use launcher arguments `--mode smoke --seconds 5400 --label smoke_v1 --plan /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_floor_preparation/deployment/smoke_plan.json`. Allocation is 4 CPUs/96 GiB; native numerical libraries remain one thread. The 90-minute actual-start timeout, per-stage controller limits, iteration limits, and 260 actual policy-call cap stop work without automatic retry. The launcher writes start/terminal receipts and a heartbeat every five minutes; native controller checkpoints preserve completed reference/seed inputs for explicit later reuse.

The short six-date smoke verifies the changed-shock equilibrium and fresh replay. Its terminal comparison remains diagnostic; it cannot certify the full 104/128-date transition or production readiness. No model job is authorized by running these preparation scripts alone.
