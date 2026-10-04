# V8 isolated output-routing package

V8 was staged and its zero-solve native preflight failed because the new handoff hash was not propagated to `identity.reference_sha256`. No model job was submitted. The immutable failed stage and receipt are preserved in `deployment_v8/`.

- Parent v7 inventory SHA-256: `5f1ee713bd1c5ab62d2a1776667f6191d61e69794cfe9a1fd22e6a6841397e92`.
- V8 inventory: `deployment_v8/inventory.json`, SHA-256 `335ee15a617b0d1c20c483ad231e8ce0d3f4d09977aaceb97cefff3239ee44e2`.
- V8 archive: `deployment_v8/stage.tar.gz`, SHA-256 `24427fcdfd094b89bee5ada6098cb2149b25847256cb111a80dac61e2d606e09`.
- Remote root: `/scratch/td2248/projects/current_estate_transition_20261003_v8`.

The detached runtime delta is exactly the reviewed four lines in `floor_runtime.py`: retain the old inherited-evidence directory, route the copied parameter object to `folder/inherited_state_evidence`, and record `output_override.json` with `economic_change=false`. Its source SHA changes from `94914395f158ed58e0d50943708270d74525a875b880ddddf09e8b91cfe09880` to `379a3f31f183cdf01ddb0a7a9f7127662c41f43f7654cc26e7154e49b781d0bb`.

The only staged source/data differences are that file and new `plans_v8/{handoff,smoke_plan,fit_plan,panel_config}.json`. They refresh the source, handoff, plan, and panel pins; the other 704 v7 source/data hashes are unchanged. The stage-local launchers only retarget remote output paths. Fit settings retain 12 fit evaluations, 12 path iterations, endpoint 48/1,800 seconds, 20,000 calls, 21,480 seconds and six hours. The grid is unchanged; intended release is only index 4 (`psi=0.164607063007583`).

## Local verification

```sh
python3 code/cluster/estate_birth_transition/prepare_output_fix.py verify
output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python output/model/transition_readiness_v1/current_baseline_20261003/work/failure_output_fix_v1_mock.py
git diff --no-index -- output/model/transition_readiness_v1/current_baseline_20261003/deployment_v7/source/code/model/experiments/transition_readiness/floor_runtime.py output/model/transition_readiness_v1/current_baseline_20261003/deployment_v8/source/code/model/experiments/transition_readiness/floor_runtime.py
```

The verifier authenticates v7 before hardlink-copying, confirms the target was unlinked before write, checks every v8 hash plus all transitive pins, and rejects numerical/economic plan drift. The five zero-solve mock cases passed: writable redirected evidence for retained-below-tolerance and rejected cases (positive), and stale read-only routing reproducing errno 30 (negative). Its active baseline hash equals v7 and the v8 source diff is exact. Identity-negative check:

```sh
python3 - <<'PY'
import copy,json
from pathlib import Path
d=Path('output/model/transition_readiness_v1/current_baseline_20261003/deployment_v8')
i=json.loads((d/'inventory.json').read_text()); b=copy.deepcopy(i)
b['identity']['source_pins']['code/model/experiments/transition_readiness/floor_runtime.py']='0'*64
assert b['identity'] != i['identity']
p=json.loads(Path('output/model/transition_readiness_v1/current_baseline_20261003/plans_v8/fit_plan.json').read_text()); q=copy.deepcopy(p); q['handoff']['sha256']='0'*64
assert q['handoff'] != p['handoff']; print('PASS identity-negative mutations detected')
PY
```

## Lead commands (not run)

```sh
stage=output/model/transition_readiness_v1/current_baseline_20261003/deployment_v8
remote=/scratch/td2248/projects/current_estate_transition_20261003_v8
ssh -o BatchMode=yes torch "test ! -e '$remote' && mkdir -p '$remote/logs' && cp -a --reflink=auto /scratch/td2248/projects/current_estate_transition_20261003_v7/source '$remote/source'"
scp "$stage/stage.tar.gz" "torch:$remote/stage.tar.gz"
ssh -o BatchMode=yes torch "cd '$remote' && tar -xzf stage.tar.gz && /share/apps/anaconda3/2025.06/bin/python deploy.py verify --stage '$remote'"
ssh -o BatchMode=yes torch "cd '$remote' && sbatch --parsable --partition=cs --time=01:00:00 --nodes=1 --ntasks=1 --export=ALL,TRANSITION_MODE=smoke,TRANSITION_WALL_SECONDS=3600 panel_launch_torch_v7.sh"
```

The native exact-loop smoke is mandatory; v5/v7 smoke cannot bridge the changed runtime identity. Lead must inspect the v8 smoke receipt (PASS, v8 identity, fresh seed, native setup/root/replay/accounting gates) and verify that `output_override.json` is candidate-owned. Then create and copy a fresh `fit_review_gate.json` binding the v8 `identity`, `plans`, inventory SHA and smoke-receipt SHA, with `reviewed_by_lead: true` and a reviewer. Only then release index 4:

```sh
ssh -o BatchMode=yes torch "cd '$remote' && sbatch --parsable --partition=cs --time=06:00:00 --nodes=1 --ntasks=1 --array=4%1 --export=ALL,TRANSITION_MODE=panel,TRANSITION_WALL_SECONDS=21600 panel_launch_torch_v7.sh"
```

The completed first candidate is preserved, but there is no supported serialized-controller resume argument. The later exception occurred before accounting a completed record. Thus exact tested continuation is unavailable; only a fresh post-smoke index-4 restart is documented.

## Corrected v9 metadata

The verified correction is in `deployment_v9/` and `plans_v9/`, built by `code/cluster/estate_birth_transition/prepare_output_fix.py`. Its inventory SHA-256 is `15cc5f6036d57b985ea2ea9098f015dd527d2ce11ec94bd458a49913f017d977`. The four-line runtime patch is byte-identical to v8; the complete smoke and fit plans differ only in refreshed reference identity and handoff path. The new reference hash is the actual handoff SHA-256. `deployment_v9/lead_verification.json` records the lead check. Actual remote preflight and fresh smoke are still required before only index 4 is released. V8 commands above are historical preparation evidence, not commands to rerun.

V9 was staged and passed host/container/native-constructor preflight (zero native calls, eight threads, exact identity). Smoke job **19146797** was submitted with a one-hour wall time, one `cs` node, eight CPUs and 48 GiB at 08:00 UTC; it was pending at the first scheduler check. Authoritative receipts are in `deployment_v9/`; never resubmit after this receipt. The copied deployer's `review-smoke` only implements the old cross-stage cap-change bridge, and its panel submission would select the old broad subset: neither is the v9 release route. After independently validating the actual v9 smoke and unchanged panel tests, record a fresh same-stage review gate (exact identity/plans/inventory/smoke receipt hashes), verify it locally and remotely with `verify-gate`, and use an atomic submission guard for **array index 4 only**. Gate acceptance must compare the actual remote smoke receipt. Preserve old arrays and do not restore any old numerical state.
