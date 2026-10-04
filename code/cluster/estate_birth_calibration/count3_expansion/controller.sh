#!/usr/bin/env bash
# Requires --dependency=afterok:19136605; one native smoke before five production tasks.
#SBATCH --job-name=estatec3ctl
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/estate_birth_count3_expansion_20261004_v1/logs/%x-%j.out
set -euo pipefail
stage=/scratch/td2248/projects/estate_birth_count3_expansion_20261004_v1
parent=/scratch/td2248/projects/estate_birth_calibration_20261003_v3
base=/scratch/td2248/projects/estate_birth_continuation_20261004_v1
python=/share/apps/anaconda3/2025.06/bin/python
cutoff=1791122400
cd "$stage"
mkdir -p control
phase=initialize
terminal() {
  status=$?
  trap - EXIT
  "$python" - "$status" "$stage" "$phase" <<'PY'
import json,os,sys,time
from pathlib import Path
status,stage,phase=sys.argv[1:]
path=Path(stage)/'control/controller_terminal.json'
path.write_text(json.dumps(dict(exit_code=int(status),controller_job_id=os.getenv('SLURM_JOB_ID'),
    base_controller_job_id='19136605',phase=phase,epoch=time.time(),no_auto_retry=True),indent=2)+'\n')
PY
  exit "$status"
}
trap terminal EXIT
trap 'exit 143' TERM
[[ ! -e control/start_plan_receipt.json && ! -e control/production_submission.json ]] || {
  echo 'Refusing duplicate expansion controller'; exit 2;
}
phase=verify_stage
"$python" verify_stage.py --host > control/stage_verification.json
phase=verify_base_smokes
"$python" "$base/verify_stage.py" --host > control/base_stage_verification.json
"$python" "$base/verify_smoke_gate.py" "$base" > control/base_smoke_gate_recheck.json
phase=prepare_starts
"$python" prepare_expansion.py --parent "$parent" --base "$base" --stage "$stage" > control/prepare_stdout.json
phase=verify_starts
"$python" verify_plans.py "$stage" > control/plan_verification.json
phase=smoke_count3
SLURM_ARRAY_TASK_ID=0 ESTATE_RUN_MODE=smoke bash launch_torch.sh
phase=smoke_gate
"$python" verify_smoke_gate.py "$stage" > control/smoke_gate.json
phase=release_budget
now=$(date +%s)
(( now + 2700 < cutoff )) || { echo 'Insufficient time for search plus native reserve'; exit 124; }
wall_minutes=$(((cutoff-now+59)/60))
"$python" - "$stage" <<'PY'
import os,sys
s=os.statvfs(sys.argv[1])
free=s.f_bavail*s.f_frsize
assert free>=7800*1024**3,'Insufficient combined old+new scratch reserve (7800 GiB)'
PY
phase=submit_production
job=$(sbatch --parsable --dependency="afterok:${SLURM_JOB_ID:?}" --array=0-4%5 \
  --time="$wall_minutes" --export=ALL,ESTATE_RUN_MODE=production launch_torch.sh)
job=${job%%;*}
phase=write_submission_receipt
"$python" - "$job" "$stage" "$wall_minutes" <<'PY'
import hashlib,json,os,sys,time
from pathlib import Path
job,stage,wall=sys.argv[1:]
root=Path(stage)
plan=json.loads((root/'control/start_plan_receipt.json').read_text())
receipt=dict(status='conditionally_submitted_after_expansion_controller_success',job_id=job,
  parent_array_job_id='19127370',base_controller_job_id='19136605',controller_job_id=os.getenv('SLURM_JOB_ID'),
  dependency='afterok:'+os.getenv('SLURM_JOB_ID'),array='0-4%5',additional_cores_max=5,
  prior_production_cores_max=10,total_estate_production_cores_max=15,memory_GiB_per_task=24,
  maximum_objective_calls_per_task=500,final_native_reserve_seconds=1800,
  hard_stop_epoch=1791122400,requested_wall_minutes=int(wall),
  stage_inventory_sha256=hashlib.sha256((root/'inventory.json').read_bytes()).hexdigest(),
  start_plan_sha256=plan['plan_sha256'],
  smoke_gate_sha256=hashlib.sha256((root/'control/smoke_gate.json').read_bytes()).hexdigest(),
  submission_epoch=time.time(),no_auto_retry=True,no_scientific_adoption=True)
(root/'control/production_submission.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
PY
phase=complete
