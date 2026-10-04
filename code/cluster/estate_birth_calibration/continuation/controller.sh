#!/usr/bin/env bash
# Submit ONLY with --dependency=afterok:19127370; it gates and conditionally releases production.
#SBATCH --job-name=estatecontctl
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=03:15:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/estate_birth_continuation_20261004_v1/logs/%x-%j.out
set -euo pipefail
stage=/scratch/td2248/projects/estate_birth_continuation_20261004_v1
parent=/scratch/td2248/projects/estate_birth_calibration_20261003_v3
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
    parent_array_job_id='19127370',phase=phase,epoch=time.time(),no_auto_retry=True),indent=2)+'\n')
PY
  exit "$status"
}
trap terminal EXIT
trap 'exit 143' TERM
[[ ! -e control/gate/parent_gate.json && ! -e control/production_submission.json ]] || {
  echo 'Refusing duplicate controller execution'; exit 2;
}
phase=verify_stage
"$python" verify_stage.py --host
phase=parent_accounting
sacct -n -P -j 19127370 --format=JobID,JobIDRaw,State,ExitCode > control/parent_sacct.txt
phase=parent_gate
"$python" prepare_continuation.py --parent "$parent" --stage "$stage" \
  --receipt-root "$stage/control/gate" --sacct control/parent_sacct.txt > control/prepare_stdout.json
phase=verify_plans
"$python" verify_plans.py "$stage" > control/plan_verification.json
for task in 0 5; do
  phase="smoke_task_${task}"
  SLURM_ARRAY_TASK_ID="$task" ESTATE_RUN_MODE=smoke bash launch_torch.sh
done
phase=smoke_gate
"$python" verify_smoke_gate.py "$stage" > control/smoke_gate.json
phase=release_budget
now=$(date +%s)
(( now + 2700 < cutoff )) || { echo 'Insufficient time for search plus native reserve'; exit 124; }
wall_minutes=$(((cutoff-now+59)/60))
"$python" - "$stage" <<'PY'
import os,sys
s=os.statvfs(sys.argv[1]);assert s.f_bavail*s.f_frsize>=5200*1024**3,'Insufficient shared scratch plan reserve'
PY
# The new array is eligible only after this controller exits successfully.
phase=submit_production
job=$(sbatch --parsable --dependency="afterok:${SLURM_JOB_ID:?}" --array=0-9%10 \
  --time="$wall_minutes" --export=ALL,ESTATE_RUN_MODE=production launch_torch.sh)
job=${job%%;*}
phase=write_submission_receipt
"$python" - "$job" "$stage" <<'PY'
import hashlib,json,os,sys,time
from pathlib import Path
job,stage=sys.argv[1:]
root=Path(stage)
gate=json.loads((root/'control/gate/parent_gate.json').read_text())
receipt=dict(status='conditionally_submitted_after_controller_success',job_id=job,
  parent_array_job_id='19127370',controller_job_id=os.getenv('SLURM_JOB_ID'),
  dependency='afterok:'+os.getenv('SLURM_JOB_ID'),array='0-9%10',cores_max=10,
  memory_GiB_per_task=24,maximum_objective_calls_per_task=500,
  final_native_reserve_seconds=1800,hard_stop_epoch=1791122400,
  stage_inventory_sha256=hashlib.sha256((root/'inventory.json').read_bytes()).hexdigest(),
  parent_inventory_sha256=gate['parent_inventory_sha256'],
  parent_gate_sha256=hashlib.sha256((root/'control/gate/parent_gate.json').read_bytes()).hexdigest(),
  smoke_gate_sha256=hashlib.sha256((root/'control/smoke_gate.json').read_bytes()).hexdigest(),
  submission_epoch=time.time(),no_auto_retry=True,no_scientific_adoption=True)
(root/'control/production_submission.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
PY
phase=complete
