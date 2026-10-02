#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/transition_readiness_v1/normalized_resumed_fit_v2
/share/apps/anaconda3/2025.06/bin/python - "$remote" <<'PY'
import json,hashlib,sys
from pathlib import Path
r=Path(sys.argv[1]);i=json.loads((r/'inventory.json').read_text());assert hashlib.sha256((r/'inventory.json').read_bytes()).hexdigest()=='4c99010afa74f22b4af20c60a7012fe79727183f3e6f76660ae1c358729e8494'
for rel,h in i['files'].items():assert hashlib.sha256((r/'source'/rel).read_bytes()).hexdigest()==h,rel
assert hashlib.sha256((r/'floor_launch.sh').read_bytes()).hexdigest()==i['files']['code/model/experiments/transition_readiness/floor_launch.sh']
p=r/'source/output/model/transition_readiness_v1/normalized_restart_v1/resume_preparation/deployment/v2/resume_plan.json';assert hashlib.sha256(p.read_bytes()).hexdigest()=='12b7088ba4d6ffef560e302a698a894bb4997c9fed808918b171d47fe7c1b2ff';plan=json.loads(p.read_text())
a=json.loads((r/'results/restore_smoke_v2/restore_receipt.json').read_text());assert a['status']=='actual_reference_and_measured_J_restored' and a['native_calls']==0 and a['matrix_shape']==[24,24] and a['reference_verified'];assert a['identity']==plan['identity'] and a['controller']==plan['source_files']['controller'] and a['runtime']==plan['source_files']['runtime']
print('PASS frozen metadata restore inventory',len(i['files']))
PY
deadline=1790911128;remaining=$((deadline-$(date +%s)))
[[ "$remaining" -gt 120 ]] || { echo 'Original resume deadline exhausted'; exit 124; }
mkdir "$remote/submission.lock"
slurm_minutes=$(((remaining+59)/60+1))
job=$(sbatch --parsable --time="$slurm_minutes" --output="$remote/results/normalized_resumed_fit_v2_slurm.out" "$remote/floor_launch.sh" --mode fit --seconds "$remaining" --deadline-epoch "$deadline" --label normalized_resumed_fit_v2 --plan /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/normalized_restart_v1/resume_preparation/deployment/v2/resume_plan.json)
/share/apps/anaconda3/2025.06/bin/python - "$remote" "$job" "$deadline" <<'PY'
import json,os,sys,time
from pathlib import Path
p=Path(sys.argv[1])/'submission_receipt.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(dict(job_id=sys.argv[2],deadline_epoch=int(sys.argv[3]),submitted_epoch=time.time(),sole_submitter='lead',maximum_policy_calls=1214),indent=2)+'\n');os.replace(t,p)
PY
printf '%s\n' "$job"
