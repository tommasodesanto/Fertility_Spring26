#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/utility_floor_round2_v1
cd "$remote"
[[ ! -e submission_receipt.json && ! -e submission_rows.tsv ]] || { echo 'Refusing duplicate submission'; exit 2; }
[[ $(date +%s) -lt 1790833016 ]] || { echo 'Absolute deadline reached'; exit 124; }
: > submission_rows.tsv
for arm in floor; do
  smoke=$(sbatch --parsable --time=00:30:00 --export="ALL,UTILITY_STAGE=smoke,UTILITY_ARM=$arm" launch_torch.sh)
  smoke=${smoke%%;*}
  printf '%s\tsmoke\t%s\n' "$arm" "$smoke" >> submission_rows.tsv
  parent="$remote/results/smoke/${arm}_s0/smoke"
  search=$(sbatch --parsable --array=0-7 --time=04:00:00 --dependency="afterok:$smoke" --kill-on-invalid-dep=yes --export="ALL,UTILITY_STAGE=search,UTILITY_ARM=$arm,SMOKE_PARENT=$parent" launch_torch.sh)
  search=${search%%;*}
  printf '%s\tsearch\t%s\n' "$arm" "$search" >> submission_rows.tsv
  echo "$arm smoke=$smoke search=$search"
done
python - <<'PY'
import csv,json,time
from pathlib import Path
rows=[]
with open('submission_rows.tsv') as stream:
 for arm,stage,job in csv.reader(stream,delimiter='\t'):
  rows.append(dict(arm=arm,stage=stage,job_id=job))
r=dict(status='submitted',submission_epoch=time.time(),common_deadline_epoch=1790833016,individual_smokes=1,dependent_search_arrays=1,total_search_lanes=8,rows=rows,automatic_restarts=False)
p=Path('submission_receipt.json');t=p.with_suffix('.tmp');t.write_text(json.dumps(r,indent=2)+'\n');t.replace(p)
PY
