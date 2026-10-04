#!/usr/bin/env bash
# Explicit lead-only production release; fails closed until chain-0 smoke passes.
set -euo pipefail
remote=/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v3; cd "$remote"
[[ ! -e submission_receipt.json ]] || { echo "duplicate production submission refused"; exit 2; }
/share/apps/anaconda3/2025.06/bin/python verify_stage.py --host
/share/apps/anaconda3/2025.06/bin/python verify_smoke_gate.py "$remote"
available_gib=$(df -BG "$remote" | awk 'NR==2 {gsub(/G/,"",$4); print $4}')
[[ "$available_gib" -ge 2200 ]] || { echo "disk reserve below 2200 GiB"; exit 2; }
job=$(sbatch --parsable --array=0-3%4 --time=06:00:00 --export=ALL,CES_RUN_MODE=production launch_torch.sh); job=${job%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$job" <<'PY'
import hashlib,json,sys,time
from pathlib import Path
Path('submission_receipt.json').write_text(json.dumps(dict(status='production_submitted',job_id=sys.argv[1],array='0-3%4',account='torch_pr_570_general',chains=4,cpus_per_chain=1,memory_GiB_per_chain=24,wall_seconds_per_chain=21600,maximum_objective_calls_per_chain=500,final_native_reserve_seconds=1800,inventory_sha256=hashlib.sha256(Path('inventory.json').read_bytes()).hexdigest(),submission_epoch=time.time(),no_auto_retry=True),indent=2)+'\n')
PY
