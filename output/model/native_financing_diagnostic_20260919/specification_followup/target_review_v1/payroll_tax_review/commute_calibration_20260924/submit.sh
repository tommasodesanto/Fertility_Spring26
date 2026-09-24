#!/usr/bin/env bash
set -euo pipefail
bundle=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/commute_calibration_20260924_v1
run_root="$bundle/results/run_001"
python3 - "$bundle/manifest.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
manifest=json.load(open(sys.argv[1]))
assert manifest['status']=='ready_for_reviewed_launch'
for name,item in manifest['files'].items():
    p=Path(item['path']); assert p.is_file(),name
    assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'],name
assert manifest['native_source_manifest_sha256']=='76406fcc10206d9e30bcc29d4e18accdf7fffdd9219e50e6e11c1c3360f01336'
assert manifest['checkpoint_sha256']=='83a28e46b36e2fbe30338d366611f3ec209f0c5a68309ee4ee9fa8523b66adee'
assert manifest['objective_sha256']=='c28e5d620d463dc3a592c038ca2771b47bf37ed1a9cf5ba31f7dc793b58b61db'
PY
test ! -e "$run_root"
if squeue -u td2248 -h -o '%j' | grep -q '^e5f-commute-'; then
  echo "Existing commute Slurm job found; refusing duplicate submission" >&2
  exit 2
fi
mkdir -p "$run_root"
smoke=$(sbatch --parsable --account=torch_pr_570_general --partition=cs \
  --ntasks=1 --cpus-per-task=1 --mem=16G --time=01:00:00 \
  --job-name=e5f-commute-smoke --output="$run_root/slurm_smoke_%j.out" \
  --error="$run_root/slurm_smoke_%j.err" \
  "$bundle/run_task.sh" smoke "$run_root")
production=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$smoke" \
  --account=torch_pr_570_general --partition=cs --array=1-8%8 \
  --ntasks=1 --cpus-per-task=1 --mem=16G --time=01:00:00 \
  --job-name=e5f-commute-worker --output="$run_root/slurm_worker_%A_%a.out" \
  --error="$run_root/slurm_worker_%A_%a.err" \
  "$bundle/run_task.sh" worker "$run_root")
export_job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterany:$production" \
  --account=torch_pr_570_general --partition=cs \
  --ntasks=1 --cpus-per-task=1 --mem=16G --time=00:10:00 \
  --job-name=e5f-commute-export --output="$run_root/slurm_export_%j.out" \
  --error="$run_root/slurm_export_%j.err" \
  "$bundle/run_task.sh" export "$run_root")
printf 'smoke=%s\nproduction=%s\nexport=%s\nrun_root=%s\n' \
  "$smoke" "$production" "$export_job" "$run_root"
