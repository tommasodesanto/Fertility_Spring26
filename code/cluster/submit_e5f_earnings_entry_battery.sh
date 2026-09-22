#!/usr/bin/env bash
# Stage-check or submit four smokes followed by cell-matched 10-worker arrays.
set -euo pipefail

project_root="$(cd "$(dirname "$0")/../.." && pwd)"
local_bundle="${E5F_EB_LOCAL_BUNDLE:?set E5F_EB_LOCAL_BUNDLE to the local prepared bundle directory}"
remote_bundle="${E5F_EB_REMOTE_BUNDLE:?set E5F_EB_REMOTE_BUNDLE to its absolute shared-filesystem path on Torch}"
local_manifest="$local_bundle/manifest.json"
remote_manifest="$remote_bundle/manifest.json"
local_driver="$local_bundle/tools/run_e5f_earnings_entry_battery.py"
remote_driver="$remote_bundle/tools/run_e5f_earnings_entry_battery.py"
local_launcher="$project_root/code/cluster/submit_e5f_earnings_entry_battery.sh"
remote_launcher="$remote_bundle/tools/submit_e5f_earnings_entry_battery.sh"
python_path="/share/apps/anaconda3/2025.06/bin/python"
mode="${1:-stage}"
output_root="${E5F_EB_OUTPUT_ROOT:?set E5F_EB_OUTPUT_ROOT to a new absolute output directory}"
manifest_sha="${E5F_EB_MANIFEST_SHA256:?set E5F_EB_MANIFEST_SHA256 to the pinned manifest SHA-256}"
driver_sha="${E5F_EB_DRIVER_SHA256:?set E5F_EB_DRIVER_SHA256 to the pinned worker SHA-256}"
launcher_sha="${E5F_EB_LAUNCHER_SHA256:?set E5F_EB_LAUNCHER_SHA256 to the pinned launcher SHA-256}"

validate_local() {
  python3 - "$local_manifest" "$manifest_sha" "$local_driver" "$driver_sha" "$local_launcher" "$launcher_sha" <<'PY'
import hashlib, json, pathlib, sys
mp, mh, dp, dh, lp, lh = sys.argv[1:]
def sha(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
assert sha(mp) == mh, 'manifest hash mismatch'
assert sha(dp) == dh, 'worker hash mismatch'
assert sha(lp) == lh, 'launcher hash mismatch'
m = json.load(open(mp))
assert m.get('schema') == 'e5f_earnings_entry_battery_manifest_v1' and m.get('status') == 'ready'
assert len(m.get('cells', {})) == 4 and len(m.get('smoke', {}).get('cases', [])) == 4
assert len(m.get('production', {}).get('cases', [])) == 40
for cell in m['cells']:
    assert sum(x['cell_id'] == cell for x in m['smoke']['cases']) == 1
    assert sum(x['cell_id'] == cell for x in m['production']['cases']) == 10
print('\n'.join(m['cells']))
PY
}

remote_check() {
  ssh torch "test -x '$python_path' && test -r '$remote_manifest' && test -r '$remote_driver' && test -r '$remote_launcher' && test \"\$(sha256sum '$remote_manifest' | awk '{print \$1}')\" = '$manifest_sha' && test \"\$(sha256sum '$remote_driver' | awk '{print \$1}')\" = '$driver_sha' && test \"\$(sha256sum '$remote_launcher' | awk '{print \$1}')\" = '$launcher_sha' && if [ -d '$output_root' ]; then test -z \"\$(find '$output_root' -mindepth 1 -maxdepth 1 -print -quit)\"; else mkdir -p '$output_root'; fi"
}

submit() {
  local cell smoke production
  while IFS= read -r cell; do
    smoke=$(ssh torch sbatch --parsable --account=torch_pr_570_general --partition=cpu_short \
      --ntasks=1 --cpus-per-task=1 --mem=32G --time=01:00:00 --job-name="e5feb-smoke" \
      --output="$output_root/slurm_smoke_${cell}_%j.out" --error="$output_root/slurm_smoke_${cell}_%j.err" <<SBATCH
#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
test "\$(sha256sum '$remote_manifest' | awk '{print \$1}')" = '$manifest_sha'
test "\$(sha256sum '$remote_driver' | awk '{print \$1}')" = '$driver_sha'
test "\$(sha256sum '$remote_launcher' | awk '{print \$1}')" = '$launcher_sha'
exec '$python_path' '$remote_driver' --manifest '$remote_manifest' --stage smoke --cell-id '$cell' --task-id 1 --output-root '$output_root'
SBATCH
    )
    production=$(ssh torch sbatch --parsable --dependency="afterok:$smoke" --account=torch_pr_570_general \
      --partition=cpu_short --array=1-10%10 --ntasks=1 --cpus-per-task=1 --mem=32G --time=01:00:00 \
      --job-name="e5feb-prod" --output="$output_root/slurm_production_${cell}_%A_%a.out" \
      --error="$output_root/slurm_production_${cell}_%A_%a.err" <<SBATCH
#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
test "\$(sha256sum '$remote_manifest' | awk '{print \$1}')" = '$manifest_sha'
test "\$(sha256sum '$remote_driver' | awk '{print \$1}')" = '$driver_sha'
test "\$(sha256sum '$remote_launcher' | awk '{print \$1}')" = '$launcher_sha'
exec '$python_path' '$remote_driver' --manifest '$remote_manifest' --stage production --cell-id '$cell' --task-id "\$SLURM_ARRAY_TASK_ID" --output-root '$output_root'
SBATCH
    )
    printf 'cell=%s smoke_job=%s production_array=%s\n' "$cell" "$smoke" "$production"
  done < <(validate_local)
}

if [[ "$mode" != stage && "$mode" != submit ]]; then
  echo "usage: $0 [stage|submit]" >&2
  exit 2
fi
validate_local >/dev/null
remote_check
if [[ "$mode" == submit ]]; then submit; fi
