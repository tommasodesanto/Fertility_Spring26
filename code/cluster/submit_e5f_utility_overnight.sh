#!/usr/bin/env bash
# Stage/check or submit exact-loop smokes, forty persistent workers, four verifiers.
set -euo pipefail

mode="${1:-stage}"
manifest="${E5F_UTILITY_MANIFEST:?absolute local manifest required}"
manifest_sha="${E5F_UTILITY_MANIFEST_SHA256:?manifest hash required}"
output_root="${E5F_UTILITY_OUTPUT_ROOT:?new absolute shared output root required}"
root="$(cd "$(dirname "$0")/../.." && pwd)"
runner="$root/code/model/tools/run_e5f_utility_overnight.py"
launcher="$root/code/cluster/submit_e5f_utility_overnight.sh"
runner_sha="${E5F_UTILITY_RUNNER_SHA256:?runner hash required}"
launcher_sha="${E5F_UTILITY_LAUNCHER_SHA256:?launcher hash required}"
remote_runner="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["runner"]["path"])' "$manifest")"
remote_launcher="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["launcher"]["path"])' "$manifest")"
remote_manifest="$(dirname "$(dirname "$remote_runner")")/manifest.json"

validate_local() {
  python3 - "$manifest" "$manifest_sha" "$runner" "$runner_sha" "$launcher" "$launcher_sha" <<'PY'
import hashlib,json,pathlib,sys
manifest, mh, runner, rh, launcher, lh = sys.argv[1:]
def sha(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
assert sha(manifest)==mh, 'manifest hash mismatch'
assert sha(runner)==rh, 'runner hash mismatch'
assert sha(launcher)==lh, 'launcher hash mismatch'
m=json.load(open(manifest))
assert m['schema']=='e5f_utility_overnight_manifest_v1' and m['status']=='ready'
assert m['total_workers']==40 and m['entry_rule']=='fixed_reference_marginal'
assert m['runner']['sha256']==rh and m['launcher']['sha256']==lh
assert m['production_stop_utc']=='2026-09-23T12:00:00Z'
assert m['verification_stop_utc']=='2026-09-23T13:00:00Z'
PY
}

remote_check() {
  ssh torch "test -r '$remote_manifest' && test -r '$remote_runner' && test -r '$remote_launcher' && test \"\$(sha256sum '$remote_manifest' | awk '{print \$1}')\" = '$manifest_sha' && test \"\$(sha256sum '$remote_runner' | awk '{print \$1}')\" = '$runner_sha' && test \"\$(sha256sum '$remote_launcher' | awk '{print \$1}')\" = '$launcher_sha' && if [ -d '$output_root' ]; then test -z \"\$(find '$output_root' -mindepth 1 -maxdepth 1 -print -quit)\"; else mkdir -p '$output_root'; fi"
}

submit() {
  local cell smoke production verification finish prod_range verify_range
  for cell in 1 2 3 4; do
    case "$cell" in
      1) prod_range=1-18; verify_range=1-2 ;;
      2) prod_range=19-36; verify_range=3-4 ;;
      3) prod_range=37-38; verify_range=5-6 ;;
      4) prod_range=39-40; verify_range=7-8 ;;
    esac
    smoke=$(ssh torch sbatch --parsable --account=torch_pr_570_general --partition=cs \
    --array="$cell" --ntasks=1 --cpus-per-task=1 --mem=16G --time=03:00:00 \
    --job-name="e5futil-smoke-$cell" --output="$output_root/slurm_smoke_%A_%a.out" \
    --error="$output_root/slurm_smoke_%A_%a.err" <<SBATCH
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
exec python3 '$remote_runner' --manifest '$remote_manifest' --manifest-sha256 '$manifest_sha' --stage smoke --task-id "\$SLURM_ARRAY_TASK_ID" --output-root '$output_root'
SBATCH
    )
    production=$(ssh torch sbatch --parsable --dependency="afterok:$smoke" --account=torch_pr_570_general \
    --partition=cs --array="$prod_range" --ntasks=1 --cpus-per-task=1 --mem=16G --time=10:00:00 \
    --job-name="e5futil-prod-$cell" --output="$output_root/slurm_production_%A_%a.out" \
    --error="$output_root/slurm_production_%A_%a.err" <<SBATCH
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
exec python3 '$remote_runner' --manifest '$remote_manifest' --manifest-sha256 '$manifest_sha' --stage production --task-id "\$SLURM_ARRAY_TASK_ID" --output-root '$output_root'
SBATCH
    )
    verification=$(ssh torch sbatch --parsable --dependency="afterok:$production" --account=torch_pr_570_general \
    --partition=cs --array="$verify_range" --ntasks=1 --cpus-per-task=1 --mem=16G --time=01:00:00 \
    --job-name="e5futil-verify-$cell" --output="$output_root/slurm_verification_%A_%a.out" \
    --error="$output_root/slurm_verification_%A_%a.err" <<SBATCH
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
exec python3 '$remote_runner' --manifest '$remote_manifest' --manifest-sha256 '$manifest_sha' --stage verification --task-id "\$SLURM_ARRAY_TASK_ID" --output-root '$output_root'
SBATCH
    )
    finish=$(ssh torch sbatch --parsable --dependency="afterok:$verification" --account=torch_pr_570_general \
    --partition=cs --array="$cell" --ntasks=1 --cpus-per-task=1 --mem=16G --time=00:20:00 \
    --job-name="e5futil-final-$cell" --output="$output_root/slurm_final_%A_%a.out" \
    --error="$output_root/slurm_final_%A_%a.err" <<SBATCH
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
exec python3 '$remote_runner' --manifest '$remote_manifest' --manifest-sha256 '$manifest_sha' --stage finalize --task-id "\$SLURM_ARRAY_TASK_ID" --output-root '$output_root'
SBATCH
    )
    printf 'cell=%s smoke=%s production=%s verification=%s finalize=%s\n' "$cell" "$smoke" "$production" "$verification" "$finish"
  done
  printf 'output_root=%s\n' "$output_root"
}

validate_local
case "$mode" in
  stage) remote_check ;;
  submit) remote_check; submit ;;
  *) echo "usage: $0 [stage|submit]" >&2; exit 2 ;;
esac
