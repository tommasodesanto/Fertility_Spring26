#!/usr/bin/env bash
# Local launcher (experiment only): stage and submit the smoothed-transition batch on Torch.
# Usage: code/cluster/submit_e5f_ssj_smoothed_transition.sh <kappa> [tag]
# Optional env: E5F_SMOOTH_STATIONARY_START="price,pension,rebate" (warm start for the probe-scale
# stationary solve), E5F_SMOOTH_STATIONARY_EVALS=16|24 (frozen terminal solver budget).
# Example: code/cluster/submit_e5f_ssj_smoothed_transition.sh 0.05 20260916a
set -euo pipefail
KAPPA="${1:?kappa required, e.g. 0.05}"; TAG="${2:-$(date -u +%Y%m%d)a}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
R=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches
B="$R/announced_original_queue_20260913c_ssj_smoothed_k${KAPPA}_${TAG}"
ssh -o BatchMode=yes torch "mkdir -p $B/source"
scp -q "$ROOT/code/cluster/run_e5f_ssj_smoothed_transition.py" "$ROOT/code/cluster/run_e5f_ssj_announced_rescue.py" \
       "$ROOT/code/model/tools/e5f_ssj_toeplitz_jacobian.py" "$ROOT/code/model/tools/e5f_ssj_scaled_step_root.py" \
       "torch:$B/source/"
scp -q "$ROOT/code/cluster/prepare_e5f_ssj_smoothed_transition.py" "torch:$B/"
ssh -o BatchMode=yes torch "cd $B && E5F_SMOOTH_STATIONARY_START='${E5F_SMOOTH_STATIONARY_START:-}' E5F_SMOOTH_STATIONARY_EVALS='${E5F_SMOOTH_STATIONARY_EVALS:-16}' /share/apps/anaconda3/2025.06/bin/python prepare_e5f_ssj_smoothed_transition.py $B $KAPPA && sbatch run.sbatch"
LOCAL="$ROOT/output/model/e5f_sequence_space_prototype_20260913/smoothed_transition_k${KAPPA}_${TAG}"
mkdir -p "$LOCAL"; scp -q "torch:$B/manifest.json" "torch:$B/run.sbatch" "$LOCAL/"
echo "batch: $B"; echo "local mirror: $LOCAL"
