#!/usr/bin/env bash
# Isolated fixed-price joint-choice diagnostic, never a calibration search.
#SBATCH --job-name=e5fjoint
#SBATCH --partition=cpu_short
#SBATCH --time=00:35:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --account=torch_pr_570_general
#SBATCH --output=logs/slurm_e5fjoint_%j.out
#SBATCH --error=logs/slurm_e5fjoint_%j.err
set -euo pipefail
: "${E5F_JOINT_STAGE:?smoke or panel required}"
: "${E5F_JOINT_CONTRACT_SHA256:?required contract hash}"
: "${E5F_JOINT_CHECKPOINT:?required reference checkpoint}"
SNAPSHOT="$(cd "${SLURM_SUBMIT_DIR:?submit from frozen code/cluster}/../.." && pwd)"
module load anaconda3/2025.06
export PYTHONPATH="$SNAPSHOT/code/model:$SNAPSHOT/code/model/tools"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/e5f_joint_mpl_${SLURM_JOB_ID}"
cd "$SNAPSHOT"
ARGS=(--stage "$E5F_JOINT_STAGE" --contract "$SNAPSHOT/contract.json"
      --contract-sha256 "$E5F_JOINT_CONTRACT_SHA256"
      --checkpoint "$E5F_JOINT_CHECKPOINT" --outdir "$SNAPSHOT/output/$E5F_JOINT_STAGE")
if [[ "$E5F_JOINT_STAGE" = panel ]]; then
  : "${E5F_JOINT_SMOKE_SHA256:?required completed exact-loop smoke hash}"
  ARGS+=(--smoke "$SNAPSHOT/output/smoke" --smoke-sha256 "$E5F_JOINT_SMOKE_SHA256")
elif [[ "$E5F_JOINT_STAGE" != smoke ]]; then
  echo 'Unsupported stage' >&2; exit 2
fi
python3 -u code/model/tools/run_e5f_joint_nested_experiment.py "${ARGS[@]}"
