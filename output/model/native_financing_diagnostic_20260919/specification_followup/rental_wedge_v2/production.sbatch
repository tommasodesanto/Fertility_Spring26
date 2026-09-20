#!/usr/bin/env bash
#SBATCH --job-name=rw_production
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --chdir=/scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/logs/%x_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/logs/%x_%j.err
#SBATCH --kill-on-invalid-dep=yes
set -euo pipefail
module load anaconda3/2025.06
unset PYTHONPATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0
export NUMBA_CACHE_DIR="/scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/numba_cache" MPLBACKEND=Agg PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
cd "/scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2"
echo "job ${SLURM_JOB_ID:-none} stage production start $(date -u +%FT%TZ) host $(hostname)"
python -B /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/code/model/tools/run_e5f_isolated_rental_wedge.py --mode verify --stage production --phase before --launch-manifest /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/launch_manifest.json --results /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/results
python -B /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/code/model/tools/run_e5f_isolated_rental_wedge.py --mode check-smoke --launch-manifest /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/launch_manifest.json --results /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/results
DEADLINE=$(( $(date +%s) + 3600 ))
timeout --signal=TERM --kill-after=30s 600s python -B /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/code/model/tools/run_e5f_isolated_rental_wedge.py --mode case --case cap10zero --launch-manifest /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/launch_manifest.json --results /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/results --deadline-epoch $DEADLINE
timeout --signal=TERM --kill-after=30s 900s python -B /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/code/model/tools/run_e5f_isolated_rental_wedge.py --mode case --case cap10s005 --launch-manifest /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/launch_manifest.json --results /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/results --deadline-epoch $DEADLINE
timeout --signal=TERM --kill-after=30s 900s python -B /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/code/model/tools/run_e5f_isolated_rental_wedge.py --mode case --case cap10s1 --launch-manifest /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/launch_manifest.json --results /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/results --deadline-epoch $DEADLINE
python -B /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/code/model/tools/run_e5f_isolated_rental_wedge.py --mode combine --launch-manifest /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/launch_manifest.json --results /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/results
python -B /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/code/model/tools/run_e5f_isolated_rental_wedge.py --mode verify --stage production --phase after --launch-manifest /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/launch_manifest.json --results /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_v2/results
echo "production stage complete $(date -u +%FT%TZ)"
