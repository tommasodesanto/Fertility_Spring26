#!/bin/bash
#SBATCH --job-name=fertility_laptop_compare
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cpu_short
#SBATCH --time=00:17:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
module load anaconda3/2025.06
cd /scratch/td2248/projects/Fertility_Spring26_laptop_benchmark_20260919
hostname
lscpu | head -18
python -u -c 'import runpy,resource,sys; sys.argv=sys.argv[1:]; runpy.run_path(sys.argv[0],run_name="__main__"); print("PEAK_RSS_KIB",resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)' benchmark_saved_stationary.py \
 --source-root /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/corrected_initial_source_v2 \
 --checkpoint initial_state.pkl.gz --manifest manifest.json --output torch_run
