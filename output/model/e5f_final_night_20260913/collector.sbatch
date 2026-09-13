#!/bin/bash
#SBATCH --account=torch_pr_570_general
set -euo pipefail
module load anaconda3/2025.06
python -B /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/collect_e5f_final_night_v2.py /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913 --output /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/monitor_summary.json --jobs 17595967,17597260,17603133,17605279 --watch-seconds 43200 --interval 300
