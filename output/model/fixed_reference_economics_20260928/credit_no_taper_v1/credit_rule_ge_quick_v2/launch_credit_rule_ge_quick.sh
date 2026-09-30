#!/usr/bin/env bash
# Prepared only; this launcher never submits itself.
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:16:00
#SBATCH --array=0-1
set -euo pipefail
reference=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
remote=/scratch/td2248/projects/fixed_reference_credit_rule_ge_quick_20260929
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
source_root=$remote/source_packet_v2
source=$source_root/credit_rule_ge_quick_v2
case=(ours author); selected=${case[$SLURM_ARRAY_TASK_ID]}
[[ $(date -u +%s) -lt 1790730025 ]] || { echo 'absolute deadline expired' >&2; exit 2; }
[[ -f "$source/run_credit_rule_ge_quick.py" && -f "$source/rule_interface.py" && -f "$source/run_ge_accounting_source.py" && -f "$source/plan.json" && -f "$source_root/credit_rule_quick_v1/run_credit_rule_quick.py" && -f "$source_root/credit_rule_quick_v1/strict_tenure.py" ]] || { echo 'stage both pinned sibling packets first' >&2; exit 2; }
[[ $(sha256sum "$source/run_credit_rule_ge_quick.py" | awk '{print $1}') == f5c0c5e3a099f5f398883a0dd17ff5b3c16561193409e90dce0031f8b466c7b9 ]] || exit 2
[[ $(sha256sum "$source/rule_interface.py" | awk '{print $1}') == 90e1e46b5f894dcf3e83b75842eb0ab44c0c549a848fbae00730218264846bdf ]] || exit 2
[[ $(sha256sum "$source/run_ge_accounting_source.py" | awk '{print $1}') == 59fa315ef6d50aaa8be2abcfb29f01cf406f3756e363a0bd088e473a424311b8 ]] || exit 2
[[ $(sha256sum "$source_root/credit_rule_quick_v1/strict_tenure.py" | awk '{print $1}') == c64f98bf72484facc36be5eb4d38d581c1f813b419feae9ab28eee3e96117677 ]] || exit 2
mkdir -p "$remote/results_v2" "$remote/cache"; [[ ! -e "$remote/results_v2/$selected" ]] || { echo 'refusing repeat' >&2; exit 2; }
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_CACHE_DIR=/work/numba_cache MPLBACKEND=Agg
apptainer exec --bind "$reference:$original:ro,$source_root:/work/source_packet:ro,$remote/results_v2:/work/results,$remote/cache:/work/numba_cache" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python /work/source_packet/credit_rule_ge_quick_v2/run_credit_rule_ge_quick.py --plan /work/source_packet/credit_rule_ge_quick_v2/plan.json --case "$selected" --output "/work/results/$selected"
