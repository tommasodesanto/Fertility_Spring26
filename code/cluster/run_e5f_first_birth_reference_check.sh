#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
set -euo pipefail
umask 077
module load stata/19.0
phase=${1:?Specify toy or full}
task_root=${2:?Specify frozen task directory}
outdir="$task_root/$phase"
mkdir -p "$outdir"
cd "$outdir"
cat > entry.do <<EODO
clear all
set more off
set processors 8
version 17.0
sysdir set PLUS "$task_root/ado/"
adopath ++ "$task_root/ado/"
log using "$outdir/sa_rooms_first_birth_household_aligned_v1.log", replace text
EODO
if [[ "$phase" == full ]]; then
    printf 'use "%s/analysis_sample.dta", clear\n' "$task_root" >> entry.do
fi
printf 'do "%s/reference_check.do" %s "%s"\n' "$task_root" "$phase" "$outdir" >> entry.do
stata-mp -bq do entry.do &
stata_pid=$!
trap 'kill "$stata_pid" 2>/dev/null || true' TERM INT
while kill -0 "$stata_pid" 2>/dev/null; do
    date -u '+%Y-%m-%dT%H:%M:%SZ' > heartbeat.txt
    sleep 10
done
wait "$stata_pid"
grep -q '^FIRST_BIRTH_REFERENCE_INVARIANCE_CHECK_COMPLETED$' sa_rooms_first_birth_household_aligned_v1.log
printf 'pass\n' > completion.txt
