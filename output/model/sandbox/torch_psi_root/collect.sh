#!/usr/bin/env bash
#
# collect.sh — pull finished sandbox_psi_root_20260918a outputs from Torch.
#
# For each spec, copies the four run_ss.py output files from the batch's
# output_runs/<spec>/ into output/model/sandbox/torch_psi_root/<spec>/.
# Run from the repository root. Safe to re-run; skips specs with no remote
# output yet.
set -uo pipefail

BATCH=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/sandbox_psi_root_20260918a
LOCAL_ROOT=output/model/sandbox/torch_psi_root
FILES="summary.md moments.csv parameters.csv graphs.pdf"

SPECS="sw_penalty_psi_root sw_mortgage_psi_root sw_wedge_psi_root sw_estate_psi_root all_switches_psi_root earnings_e6b_psi_root earnings_bgm_psi_root credit_line_modest_psi_root s1_concave_benefit kappa_h_zero"

for spec in $SPECS; do
  dest="$LOCAL_ROOT/$spec"
  mkdir -p "$dest"
  # shellcheck disable=SC2086
  if scp -o BatchMode=yes "torch:$BATCH/output_runs/$spec/"\{summary.md,moments.csv,parameters.csv,graphs.pdf\} "$dest/" 2>/dev/null; then
    echo "collected $spec"
  else
    echo "not ready: $spec"
  fi
done
# Smoke output (diagnostic only) can be pulled with:
#   scp -o BatchMode=yes "torch:$BATCH/output_runs/sw_penalty_psi_fixed/"\{summary.md,moments.csv,parameters.csv,graphs.pdf\} "$LOCAL_ROOT/smoke_sw_penalty_psi_fixed/"
