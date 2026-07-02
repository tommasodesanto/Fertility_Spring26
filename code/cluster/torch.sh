#!/usr/bin/env bash
#
# torch.sh — standard NYU Torch HPC workflow wrapper.
#
# Runs LOCALLY on your Mac and shells over `ssh torch`. It encodes the account
# rule, the remote working directory, and the exact submit / smoke / monitor /
# collect commands so the cluster ceremony is one consistent call for any agent
# (especially cheaper models driving overnight runs). Full procedure:
# docs/workflow/delegation_and_cluster_playbook.md
#
# This wrapper NEVER works around credentials. If the SSH probe fails it stops
# and asks you to refresh the Torch login.

set -uo pipefail

SSH_HOST="${TORCH_SSH_HOST:-torch}"
CLUSTER_USER="${TORCH_USER:-td2248}"

# Remote cluster dir. $SCRATCH is kept LITERAL here so it expands ON Torch, not
# on the Mac. Point TORCH_REMOTE_DIR at a dated snapshot when CALIBRATION_STATUS.md
# names one, e.g. /scratch/td2248/projects/<SNAPSHOT>/code/cluster
REMOTE_DIR='$SCRATCH/projects/Fertility_Spring26/code/cluster'
[ -n "${TORCH_REMOTE_DIR:-}" ] && REMOTE_DIR="$TORCH_REMOTE_DIR"
REMOTE_MODEL_DIR="${REMOTE_DIR%/cluster}/model"

err() { printf '%s\n' "$*" >&2; }

usage() {
  cat >&2 <<'EOF'
torch.sh — standard NYU Torch HPC workflow (runs locally, over `ssh torch`).

Usage:
  torch.sh status
  torch.sh smoke   <family> <tag>
  torch.sh submit  <family> <tag> [extra sbatch args...]
  torch.sh logs    <family> [--follow]
  torch.sh health  <family> <tag>
  torch.sh collect <family> <tag>

Families:
  intergen-twohour   submit_intergen_housing_fertility_twohour_panel.sh
  intergen-de        submit_intergen_housing_fertility_global_de.sh
  intergen-polish    submit_intergen_housing_fertility_local_polish.sh   (needs INTERGEN_SEED_THETA_JSON)
  direct             submit_python_direct_geometry_overnight.sh

Env overrides:
  TORCH_SSH_HOST    ssh alias                       (default: torch)
  TORCH_USER        cluster user for squeue         (default: td2248)
  TORCH_REMOTE_DIR  remote .../code/cluster path     (default: $SCRATCH/projects/Fertility_Spring26/code/cluster)
  TORCH_ENV         extra "KEY=VAL KEY=VAL" prefixed to a submit
  TORCH_DRYRUN=1    print the ssh command(s) without executing (local check)
  TORCH_LOG_FILES   recent log files for `logs`     (default: 4)
  TORCH_LOG_LINES   tail lines per log for `logs`   (default: 40)

Long-run discipline (see CLAUDE.md "Long-Run Search Safety"):
  1) smoke <family> <tag>   before any real run
  2) submit with an explicit run tag, walltime, array size, seed base
  3) health <family> <tag>  and confirm best.json + cases.jsonl are being written
  4) verify the active target set / snapshot path in CALIBRATION_STATUS.md
EOF
}

# family -> metadata. Sets globals: SCRIPT TAGVAR PREFIX COLLECTOR LOGGLOB
family_meta() {
  case "$1" in
    intergen-twohour)
      SCRIPT=submit_intergen_housing_fertility_twohour_panel.sh
      TAGVAR=INTERGEN_RUN_TAG; PREFIX=intergen_housing_fertility
      COLLECTOR=collect_intergen_panel_results.py; LOGGLOB=slurm_ihf_2hr ;;
    intergen-de)
      SCRIPT=submit_intergen_housing_fertility_global_de.sh
      TAGVAR=INTERGEN_RUN_TAG; PREFIX=intergen_housing_fertility
      COLLECTOR=collect_intergen_panel_results.py; LOGGLOB=slurm_ihf_de ;;
    intergen-polish)
      SCRIPT=submit_intergen_housing_fertility_local_polish.sh
      TAGVAR=INTERGEN_RUN_TAG; PREFIX=intergen_housing_fertility
      COLLECTOR=collect_intergen_panel_results.py; LOGGLOB=slurm_ihf_polish ;;
    direct)
      SCRIPT=submit_python_direct_geometry_overnight.sh
      TAGVAR=DT_DIRECT_RUN_TAG; PREFIX=python_direct_geometry
      COLLECTOR=collect_direct_geometry_results.py; LOGGLOB=slurm_py_direct ;;
    *)
      err "Unknown family: $1"
      err "Families: intergen-twohour | intergen-de | intergen-polish | direct"
      return 1 ;;
  esac
}

run_remote() {
  # $1 = human description, $2 = remote command string
  local desc="$1" cmd="$2" rc=0
  err ">> torch $desc"
  err "   ssh $SSH_HOST \"$cmd\""
  if [ -n "${TORCH_DRYRUN:-}" ]; then
    err "   (TORCH_DRYRUN set; not executed)"
    return 0
  fi
  ssh -o BatchMode=yes "$SSH_HOST" "$cmd" || rc=$?
  if [ "$rc" -eq 255 ]; then
    err "!! Torch SSH failed (auth/connection, exit 255)."
    err "   Refresh your Torch login and retry. Do not work around credentials."
    exit 2
  fi
  return "$rc"
}

cmd_status() {
  err ">> torch status  (auth probe + queue)"
  if [ -n "${TORCH_DRYRUN:-}" ]; then
    err "   ssh -o BatchMode=yes $SSH_HOST 'hostname; echo SCRATCH=\$SCRATCH; squeue -u $CLUSTER_USER'"
    err "   (TORCH_DRYRUN set; not executed)"
    return 0
  fi
  local rc=0
  ssh -o BatchMode=yes "$SSH_HOST" 'hostname; echo "SCRATCH=$SCRATCH"; squeue -u '"$CLUSTER_USER" || rc=$?
  if [ "$rc" -ne 0 ]; then
    err "!! Torch SSH probe failed (exit $rc)."
    err "   Refresh your Torch login and retry. Do not work around credentials."
    exit 2
  fi
}

cmd_smoke() {
  local family="${1:?family required}" tag="${2:?tag required}"
  family_meta "$family" || exit 1
  local pre
  case "$family" in
    direct)
      pre="DT_DIRECT_RUN_TAG=$tag DT_DIRECT_SETUP=fast DT_DIRECT_MAX_ITER_EQ=4 DT_DIRECT_MAX_EVALS=2 DT_DIRECT_BUDGET_SEC=1800"
      run_remote "smoke $family $tag" "cd $REMOTE_DIR && $pre sbatch --array=1-2%2 $SCRIPT" ;;
    intergen-de)
      pre="INTERGEN_RUN_TAG=$tag INTERGEN_MINUTES=8 INTERGEN_GLOBAL_EVALS_PER_TASK=6 INTERGEN_GLOBAL_POP_SIZE=6"
      run_remote "smoke $family $tag" "cd $REMOTE_DIR && $pre sbatch --array=1-2%2 --time=00:20:00 $SCRIPT" ;;
    intergen-twohour)
      pre="INTERGEN_RUN_TAG=$tag INTERGEN_MINUTES=8 INTERGEN_CASES_PER_TASK=3"
      run_remote "smoke $family $tag" "cd $REMOTE_DIR && $pre sbatch --array=1-2%2 --time=00:20:00 $SCRIPT" ;;
    intergen-polish)
      : "${INTERGEN_SEED_THETA_JSON:?intergen-polish needs INTERGEN_SEED_THETA_JSON (a seed-theta path on Torch)}"
      pre="INTERGEN_RUN_TAG=$tag INTERGEN_MINUTES=8 INTERGEN_LOCAL_MAX_EVALS=4 INTERGEN_SEED_THETA_JSON=$INTERGEN_SEED_THETA_JSON"
      run_remote "smoke $family $tag" "cd $REMOTE_DIR && $pre sbatch --array=1-2%2 --time=00:20:00 $SCRIPT" ;;
  esac
}

cmd_submit() {
  local family="${1:?family required}" tag="${2:?tag required}"; shift 2 || true
  family_meta "$family" || exit 1
  local sbatch_args="$*"
  local env_prefix="${TORCH_ENV:-}"
  if [ "$family" = intergen-polish ]; then
    : "${INTERGEN_SEED_THETA_JSON:?intergen-polish needs INTERGEN_SEED_THETA_JSON (a seed-theta path on Torch)}"
    env_prefix="$env_prefix INTERGEN_SEED_THETA_JSON=$INTERGEN_SEED_THETA_JSON"
  fi
  run_remote "submit $family $tag" "cd $REMOTE_DIR && ${env_prefix:+$env_prefix }$TAGVAR=$tag sbatch $sbatch_args $SCRIPT"
}

cmd_logs() {
  local family="${1:?family required}"; shift || true
  family_meta "$family" || exit 1
  if [ "${1:-}" = "--follow" ]; then
    run_remote "logs $family --follow" "cd $REMOTE_DIR && tail -f \$(ls -t logs/${LOGGLOB}_*.out 2>/dev/null | head -n 1)"
  else
    run_remote "logs $family" "cd $REMOTE_DIR && for f in \$(ls -t logs/${LOGGLOB}_*.out 2>/dev/null | head -n ${TORCH_LOG_FILES:-4}); do echo \"== \$f ==\"; tail -n ${TORCH_LOG_LINES:-40} \"\$f\"; done"
  fi
}

cmd_health() {
  local family="${1:?family required}" tag="${2:?tag required}"
  family_meta "$family" || exit 1
  local rdir="results_${PREFIX}_$tag"
  run_remote "health $family $tag" "cd $REMOTE_DIR && echo '== best.json ==' && find $rdir -name best.json -printf '%TY-%Tm-%Td %TH:%TM  %p\n' 2>/dev/null | sort | tail -n 8; echo '== jsonl (cases) ==' && find $rdir -name '*.jsonl' -printf '%TY-%Tm-%Td %TH:%TM  %p\n' 2>/dev/null | sort | tail -n 8"
}

cmd_collect() {
  local family="${1:?family required}" tag="${2:?tag required}"
  family_meta "$family" || exit 1
  run_remote "collect $family $tag" "cd $REMOTE_MODEL_DIR && python tools/$COLLECTOR --results-dir ../cluster/results_${PREFIX}_$tag"
}

main() {
  local sub="${1:-}"; shift || true
  case "$sub" in
    status|probe)    cmd_status "$@" ;;
    smoke)           cmd_smoke "$@" ;;
    submit)          cmd_submit "$@" ;;
    logs)            cmd_logs "$@" ;;
    health)          cmd_health "$@" ;;
    collect)         cmd_collect "$@" ;;
    help|-h|--help)  usage ;;
    "")              usage; exit 1 ;;
    *)               err "Unknown subcommand: $sub"; err ""; usage; exit 1 ;;
  esac
}

main "$@"
