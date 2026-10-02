#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:30:00
set -euo pipefail
umask 077
python_bin=${ROOMS_PYTHON_BIN:-python3}
phase=${1:?Specify toy or full}
task_root=${2:?Specify frozen task directory}
[[ "$phase" == toy || "$phase" == full ]] || { echo "phase must be toy or full" >&2; exit 64; }
driver="$task_root/code/data/psid_followup_mar2026/rooms_income_sensitivity.py"
cd "$task_root"; started=$(date +%s); status=failed; stata_pid=""
receipt() { "$python_bin" - "$task_root" "$phase" "$started" "$1" "$status" <<'PY'
import json,pathlib,sys,time
p,phase,start,code,status=sys.argv[1:]; pathlib.Path(p,"run_receipt.json").write_text(json.dumps({"phase":phase,"status":status,"exit_code":int(code),"elapsed_seconds":time.time()-int(start),"processors":8},indent=2)+"\n")
PY
}
finish() { code=$?; trap - EXIT TERM INT; receipt "$code"; exit "$code"; }
terminate() { [[ -n "$stata_pid" ]] && kill "$stata_pid" 2>/dev/null || true; status=signal_143; receipt 143; trap - EXIT TERM INT; exit 143; }
trap finish EXIT; trap terminate TERM INT
module load stata/19.0
if [[ "$phase" == full ]]; then
  [[ -f entry.do && -f run_config.json && -f SHA256SUMS ]] || { echo "stage is incomplete" >&2; exit 66; }
  sha256sum -c SHA256SUMS
  "$python_bin" "$driver" validate-stage "$task_root"
  stata-mp -bq do entry.do & stata_pid=$!
  while kill -0 "$stata_pid" 2>/dev/null; do date -u +%Y-%m-%dT%H:%M:%SZ > heartbeat.txt; sleep 30; done
  wait "$stata_pid"; stata_pid=""
  for fit in baseline_full common_no_income common_income; do grep -Fq "INCOME_FIT_PASS $fit" stata_run.log || exit 67; done
  grep -Fq ROOMS_INCOME_THREE_FIT_PASS stata_run.log || exit 67
  "$python_bin" "$driver" finalize-private "$task_root"
else
  "$python_bin" "$driver" toy "$task_root" --reference-estimator "$task_root/reference_estimator.do"
  [[ -f TOY_PASS ]] || exit 69
fi
for fit in baseline_full common_no_income common_income; do [[ -f "results/$fit/completion.txt" ]] || exit 68; done
status=pass
