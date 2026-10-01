#!/bin/zsh
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1
here=${0:A:h}
out="$here/results_$(date -u +%Y%m%dT%H%M%SZ)"
deadline=${1:-$(( $(date +%s) + 600 ))}
print -r -- "OUTPUT=$out DEADLINE_EPOCH=$deadline"
"${here:h:h:h:h:h:h:h:h}/code/model/.venv/bin/python" "$here/run_renter_only.py" --out "$out" --deadline-epoch "$deadline" > "$here/current_stdout.txt" 2> "$here/current_stderr.txt" &
pid=$!
(
  while kill -0 "$pid" 2>/dev/null; do
    rss=$(ps -o rss= -p "$pid" | tr -d ' ')
    [[ -z "$rss" ]] && break
    if (( rss > 25165824 )); then
      print -r -- "RSS_EXCEEDED_KB=$rss PID=$pid" > "$here/memory_stop.txt"
      kill -TERM "$pid" 2>/dev/null || true
      break
    fi
    sleep 1
  done
) &
watch=$!
set +e
wait "$pid"
run_exit=$?
kill "$watch" 2>/dev/null || true
wait "$watch" 2>/dev/null || true
print -r -- "EXIT_STATUS=$run_exit OUTPUT=$out"
exit "$run_exit"
