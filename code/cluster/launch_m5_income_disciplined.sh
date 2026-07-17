#!/usr/bin/env bash
# One-shot M5 launch: deploy -> smoke -> (on pass) nested ref + 8 chains + collector.
# Prereq: a working `ssh torch` login. Run from repo root:
#   bash code/cluster/launch_m5_income_disciplined.sh
set -uo pipefail
cd "$(dirname "$0")/../.."

FILES=(
  code/model/intergen_housing_fertility/solver.py
  code/model/intergen_housing_fertility/calibration.py
  code/model/intergen_housing_fertility/tests/test_bequest_target_moments.py
  code/model/tools/run_intergen_bequest_exit_chain.py
  code/model/tools/collect_intergen_internal_bequest_recalibration.py
  code/cluster/submit_intergen_income_disciplined_recalibration.sh
  code/cluster/submit_intergen_income_disciplined_recalibration_smoke.sh
  code/cluster/submit_intergen_income_disciplined_recalibration_collector.sh
  code/cluster/submit_intergen_income_disciplined_nested_reference.sh
)
REMOTE='$SCRATCH/projects/Fertility_Spring26'

echo ">> auth probe"
ssh -o BatchMode=yes torch 'echo OK' || { echo "!! Torch login not working. Refresh it first."; exit 2; }

echo ">> deploy ${#FILES[@]} files"
for f in "${FILES[@]}"; do
  rsync -c "$f" "torch:$REMOTE/$f" || { echo "!! deploy failed: $f"; exit 2; }
done

echo ">> submit smoke"
SMOKE=$(ssh torch "cd $REMOTE/code/cluster && sbatch --parsable submit_intergen_income_disciplined_recalibration_smoke.sh") || exit 2
echo "   smoke job $SMOKE; polling (max 45 min)"
for i in $(seq 1 45); do
  sleep 60
  STATE=$(ssh -o BatchMode=yes torch "squeue -j $SMOKE -h -o %T 2>/dev/null | head -1")
  [ -z "$STATE" ] && break
  echo "   [$i min] smoke: $STATE"
done

echo ">> smoke verdict"
DEAD=$(ssh torch "grep -l 'InfeasibleThetaError\|Traceback' \$(ls -t $REMOTE/code/cluster/logs/*income_disciplined*smoke* 2>/dev/null | head -4) 2>/dev/null" || true)
OK=$(ssh torch "grep -h 'status=ok' \$(ls -t $REMOTE/code/cluster/logs/*income_disciplined*smoke* 2>/dev/null | head -4) 2>/dev/null | wc -l" || echo 0)
echo "   ok-evals: $OK; error files: ${DEAD:-none}"
if [ -n "$DEAD" ] || [ "${OK:-0}" -eq 0 ]; then
  echo "!! smoke unhealthy -- NOT launching production. Inspect logs above."
  exit 3
fi

echo ">> launch production chain"
NREF=$(ssh torch "cd $REMOTE/code/cluster && sbatch --parsable submit_intergen_income_disciplined_nested_reference.sh")
MAIN=$(ssh torch "cd $REMOTE/code/cluster && sbatch --parsable submit_intergen_income_disciplined_recalibration.sh")
COLL=$(ssh torch "cd $REMOTE/code/cluster && sbatch --parsable --dependency=afterany:$NREF:$MAIN submit_intergen_income_disciplined_recalibration_collector.sh")
echo "LAUNCHED: nested=$NREF main=$MAIN collector=$COLL"
echo "$(date '+%F %T') nested=$NREF main=$MAIN collector=$COLL smoke=$SMOKE" >> output/model/intergen_income_disciplined_recalibration_20260716/JOB_IDS.txt
echo ">> morning: results under \$SCRATCH .../output/model/intergen_income_disciplined_recalibration_20260716/ ; pull final_report + collector outputs."
