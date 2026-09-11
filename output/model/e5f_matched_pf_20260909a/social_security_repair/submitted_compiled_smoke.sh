#!/bin/bash
#SBATCH --job-name=e5f_ss_smoke
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_paygo_audit_20260911a/smoke_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_paygo_audit_20260911a/smoke_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_paygo_smoke_a654219c
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/output/cache/mpl"
export NUMBA_CACHE_DIR="$PWD/output/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python - <<'PY'
import json, os, time, unittest
from pathlib import Path
import numba
assert not numba.config.DISABLE_JIT
started=time.monotonic()
print('Compiled Social Security smoke: six two-date paths, 24 dated Bellman calls, plus tiny fixture setup.',flush=True)
suite=unittest.defaultTestLoader.loadTestsFromName('test_e5f_social_security_compiled')
result=unittest.TextTestRunner(verbosity=2).run(suite)
packet=dict(status='passed' if result.wasSuccessful() else 'failed',
    source_commit='a654219c',job_id=os.environ['SLURM_JOB_ID'],
    tests_run=result.testsRun,errors=len(result.errors),failures=len(result.failures),
    elapsed_seconds=time.monotonic()-started,conditional_paths_only=True,
    equilibrium_certified=False,calibration_launched=False)
Path('output/summary.json').write_text(json.dumps(packet,indent=2)+'\n')
print(json.dumps(packet),flush=True)
raise SystemExit(0 if result.wasSuccessful() else 1)
PY
