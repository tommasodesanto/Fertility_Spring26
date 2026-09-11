#!/bin/bash
#SBATCH --job-name=e5f_utility_fiscal
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_paygo_parenthood_1a5e88b7/smoke_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_paygo_parenthood_1a5e88b7/smoke_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_paygo_parenthood_1a5e88b7
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python - <<'PY'
import hashlib,json,os,time,unittest
from pathlib import Path
import numba
assert not numba.config.DISABLE_JIT
for name,pin in json.loads(Path('source_pins.json').read_text()).items():
    assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==pin,name
started=time.monotonic()
print('Compiled smoke: 11 tiny utility tests + 4 arithmetic tests + six two-date fiscal paths / 24 dated Bellman calls, plus tiny fixture setup.',flush=True)
suite=unittest.defaultTestLoader.loadTestsFromNames(['test_e5f_parenthood_utility','test_e5f_stationary_paygo','test_e5f_social_security_compiled'])
result=unittest.TextTestRunner(verbosity=2).run(suite)
packet=dict(status='passed' if result.wasSuccessful() else 'failed',source_commit='1a5e88b7',job_id=os.environ['SLURM_JOB_ID'],tests_run=result.testsRun,errors=len(result.errors),failures=len(result.failures),elapsed_seconds=time.monotonic()-started,conditional_paths_only=True,equilibrium_certified=False,calibration_launched=False)
Path('output/summary.json').write_text(json.dumps(packet,indent=2)+'\n')
print(json.dumps(packet),flush=True)
raise SystemExit(0 if result.wasSuccessful() else 1)
PY
