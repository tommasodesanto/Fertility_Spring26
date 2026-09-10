#!/bin/bash
#SBATCH --job-name=e5f_h100_continue
#SBATCH --output=logs/h100_continue_%j.out
#SBATCH --error=logs/h100_continue_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=03:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --nodelist=cs693
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_root_h100_20260910a
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR=$PWD/output/cache/mpl_h100_continuation
export NUMBA_CACHE_DIR=$PWD/output/cache/numba_h100_continuation
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
NUMBA_DISABLE_JIT=1 python -m unittest test_run_e5f_matched_pf_historical_root test_e5f_matched_pf_history test_e5f_matched_pf_path_root test_collect_e5f_matched_pf_price_jacobian
python - <<'PYCONT'
from pathlib import Path
import hashlib,json,subprocess,sys
import run_e5f_matched_pf_historical_root as driver
r=Path.cwd()
raw=r/'horizon100_root_contract.json'
assert hashlib.sha256(raw.read_bytes()).hexdigest()=='dff9ae51142c051ac02614bdace5fc21e4d938ac92f5cac2305637319234510e'
c=json.loads(raw.read_text())
parent=r/'output/historical_root_h100_01/sequential'
parent_contract=json.loads((parent/'contract.json').read_text())
assert parent_contract['contract_sha256']=='dff9ae51142c051ac02614bdace5fc21e4d938ac92f5cac2305637319234510e'
assert parent_contract['source_sha256']==c['source_sha256']
summary=json.loads((parent/'summary.json').read_text())
if summary.get('finite_horizon_market_converged') is True:
    print('Continuation not needed: parent already converged. No model solve.')
    raise SystemExit(0)
c.update(seconds=10800,maximum_path_evaluations=3)
for name,filename in [('restart_contract','contract.json'),('restart_history','root_history.json'),('restart_summary','summary.json')]:
    path=parent/filename
    c[name]=str(path)
    c[name+'_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
packet_path=Path(c['jacobian_packet'])
driver.primitive.verify(packet_path,c['jacobian_packet_sha256'])
packet=json.loads(packet_path.read_text())
driver.validate_jacobian_packet(packet,c,'sequential')
driver.load_restart(c,packet,'sequential')
contract=r/'horizon100_continuation_contract.json'
output=r/'output/historical_root_h100_continue_01/sequential'
if contract.exists() or output.exists():
    raise FileExistsError('Continuation already exists; refusing duplicate execution')
contract.write_text(json.dumps(c,indent=2)+'\n')
pin=hashlib.sha256(contract.read_bytes()).hexdigest()
driver.joined.load_smoke_contract(contract,pin,'sequential',maximum_seconds=driver.MAXIMUM_ROOT_SECONDS)
subprocess.run([sys.executable,'code/model/tools/run_e5f_matched_pf_historical_root.py',
    '--contract',str(contract),'--contract-sha256',pin,'--arm','sequential','--output',str(output)],check=True)
PYCONT
