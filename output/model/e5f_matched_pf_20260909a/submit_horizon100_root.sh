#!/bin/bash
#SBATCH --job-name=e5f_h100_root
#SBATCH --output=logs/h100_root_%j.out
#SBATCH --error=logs/h100_root_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_root_h100_20260910a
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_root_h100_20260910a/output/cache/mpl_root
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_root_h100_20260910a/output/cache/numba_root
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
NUMBA_DISABLE_JIT=1 python -m unittest test_run_e5f_matched_pf_historical_root test_e5f_matched_pf_history test_e5f_matched_pf_path_root test_collect_e5f_matched_pf_price_jacobian
python - <<'PYROOT'
import hashlib,json,subprocess,sys
from pathlib import Path
import run_e5f_matched_pf_historical_root as driver
r=Path('/scratch/td2248/projects/Fertility_Spring26_matched_pf_root_h100_20260910a')
e=Path('/scratch/td2248/projects/Fertility_Spring26_matched_pf_horizon_20260910a')
base=e/'horizon100_sequential_contract.json'
assert hashlib.sha256(base.read_bytes()).hexdigest()=='cd1fbfd7264f9354f3f4596fb2144d9f40a8c962eb9b77f0e0cf378ad4224ff3'
c=json.loads(base.read_text())
overrides={'code/model/tools/run_e5f_matched_pf_historical_root.py': '5c861bab36aaef23a7ad64e9a8255def81d25f9aa893e90de30f863a71e7859b', 'code/model/tools/run_e5f_matched_pf_history.py': 'b19512aa585bd61bcf50b29f6c145e4c35b71b4c6eeb39b391407bb6902551d3', 'code/model/tools/test_run_e5f_matched_pf_historical_root.py': 'bd27139c743875d9ac5804bb0b541e074beaba5072ce9c44a6b528bde477e606'}
old={'code/model/tools/run_e5f_matched_pf_historical_root.py': '1d86e998660e61eb2da571d913a7956b4ed1ec65b469c56e733af22b1ba245fb', 'code/model/tools/run_e5f_matched_pf_history.py': 'e394ade1e8abeb6291e1e749b029fe81b83cb5aba8eaef82e1c1071d711cd360', 'code/model/tools/test_run_e5f_matched_pf_historical_root.py': 'c77e42d505812c06191cf3c5d720fe406e98b07384e84a12c04cbc9267b6bb71'}
for name, expected in old.items():
    assert c['source_sha256'][name]==expected
    assert hashlib.sha256((e/name).read_bytes()).hexdigest()==expected
for name, expected in overrides.items():
    assert hashlib.sha256((r/name).read_bytes()).hexdigest()==expected
packet_path=e/'output/jacobian_horizon100_sequential.json'
packet=json.loads(packet_path.read_text())
c.update(source_root=str(r),experiment='normalized_historical_path_root',
         maximum_path_evaluations=6,seconds=21000,
         jacobian_packet=str(packet_path),
         jacobian_packet_sha256=hashlib.sha256(packet_path.read_bytes()).hexdigest())
c['source_sha256'].update(overrides)
evaluator='code/model/tools/run_e5f_matched_pf_baseline.py'
c['reviewed_evaluator_driver_change']=dict(path=evaluator,
    from_sha256=packet['shared_contract']['source_sha256'][evaluator],
    to_sha256=c['source_sha256'][evaluator],
    scope='explicit supplied-price hook and optional dated-state checkpoint; unchanged economic evaluator')
c['reviewed_root_runtime_changes']=dict(
    scope='root runtime ceiling and explicit provenance checks only; unchanged economic evaluator and numerical root',
    files={name:dict(from_sha256=old[name],to_sha256=overrides[name]) for name in sorted(overrides)})
driver.validate_jacobian_packet(packet,c,'sequential')
contract=r/'horizon100_root_contract.json'
if contract.exists():
    raise FileExistsError(contract)
contract.write_text(json.dumps(c,indent=2)+'\n')
pin=hashlib.sha256(contract.read_bytes()).hexdigest()
driver.joined.load_smoke_contract(contract,pin,'sequential',maximum_seconds=driver.MAXIMUM_ROOT_SECONDS)
subprocess.run([sys.executable,'code/model/tools/run_e5f_matched_pf_historical_root.py',
    '--contract',str(contract),'--contract-sha256',pin,'--arm','sequential',
    '--output',str(r/'output/historical_root_h100_01/sequential')],check=True)
PYROOT
