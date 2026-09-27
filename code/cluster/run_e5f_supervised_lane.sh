#!/usr/bin/env bash
# Run a pinned smoke, then wait for the lead's separate acceptance/promotion.
# No source/contract modification and no scientific decision happens here.
set -euo pipefail
contract="${1:?smoke contract}"
results="${2:?result root}"
expected="${3:?reviewed contract SHA256}"
mode="${4:-smoke-only}"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
fi
python_bin="${E5F_PYTHON:-python3}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 NUMBA_DISABLE_JIT=0
export EXPECTED_UTILITY_OVERNIGHT_SHA256="$expected"
driver="$("$python_bin" - "$contract" "$expected" <<'PY'
import hashlib,json,sys
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
assert sha(sys.argv[1])==sys.argv[2],'Launch contract changed'
c=json.load(open(sys.argv[1]));r=c['files']['driver']
assert sha(r['path'])==r['sha256'],'Launch controller changed'
print(r['path'])
PY
)"
mkdir -p "$results"
"$python_bin" "$driver" --stage preflight --contract "$contract" --output "$results/preflight"
"$python_bin" "$driver" --stage smoke --contract "$contract" --output "$results/smoke"
if [[ "$mode" == "smoke-only" ]]; then exit 0; fi
if [[ "$mode" != "await-acceptance" ]]; then exit 2; fi
# Approval is a pinned JSON receipt written by the coordinating lead after all
# acceptance comparisons. Stop waiting at the global search cutoff.
"$python_bin" - "$contract" "$results" <<'PY'
import json, pathlib, sys, time
c=json.load(open(sys.argv[1]));out=pathlib.Path(sys.argv[2])
deadline=c['budget']['absolute_end_epoch']-c['budget']['repeat_seconds']-c['budget']['export_seconds']
while not (out/'acceptance.json').exists():
    if time.time()>=deadline: raise SystemExit('No acceptance before search cutoff')
    tmp=out/'waiting_for_acceptance.tmp'
    tmp.write_text(json.dumps({'status':'smoke_passed_waiting_for_lead','epoch':time.time()}))
    tmp.replace(out/'waiting_for_acceptance.json')
    time.sleep(30)
PY
exec "$python_bin" - "$results/acceptance.json" "$contract" "$expected" "$results" <<'PY'
import hashlib,json,os,sys,pathlib
r=json.load(open(sys.argv[1]));assert r['status']=='accepted_for_bounded_search'
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
assert sha(r['contract_path'])==r['contract_sha256']
assert sha(r['smoke_receipt'])==r['smoke_sha256']
original_path=pathlib.Path(sys.argv[2]).resolve();results=pathlib.Path(sys.argv[4]).resolve()
assert sha(original_path)==sys.argv[3]
assert pathlib.Path(r['smoke_receipt']).resolve()==results/'smoke/complete.json'
assert pathlib.Path(r['search_output']).resolve()==results/'search'
smoke=json.load(open(r['smoke_receipt']))
assert pathlib.Path(smoke['contract_path']).resolve()==original_path
assert smoke['contract_sha256']==sys.argv[3]
original=json.load(open(original_path));c=json.load(open(r['contract_path']))
def science(x):
    return {k:v for k,v in x.items() if k not in ('status','approval','verified_smoke','production_blockers')}
assert science(original)==science(c),'Promoted contract belongs to another lane'
assert c['approval']['production_authorized']
os.environ['EXPECTED_UTILITY_OVERNIGHT_SHA256']=r['contract_sha256']
os.execv(sys.executable,[sys.executable,c['files']['driver']['path'],'--stage','search',
    '--contract',r['contract_path'],'--output',r['search_output'],
    '--smoke-receipt',r['smoke_receipt'],'--smoke-sha256',r['smoke_sha256']])
PY
