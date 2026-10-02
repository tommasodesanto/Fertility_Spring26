"""Explicit read-only local file overlay, equivalent to two frozen source binds."""
import builtins,hashlib,importlib.machinery,io,json,os,runpy,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;PACKET=HERE.parent;ROOT=PACKET.parents[3]
DIGESTS={'e5f_exact_policy_cache.py':'d51bbd13026288db0194b44b420bb49ff6970588f40a906242978d49038a5d6f','test_e5f_exact_policy_cache.py':'50784af81c5e6fd64d345d504c8c7fa81209b05a6a655789a224022485208ae5'}
MAPPING={str(ROOT/'code/model/tools'/name):str(HERE/'frozen_sources'/name) for name in DIGESTS}
for name,digest in DIGESTS.items():
 assert hashlib.sha256((HERE/'frozen_sources'/name).read_bytes()).hexdigest()==digest,name

def mapped(path):
 try:return MAPPING.get(os.path.abspath(os.fspath(path)),path)
 except TypeError:return path
# Restrict redirects to reads. Writes to original or frozen source are forbidden.
old_io_open=io.open;old_open=builtins.open

def redirected_open(path,mode='r',*args,**kwargs):
 target=mapped(path)
 if target!=path and any(k in mode for k in 'wax+'):raise RuntimeError('Frozen overlay write forbidden')
 return old_io_open(target,mode,*args,**kwargs)

def redirected_builtin(path,mode='r',*args,**kwargs):
 target=mapped(path)
 if target!=path and any(k in mode for k in 'wax+'):raise RuntimeError('Frozen overlay write forbidden')
 return old_open(target,mode,*args,**kwargs)
io.open=redirected_open;builtins.open=redirected_builtin
Path.open=lambda self,mode="r",*args,**kwargs:redirected_open(self,mode,*args,**kwargs)
old_code=importlib.machinery.SourceFileLoader.get_code

def source_code(loader,fullname):
 target=mapped(loader.path)
 if target!=loader.path:
  # Ignore any current-main bytecode and compile the authenticated frozen source.
  return loader.source_to_code(old_io_open(target,'rb').read(),loader.path)
 return old_code(loader,fullname)
importlib.machinery.SourceFileLoader.get_code=source_code
import importlib,numpy as np
# Existing setup-guide NumPy2 checkpoint namespace compatibility on NumPy1.
import pathlib
sys.modules.setdefault('pathlib._local',pathlib)
sys.modules.setdefault('numpy._core',np.core)
sys.modules.setdefault('numpy._core.multiarray',importlib.import_module('numpy.core.multiarray'))
sys.modules.setdefault('numpy._core.numeric',importlib.import_module('numpy.core.numeric'))
sys.path.insert(0,str(PACKET))
if '--preflight-context' in sys.argv:
    sys.argv.remove('--preflight-context')
    runpy.run_path(str(HERE/'preflight_context.py'),run_name='__main__')
else:
    assert '--chain' in sys.argv, 'Explicit local chain required'
    chain=int(sys.argv[sys.argv.index('--chain')+1])
    assert 48 <= chain <= 57, 'Only local chain IDs 48..57 may launch here'
    pins=json.loads((HERE/'local_manifest.json').read_text())
    for name,digest in pins.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    local_plan=json.loads((HERE/'local_plan.json').read_text())
    parent_plan=json.loads((PACKET/'plan.json').read_text())
    assert local_plan['base_target_contract']==parent_plan['base_target_contract']
    assert local_plan['profiles']==parent_plan['profiles']
    assert local_plan['target_fingerprint']==parent_plan['target_fingerprint']
    assert local_plan['weight_fingerprint']==parent_plan['weight_fingerprint']
    assert local_plan['nearby_starts'][:48]==parent_plan['nearby_starts']
    runpy.run_path(str(HERE/'run_local_psi.py'),run_name='__main__')
