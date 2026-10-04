import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(sys.argv[1]).parent))
import two_shock as d
import two_shock_runtime as r
from numba import get_num_threads
assert get_num_threads()==1
plan=json.load(open(sys.argv[2]));d.preflight(plan)
out=Path(sys.argv[3]);runtime=r.build_runtime(plan=plan,output=out,smoke=True)
assert runtime.rt.total_native_calls==0
assert runtime.rt.identity()==plan['identity']
d.write(out/'native_status.json',dict(status='PASS_ZERO_SOLVES_NATIVE_CONSTRUCTOR',native_calls=0,numba_threads=1,identity=runtime.rt.identity(),scientific_validation=False))
