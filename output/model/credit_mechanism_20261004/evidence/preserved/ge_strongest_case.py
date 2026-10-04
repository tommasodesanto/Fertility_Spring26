"""Stationary GE (birth-renewal price root, fixed H0 closure) for the strongest credit case:
renter cap 4, +1 room per child at home, psi_child re-normalized (0.322215), LTV 0.80 vs 0.95,
plus an unchanged-default validation run. All outputs under this scratch folder."""
import os, sys, json, time, traceback
for v in ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
OUT=os.path.dirname(os.path.abspath(__file__)); os.environ["NUMBA_CACHE_DIR"]=os.path.join(OUT,"numba_cache")
sys.dont_write_bytecode=True
import numpy as np
ROOT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"; sys.path.insert(0,ROOT+"/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE
from production.equilibrium import solve_stationary_ge
PSI=json.load(open(os.path.join(OUT,"credit_renormalized_psi.json")))["cap4"]["psi_star"]
CASES=[("ge_default",{},{},0.0),
       ("ge_cap4_hb1_phi80",{"psi_child":PSI},{"phi":[0.80]*4,"hR_max":4.0},1.0),
       ("ge_cap4_hb1_phi95",{"psi_child":PSI},{"phi":[0.95]*4,"hR_max":4.0},1.0)]
summary={}
for name,pars,ext,hb in CASES:
    t0=time.time(); out=os.path.join(OUT,"ge_runs",name)
    try:
        P,grid=load_inputs(parameters=pars,external_inputs=ext); P.hbar_child_rooms=float(hb)
        res=solve_stationary_ge(P,grid,out=out,price_start=DEFAULT_PRICE,budget_seconds=2400,max_lifecycle=32)
        keep={k:v for k,v in (res.items() if isinstance(res,dict) else []) if isinstance(v,(int,float,str,bool))}
        summary[name]=dict(status="ok",seconds=round(time.time()-t0,1),result_keys=sorted(res.keys()) if isinstance(res,dict) else str(type(res)),scalars=keep)
    except Exception as e:
        summary[name]=dict(status="failed",seconds=round(time.time()-t0,1),error=repr(e),tb=traceback.format_exc()[-2000:])
    json.dump(summary,open(os.path.join(OUT,"ge_strongest_case_summary.json"),"w"),indent=1,default=str)
    print(name,summary[name]["status"],summary[name]["seconds"],flush=True)
