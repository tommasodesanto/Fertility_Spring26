"""Fixed-price credit x child earnings penalty (after-tax earnings x (1-pen) while >=1 child at home,
working ages; engine option child_earnings_penalty). Chain-13 inputs; no market clearing/recalibration."""
import os, sys, json
for v in ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
OUT=os.path.dirname(os.path.abspath(__file__)); os.environ["NUMBA_CACHE_DIR"]=os.path.join(OUT,"numba_cache")
sys.dont_write_bytecode=True
import numpy as np
ROOT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"; sys.path.insert(0,ROOT+"/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE
from production.equilibrium import solve_at_price
from production.engine.parameters import child_earnings_penalty_active

def run(phi,hR,hb,pen):
    P,grid=load_inputs(external_inputs={"phi":[phi]*4,"hR_max":float(hR)})
    P.hbar_child_rooms=float(hb); P.child_earnings_penalty=np.array([0.0,pen,pen,pen])
    out=solve_at_price(P,grid,DEFAULT_PRICE); sol,Q=out["solution"],out["P"]
    assert child_earnings_penalty_active(Q)==(pen>0)
    post=np.asarray(sol.g_beginning_distribution); gch=np.asarray(sol.g)
    F1,F2,F3=(np.asarray(getattr(Q,k),float) for k in ("_first_births_by_age","_second_births_by_age","_third_births_by_age"))
    pp=post.sum(axis=(0,1,2,4,6)); pre=pp.copy(); pre[:,0]+=F1; pre[:,1]+=-F1+F2; pre[:,2]+=-F2+F3; pre[:,3]+=-F3
    s=0.125*pre[1]/pre[1].sum()+0.875*pp[1]/pp[1].sum(); endc=pp[8]/pp[8].sum()
    return dict(phi=phi,hR=hR,hb=hb,pen=pen,births=float(F1.sum()+F2.sum()+F3.sum()),first=float(F1.sum()),
        second=float(F2.sum()),third=float(F3.sum()),age25=float(s@np.arange(4)),end_ceb=float(endc@np.arange(4)),
        own_22=float(gch[:,1:,:,1].sum()/gch[:,:,:,1].sum()),own_30=float(gch[:,1:,:,3].sum()/gch[:,:,:,3].sum()),
        nan=bool(np.isnan(gch).any()),entry_censored=float(Q._entry_censored_mass))

if __name__=="__main__":
    R=[]
    for (hb,hR) in ((0.0,6.0),(1.0,4.0)):
        for pen in (0.10,0.20):
            for phi in (0.80,0.95):
                r=run(phi,hR,hb,pen); R.append(r)
                print(json.dumps({k:(round(v,5) if isinstance(v,float) else v) for k,v in r.items()}),flush=True)
                json.dump(R,open(os.path.join(OUT,"credit_earnings_penalty.json"),"w"),indent=1)
