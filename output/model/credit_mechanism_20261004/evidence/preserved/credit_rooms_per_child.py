"""Fixed-price credit x rental cap x per-child room requirement, chain-13 inputs.
Same solver as credit_space_2x2.py; hbar_child_rooms set on P before the solve
(precompute_shared rebuilds h_bar). No market clearing, no recalibration."""
import os, sys, json
for v in ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
OUT=os.path.dirname(os.path.abspath(__file__)); os.environ["NUMBA_CACHE_DIR"]=os.path.join(OUT,"numba_cache")
sys.dont_write_bytecode=True
import numpy as np
ROOT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"; sys.path.insert(0,ROOT+"/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE
from production.equilibrium import solve_at_price
from production.engine.parameters import get_fecundity_by_age, independent_child_maturation_active

def run(phi,hR,hb):
    P,grid=load_inputs(external_inputs={"phi":[phi]*4,"hR_max":float(hR)})
    P.hbar_child_rooms=float(hb)
    out=solve_at_price(P,grid,DEFAULT_PRICE); sol,Q=out["solution"],out["P"]
    post=np.asarray(sol.g_beginning_distribution); gch=np.asarray(sol.g)
    F1,F2,F3=(np.asarray(getattr(Q,k),float) for k in ("_first_births_by_age","_second_births_by_age","_third_births_by_age"))
    pp=post.sum(axis=(0,1,2,4,6)); pre=pp.copy(); pre[:,0]+=F1; pre[:,1]+=-F1+F2; pre[:,2]+=-F2+F3; pre[:,3]+=-F3
    s=0.125*pre[1]/pre[1].sum()+0.875*pp[1]/pp[1].sum()
    endc=pp[8]/pp[8].sum()
    own_par=lambda j,n: float(gch[:,1:,:,j,:,n].sum()/max(gch[:,:,:,j,:,n].sum(),1e-300))
    bad=[a for a in dir(sol) if any(t in a.lower() for t in ("infeas","dead","censor"))]
    return dict(phi=phi,hR=hR,hb=hb,
        hbar_by_n=[float(x) for x in np.asarray(Q._h_bar if hasattr(Q,"_h_bar") else []).ravel()[:0]],
        births=float(F1.sum()+F2.sum()+F3.sum()),first=float(F1.sum()),second=float(F2.sum()),third=float(F3.sum()),
        age25=float(s@np.arange(4)),end_dist=endc.tolist(),end_ceb=float(endc@np.arange(4)),
        second_22=float(F2[1]),second_26=float(F2[2]),
        own_22=float(gch[:,1:,:,1].sum()/gch[:,:,:,1].sum()),own_30=float(gch[:,1:,:,3].sum()/gch[:,:,:,3].sum()),
        own_par30=[own_par(3,n) for n in range(4)],
        mass=float(gch.sum()),nan=bool(np.isnan(gch).any() or np.isnan(np.asarray(sol.fert_probs)).any()),
        entry_censored=float(getattr(Q,"_entry_censored_mass",np.nan)),diag_attrs=bad,
        indep_maturation=bool(independent_child_maturation_active(Q)))

if __name__=="__main__":
    R=[]
    for hb in (0.0,0.5,1.0):
        for hR in (6.0,4.0):
            for phi in (0.80,0.95):
                r=run(phi,hR,hb); R.append(r)
                print(json.dumps({k:(round(v,5) if isinstance(v,float) else v) for k,v in r.items() if k not in("end_dist","hbar_by_n")}),flush=True)
                json.dump(R,open(os.path.join(OUT,"credit_rooms_per_child.json"),"w"),indent=1)
