"""Mass of childless renters who would have a first birth if space/finance constraints were relaxed.
Fixed price (chain-13 GE price). Weights: reconstructed pre-fertility childless inherited-renter mass."""
import os, sys, json
for v in ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
OUT=os.path.dirname(os.path.abspath(__file__)); os.environ["NUMBA_CACHE_DIR"]=os.path.join(OUT,"numba_cache")
sys.dont_write_bytecode=True
import numpy as np
ROOT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"; sys.path.insert(0,ROOT+"/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE
from production.equilibrium import solve_at_price
from production.engine.parameters import get_fecundity_by_age
def run(phi,hR):
    P,g=load_inputs(external_inputs={"phi":[phi]*4,"hR_max":hR}); o=solve_at_price(P,g,DEFAULT_PRICE); s,Q=o["solution"],o["P"]
    fec=np.asarray(get_fecundity_by_age(Q),float); fp=np.asarray(s.fert_probs)[...,1]
    post0=np.asarray(s.g_beginning_distribution)[...,0,0]; pre0=post0/np.clip(1-fec[None,None,None,:,None]*fp,1e-12,None)
    w=pre0[:,0,0]; p=fp[:,0,0]                       # childless inherited renters (b,J,z)
    births=(w*p*fec[None,:,None])                   # expected first births
    F=[np.asarray(getattr(Q,k),float).sum() for k in ("_first_births_by_age","_second_births_by_age","_third_births_by_age")]
    return dict(phi=phi,hR=hR,w=w,births=births,total_births=float(sum(F)),
                renter_first_births_18_33=float(births[:,:4].sum()),childless_renter_mass_18_33=float(w[:,:4].sum()))
base=run(0.80,6.0); cases={"finance (99% LTV)":run(0.99,6.0),"space (rental cap 12)":run(0.80,12.0),"both":run(0.99,12.0)}
print(f"baseline: childless renters 18-33 mass {base['childless_renter_mass_18_33']:.4f}, their first births {base['renter_first_births_18_33']:.5f}, all births {base['total_births']:.5f}")
res={}
for k,c in cases.items():
    d1=c["renter_first_births_18_33"]-base["renter_first_births_18_33"]
    # by income state at ages 22-25, using each case's own weights
    bz_b=base["births"][:,1].sum(0); bz_c=c["births"][:,1].sum(0); wz=base["w"][:,1].sum(0)
    res[k]=dict(extra_renter_first_births_18_33=d1,pct=100*d1/base["renter_first_births_18_33"],
               all_births_pct=100*(c["total_births"]/base["total_births"]-1),
               rate22_by_z_base=(bz_b/np.maximum(wz,1e-300)).tolist(),rate22_by_z_case=(bz_c/np.maximum(c["w"][:,1].sum(0),1e-300)).tolist())
    print(f"{k:24s} extra first births among childless renters 18-33: {d1:+.5f} ({res[k]['pct']:+.2f}%) | all births {res[k]['all_births_pct']:+.2f}%")
    print("   first-birth rate 22-25 by income state, base:",[round(x,3) for x in res[k]['rate22_by_z_base']])
    print("   first-birth rate 22-25 by income state, case:",[round(x,3) for x in res[k]['rate22_by_z_case']])
json.dump(res,open(os.path.join(OUT,"would_be_parents.json"),"w"),indent=1)
