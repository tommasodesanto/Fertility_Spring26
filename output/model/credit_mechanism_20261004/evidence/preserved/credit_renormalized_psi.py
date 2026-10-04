"""Credit test after re-normalizing the child benefit psi so births return to the chain-13 level,
with +1 room per child at home. Fixed price (chain-13 GE price), no market clearing, no full
recalibration: only psi_child moves (diagnostic; values above the 0.5 search bound are labeled)."""
import os, sys, json
for v in ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
OUT=os.path.dirname(os.path.abspath(__file__)); os.environ["NUMBA_CACHE_DIR"]=os.path.join(OUT,"numba_cache")
sys.dont_write_bytecode=True
import numpy as np
ROOT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"; sys.path.insert(0,ROOT+"/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE, DEFAULT_PARAMETERS
from production.equilibrium import solve_at_price

def run(phi,hR,hb,psi=None):
    pars={} if psi is None else {"psi_child":float(psi)}
    P,grid=load_inputs(parameters=pars,external_inputs={"phi":[phi]*4,"hR_max":float(hR)})
    P.hbar_child_rooms=float(hb)
    out=solve_at_price(P,grid,DEFAULT_PRICE); sol,Q=out["solution"],out["P"]
    post=np.asarray(sol.g_beginning_distribution); gch=np.asarray(sol.g)
    F1,F2,F3=(np.asarray(getattr(Q,k),float) for k in ("_first_births_by_age","_second_births_by_age","_third_births_by_age"))
    pp=post.sum(axis=(0,1,2,4,6)); pre=pp.copy(); pre[:,0]+=F1; pre[:,1]+=-F1+F2; pre[:,2]+=-F2+F3; pre[:,3]+=-F3
    s=0.125*pre[1]/pre[1].sum()+0.875*pp[1]/pp[1].sum(); endc=pp[8]/pp[8].sum()
    return dict(phi=phi,hR=hR,hb=hb,psi=float(Q.psi_child),births=float(F1.sum()+F2.sum()+F3.sum()),
        first=float(F1.sum()),second=float(F2.sum()),third=float(F3.sum()),age25=float(s@np.arange(4)),
        end_ceb=float(endc@np.arange(4)),childless_end=float(endc[0]),
        own_22=float(gch[:,1:,:,1].sum()/gch[:,:,:,1].sum()),own_30=float(gch[:,1:,:,3].sum()/gch[:,:,:,3].sum()),
        nan=bool(np.isnan(gch).any()),entry_censored=float(Q._entry_censored_mass))

def show(r): print(json.dumps({k:(round(v,6) if isinstance(v,float) else v) for k,v in r.items()}),flush=True)

if __name__=="__main__":
    R={"baseline":run(0.80,6.0,0.0)}; show(R["baseline"])
    assert abs(R["baseline"]["births"]-0.11528810162653443)<1e-12 and abs(R["baseline"]["age25"]-0.5338046821453732)<1e-12, "baseline drift"
    target=R["baseline"]["births"]; psi0=float(DEFAULT_PARAMETERS["psi_child"])
    for hR in (4.0,6.0):
        lo,hi=psi0,0.5; rl=run(0.80,hR,1.0,lo)
        rh=run(0.80,hR,1.0,hi)
        while rh["births"]<target and hi<3.0: lo,rl=hi,rh; hi*=1.6; rh=run(0.80,hR,1.0,hi)
        for _ in range(30):
            mid=0.5*(lo+hi); rm=run(0.80,hR,1.0,mid)
            if rm["births"]<target: lo=mid
            else: hi=mid
            if abs(rm["births"]/target-1)<2e-4: break
        psi_star=mid; print(f"cap {hR}: psi* = {psi_star:.5f} (births {rm['births']:.6f} vs target {target:.6f})",flush=True)
        a=run(0.80,hR,1.0,psi_star); b=run(0.95,hR,1.0,psi_star); show(a); show(b)
        R[f"cap{int(hR)}"]={"psi_star":psi_star,"phi80":a,"phi95":b,"credit_effect_births":b["births"]/a["births"]-1}
        print(f"cap {hR}: credit effect on births {100*(b['births']/a['births']-1):+.3f}%",flush=True)
        json.dump(R,open(os.path.join(OUT,"credit_renormalized_psi.json"),"w"),indent=1)
