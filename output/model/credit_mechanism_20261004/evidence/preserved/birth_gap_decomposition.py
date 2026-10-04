"""Decompose how credit vs family-space shocks shift the net value of a first birth,
D = V_child - cost - V_wait = kappa*logit(p_attempt)/pi, at native nodes for childless
inherited renters, in wealth units (divide by dV/db of fert_value at the reference case).
Fixed price (chain-13 GE price). Weights: reconstructed pre-fertility mass, reference case."""
import os, sys, json
for v in ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
OUT=os.path.dirname(os.path.abspath(__file__)); os.environ["NUMBA_CACHE_DIR"]=os.path.join(OUT,"numba_cache")
sys.dont_write_bytecode=True
import numpy as np
ROOT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"; sys.path.insert(0,ROOT+"/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE, DEFAULT_PARAMETERS
from production.equilibrium import solve_at_price
from production.engine.parameters import get_fecundity_by_age
PSI=json.load(open(os.path.join(OUT,"credit_renormalized_psi.json")))["cap4"]["psi_star"]

def solve(phi,hR,hb,psi=None):
    P,grid=load_inputs(parameters={} if psi is None else {"psi_child":psi},external_inputs={"phi":[phi]*4,"hR_max":hR})
    P.hbar_child_rooms=hb; o=solve_at_price(P,grid,DEFAULT_PRICE); s,Q=o["solution"],o["P"]
    fec=np.asarray(get_fecundity_by_age(Q),float); k=float(Q.kappa_fert)
    fv=np.asarray(s.fert_value)[:,0,0]                     # (Nb,J,Nz) inherited renters
    fp=np.asarray(s.fert_probs)[:,0,0]                     # (Nb,J,Nz,4)
    p1=fp[...,1]; p0=fp[...,0]
    with np.errstate(divide="ignore",invalid="ignore"):
        Vwait=fv+k*np.log(p0); D=k*(np.log(p1)-np.log(p0))/fec[None,:,None]
    post0=np.asarray(s.g_beginning_distribution)[:,0,0,:,:,0,0]
    pre0=post0/np.clip(1-fec[None,:,None]*p1,1e-12,None)
    return dict(b=np.asarray(o["b_grid"]),fv=fv,Vwait=Vwait,D=D,p1=p1,w=pre0)

cases={"ref":solve(0.80,6.0,0.0),"credit":solve(0.95,6.0,0.0),
       "space_ref":solve(0.80,4.0,1.0,PSI),"space_credit":solve(0.95,4.0,1.0,PSI),
       "space_shock":solve(0.80,4.0,1.0)}
b=cases["ref"]["b"]
def slope(fv):  # forward difference of value in b, wealth units
    s=np.full_like(fv,np.nan); s[:-1]=(fv[1:]-fv[:-1])/(b[1:]-b[:-1])[:,None,None]; return s
rows=[]
for j,age in ((0,"18-21"),(1,"22-25"),(2,"26-29")):
    for z in (3,4,5):                                      # income states 4,5,6 (1-based)
        for refname,altname,label in (("ref","credit","credit, current model"),
                                      ("space_ref","space_credit","credit, 4-room cap +1 room/child"),
                                      ("ref","space_shock","space: 4-room cap +1 room/child")):
            R,A=cases[refname],cases[altname]
            Vb=slope(R["fv"])[:,j,z]; w=R["w"][:,j,z]
            ok=(w>1e-14)&np.isfinite(Vb)&(Vb>0)&np.isfinite(R["D"][:,j,z])&np.isfinite(A["D"][:,j,z])&np.isfinite(A["Vwait"][:,j,z])&np.isfinite(R["Vwait"][:,j,z])
            if w[ok].sum()<=0: continue
            W=w[ok]/w[ok].sum()
            dVwait=(A["Vwait"][:,j,z]-R["Vwait"][:,j,z])[ok]/Vb[ok]
            dD=(A["D"][:,j,z]-R["D"][:,j,z])[ok]/Vb[ok]
            dVchild=dVwait+dD
            scale=(R["D"][:,j,z][ok]*0+cases["ref"]["D"][:,j,z][ok]*0)  # placeholder
            kap=float(0.11652185618155607)
            rows.append(dict(age=age,z=z+1,comparison=label,mass_share=float(w[ok].sum()/R["w"][:,j,:].sum()),
                p_ref=float(W@R["p1"][:,j,z][ok]),p_alt=float(W@A["p1"][:,j,z][ok]),
                value_without_child=float(W@dVwait),value_with_child=float(W@dVchild),shift_in_birth_value=float(W@dD),
                logit_unit_in_wealth=float(W@(kap/Vb[ok])),D_ref_wealth=float(W@(R["D"][:,j,z][ok]/Vb[ok]))))
json.dump(rows,open(os.path.join(OUT,"birth_gap_decomposition.json"),"w"),indent=1)
for r in rows: print(json.dumps({k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items()}))
