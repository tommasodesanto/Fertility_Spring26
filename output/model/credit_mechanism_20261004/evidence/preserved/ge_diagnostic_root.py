"""Uncertified diagnostic stationary GE: secant on price for the solver's birth-renewal residual
(entry - births_adj/2.1 = 0) using production solve_at_price; fixed-H0 population scale
N = H0*(ucr*p/r_bar)^xi / housing_demand (same formula as production/equilibrium.py:60-62).
Skips the frozen reporting/budget audit gates (which failed on ~1e-7 infeasible mass)."""
import os, sys, json
for v in ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
OUT=os.path.dirname(os.path.abspath(__file__)); os.environ["NUMBA_CACHE_DIR"]=os.path.join(OUT,"numba_cache")
sys.dont_write_bytecode=True
import numpy as np
ROOT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"; sys.path.insert(0,ROOT+"/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE
from production.equilibrium import solve_at_price
PSI=json.load(open(os.path.join(OUT,"credit_renormalized_psi.json")))["cap4"]["psi_star"]

def solve(pars,ext,hb,p):
    P,grid=load_inputs(parameters=pars,external_inputs=ext); P.hbar_child_rooms=float(hb)
    o=solve_at_price(P,grid,p); return o["solution"],o["P"]

def stats(sol,Q,p):
    post=np.asarray(sol.g_beginning_distribution); gch=np.asarray(sol.g)
    F=[np.asarray(getattr(Q,k),float) for k in ("_first_births_by_age","_second_births_by_age","_third_births_by_age")]
    pp=post.sum(axis=(0,1,2,4,6)); pre=pp.copy(); pre[:,0]+=F[0]; pre[:,1]+=-F[0]+F[1]; pre[:,2]+=-F[1]+F[2]; pre[:,3]+=-F[2]
    s=0.125*pre[1]/pre[1].sum()+0.875*pp[1]/pp[1].sum(); endc=pp[8]/pp[8].sum()
    fac=(Q.user_cost_rate*p/Q.r_bar[0])**Q.xi_supply[0]
    return dict(price=p,residual=float(sol.adult_entry_stationary_residual),entry=float(sol.entry_rate),
        population_scale=float(Q.H0[0]*fac/float(np.asarray(sol.housing_demand).ravel()[0])),
        rooms_per_household=float(np.asarray(sol.housing_demand).ravel()[0]),
        births=float(sum(f.sum() for f in F)),age25=float(s@np.arange(4)),end_ceb=float(endc@np.arange(4)),
        childless_end=float(endc[0]),own_22=float(gch[:,1:,:,1].sum()/gch[:,:,:,1].sum()),
        own_30=float(gch[:,1:,:,3].sum()/gch[:,:,:,3].sum()),own_all=float(gch[:,1:].sum()/gch.sum()))

def root(pars,ext,hb):
    p0,p1=DEFAULT_PRICE,DEFAULT_PRICE*1.01
    s0,Q0=solve(pars,ext,hb,p0); r0=float(s0.adult_entry_stationary_residual)
    if abs(r0)<1e-11: return stats(s0,Q0,p0)|{"iterations":1}
    s1,Q1=solve(pars,ext,hb,p1); r1=float(s1.adult_entry_stationary_residual)
    for it in range(15):
        p2=p1-r1*(p1-p0)/(r1-r0); p0,r0=p1,r1; p1=p2
        s1,Q1=solve(pars,ext,hb,p1); r1=float(s1.adult_entry_stationary_residual)
        if abs(r1)<1e-11: break
    return stats(s1,Q1,p1)|{"iterations":it+3}

CASES=[("default_phi80",{}, {"phi":[0.80]*4},0.0),("default_phi95",{},{"phi":[0.95]*4},0.0),
       ("cap4_hb1_phi80",{"psi_child":PSI},{"phi":[0.80]*4,"hR_max":4.0},1.0),
       ("cap4_hb1_phi95",{"psi_child":PSI},{"phi":[0.95]*4,"hR_max":4.0},1.0)]
R={}
for name,pars,ext,hb in CASES:
    R[name]=root(pars,ext,hb); print(name,json.dumps({k:round(v,6) for k,v in R[name].items()}),flush=True)
    json.dump(R,open(os.path.join(OUT,"ge_diagnostic_root.json"),"w"),indent=1)
