"""Stationary comparison of a doubled property tax (annual 1.06% -> 2.12%), chain-13 inputs.
Price solves the birth condition (solver renewal residual = 0); population solves the housing market:
N = H0*(ucr*p/rbar)^xi / rooms_per_household. B rebates the extra revenue lump-sum (iterated to balance).
Diagnostic (uncertified) root, validated on the unchanged baseline."""
import os, sys, json
for v in ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
OUT=os.path.dirname(os.path.abspath(__file__)); os.environ["NUMBA_CACHE_DIR"]=os.path.join(OUT,"numba_cache")
sys.dont_write_bytecode=True
import numpy as np
ROOT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"; sys.path.insert(0,ROOT+"/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE
from production.equilibrium import solve_at_price
TAU0=0.0423934430954904; TAU1=2*TAU0
def solve(tau,reb,p):
    P,g=load_inputs(external_inputs={"tau_H":tau,"property_tax_lump_sum_transfer":reb}); o=solve_at_price(P,g,p); return o["solution"],o["P"]
def stats(s,Q,p,tau):
    d=float(np.asarray(s.housing_demand).ravel()[0]); gch=np.asarray(s.g)
    F=sum(np.asarray(getattr(Q,k),float).sum() for k in ("_first_births_by_age","_second_births_by_age","_third_births_by_age"))
    return dict(price=p,ucr=float(Q.user_cost_rate),rent_per_room=float(Q.user_cost_rate*p),rooms=d,
                N=float(Q.H0[0]*(Q.user_cost_rate*p/Q.r_bar[0])**Q.xi_supply[0]/d),births_per_hh=float(F),
                own_all=float(gch[:,1:].sum()/gch.sum()),own_22=float(gch[:,1:,:,1].sum()/gch[:,:,:,1].sum()),
                revenue_per_hh=float(tau*p*d))
def root(tau,reb,p0=DEFAULT_PRICE):
    s0,Q0=solve(tau,reb,p0); r0=float(s0.adult_entry_stationary_residual)
    if abs(r0)<1e-11: return stats(s0,Q0,p0,tau)
    p1=p0*0.97; s1,Q1=solve(tau,reb,p1); r1=float(s1.adult_entry_stationary_residual)
    for _ in range(20):
        p2=p1-r1*(p1-p0)/(r1-r0); p0,r0=p1,r1; p1=p2; s1,Q1=solve(tau,reb,p1); r1=float(s1.adult_entry_stationary_residual)
        if abs(r1)<1e-11: break
    return stats(s1,Q1,p1,tau)
R={}
R["baseline"]=b=root(TAU0,0.0); print("baseline",json.dumps({k:round(v,6) for k,v in b.items()}),flush=True)
assert abs(b["price"]-0.7760569760205563)<1e-9 and abs(b["N"]-1)<1e-9
# fixed-price incentive effect (old price)
fa,Qa=solve(TAU1,0.0,b["price"]); fpa=stats(fa,Qa,b["price"],TAU1)
R["A_not_rebated"]=a=root(TAU1,0.0); print("A",json.dumps({k:round(v,6) for k,v in a.items()}),flush=True)
reb=(TAU1-TAU0)*b["price"]*b["rooms"]
for it in range(6):
    c=root(TAU1,reb,a["price"]); new=(TAU1-TAU0)*c["price"]*c["rooms"]
    print(f"rebate iter {it}: rebate {reb:.6f} -> implied {new:.6f}",flush=True)
    if abs(new-reb)<1e-6: break
    reb=new
c["rebate_per_hh"]=reb; R["B_rebated"]=c; print("B",json.dumps({k:round(v,6) for k,v in c.items()}),flush=True)
fb,Qb=solve(TAU1,reb,b["price"]); fpb=stats(fb,Qb,b["price"],TAU1)
R["fixed_price_births_change_pct"]={"A":100*(fpa["births_per_hh"]/b["births_per_hh"]-1),"B":100*(fpb["births_per_hh"]/b["births_per_hh"]-1)}
json.dump(R,open(os.path.join(OUT,"property_tax_ss.json"),"w"),indent=1)
for k in ("A_not_rebated","B_rebated"):
    r=R[k]; xi=0.63
    print(f"{k}: price {100*(r['price']/b['price']-1):+.2f}% | rent/room {100*(r['rent_per_room']/b['rent_per_room']-1):+.2f}% | supply {100*xi*np.log(r['rent_per_room']/b['rent_per_room']):+.2f}% | rooms/hh {100*np.log(r['rooms']/b['rooms']):+.2f}% | POPULATION {100*(r['N']-1):+.2f}% | own all {b['own_all']:.3f}->{r['own_all']:.3f} | births at old price {R['fixed_price_births_change_pct'][k[0]]:+.2f}%")
