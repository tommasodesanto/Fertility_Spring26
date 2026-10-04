"""Credit/space relaxations in the 'benefit world': transfer floor G0=0.75, Gn=0.6 per child at home,
psi re-normalized to restore chain-13 births. Fixed price (chain-13 GE price); transfers unfinanced."""
import os, sys, json
for v in ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
OUT=os.path.dirname(os.path.abspath(__file__)); os.environ["NUMBA_CACHE_DIR"]=os.path.join(OUT,"numba_cache")
sys.dont_write_bytecode=True
import numpy as np
ROOT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"; sys.path.insert(0,ROOT+"/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE
from production.equilibrium import solve_at_price
PSI=[r for r in json.load(open(os.path.join(OUT,"gradient_fixes_renorm.json"))) if r["Gn"]==0.6][0]["psi"]
GROUPS={"low":[0,1,2,3],"mid":[4],"high":[5,6,7,8]}
def run(ext):
    P,grid=load_inputs(parameters={"psi_child":PSI},external_inputs=ext); P.transfer_floor_G0=0.75; P.transfer_floor_Gn=0.6
    o=solve_at_price(P,grid,DEFAULT_PRICE); s,Q=o["solution"],o["P"]
    post=np.asarray(s.g_beginning_distribution); gch=np.asarray(s.g)
    F=[np.asarray(getattr(Q,k),float) for k in ("_first_births_by_age","_second_births_by_age","_third_births_by_age")]
    pp=post.sum(axis=(0,1,2,4,6)); pre=pp.copy(); pre[:,0]+=F[0]; pre[:,1]+=-F[0]+F[1]; pre[:,2]+=-F[1]+F[2]; pre[:,3]+=-F[2]
    sh=0.125*pre[1]/pre[1].sum()+0.875*pp[1]/pp[1].sum(); endc=pp[8]/pp[8].sum()
    pz=post[:,:,:,1].sum(axis=(0,1,2,5)); at26={k:float(pz[ix].sum(0)@np.arange(4)/pz[ix].sum()) for k,ix in GROUPS.items()}
    oz=gch[:,1:,:,1].sum(axis=(0,1,3,4)); mz=gch[:,:,:,1].sum(axis=(0,1,2,4,5)); 
    own22={k:float(gch[:,1:,:,1][...,ix,:,:].sum()/gch[:,:,:,1][...,ix,:,:].sum()) for k,ix in GROUPS.items()}
    return dict(births=float(sum(f.sum() for f in F)),first=float(F[0].sum()),childless_end=float(endc[0]),end_ceb=float(endc@np.arange(4)),
                age25=float(sh@np.arange(4)),at26=at26,own22=own22,own22_all=float(gch[:,1:,:,1].sum()/gch[:,:,:,1].sum()))
cases=[("benefit world (LTV 80%)",{}),("loans to 95%",{"phi":[0.95]*4}),("loans to 99%",{"phi":[0.99]*4}),
       ("renter borrowing 1 yr earnings",{"unsecured_credit_limit":0.74}),("rental cap 12 rooms",{"hR_max":12.0})]
R={}; base=None
for name,ext in cases:
    r=run(ext); R[name]=r; base=base or r
    print(f"{name:32s} births {100*(r['births']/base['births']-1):+6.2f}% | childless {100*r['childless_end']:.1f}% | completed {r['end_ceb']:.3f} | age25 {r['age25']:.3f} | CEB@26 L/M/H "+"/".join(f"{r['at26'][k]:.2f}" for k in GROUPS)+f" | own22 {r['own22_all']:.2f} (L/M/H "+"/".join(f"{r['own22'][k]:.2f}" for k in GROUPS)+")",flush=True)
    json.dump(dict(psi=PSI,results=R),open(os.path.join(OUT,"credit_in_benefit_world.json"),"w"),indent=1)
