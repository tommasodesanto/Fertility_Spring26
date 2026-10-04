"""Fixed-price tests of two in-engine levers for the income-fertility gradient (chain-13 inputs):
(a) means-tested floor rising with children at home: resources topped up to G0 + Gn*m (kernels.py:676-681);
(b) larger taste noise in birth choices (kappa_fert, kappa_fert_continuation scaled).
Income groups fixed by state: low = states 1-4, mid = 5, high = 6-9. No GE, no recalibration, no fiscal financing."""
import os, sys, json
for v in ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
OUT=os.path.dirname(os.path.abspath(__file__)); os.environ["NUMBA_CACHE_DIR"]=os.path.join(OUT,"numba_cache")
sys.dont_write_bytecode=True
import numpy as np
ROOT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"; sys.path.insert(0,ROOT+"/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE, DEFAULT_PARAMETERS
from production.equilibrium import solve_at_price
from production.engine.household import income_at_state
from production.engine.parameters import children_at_home_count
GROUPS={"low":[0,1,2,3],"mid":[4],"high":[5,6,7,8]}
K0,K1=float(DEFAULT_PARAMETERS["kappa_fert"]),float(DEFAULT_PARAMETERS["kappa_fert_continuation"])
def run(G0=0.0,Gn=0.0,kscale=1.0):
    pars={} if kscale==1.0 else {"kappa_fert":K0*kscale,"kappa_fert_continuation":K1*kscale}
    P,grid=load_inputs(parameters=pars); P.transfer_floor_G0=G0; P.transfer_floor_Gn=Gn
    o=solve_at_price(P,grid,DEFAULT_PRICE); s,Q=o["solution"],o["P"]; b=np.asarray(o["b_grid"])
    post=np.asarray(s.g_beginning_distribution); gch=np.asarray(s.g)
    F=[np.asarray(getattr(Q,k),float) for k in ("_first_births_by_age","_second_births_by_age","_third_births_by_age")]
    pp=post.sum(axis=(0,1,2,4,6)); pre=pp.copy(); pre[:,0]+=F[0]; pre[:,1]+=-F[0]+F[1]; pre[:,2]+=-F[1]+F[2]; pre[:,3]+=-F[2]
    sh=0.125*pre[1]/pre[1].sum()+0.875*pp[1]/pp[1].sum()
    def by_group(j):
        pz=post[:,:,:,j].sum(axis=(0,1,2,5)); out={}
        for k,ix in GROUPS.items():
            m=pz[ix].sum(0); out[k]=dict(ceb=float(m@np.arange(4)/m.sum()),childless=float(m[0]/m.sum()))
        return out
    # transfer outlays per household per period (post-tenure distribution, debt-blind means test)
    z=np.asarray(Q.z_grid if hasattr(Q,"z_grid") else Q.z_vals).ravel(); R=float(Q.R_gross); T=0.0
    if G0>0 or Gn>0:
        npar,ncs=gch.shape[5],gch.shape[6]
        mk=np.array([[children_at_home_count(n,c,Q) for c in range(ncs)] for n in range(npar)],float)
        for j in range(gch.shape[3]):
            for zi in range(gch.shape[4]):
                y=income_at_state(Q,0,j,float(z[zi])); res=R*np.maximum(b,0)+y
                gfl=G0+Gn*mk                                       # (npar,ncs)
                tr=np.clip(gfl[None,:,:]-res[:,None,None],0,None); tr=np.minimum(tr,gfl[None,:,:])
                T+=float((gch[:,:,:,j,zi].sum(axis=(1,2))*tr).sum())
    return dict(G0=G0,Gn=Gn,kscale=kscale,births=float(sum(f.sum() for f in F)),age25=float(sh@np.arange(4)),
                at26=by_group(1),at42_45=by_group(6),transfer_per_household_period=T,
                transfer_pct_mean_period_earnings=100*T/(4*0.740801))
if __name__=="__main__":
    cases=[dict(),dict(Gn=0.3),dict(G0=0.75),dict(G0=0.75,Gn=0.3),dict(G0=0.75,Gn=0.6),dict(kscale=2.0),dict(kscale=4.0)]
    R=[]
    for c in cases:
        try: r=run(**c)
        except Exception as e: r=dict(case=c,error=repr(e))
        R.append(r); json.dump(R,open(os.path.join(OUT,"gradient_fixes.json"),"w"),indent=1)
        if "error" in r: print(c,"ERROR",r["error"]); continue
        f=lambda d:"/".join(f"{d[k]['ceb']:.2f}" for k in GROUPS); g=lambda d:"/".join(f"{100*d[k]['childless']:.0f}%" for k in GROUPS)
        print(f"G0={r['G0']:.2f} Gn={r['Gn']:.2f} kx{r['kscale']:.0f} | births {r['births']:.4f} age25 {r['age25']:.3f} | CEB@26 L/M/H {f(r['at26'])} | CEB@42-45 {f(r['at42_45'])} | childless@42-45 {g(r['at42_45'])} | transfers {r['transfer_pct_mean_period_earnings']:.2f}% of mean earnings",flush=True)
