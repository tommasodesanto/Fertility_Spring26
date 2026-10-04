"""Fixed-price test: unsecured borrowing for renters (b' >= -d_bar), current chain-13 model."""
import os, sys, json
for v in ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
OUT=os.path.dirname(os.path.abspath(__file__)); os.environ["NUMBA_CACHE_DIR"]=os.path.join(OUT,"numba_cache")
sys.dont_write_bytecode=True
import numpy as np
ROOT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"; sys.path.insert(0,ROOT+"/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE
from production.equilibrium import solve_at_price
from production.engine.parameters import get_fecundity_by_age
R=[]
for d in (0.0,0.37,0.74,1.48):
    P,grid=load_inputs(external_inputs={"unsecured_credit_limit":d})
    o=solve_at_price(P,grid,DEFAULT_PRICE); s,Q=o["solution"],o["P"]; b=np.asarray(o["b_grid"])
    post=np.asarray(s.g_beginning_distribution); gch=np.asarray(s.g); fec=np.asarray(get_fecundity_by_age(Q),float)
    F=[np.asarray(getattr(Q,k),float) for k in ("_first_births_by_age","_second_births_by_age","_third_births_by_age")]
    pp=post.sum(axis=(0,1,2,4,6)); pre=pp.copy(); pre[:,0]+=F[0]; pre[:,1]+=-F[0]+F[1]; pre[:,2]+=-F[1]+F[2]; pre[:,3]+=-F[2]
    sh=0.125*pre[1]/pre[1].sum()+0.875*pp[1]/pp[1].sum(); endc=pp[8]/pp[8].sum()
    fp1=np.asarray(s.fert_probs)[...,1]; post0=post[...,0,0]; pre0=post0/np.clip(1-fec[None,None,None,:,None]*fp1,1e-12,None)
    w=pre0[:,0,0]                                   # childless inherited renters (Nb,J,Nz)
    att22=((w[:,1]*fp1[:,0,0,1]).sum(0)/np.maximum(w[:,1].sum(0),1e-300)).tolist()
    # use of credit: renter mass (post-tenure) at negative beginning wealth, by age
    gr=gch[:,0]; neg=[float(gr[b<-1e-9][:,:,j].sum()/max(gr[:,:,j].sum(),1e-300)) for j in range(8)]
    bp=np.asarray(s.bp_pol)[:,0,0]; minbp22=float(np.nanmin(np.where(gch[:,0,0,1]>1e-12,bp[:,1],np.nan)))
    c=np.asarray(s.c_pol)[:,0,0,1]; c22=float((gch[:,0,0,1]*c).sum()/gch[:,0,0,1].sum())
    r=dict(d_bar=d,births=float(sum(f.sum() for f in F)),first=float(F[0].sum()),second=float(F[1].sum()),
           age25=float(sh@np.arange(4)),any25=float(1-sh[0]),end_ceb=float(endc@np.arange(4)),childless_end=float(endc[0]),
           mean_age_fb=float((F[0]/F[0].sum())@(20+4*np.arange(Q.J))),
           own_22=float(gch[:,1:,:,1].sum()/gch[:,:,:,1].sum()),own_30=float(gch[:,1:,:,3].sum()/gch[:,:,:,3].sum()),
           renter_share_negative_b_by_age=neg,min_renter_bprime_22=minbp22,renter_c_22=c22,attempt22_by_z=att22)
    R.append(r); print(json.dumps({k:(round(v,4) if isinstance(v,float) else ([round(x,3) for x in v] if isinstance(v,list) else v)) for k,v in r.items()}),flush=True)
    json.dump(R,open(os.path.join(OUT,"unsecured_credit.json"),"w"),indent=1)
