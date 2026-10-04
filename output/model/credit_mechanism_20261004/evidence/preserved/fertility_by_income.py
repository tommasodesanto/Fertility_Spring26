"""Children ever born by family income and education: CPS June 2024 (local public file)
vs chain-13 model by income state (fixed price = chain-13 GE price). Read-only."""
import os, sys, csv, json
for v in ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
OUT=os.path.dirname(os.path.abspath(__file__)); os.environ["NUMBA_CACHE_DIR"]=os.path.join(OUT,"numba_cache")
sys.dont_write_bytecode=True
import numpy as np
ROOT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
rows=[]
with open(ROOT+"/code/data/cps_fertility/cache/jun24pub.csv",newline="") as f:
    r=csv.DictReader(f)
    for d in r:
        try:
            if d["PESEX"]!="2": continue
            age=int(d["PRTAGE"]); n=int(d["PTSF1"]); w=float(d["PWSSWGT"]); inc=int(d["HEFAMINC"]); ed=int(d["PEEDUCA"]); rel=int(d["PRFAMREL"])
        except (ValueError,KeyError): continue
        if not (0<=n<=5) or w<=0: continue
        if 24<=age<=26 or 40<=age<=44: rows.append((age,n,w,inc,ed,rel))
A=np.array(rows,float)
def grp_inc(c): return np.where(c<1,-1,np.where(c<=10,0,np.where(c<=14,1,2)))   # <40k, 40-100k, >=100k
def grp_ed(c):  return np.where(c<=38,0,np.where(c==39,1,np.where(c<=42,2,np.where(c==43,3,4))))
INC=["<$40k","$40-100k",">=$100k"]; ED=["<HS","HS","some college","BA","grad"]
def tab(mask,g,labels):
    out=[]
    for k,lab in enumerate(labels):
        m=mask&(g==k); w=A[m,2]; n=np.minimum(A[m,1],3)
        if w.sum()==0: continue
        out.append(dict(group=lab,n=int(m.sum()),pop_share=float(w.sum()/A[mask,2].sum()),
            mean_ceb_cap3=float(n@w/w.sum()),childless=float(w[A[m,1]==0].sum()/w.sum()),two_plus=float(w[A[m,1]>=2].sum()/w.sum())))
    return out
young=(A[:,0]>=24)&(A[:,0]<=26); old=(A[:,0]>=40)
own_fam=np.isin(A[:,5],[1,2])   # family reference person or spouse
data=dict(age24_26_income=tab(young,grp_inc(A[:,3]),INC),age24_26_income_ownfamily=tab(young&own_fam,grp_inc(A[:,3]),INC),
          age24_26_educ=tab(young,grp_ed(A[:,4]),ED),age40_44_income=tab(old,grp_inc(A[:,3]),INC),age40_44_educ=tab(old,grp_ed(A[:,4]),ED))
# model
sys.path.insert(0,ROOT+"/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE
from production.equilibrium import solve_at_price
P,g=load_inputs(); o=solve_at_price(P,g,DEFAULT_PRICE); s=o["solution"]
post=np.asarray(s.g_beginning_distribution)   # (b,ten,loc,J,z,n,cs)
def model_at(j):
    pz=post[:,:,:,j].sum(axis=(0,1,2,5))       # (z,n)
    m=pz.sum(1); res=[]
    for z in range(pz.shape[0]):
        if m[z]<=0: continue
        res.append(dict(z=z+1,mass_share=float(m[z]/m.sum()),mean_ceb_cap3=float(pz[z]@np.arange(4)/m[z]),childless=float(pz[z,0]/m[z]),two_plus=float(pz[z,2:].sum()/m[z])))
    # terciles of z by cumulative mass (states assigned by midpoint of cumulative share)
    cum=np.cumsum(m)/m.sum(); mid=cum-m/m.sum()/2; t=np.digitize(mid,[1/3,2/3])
    terc=[dict(tercile=k+1,mass_share=float(m[t==k].sum()/m.sum()),mean_ceb_cap3=float((pz[t==k].sum(0)@np.arange(4))/m[t==k].sum()),
               childless=float(pz[t==k,0].sum()/m[t==k].sum())) for k in range(3)]
    return dict(by_state=res,by_tercile=terc)
model=dict(end_of_22_25_cell=model_at(1),cell_42_45=model_at(6))
json.dump(dict(data_cps_june2024=data,model_chain13=model),open(os.path.join(OUT,"fertility_by_income.json"),"w"),indent=1)
fmt=lambda L:[{k:(round(v,3) if isinstance(v,float) else v) for k,v in d.items()} for d in L]
for k,v in data.items(): print("DATA",k); [print("  ",x) for x in fmt(v)]
for k,v in model.items(): print("MODEL",k,"terciles"); [print("  ",x) for x in fmt(v["by_tercile"])]; print("  by state ceb:",[round(x["mean_ceb_cap3"],2) for x in v["by_state"]])
