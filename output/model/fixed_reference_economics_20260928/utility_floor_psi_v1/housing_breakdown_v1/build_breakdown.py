"""Read-only saved-array housing decomposition; no lifecycle/model solve."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='1'
import csv,json,hashlib,pathlib
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[5]
BASE=pathlib.Path(__file__).resolve().parent.parent
SRC=BASE/'chain_11/postcheck_results/selected_postcheck'
OUT=pathlib.Path(__file__).resolve().parent
ARR=SRC/'phase_b_ge/selected_repeat/stage/solution_arrays.npz'
with np.load(ARR,allow_pickle=False) as a:
    g=a['g'];h=a['hR_pol'];saved_own=float(g[:,1:].sum()/g.sum());saved_hd=float(a['housing_demand'].sum());houses=np.array([2.,4.,6.,8.,10.])
assert g.shape==h.shape==(120,6,1,17,9,4,4)
assert np.isfinite(g).all() and np.min(g)>=0
obs=json.loads((SRC/'phase_b_ge/selected_root/observers.json').read_text())['housing_wealth']
fit=list(csv.DictReader((SRC/'phase_b_ge/selected_root/target_fit.csv').open()))
finalfit=list(csv.DictReader((SRC/'phase_b_ge/selected_repeat_final/target_fit.csv').open()))
assert fit==finalfit
rows=[]
def measure(label,ages,child=None,ever=None):
    # Saved g is post-housing g_current; only its actual renter slice uses hR.
    mask=np.ones((4,4),bool)
    for n in range(4):
      for cs in range(4):
        m=min(n,cs,3)
        mask[n,cs]=(child is None or m==child) and (ever is None or n==ever)
    age=np.asarray(ages)[None,None,:,None,None,None]
    rr=g[:,0]*mask*age
    hr=np.where(rr>0,h[:,0],0)
    rm=float(rr.sum());rs=float((rr*hr).sum());cap=float((rr*(np.abs(hr-6)<=1e-8)).sum())
    om=ors=0.
    for t,size in enumerate(houses,1):
      mass=float((g[:,t]*mask*age).sum());om+=mass;ors+=mass*size
    mass=rm+om
    return dict(group=label,population_mass=mass,population_share=mass/float(g.sum()),ownership=om/mass if mass else None,renter_mean_rooms=rs/rm if rm else None,owner_mean_rooms=ors/om if om else None,mean_rooms=(rs+ors)/mass if mass else None,renter_cap_fraction=cap/rm if rm else None,rental_room_contribution_to_aggregate=rs/float(g.sum()),owner_room_contribution_to_aggregate=ors/float(g.sum()),total_room_contribution_to_aggregate=(rs+ors)/float(g.sum()))
rows.append(measure('All ages 18-85',np.ones(17)))
for m in range(4):rows.append(measure('All ages; current children '+str(m),np.ones(17),m))
for lo,hi in [(18,30),(30,46),(46,66),(66,86),(25,35),(30,56)]:
 age=np.maximum(0,np.minimum(18+np.arange(17)*4+4,hi)-np.maximum(18+np.arange(17)*4,lo))/4
 rows.append(measure(f'Ages {lo}-{hi-1}',age))
 for m in range(4):rows.append(measure(f'Ages {lo}-{hi-1}; current children {m}',age,m))
 if lo in [18,25,30]:rows.append(measure(f'Ages {lo}-{hi-1}; never had children',age,None,0))
for j in range(17):
 age=np.zeros(17);age[j]=1
 for m in range(4):rows.append(measure(f'Age cell {18+4*j}-{21+4*j}; current children {m}',age,m))
allrow=rows[0]
assert abs(allrow['mean_rooms']-obs['moments']['aggregate_mean_occupied_rooms_ahs_uncapped_18_85'])<1e-12
assert abs(allrow['ownership']-saved_own)<1e-12
assert abs(allrow['mean_rooms']-saved_hd)<1e-12
assert abs(sum(r['total_room_contribution_to_aggregate'] for r in rows[1:5])-allrow['mean_rooms'])<1e-12
with (OUT/'aggregate_housing.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
selected=[r for r in rows if r['group'] in ['All ages 18-85','All ages; current children 0','All ages; current children 1','All ages; current children 2','All ages; current children 3','Ages 18-29','Ages 18-29; never had children','Ages 25-34; never had children','Ages 30-45','Ages 46-65','Ages 66-85']]
report=['# Housing decomposition of verified local chain 11 postcheck','','Original-weight loss '+str(sum(float(r['loss_contribution']or 0) for r in fit))+'. Experimental candidate; not adopted. q=0.689763713921058.','Saved g is the realized post-housing distribution. Renter housing uses only g[:,0] and hR[:,0]; owners use physical room rungs [2,4,6,8,10]. Current children m=min(n,cs,3) under pinned independent_count mode. Age intervals use uniform four-year-cell overlap.','','| Group | Mass share | Ownership | Renter rooms | Owner rooms | Renter cap % | Aggregate room contribution |','|---|---:|---:|---:|---:|---:|---:|']
for r in selected:
 report.append('|'+r['group']+'|'+ '|'.join(f'{100*r[k]:.2f}%' if k in ['population_share','ownership','renter_cap_fraction'] else f'{r[k]:.4f}' for k in ['population_share','ownership','renter_mean_rooms','owner_mean_rooms','renter_cap_fraction','total_room_contribution_to_aggregate'])+'|')
branch=next(r['branch']for r in obs['rows'] if 'branch'in r)
receipt=dict(source=str(ARR),sha256=hashlib.sha256(ARR.read_bytes()).hexdigest(),repeat_receipt=str(SRC/'native_selected_repeat.json'),fit_identity_exact=True,aggregate_rooms_matches_saved_observer=True,aggregate_ownership_matches_saved_g=True,aggregate_rooms_matches_native_demand=True,current_children_contributions_add_exact=True,original_loss=sum(float(r['loss_contribution']or 0)for r in fit),branch_summary=branch,no_model_solves=True)
(OUT/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
report+=['','The existing matched first-birth observer gives control '+str(branch['control_mean_housing'])+' rooms and treated '+str(branch['treated_mean_housing'])+' rooms, a '+str(branch['housing_response'])+' room response one four-year period later. These include renters and owners. Compact observer JSON retains branch levels but not cohort distributions; its mean near six does not itself establish a renter cap.','','Complete target fit and all 31 parameters remain in '+str(SRC/'phase_b_ge/selected_root')+'.','']
(OUT/'README.md').write_text('\n'.join(report))
for r in selected:print(r)
print('OUT',OUT)
