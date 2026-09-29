"""Torch-only saved-result fit plots; zero model solves; supplemental figures."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'): os.environ[k]='1'
os.environ['MPLBACKEND']='Agg'
import sys,json,csv,gzip,pickle,hashlib
from pathlib import Path
assert sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit()
ROOT=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
B=ROOT/'output/model/fertility_identification_20260928';out=B/'numerical_pair_v1'
EMP=ROOT/'output/model/overnight_calibration_20260928/morning_fits/empirical_profiles.csv'
sys.path.insert(0,str(ROOT/'code/model/tools'))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import e5f_evening_calibration_runtime as runtime
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
c=json.loads((B/'contract_v1/contract.json').read_text())
assert sha(B/'contract_v1/contract.json')=='68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf'
source=out/'run_v1/retained_start/case'; receipt=json.loads((source/'receipt.json').read_text())
assert sha(source/'receipt.json')=='6c051eea4f7899fb4eee0a320fc0402972dc62ff0789d91858f4ff0459b0fa1a'
assert receipt['case_checkpoint_sha256']=='9d36f69a34dae685ed8a43abc0865896040eafb682d1794cecad569441a45881'
assert sha(source/'initial_state.pkl.gz')==receipt['case_checkpoint_sha256']
obj=json.loads(Path(c['lanes']['primary']['objective']['path']).read_text())
evaluator=runtime.setup(dict(c,objective=c['lanes']['primary']['objective']),obj,out/'lifecycle_fit_runtime')
with gzip.open(source/'initial_state.pkl.gz','rb') as f: packet=pickle.load(f)
P=packet['parameters']; e=packet['evaluation']; pol=e.policy;g=e.g_current;bg=packet['b_grid'];model=evaluator.rt['model'];zv,_,_=model.income_transition_values(P)
rows=[];en=0.;em=0.
for j in range(P.J):
 a=P.age_start+P.da*j;gj=g[:,:,:,j];m=gj.sum();tm=gj.sum(axis=(0,2,3,4,5));zm=gj.sum(axis=(0,1,2,4,5))
 asset=e.g_post_fertility[:,:,:,j];ab=asset.sum(axis=(1,2,3,4,5));at=asset.sum(axis=(0,2,3,4,5))
 nw=(np.dot(ab,bg)+np.dot(at[1:],np.asarray(P.H_own)*pol.price[0]))/asset.sum()
 room=((gj[:,0]*np.minimum(pol.hR_pol[:,0,:,j],9)).sum()+np.dot(tm[1:],np.minimum(P.H_own,9)))/m
 def nc(gg):return np.dot(np.arange(gg.shape[4]),gg.sum(axis=(0,1,2,3,5)))/gg.sum()
 pre=nc(e.g_pre[:,:,:,j]);post=nc(gj)
 rows.append(dict(age_lower=a,age_upper=a+3,ownership=tm[1:].sum()/m,rooms=room,wealth=nw,children=(pre+post)/2,children_pre=pre,children_post=post,mass=m))
 if a<=65: en+=np.dot(zm,[model.annual_gross_income_at_state(P,0,j,z) for z in zv]);em+=m
for r in rows:r['wealth']/=en/em
profiles=pd.DataFrame(rows);profiles.to_csv(out/'model_profiles.csv',index=False)
fit=pd.read_csv(source/'target_fit.csv').set_index('moment')
y=next(r for r in rows if r['age_lower']==22);assert abs(.125*y['children_pre']+.875*y['children_post']-fit.loc['early_fertility','model'])<1e-10
assert abs(sum(r['wealth']*r['mass'] for r in rows)/sum(r['mass'] for r in rows if r['age_lower']<=65)-fit.loc['wealth_earnings','model'])<1e-9
assert sha(EMP)=='50f434a81cbd159de0184674bd2a95719b785777c48c7abc7afeb9bd54910351'
emp=pd.read_csv(EMP);assert not (emp.series=='Model').any()
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
 fig.savefig(out/(name+'.png'),dpi=150);plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(12,8))
for ax,metric,title,label in zip(axs.flat,['children','ownership','rooms','wealth'],['Children ever born','Homeownership','Housing size','Net worth'],['Mean children, capped at 3','Share owning','Mean rooms, capped at 9','Mean / mean working-age earnings']):
 d=emp[emp.metric==metric]; mm=profiles[profiles.age_lower<=42] if metric=='children' else profiles
 ax.plot(mm.age_lower+1.5,mm[metric],'o-',label='Model',color='#166c8a')
 ax.plot(d.age_lower+1.5,d.value,'s--',label=d.series.iloc[0],color='#c56a29')
 if metric=='children':
  ax.scatter([25],[fit.loc['early_fertility','target']],marker='*',s=140,color='#c56a29',label='Target at exact age 25')
  ax.scatter([25],[fit.loc['early_fertility','model']],marker='*',s=140,color='#166c8a',label='Model at exact age 25')
  ax.set_xlim(18,45)
 if metric=='ownership':ax.set_ylim(0,1)
 ax.set(title=title,xlabel='Age (four-year cell midpoint)',ylabel=label);ax.grid(alpha=.2);ax.legend(fontsize=8,frameon=False)
fig.suptitle('Lifecycle comparison: model versus data',fontsize=16)
fig.text(.5,.025,'2007 stationary approximation; experimental calibration matching the table. Cross-sectional profiles, not a transition.\nACS housing samples differ from aggregate calibration samples. Children use literal cap 3; fertility normalization uses a top-bin adjustment.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.08,1,.94]);save(fig,'lifecycle_fit')
qa=dict(status='passed',model_solves=0,source=str(source),checkpoint_sha256=receipt['case_checkpoint_sha256'],receipt_sha256=sha(source/'receipt.json'),early_target_replay=True,wealth_target_replay=True,empirical_profiles_sha256=sha(EMP),output_sha256={n:sha(out/n) for n in ['lifecycle_fit.png','model_profiles.csv']})
(out/'lifecycle_fit_qa.json').write_text(json.dumps(qa,indent=2)+'\n');print(json.dumps(qa))
