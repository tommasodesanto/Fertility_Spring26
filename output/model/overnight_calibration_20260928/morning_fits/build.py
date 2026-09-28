"""Torch-only saved-result fit plots; zero model solves; supplemental figures."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'): os.environ[k]='1'
os.environ['MPLBACKEND']='Agg'
import sys,json,csv,gzip,pickle,hashlib
from pathlib import Path
assert sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit()
ROOT=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
B=ROOT/'output/model/overnight_calibration_20260928';out=B/'morning_fits'
sys.path.insert(0,str(ROOT/'code/model/tools'))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import e5f_evening_calibration_runtime as runtime
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
c=json.loads((B/'contract_v1/contract.json').read_text())
assert sha(B/'contract_v1/contract.json')=='c83aaff1a90b1ba5bb0919151840e745e6816cb0a3a023ce746e750f77226c5d'
source=B/'gated_v1/search/selected_export/block'; receipt=json.loads((source/'receipt.json').read_text())
assert sha(source/'initial_state.pkl.gz')==receipt['case_checkpoint_sha256']
obj=json.loads(Path(c['lanes']['block']['objective']['path']).read_text())
evaluator=runtime.setup(dict(c,objective=c['lanes']['block']['objective']),obj,out/'runtime')
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
fit=pd.read_csv(out/'target_fit_primary_rescore.csv').set_index('moment')
y=next(r for r in rows if r['age_lower']==22);assert abs(.125*y['children_pre']+.875*y['children_post']-fit.loc['early_fertility','model'])<1e-10
assert abs(sum(r['wealth']*r['mass'] for r in rows)/sum(r['mass'] for r in rows if r['age_lower']<=65)-fit.loc['wealth_earnings','model'])<1e-9
emp=pd.read_csv(out/'empirical_profiles.csv');assert not (emp.series=='Model').any()
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
pdf=PdfPages(out/'fit_plots.pdf')
def save(fig,name):
 fig.savefig(out/(name+'.png'),dpi=150);pdf.savefig(fig);plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(12,8))
for ax,metric,title,label in zip(axs.flat,['children','ownership','rooms','wealth'],['Children ever born','Homeownership','Housing size','Net worth'],['Mean children, capped at 3','Share owning','Mean rooms, capped at 9','Mean / mean working-age earnings']):
 d=emp[emp.metric==metric]; mm=profiles[profiles.age_lower<=42] if metric=='children' else profiles
 ax.plot(mm.age_lower+1.5,mm[metric],'o-',label='Selected model',color='#166c8a')
 ax.plot(d.age_lower+1.5,d.value,'s--',label=d.series.iloc[0],color='#c56a29')
 if metric=='children':
  ax.scatter([25],[fit.loc['early_fertility','target']],marker='*',s=140,color='#c56a29',label='Target at exact age 25')
  ax.scatter([25],[fit.loc['early_fertility','model']],marker='*',s=140,color='#166c8a',label='Model at exact age 25')
  ax.set_xlim(18,45)
 if metric=='ownership':ax.set_ylim(0,1)
 ax.set(title=title,xlabel='Age (four-year cell midpoint)',ylabel=label);ax.grid(alpha=.2);ax.legend(fontsize=8,frameon=False)
fig.suptitle('Lifecycle comparison: selected September 28 calibration',fontsize=16)
fig.text(.5,.025,'Age profiles are untargeted diagnostics; stars show the targeted age-25 moment. Lines join four-year cells.\nACS housing samples differ from aggregate calibration samples. Children use literal cap 3; fertility normalization uses a top-bin adjustment.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.08,1,.94]);save(fig,'lifecycle_fit')
labels=['Childlessness, ages 40–44','One child among mothers','Mean first-birth age','Wealth / annual earnings','Annual bequests / wealth','Mean rooms','Ownership, ages 30–55','First-birth rooms response','Recent-parent ownership gap','Children ever born at 25']
names=fit[fit.role=='scored'].index.tolist()
fig,axs=plt.subplots(2,5,figsize=(15,6))
for ax,name,label in zip(axs.flat,names,labels):
 r=fit.loc[name]; ax.bar([0,1],[r.target,r.model],color=['#c56a29','#166c8a'],width=.6);ax.set_xticks([0,1],['Data','Model']);ax.set_title(label,fontsize=10,wrap=True);ax.grid(axis='y',alpha=.15)
 for i,v in enumerate([r.target,r.model]):ax.annotate(f'{v:.3f}',(i,v),xytext=(0,4),textcoords='offset points',ha='center',fontsize=9)
 ax.set_ylim(0,max(r.target,r.model)*1.22)
fig.suptitle('Ten targeted moments — original units, separate panel scales',fontsize=16);fig.tight_layout(rect=[0,0,1,.92]);save(fig,'targeted_fits')
fig,axs=plt.subplots(1,3,figsize=(11,4))
for ax,name,title in zip(axs,['nchs_share30','old_dispersion','family_rooms'],['First births at age 30+ (share)','Older wealth/income: p90 / median','Rooms gap: 3+ vs 1–2 children']):
 r=fit.loc[name];ax.bar([0,1],[r.target,r.model],color=['#c56a29','#166c8a'],width=.6);ax.set_xticks([0,1],['Data','Model']);ax.set_title(title,fontsize=10)
 for i,v in enumerate([r.target,r.model]):ax.annotate(f'{v:.3f}',(i,v),xytext=(0,4),textcoords='offset points',ha='center')
 ax.set_ylim(0,max(r.target,r.model)*1.2)
fig.suptitle('Untargeted checks — zero calibration weight',fontsize=15);fig.tight_layout(rect=[0,0,1,.9]);save(fig,'untargeted_fits');pdf.close()
qa=dict(status='saved_state_only',model_solves=0,source=str(source),checkpoint_sha256=receipt['case_checkpoint_sha256'],early_target_replay=True,wealth_target_replay=True,figures=['lifecycle_fit','targeted_fits','untargeted_fits'],empirical_profiles_sha256=sha(out/'empirical_profiles.csv'))
(out/'qa.json').write_text(json.dumps(qa,indent=2)+'\n');print(json.dumps(qa))
