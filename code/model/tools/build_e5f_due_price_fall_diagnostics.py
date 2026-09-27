"""Supplemental origin-specific policies from saved dated experiments; no solves.
Run with all numerical-library threads set to one and --root pointing at the
existing due_price_fall folder. Original prechoice owner mass weights conditional
keep-house policies; realized DUE stayer mass is separately labelled.
"""
import argparse,csv,gzip,hashlib,json,pickle
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
class SavedObjects(pickle.Unpickler):
 def find_class(self,module,name):
  if module.startswith(('run_','intergen_','e5f_','demographic_')):return SimpleNamespace
  return super().find_class(module,name)
def load(p):
 with gzip.open(p,'rb') as f:return SavedObjects(f).load()
def main(root):
 out=root/'supplement';out.mkdir(exist_ok=False)
 arms={};hashes={}
 for arm in ('baseline','due'):
  folder=root/(arm+'_v1');accounts=json.loads((folder/'dated_accounts.json').read_text())
  assert all(accounts['gates'].values()) and np.isfinite(accounts['estate']['estate']['totals']['net_negative']) and accounts['estate']['estate']['totals']['net_negative']<=1e-10
  pp=folder/'policy_and_original_state.pkl.gz';op=folder/'dated_operator.pkl.gz'
  hashes[str(pp)]=sha(pp);hashes[str(op)]=sha(op)
  assert hashes[str(pp)]==accounts['checkpoint_sha256']
  a=load(pp);d=load(op);np.testing.assert_array_equal(a['g_pre'],d['evaluation'].g_pre)
  for key in ('bp_pol','c_pol'):
   np.testing.assert_array_equal(getattr(a['policy'],key),getattr(d['evaluation'].policy,key))
  arms[arm]=(a,d,accounts)
 a,d,_=arms['due'];base=arms['baseline'][0];P=a['parameters'];g=a['g_pre'];bg=a['b_grid'];policy=a['policy']
 np.testing.assert_array_equal(g,base['g_pre']);np.testing.assert_array_equal(bg,base['b_grid'])
 J=P.J;price=float(a['price']);phi=a['shared'].phi_choice
 L=np.zeros_like(g);F=np.zeros_like(g);deathfloor=np.zeros_like(g);elig=np.zeros_like(g,dtype=bool)
 for t,h in enumerate(P.H_own,1):
  for j in range(J):
   l=-np.asarray(phi[0,t])*price*h;floor=np.minimum(bg[:,None,None,None],l[None,None,:,:])
   mortality=j==J-1 or (P.use_age_survival and P.survival_probs[j]<1)
   df=-(1-P.psi)*price*h if mortality else -np.inf
   L[:,t,0,j]=l;F[:,t,0,j]=np.maximum(floor,df);deathfloor[:,t,0,j]=df
   elig[:,t,0,j]=bg[:,None,None,None]<l
 own=g.copy();own[:,0]=0;stay=d['evaluation'].g_stay_distribution
 relaxed=policy.bp_pol_stay<L-1e-9;changed=abs(policy.bp_pol_stay-base['policy'].bp_pol)>1e-8
 binding=abs(policy.bp_pol_stay-F)<=1e-8
 deathbind=np.isfinite(deathfloor)&(abs(policy.bp_pol_stay-deathfloor)<=1e-8)
 summary=dict(original_owner_mass=float(own.sum()),original_owner_mass_eligible_below_current_ltv=float(own[elig].sum()),
  original_owner_mass_with_different_conditional_saving=float(own[changed].sum()),
  original_owner_mass_with_conditional_saving_below_current_ltv=float(own[relaxed].sum()),
  realized_due_stayer_mass=float(stay.sum()),realized_due_stayer_mass_below_current_ltv=float(stay[relaxed].sum()),
  realized_due_stayer_mass_at_total_floor=float(stay[binding].sum()),realized_due_stayer_mass_at_death_floor=float(stay[deathbind].sum()),
  negative_estates={k:arms[k][2]['estate']['estate']['totals']['net_negative'] for k in arms},
  classification='One dated partial-equilibrium operator, permanent lower-price expectations; no market-clearing claim',
  conditional_caveat='Original owner mass weights keep-the-same-house policies, not actual household choices. Realized stayer statistics use separate post-fertility origin mass. Changes can reflect continuation values as well as the immediate borrowing floor.')
 ages=np.array([P.age_start+j*P.da for j in range(J)])
 # Median persistent-income support is selected by its supplied type weights if available.
 zw=np.asarray(getattr(P,'z_weights',np.ones(len(P.z_grid))/len(P.z_grid)))
 z=int(np.searchsorted(np.cumsum(zw/zw.sum()),.5)); z=min(z,len(P.z_grid)-1)
 summary['income_index']=z;summary['income_state']=float(P.z_grid[z]);summary['income_selection']='weighted median if z_weights supplied; otherwise median support index'
 fig,axs=plt.subplots(3,2,figsize=(11,10),layout='constrained');rows=[];selections=[]
 for rr,age in enumerate((30,46,70)):
  j=int(np.argmin(abs(ages-age)));weights=own[:,:,0,j,z,:,:].sum(axis=0);t,n,m=np.unravel_index(np.argmax(weights),weights.shape)
  mass=own[:,t,0,j,z,n,m];ix=(slice(None),t,0,j,z,n,m)
  if mass.sum()<=0:raise ValueError('No original owner mass in selected income/age slice')
  cdf=np.cumsum(mass)/mass.sum();lo=max(0,int(np.searchsorted(cdf,.002))-1);hi=min(len(bg)-1,int(np.searchsorted(cdf,.995))+1)
  xx=bg[lo:hi+1];sel=slice(lo,hi+1);b0=base['policy'].bp_pol[ix];b1=policy.bp_pol_stay[ix];c0=base['policy'].c_pol[ix];c1=policy.c_pol_stay[ix]
  selection=dict(age=float(ages[j]),income_index=z,tenure=int(t),rooms=float(P.H_own[t-1]),children_ever_born=int(n),children_home_state=int(m),original_slice_mass=float(mass.sum()),plotted_mass_share=float(mass[sel].sum()/mass.sum()));selections.append(selection)
  for k,(v0,v1,label) in enumerate(((b0,b1,'End-period financial assets'),(c0,c1,'Nonhousing consumption'))):
   ax=axs[rr,k];ax.plot(xx,v0[sel],label='Baseline keep-house policy',color='#4361a8');ax.plot(xx,v1[sel],label='DUE keep-house policy',color='#ce6b24')
   ax.scatter(xx,v1[sel],s=4+60*mass[sel]/mass.max(),color='#ce6b24',alpha=.45)
   if k==0:
    ax.plot(xx,L[ix][sel],color='gray',ls='--',label='Purchase collateral floor');ax.plot(xx,F[ix][sel],color='black',ls=':',label='DUE floor incl. death solvency')
   ax.set_title(f'Age {ages[j]:g}, {P.H_own[t-1]:g} rooms, children {n}/{m} (ever/home state)',fontsize=10)
   ax.set_xlabel('Beginning financial assets');ax.set_ylabel(label);ax.grid(alpha=.2)
   if rr==0:ax.legend(fontsize=7)
  for bi in range(len(bg)):
   rows.append(dict(**selection,wealth=float(bg[bi]),original_owner_mass=float(mass[bi]),baseline_saving=float(b0[bi]),due_stayer_saving=float(b1[bi]),baseline_consumption=float(c0[bi]),due_stayer_consumption=float(c1[bi]),collateral_floor=float(L[ix][bi]),due_total_floor=float(F[ix][bi]),death_floor=float(deathfloor[ix][bi])))
 fig.suptitle('Keeping the same house after the price fall: baseline vs DUE\nConditional policies; marker size shows original owner mass',fontsize=13)
 fig.savefig(out/'owner_conditional_policies.png',dpi=150);plt.close(fig)
 with (out/'owner_conditional_policies.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 summary['slices']=selections;summary['source_sha256']=sha(__file__);summary['input_sha256']=hashes
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 (out/'README.md').write_text('Supplemental origin-specific DUE diagnostic.\n\n'+summary['conditional_caveat']+'\n\nRows select the most populated original owner product/family state at ages30,46,70 and median income support; horizontal range retains at least99.3% slice mass. Standard17plots unchanged.\n\nCommand: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python code/model/tools/build_e5f_due_price_fall_diagnostics.py --root output/model/daytime_calibration_20260927/due_price_fall\n')
 print(json.dumps({k:v for k,v in summary.items() if k not in ('input_sha256','slices')},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args().root)
