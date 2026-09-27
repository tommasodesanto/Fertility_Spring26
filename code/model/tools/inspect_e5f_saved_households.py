#!/usr/bin/env python3
"""Authenticated saved-state household audit; never solves the model."""
import argparse, gzip, hashlib, importlib.util, json, os, pickle, sys
from pathlib import Path
import numpy as np

def load(contract, case, output):
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    assert sha(contract)==os.environ['EXPECTED_UTILITY_OVERNIGHT_SHA256']
    c=json.loads(contract.read_text()); assert sha(c['files']['driver']['path'])==c['files']['driver']['sha256']
    spec=importlib.util.spec_from_file_location('household_audit_controller',c['files']['driver']['path'])
    d=importlib.util.module_from_spec(spec); sys.modules[spec.name]=d; spec.loader.exec_module(d)
    c,obj=d.verify(contract); d.verify_execution(c)
    r=d.read(case/'receipt.json'); assert r['target_system_sha256']==c['objective']['sha256']
    assert sha(case/'initial_state.pkl.gz')==r['case_checkpoint_sha256']
    output.mkdir(parents=True,exist_ok=True); (output/'runtime').mkdir(exist_ok=True)
    _,_,_,rt,_,_,_,_=d.setup(c,obj,r['point'],output/'runtime')
    with gzip.open(case/'initial_state.pkl.gz','rb') as f: packet=pickle.load(f)
    (output/'source_receipt.json').write_text(json.dumps(dict(case=str(case.resolve()),contract=str(contract.resolve()),contract_sha256=sha(contract),checkpoint_sha256=r['case_checkpoint_sha256'],target_system_sha256=r['target_system_sha256'],native_solves=0),indent=2)+'\n')
    return packet,rt,r

def weighted_quantile(x,w,qs):
    order=np.argsort(x); x=np.asarray(x)[order]; w=np.asarray(w)[order]
    return x[np.minimum(np.searchsorted(np.cumsum(w)/w.sum(),qs,side='left'),len(x)-1)].tolist()

def audit(packet,rt,out):
    import csv
    import matplotlib.pyplot as plt
    P=packet['parameters']; e=packet['evaluation']; pol=e.policy; bg=packet['b_grid']; g=e.g_current; total=g.sum()
    mass=g.sum(axis=tuple(range(1,7))); ages=P.age_start+P.da*np.arange(P.J)
    wealth=np.broadcast_to(bg[:,None,None,None,None,None,None],g.shape).copy()
    for t in range(1,g.shape[1]): wealth[:,t]+=pol.price[0]*P.H_own[t-1]
    # Net worth adds physical housing to financial assets at realized current timing.
    q=[.01,.1,.25,.5,.75,.9,.95,.99]
    nz=g>0
    result=dict(native_solves=0,current_mass=float(total),negative_mass=float(g[g<0].sum()),liquid_grid=[float(bg[0]),float(bg[-1])],quantiles=q,liquid_quantiles=weighted_quantile(bg,mass,q),current_net_worth_quantiles=weighted_quantile(wealth[nz],g[nz],q),bottom_liquid_grid_mass=float(mass[0]/total),top_liquid_grid_mass=float(mass[-1]/total),negative_liquid_mass=float(mass[bg<0].sum()/total),tenure_shares=(g.sum(axis=(0,2,3,4,5,6))/total).tolist(),children_capped_three_shares=(g.sum(axis=(0,1,2,3,4,6))/total).tolist(),grid_clipping_saving_mass=float(g[(pol.bp_pol<bg[0])|(pol.bp_pol>bg[-1])].sum()/total),nonpositive_consumption_mass=float(g[pol.c_pol<=0].sum()/total),nonfinite_occupied_consumption=int(np.sum(~np.isfinite(pol.c_pol[nz]))))
    result['fertility_top_bin_weight']=float(P.tfr_top_bin_weight)
    result['wealth_quantile_definition']='Discrete inverse CDF of saved grid distribution; no interpolation across atom at zero.'
    result['occupied_value_decreasing_steps']=int(np.sum((np.diff(pol.V,axis=0)<-1e-9)&((e.g_pre[:-1]>1e-12))))
    result['value_decrease_adjacent_mass']=float(np.sum(e.g_pre[:-1]*(np.diff(pol.V,axis=0)<-1e-9))/total)
    rows=[]
    for j,age in enumerate(ages):
        gj=g[:,:,:,j]; m=gj.sum(); bm=gj.sum(axis=(1,2,3,4,5))
        rows.append(dict(age=float(age),mass=float(m/total),ownership=float(gj[:,1:].sum()/m),mean_liquid=float(np.dot(bg,bm)/m),bottom_grid_share=float(bm[0]/m),top_grid_share=float(bm[-1]/m),mean_children_topbin_adjusted=float(np.dot(np.array([0.,1.,2.,float(P.tfr_top_bin_weight)]),gj.sum(axis=(0,1,2,3,5)))/m),mean_children_capped_three=float(np.dot(np.arange(g.shape[5]),gj.sum(axis=(0,1,2,3,5)))/m)))
    with (out/'distribution_by_age.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    (out/'distribution_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axs=plt.subplots(2,2,figsize=(11,8)); cutoff=weighted_quantile(bg,mass,[.99])[0]
    axs[0,0].plot(bg,np.cumsum(mass)/total);axs[0,0].set_xlim(min(bg[0],-1),cutoff);axs[0,0].set_title('Liquid wealth CDF (through 99th percentile)')
    axs[0,1].plot(ages,[r['ownership'] for r in rows]);axs[0,1].set_title('Ownership by age');axs[0,1].set_ylim(0,1)
    axs[1,0].bar(np.arange(len(result['tenure_shares'])),result['tenure_shares']);axs[1,0].set_title('Tenure shares: 0 renter, 1–5 owner products')
    for j in [3,6,10]:
        weights=g[:,0,0,j].sum(axis=(1,2,3));z,n,cs=np.unravel_index(np.argmax(g[:,0,0,j].sum(axis=0)),g.shape[4:])
        axs[1,1].plot(bg,pol.bp_pol[:,0,0,j,z,n,cs],label=f'age {ages[j]:g}, income {z}, children {n}')
    axs[1,1].set_xlim(0,cutoff);axs[1,1].set_ylim(0,cutoff);axs[1,1].legend(fontsize=7);axs[1,1].set_title('Savings policies: modal renter states (zoom)')
    fig.tight_layout();fig.savefig(out/'supplemental_distribution.png',dpi=150);plt.close(fig)
    return result

def sample_index(prob,rng):
    flat=np.asarray(prob).ravel(); assert np.min(flat)>=-1e-14 and abs(flat.sum()-1)<1e-9
    return np.unravel_index(rng.choice(flat.size,p=np.maximum(flat,0)/flat.sum()),prob.shape)

def simulate(packet,rt,out):
    import csv, inspect
    import matplotlib.pyplot as plt
    P=packet['parameters'];e=packet['evaluation'];p=e.policy;bg=packet['b_grid'];sh=packet['shared'];model=rt['model'];cal=rt['primitive'].pf.calendar
    assert P.I==1 and p.joint_choice is None
    rng=np.random.default_rng(20260927); shape=e.g_pre[:,:,:,0].shape; zero=np.zeros(shape)
    identity_choice=np.broadcast_to(np.arange(shape[1])[None,:,None,None,None,None,None],p.tenure_choice.shape)
    identity_idx=np.broadcast_to(np.minimum(np.arange(len(bg)),len(bg)-2),p.maps.tmx_idx.shape)
    identity_wt=np.broadcast_to((np.arange(len(bg))==len(bg)-1).astype(float),p.maps.tmx_wt.shape)
    zv,_,Pi=model.income_transition_values(P)
    def advance(cohort,j,identity=False):
        return model.advance_cohort_one_period_markov_income(cohort,j,p.loc_probs,identity_choice if identity else p.tenure_choice,None if identity else p.tenure_probs,p.bp_pol,P,bg,sh,p.maps.lmm_idx,p.maps.lmm_wt,identity_idx if identity else p.maps.tmx_idx,identity_wt if identity else p.maps.tmx_wt,bool(P.use_stochastic_aging),P.Pi_child,Pi)
    def current(cohort,j):
        return model.realize_current_choices_markov_income(cohort,j,p.loc_probs,p.tenure_choice,p.tenure_probs,p.maps.lmm_idx,p.maps.lmm_wt,p.maps.tmx_idx,p.maps.tmx_wt,use_compiled_scatter=bool(P.use_numba_scatter))
    # Composition equality verifies exact joint law for current-choice then saving/aging.
    errors=[]
    for j in [0,3,8]:
        cohort=zero.copy();idx=np.unravel_index(np.argmax(e.g_post_fertility[:,:,:,j]),shape);cohort[idx]=1
        errors.append(float(np.abs(advance(cohort,j)-advance(current(cohort,j),j,True)).sum()))
    assert max(errors)<1e-11, errors
    entry=e.g_pre[:,:,:,0]; bmass=entry.sum(axis=(1,2,3,4,5)); cdf=np.cumsum(bmass)/bmass.sum(); rows=[]
    for agent,q in enumerate([.1,.5,.9],1):
        b=int(np.searchsorted(cdf,q)); conditional=entry[b]/entry[b].sum(); state=(b,)+sample_index(conditional,rng)
        for j in range(P.J):
            full=np.zeros_like(e.g_pre); full[state[:3]+(j,)+state[3:]]=1
            born,_,_=cal.apply_fertility(full,p.fert_probs,P,cal.policy_continuation_birth_probs(p,P))
            post=sample_index(born[:,:,:,j],rng); cohort=zero.copy();cohort[post]=1
            realised=sample_index(current(cohort,j),rng);b,t,loc,z,n,cs=realised;idx=realised[:3]+(j,)+realised[3:]
            row=dict(agent=agent,entry_wealth_quantile=q,age=float(P.age_start+j*P.da),liquid_wealth=float(bg[b]),net_worth=float(bg[b]+(p.price[0]*P.H_own[t-1] if t else 0)),income_state=int(z),annual_gross_income=float(model.annual_gross_income_at_state(P,loc,j,zv[z])),children_capped_three=int(n),child_state=int(cs),tenure=int(t),housing_rooms=float(P.H_own[t-1] if t else p.hR_pol[idx]),consumption=float(p.c_pol[idx]),next_assets_policy=float(p.bp_pol[idx]),survives_next=False)
            survive=j<P.J-1 and (not P.use_age_survival or rng.random()<P.survival_probs[j]);row['survives_next']=bool(survive);rows.append(row)
            if not survive:break
            cohort=zero.copy();cohort[realised]=1;state=sample_index(advance(cohort,j,True),rng)
    with (out/'three_household_lives.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    fig,axs=plt.subplots(3,2,figsize=(11,10))
    for ax,key,title in zip(axs.ravel(),['annual_gross_income','net_worth','children_capped_three','housing_rooms','tenure','consumption'],['Annual income (gross earnings; pension in retirement)','Net worth','Children ever born (capped at three)','Housing rooms','Tenure: 0 renter, 1–5 owner product','Period consumption']):
        for a in range(1,4):
            rr=[r for r in rows if r['agent']==a];ax.step([r['age'] for r in rr],[r[key] for r in rr],where='post',label=f'Entry wealth q{[10,50,90][a-1]}')
        ax.set_title(title);ax.set_xlabel('Age');ax.legend(fontsize=7)
    fig.suptitle('Three illustrative simulated lives — selected entry ranks, not representative averages');fig.tight_layout();fig.savefig(out/'three_household_lives.png',dpi=150);plt.close(fig)
    receipt=dict(seed=20260927,agents=3,native_solves=0,transition_factorization_l1=errors,fertility_function=inspect.getsourcefile(cal.apply_fertility),forward_function=inspect.getsourcefile(model.advance_cohort_one_period_markov_income),description='Native fertility draw; native current-choice realization; native saving/income/child-aging kernel conditional on realized current state; survival sampled separately. Same interpolation lotteries as KFE. Entry wealth ranks selected, remaining entry states drawn conditionally. Nonrepresentative examples; period grid retained, no interpolation to annual histories.',row_count=len(rows))
    (out/'simulation_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--contract',type=Path,required=True);p.add_argument('--case',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); packet,rt,r=load(a.contract,a.case,a.output)
    print(audit(packet,rt,a.output));print(simulate(packet,rt,a.output))
