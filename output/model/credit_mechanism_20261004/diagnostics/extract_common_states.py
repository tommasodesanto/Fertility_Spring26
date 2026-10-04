"""Cached paired birth decomposition. No solver or aggregate-policy imports.

Inverts the exact independent-count sequential birth map before tenure, checks
the inversion against native birth flows, and decomposes changes in birth flows
into common-policy and composition terms. Never clips probabilities or mass.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DEFAULT = ROOT / 'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/credit_relaxation'

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''): h.update(block)
    return h.hexdigest()

def csvwrite(name, rows):
    with (HERE / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def load(path):
    p = json.loads((path / 'executed_P.json').read_text())
    for key, value in {'sequential_births':True, 'joint_nested_choice':False,
                       'readiness_gate_enabled':False, 'child_state_mode':'independent_count'}.items():
        if p[key] != value: raise ValueError(f'Unsupported birth map: {key}={p[key]}')
    names = ('b_grid','p_eq','g_beginning_distribution','fert_probs','fert2_probs','fert_value','V','type_values')
    with np.load(path / 'solution_arrays.npz', allow_pickle=False) as z:
        a = {k: z[k].copy() for k in names}
    post = a['g_beginning_distribution']
    pre = post.copy(); rate = np.zeros_like(post)
    ages = p['age_start'] + np.arange(p['J']) * p['da']
    fec = np.clip(1 - p['fecundity_omega1'] * np.exp(p['fecundity_omega2']*(ages-p['age_start'])),0,1)
    if p['fecundity_omega1'] == 0: fec[:] = 1
    if p['fecundity_terminal_decay']:
        fec *= np.exp(-p['fecundity_terminal_decay'] * np.maximum(ages-p['fecundity_tail_start_age'],0))
    if p['fecundity_omega1'] != 0: fec[ages >= p['fecundity_terminal_age']] = 0
    for j in range(p['A_f_start']-1,p['A_f_end']):
        rate[:,:,:,j,:,0,0] = fec[j]*a['fert_probs'][:,:,:,j,:,1]
        for n in (1,2):
            for m in range(n+1):
                rate[:,:,:,j,:,n,m] = fec[j]*a['fert2_probs'][:,:,:,j,:,1,n-1,m]
        # Snapshot-based native map has inflow from the prebirth lower count.
        for n in range(4):
            for m in range(n+1):
                inflow = 0 if n == 0 or m == 0 else pre[:,:,:,j,:,n-1,m-1]*rate[:,:,:,j,:,n-1,m-1]
                den = 1-rate[:,:,:,j,:,n,m]
                if np.any(den <= 0): raise ValueError('Non-invertible birth map; needs native exposure capture')
                pre[:,:,:,j,:,n,m] = (post[:,:,:,j,:,n,m]-inflow)/den
    reconstructed = pre*(1-rate)
    for n in range(1,4):
        for m in range(1,n+1): reconstructed[...,n,m] += pre[...,n-1,m-1]*rate[...,n-1,m-1]
    flow = (pre*rate).sum(axis=(0,1,2,4,6))
    native = np.array([p['_first_births_by_age'],p['_second_births_by_age'],p['_third_births_by_age']]).T
    # Signed roundoff remains in the inversion and is explicitly reported.
    checks = {'minimum_reconstructed_prebirth_mass':float(pre.min()),
              'negative_prebirth_mass':float(-pre[pre<0].sum()),
              'forward_reconstruction_max_error':float(np.max(abs(reconstructed-post))),
              'native_flow_max_error':float(np.max(abs(flow[:,:3]-native))),
              'age_mass_max_error':float(np.max(abs(pre.sum(axis=(0,1,2,4,5,6))-post.sum(axis=(0,1,2,4,5,6)))))}
    if pre.min() < -1e-12 or checks['native_flow_max_error'] > 1e-10 or checks['forward_reconstruction_max_error']>1e-10:
        raise ValueError(f'Exposure inversion fails: {checks}')
    return p,a,pre,rate,ages,fec,checks

def main():
    global HERE
    ap=argparse.ArgumentParser(); ap.add_argument('--baseline',type=Path,default=DEFAULT/'phi_080'); ap.add_argument('--relaxed',type=Path,default=DEFAULT/'phi_100'); ap.add_argument('--output',type=Path,default=HERE); ap.add_argument('--kind',choices=['credit','price110'],default='credit'); args=ap.parse_args()
    HERE=args.output; HERE.mkdir(parents=True,exist_ok=True)
    p0,a0,g0,h0,ages,fec,check0=load(args.baseline)
    p1,a1,g1,h1,_,_,check1=load(args.relaxed)
    changes=[k for k in p0 if not k.startswith('_') and p0[k]!=p1[k]]
    pair_ok=(changes==['phi'] and np.array_equal(a0['p_eq'],a1['p_eq'])) if args.kind=='credit' else (changes==[] and np.array_equal(1.1*a0['p_eq'],a1['p_eq']))
    if not pair_ok or not np.array_equal(a0['b_grid'],a1['b_grid']):
        raise ValueError(f'Pair differs from explicit {args.kind} contract: {changes}')
    # Common fixed income boundaries use baseline childless prebirth exposure.
    earnings=np.asarray(p0['income'])[0,:,None]*a0['type_values'][None,:]/(p0['period_years']*(1-p0['tau_pay']))
    fertile=np.arange(p0['J'])<p0['A_f_end']
    income_weights=g0[...,0,0].sum(axis=(0,1,2)); income_weights[~fertile]=0
    x=earnings.ravel(); w=income_weights.ravel(); order=np.argsort(x,kind='stable'); cw=np.cumsum(w[order])
    cuts=[float(x[order[np.searchsorted(cw,q*cw[-1])]]) for q in (1/3,2/3)]
    income_bin=np.searchsorted(cuts,earnings,side='left')
    wealth_bin=np.searchsorted([0.,1.],a0['b_grid'],side='left')
    grouped=[]; totals=[]; gap_rows=[]; debt_rows=[]
    for n in range(3):
        for m in range(n+1):
            # n=0,m=0 is first birth; n>0 are next-child attempts.
            if n==0:
                pr0=a0['fert_probs'][...,:2]; pr1=a1['fert_probs'][...,:2]
                iv0=a0['fert_value']; iv1=a1['fert_value']; scale=p0['kappa_fert']
            else:
                pr0=a0['fert2_probs'][...,n-1,m]; pr1=a1['fert2_probs'][...,n-1,m]
                iv0=a0['V'][...,n,m]; iv1=a1['V'][...,n,m]; scale=p0['kappa_fert_continuation'] or p0['kappa_fert']
            valid=np.all(pr0>0,axis=-1)&np.all(pr1>0,axis=-1)&np.all(np.isfinite(pr0),axis=-1)&np.all(np.isfinite(pr1),axis=-1)
            with np.errstate(divide='ignore',invalid='ignore'):
                wait_gain=(iv1+scale*np.log(pr1[...,0]))-(iv0+scale*np.log(pr0[...,0]))
                try_gain=(iv1+scale*np.log(pr1[...,1]))-(iv0+scale*np.log(pr0[...,1]))
                gap0=scale*(np.log(pr0[...,1])-np.log(pr0[...,0])); gap1=scale*(np.log(pr1[...,1])-np.log(pr1[...,0])); dgap=gap1-gap0
            # Soft-logit inversion identity: delta try - delta wait = delta gap.
            identity=float(np.max(abs(try_gain[valid]-wait_gain[valid]-dgap[valid])))
            for j,age in enumerate(ages):
                if not fertile[j]: continue
                b0=g0[:,:,:,j,:,n,m]; b1=g1[:,:,:,j,:,n,m]; r0=h0[:,:,:,j,:,n,m]; r1=h1[:,:,:,j,:,n,m]
                policy=float(np.sum(b0*(r1-r0))); composition=float(np.sum((b1-b0)*r1)); actual=float(np.sum(b1*r1)-np.sum(b0*r0))
                factor=1+(p0['tfr_top_bin_weight']-3) if n==2 else 1
                totals.append(dict(age=int(age),children_ever_born=n,children_at_home=m,exposure0=float(b0.sum()),exposure1=float(b1.sum()),birth_flow0=float(np.sum(b0*r0)),birth_flow1=float(np.sum(b1*r1)),policy=policy,composition=composition,shapley_policy=float(np.sum(.5*(b0+b1)*(r1-r0))),shapley_composition=float(np.sum(.5*(r0+r1)*(b1-b0))),top_bin_birth_factor=factor,adjusted_birth_flow0=factor*float(np.sum(b0*r0)),adjusted_birth_flow1=factor*float(np.sum(b1*r1)),adjusted_policy=factor*policy,adjusted_composition=factor*composition,total=actual,identity_error=actual-policy-composition))
                for ten in range(b0.shape[1]):
                    for debt in (True,False):
                        mask=(a0['b_grid']<0) if debt else (a0['b_grid']>=0)
                        u0=b0[mask,ten,0,:]; u1=b1[mask,ten,0,:]; d0=r0[mask,ten,0,:]; d1=r1[mask,ten,0,:]
                        if u0.sum()==0 and u1.sum()==0: continue
                        debt_rows.append(dict(age=int(age),beginning_tenure=ten,children_ever_born=n,children_at_home=m,assets='debt' if debt else 'nonnegative',exposure0=float(u0.sum()),exposure1=float(u1.sum()),exposure_change=float(u1.sum()-u0.sum()),baseline_zero_new_exposure=float(u1[u0==0].sum()),birth_flow0=float(np.sum(u0*d0)),birth_flow1=float(np.sum(u1*d1)),policy=float(np.sum(u0*(d1-d0))),composition=float(np.sum((u1-u0)*d1)),total=float(np.sum(u1*d1)-np.sum(u0*d0))))
                for ten in range(b0.shape[1]):
                    for inc in range(3):
                        for wb in range(3):
                            mask=(wealth_bin[:,None]==wb)&(income_bin[j,None,:]==inc)
                            w0=b0[:,ten,0,:]*mask; w1=b1[:,ten,0,:]*mask; mass=float(w0.sum()); mass1=float(w1.sum())
                            if mass<=0 and mass1<=0: continue
                            v=valid[:,ten,0,j,:]&mask; support=float(w0[v].sum()); unsupported=float(w0[~v].sum())
                            def av(arr): return float(np.sum(w0[v]*arr[:,ten,0,j,:][v])/support) if support>0 else None
                            delta=av(dgap); wait=av(wait_gain)
                            row=dict(age=int(age),beginning_tenure=ten,children_ever_born=n,children_at_home=m,income_group=inc,wealth_group=wb,exposure0=mass,exposure1=mass1,birth_rate0=float(np.sum(w0*r0[:,ten,0,:])/mass) if mass>0 else None,birth_rate1_common=float(np.sum(w0*r1[:,ten,0,:])/mass) if mass>0 else None,policy=float(np.sum(w0*(r1[:,ten,0,:]-r0[:,ten,0,:]))),composition=float(np.sum((w1-w0)*r1[:,ten,0,:])),supported_gap_mass=support,unsupported_gap_mass=unsupported,gap0=av(gap0),gap1=av(gap1),delta_gap=delta,success_gap_delta=delta/fec[j] if delta is not None and fec[j]>0 else None,wait_credit_gain=wait,try_credit_gain=av(try_gain),success_credit_gain=wait+delta/fec[j] if delta is not None and wait is not None and fec[j]>0 else None)
                            grouped.append(row)
            gap_rows.append(dict(children_ever_born=n,children_at_home=m,logit_identity_error=identity))
    csvwrite('flow_decomposition.csv',totals); csvwrite('common_state_groups.csv',grouped)
    csvwrite('age_debt_tenure_decomposition.csv',debt_rows)
    compact=[]
    for inc in range(3):
        for wb in range(3):
            rr=[r for r in grouped if r['children_ever_born']==0 and r['income_group']==inc and r['wealth_group']==wb]
            mass=sum(r['exposure0'] for r in rr); support=sum(r['supported_gap_mass'] for r in rr)
            row=dict(income_group=inc,wealth_group=wb,exposure0=mass,birth_rate0=sum(r['birth_rate0']*r['exposure0'] for r in rr if r['birth_rate0'] is not None)/mass,policy_response_pp=100*sum(r['policy'] for r in rr)/mass,gap_supported_share=support/mass)
            for k in ('gap0','gap1','delta_gap','success_gap_delta','wait_credit_gain','try_credit_gain','success_credit_gain'):
                row[k]=sum(r[k]*r['supported_gap_mass'] for r in rr if r[k] is not None)/support if support>0 else None
            compact.append(row)
    csvwrite('first_birth_income_wealth_summary.csv',compact)
    meta=dict(reference='post-interest chain13 without Estate-A; diagnostic, no recalibration',kind=args.kind,pair_changes=changes,price0=float(a0['p_eq'][0]),price1=float(a1['p_eq'][0]),income_cuts=cuts,income_definition='annual gross labor earnings: P.income*z/[period_years*(1-tau_pay)]; current Markov earnings state, fixed baseline childless exposure boundaries',wealth_groups=['b<=0','0<b<=1','b>1'],checks=[check0,check1],logit_checks=gap_rows,sources=[dict(path=str(p.resolve()),arrays_sha256=sha(p/'solution_arrays.npz'),parameters_sha256=sha(p/'executed_P.json'),standard_diagnostics=str(p/'standard_diagnostics')) for p in (args.baseline,args.relaxed)],validity='No aggregate-policy consumption/saving claims. Signed inversion roundoff retained. All gap support denominators reported; no dropped mass, probability clipping or gate changes.')
    (HERE/'receipt.json').write_text(json.dumps(meta,indent=2)+'\n')
    # Supplemental figures retain all original standard graphs via source links.
    first=[r for r in grouped if r['children_ever_born']==0]
    fig,ax=plt.subplots(1,3,figsize=(13,3.7),sharey=True)
    for inc in range(3):
        for wb in range(3):
            y=[]
            for age in ages[fertile]:
                rr=[r for r in first if r['age']==age and r['income_group']==inc and r['wealth_group']==wb]
                mass=sum(r['exposure0'] for r in rr)
                y.append(100*sum(r['policy'] for r in rr)/mass if mass>0 else np.nan)
            ax[inc].plot(ages[fertile],y,'o-',label=meta['wealth_groups'][wb])
        ax[inc].axhline(0,color='gray',lw=.7); ax[inc].set_title(f'Current earnings group {inc+1}'); ax[inc].set_xlabel('Age cell starts'); ax[inc].grid(alpha=.2)
    ax[0].set_ylabel('Common-state first-birth response (pp)'); ax[-1].legend(); fig.tight_layout(); fig.savefig(HERE/'supplemental_first_birth_response.png',dpi=160); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4))
    for n in range(3):
        rr=[r for r in totals if r['children_ever_born']==n]
        ax.bar(n-.18,sum(r['policy'] for r in rr),.36,color='#1f5a99',label='Policy at baseline exposure' if n==0 else None)
        ax.bar(n+.18,sum(r['composition'] for r in rr),.36,color='#d47b24',label='Exposure/distribution' if n==0 else None)
    ax.axhline(0,color='gray',lw=.7); ax.set_xticks(range(3),['First birth','Second birth','Third birth']); ax.set_ylabel('Change in household birth flow'); ax.legend(); fig.tight_layout(); fig.savefig(HERE/'supplemental_birth_decomposition.png',dpi=160); plt.close(fig)
    print(json.dumps({'checks':meta['checks'],'rows':len(grouped),'income_cuts':cuts,'totals':{str(n):{k:sum(r[k] for r in totals if r['children_ever_born']==n) for k in ('exposure0','birth_flow0','birth_flow1','policy','composition','total')} for n in range(3)}},indent=2))

if __name__=='__main__': main()
