"""Read-only accounting of frozen-reference renter debt taper incidence.

Torch Slurm only. Authenticates the block0506/160-node checkpoint, loads its
saved evaluation, and runs no Bellman, KFE, equilibrium, or transition solves.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, os, pickle, sys, time
from pathlib import Path

ROOT=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
BASE=ROOT/'output/model/fertility_identification_20260928'
SOURCE=BASE/'resume_v1/selected_export/primary'
LABEL='2007 stationary reference — block0506, September 28 verified export'
CHECKPOINT_SHA='b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d'
CONTRACT_SHA='68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf'
BIND_TOL=1e-8
MASS_TOL=2e-10
ENTRY_SHA='a24255f77656453c7a57d2cb18188f03d6696129dc2d13e97afbbb55c8cbe237'
AUDIT_ROOT=ROOT/'output/model/fixed_reference_economics_20260928/renter_taper_audit_v1'

def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''): h.update(block)
    return h.hexdigest()

def main(out):
    if sys.platform!='linux' or not os.environ.get('SLURM_JOB_ID','').isdigit():
        raise RuntimeError('Torch Slurm execution required')
    out=Path(out).resolve()
    if AUDIT_ROOT.resolve() not in out.parents or out==AUDIT_ROOT.resolve():
        raise ValueError('Output must be a new child of renter_taper_audit_v1')
    out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    sys.path.insert(0,str(ROOT/'code/model/tools'))
    import numpy as np
    import run_e5f_fertility_identification as driver
    import e5f_evening_calibration_runtime as runtime
    from e5f_overnight_estate_audit import policy_mass_branches
    manifest=json.loads((BASE/'fixed_reference_manifest.json').read_text())
    assert manifest['label']==LABEL and manifest['checkpoint']['sha256']==CHECKPOINT_SHA
    contract=BASE/'contract_v1/contract.json'
    assert sha(contract)==CONTRACT_SHA
    c,objectives=driver.verify(contract)
    hashes=json.loads((SOURCE/'artifact_hashes.json').read_text())
    for name,digest in hashes.items(): assert sha(SOURCE/name)==digest,name
    assert sum(k.endswith('.png') for k in hashes)==17
    checkpoint=SOURCE/'initial_state.pkl.gz'
    assert sha(checkpoint)==CHECKPOINT_SHA
    # Runtime setup authenticates portable observer/model sources, but does not solve.
    evaluator=runtime.setup(dict(c,objective=c['lanes']['primary']['objective']),
                            objectives['primary'],out/'runtime_preparation')
    with gzip.open(checkpoint,'rb') as f: packet=pickle.load(f)
    P,e,bg,shared=(packet[k] for k in ('parameters','evaluation','b_grid','shared'))
    bg=np.asarray(bg,dtype=float); model=evaluator.rt['model']; policy=e.policy
    if int(P.I)!=1 or len(bg)!=160:
        raise ValueError('Audit requires the authenticated one-market 160-node checkpoint')
    if not np.allclose(np.asarray(policy.loc_probs,dtype=float),1.0,rtol=0,atol=1e-12):
        raise ValueError('One-market location stay probability is not identically one')
    if float(getattr(P,'lambda_d',float('nan'))) != 0.0:
        raise ValueError('Frozen reference lambda_d is not zero')
    branches=policy_mass_branches(e,P)
    if len(branches)!=2: raise ValueError('Expected native owner-stayer branch split')
    ages=P.age_start+P.da*np.arange(P.J)
    rows=[]
    sale_rows=[]
    price=float(e.policy.price[0])
    for j,age in enumerate(ages):
        mass_all,saving_all,_=branches[0]
        mass=np.asarray(mass_all[:,0,0,j],dtype=float)
        saving=np.asarray(saving_all[:,0,0,j],dtype=float)
        b=np.broadcast_to(bg[:,None,None,None],mass.shape)
        neg=(b<0)&(mass>0)
        floor=np.broadcast_to(model.renter_borrowing_floor(P,bg,j)[:,None,None,None],mass.shape)
        expected=np.minimum(np.asarray(P.debt_taper_weights)[j+1]*np.minimum(bg,0.),-np.asarray(P.debt_caps)[j+1])
        if not np.allclose(floor[:,0,0,0],expected,rtol=0,atol=1e-12):
            raise ValueError('Native renter floor differs from authenticated next-age schedule')
        delta=np.maximum(floor-np.minimum(b,0.),0.)
        # Binding incidence is conditioned on negative b: equality with the
        # native floor only. Zero-saving/grid boundaries are not counted.
        bind=neg&(np.abs(saving-floor)<=BIND_TOL)
        nm=float(mass[neg].sum()); allm=float(mass.sum())
        rows.append(dict(age_left=float(age),occupied_renter_mass=allm,
            negative_b_renter_mass=nm,negative_b_share=nm/allm if allm else None,
            mean_current_debt_given_negative_b=float(np.sum(mass[neg]*(-b[neg]))/nm) if nm else None,
            mean_mandatory_repayment_given_negative_b=float(np.sum(mass[neg]*delta[neg])/nm) if nm else None,
            negative_b_mass_at_native_taper_floor=float(mass[bind].sum()),
            negative_b_floor_binding_share=float(mass[bind].sum()/nm) if nm else None,
            taper_weight_next_age=float(np.asarray(P.debt_taper_weights)[j+1]),
            cap_next_age=float(np.asarray(P.debt_caps)[j+1])))
        # Exact owner-origin to renter tracing for this one-market case.
        # Match native realization: cast menu probabilities to float64,
        # normalize per state, then apply the saved transaction interpolation.
        gpost=np.asarray(e.g_post_fertility,dtype=float)
        raw_shortfall=short_mass=incoming_neg=owner_renter_mass=0.0
        incoming_nodes=np.zeros((len(bg),P.Nz,P.n_parity,P.n_child_states),dtype=float)
        reconstructed=np.zeros_like(incoming_nodes)
        for old in range(0,1+P.n_house):
            sale=0.0 if old==0 else (1-float(P.psi))*price*float(P.H_own[old-1])
            for zz in range(P.Nz):
                for n in range(P.n_parity):
                    for m in range(P.n_child_states):
                        src=gpost[:,old,0,j,zz,n,m]
                        probs=np.asarray(policy.tenure_probs[:,old,0,j,zz,n,m,:],dtype=float)
                        den=probs.sum(axis=-1)
                        prob=np.divide(probs[:,0],den,out=np.zeros_like(den),where=den>0)
                        weighted=src*prob
                        pos=bg+sale
                        if old>0:
                            owner_renter_mass+=float(weighted.sum())
                            short=pos < -1e-9
                            raw_shortfall+=float(weighted[short].sum())
                            short_mass+=float(np.sum(weighted[short]*(-pos[short])))
                        # `maps` exposes named arrays; axis order is [location, old tenure, new tenure, n, m, b].
                        idx=np.asarray(policy.maps.tmx_idx[0,old,0,n,m,:],dtype=int)
                        wt=np.asarray(policy.maps.tmx_wt[0,old,0,n,m,:],dtype=float)
                        for ib in range(len(bg)):
                            w=weighted[ib]
                            if w==0: continue
                            reconstructed[idx[ib],zz,n,m]+=w*(1-wt[ib])
                            reconstructed[idx[ib]+1,zz,n,m]+=w*wt[ib]
                            if old>0:
                                incoming_nodes[idx[ib],zz,n,m]+=w*(1-wt[ib])
                                incoming_nodes[idx[ib]+1,zz,n,m]+=w*wt[ib]
        incoming_neg=float(incoming_nodes[bg<0].sum())
        actual=e.g_current[:,0,0,j]
        if abs(float(incoming_nodes.sum())-owner_renter_mass)>MASS_TOL:
            raise ValueError('Owner-to-renter interpolation did not conserve mass')
        if np.max(np.abs(reconstructed-actual))>MASS_TOL:
            raise ValueError('Native renter-origin/owner-origin replay does not match saved current renter mass')
        sale_rows.append(dict(age_left=float(age),owner_to_renter_mass=owner_renter_mass,
            strict_sale_shortfall_mass=raw_shortfall,
            strict_sale_shortfall_share=raw_shortfall/owner_renter_mass if owner_renter_mass else None,
            mean_shortfall_conditional=float(short_mass/raw_shortfall) if raw_shortfall else 0.0,
            negative_post_sale_renter_mass_after_native_grid_interpolation=incoming_neg,
            full_current_renter_mass=float(actual.sum())))
    with (out/'owner_to_renter_sale_accounting.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(sale_rows[0])); w.writeheader(); w.writerows(sale_rows)
    with (out/'renter_taper_by_age.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    entry=None
    eg=getattr(P,'fixed_reference_entry_grid',None)
    ec=getattr(P,'fixed_reference_entry_conditional',None)
    if eg is not None and ec is not None:
        eg=np.asarray(eg,dtype=float); ec=np.asarray(ec,dtype=float)
        if ec.ndim!=2 or ec.shape[0]!=len(eg): raise ValueError('Unexpected conditional entry matrix shape')
        if ec.shape[1]!=P.Nz or sha_bytes(ec.tobytes(order='C'))!=ENTRY_SHA:
            raise ValueError('Authenticated fixed-reference entry conditional array differs')
        # Use the same native income-state weights as the earnings-entry audit.
        _,income_weights,_=model.income_transition_values(P)
        zweights=np.asarray(income_weights,dtype=float).reshape(-1)
        if len(zweights)!=ec.shape[1] or zweights.sum()<=0: raise ValueError('Entry income weights mismatch')
        zweights=zweights/zweights.sum()
        entry_weights=ec@zweights
        if abs(float(entry_weights.sum())-1)>MASS_TOL: raise ValueError('Entry weights do not sum to one')
        entry=dict(grid_nodes=int(len(eg)),conditional_shape=list(ec.shape),
            conditional_sha256=sha_bytes(ec.tobytes(order='C')),income_state_weights=zweights.tolist(),
            negative_entry_mass=float(entry_weights[eg<0].sum()),
            negative_entry_share=float(entry_weights[eg<0].sum()),
            mean_entry_net_financial_wealth=float(entry_weights@eg))
    else:
        raise ValueError('Authenticated fixed-reference entry conditional and grid are required')
    # All gates must pass before any PASS receipt is emitted.
    for name,digest in hashes.items(): assert sha(SOURCE/name)==digest,name
    if abs(float(e.g_current.sum())-1.)>MASS_TOL: raise ValueError('Saved current population mass is not one')
    negative_mass=float(sum(r['negative_b_renter_mass'] for r in rows))
    all_renters=float(sum(r['occupied_renter_mass'] for r in rows))
    exposed=[r for r in rows if r['taper_weight_next_age']<1.0]
    exposed_negative=float(sum(r['negative_b_renter_mass'] for r in exposed))
    exposed_binding=float(sum(r['negative_b_mass_at_native_taper_floor'] for r in exposed))
    all_mass=float(np.asarray(e.g_current).sum())
    incidence=dict(total_occupied_renter_mass=all_renters,negative_current_renter_mass=negative_mass,
        negative_renter_share_of_renters=negative_mass/all_renters if all_renters else None,
        negative_renter_mass_share_of_all_households=negative_mass/all_mass if all_mass else None,
        taper_exposed_negative_mass=exposed_negative,
        taper_exposed_negative_share_of_renters=exposed_negative/all_renters if all_renters else None,
        taper_exposed_negative_share_of_all_households=exposed_negative/all_mass if all_mass else None,
        taper_exposed_floor_binding_mass=exposed_binding,
        taper_exposed_floor_binding_share_of_renters=exposed_binding/all_renters if all_renters else None,
        taper_exposed_floor_binding_share_of_all_households=exposed_binding/all_mass if all_mass else None)
    receipt=dict(status='PASS',reference_label=LABEL,model_solves=0,
        source_changes=[],economic_changes=[],checkpoint_sha256=CHECKPOINT_SHA,
        contract_sha256=CONTRACT_SHA,slurm_job=os.environ['SLURM_JOB_ID'],
        script_sha256=sha(__file__),elapsed_seconds=time.monotonic()-started,
        standard_plots_unchanged=17,entry_distribution=entry,incidence=incidence,
        interpretation=dict(binding='Native renter floor equality only among b<0; not generic zero/grid binding.',
            incidence='Accounting of saved equilibrium masses and saved policies; not a causal effect or predicted reform response.',
            sale_origin='Owner-to-renter origin is traced from saved post-fertility mass, native float64-normalized tenure probabilities, and saved native transaction interpolation; raw sale-shortfall and mapped negative renter mass are reported separately.'))
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps(dict(status='PASS',output=str(out),ages=len(rows),entry=entry is not None,
                          elapsed_seconds=time.monotonic()-started),sort_keys=True),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True); main(ap.parse_args().output)
