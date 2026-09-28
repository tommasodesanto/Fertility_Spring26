"""Occupied borrowing-limit accounting at the frozen reference; no solves.

Run only under Torch Slurm in the original pinned container mount. Every mean
and binding mass uses date-owned renter/buyer/stayer mass, never a conditional
policy weighted as though all owners were buyers. Figures are supplemental.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import time

ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
BASE = ROOT/'output/model/fertility_identification_20260928'
SOURCE = BASE/'resume_v1/selected_export/primary'
LABEL = '2007 stationary reference — block0506, September 28 verified export'
CONTRACT_SHA = '68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf'
CHECKPOINT_SHA = 'b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d'
BINDING_TOL = 1e-8
VIOLATION_TOL = 1e-9
MASS_TOL = 2e-10


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def main(out):
    if sys.platform != 'linux' or not os.environ.get('SLURM_JOB_ID','').isdigit():
        raise RuntimeError('Torch Slurm execution required; never run on the Mac')
    out = out.resolve()
    if out == SOURCE.resolve() or SOURCE.resolve() in out.parents:
        raise ValueError('Reference export is read-only')
    out.mkdir(parents=True,exist_ok=False)
    started = time.monotonic()
    sys.path.insert(0,str(ROOT/'code/model/tools'))
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import run_e5f_fertility_identification as driver
    import e5f_evening_calibration_runtime as runtime
    from e5f_overnight_estate_audit import policy_mass_branches

    def write(name,obj):
        (out/name).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')

    def table(name,rows):
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

    print('Authenticating frozen reference; zero solves',flush=True)
    manifest=read(BASE/'fixed_reference_manifest.json')
    assert manifest['label']==LABEL and manifest['checkpoint']['sha256']==CHECKPOINT_SHA
    cp=BASE/'contract_v1/contract.json'
    assert sha(cp)==CONTRACT_SHA
    c,objectives=driver.verify(cp)
    hashes=read(SOURCE/'artifact_hashes.json')
    for name,digest in hashes.items():assert sha(SOURCE/name)==digest,name
    assert sum(name.endswith('.png') for name in hashes)==17
    checkpoint=SOURCE/'initial_state.pkl.gz'
    assert sha(checkpoint)==CHECKPOINT_SHA
    evaluator=runtime.setup(dict(c,objective=c['lanes']['primary']['objective']),
                            objectives['primary'],out/'runtime_preparation')
    with gzip.open(checkpoint,'rb') as f:packet=pickle.load(f)
    P,e,bg,shared=(packet[k] for k in ('parameters','evaluation','b_grid','shared'))
    bg=np.asarray(bg,dtype=float);model=evaluator.rt['model'];policy=e.policy
    assert int(P.I)==1 and P.native_purchase_income and P.native_due_stayer_credit
    assert not bool(getattr(P,'native_solvency_credit',False))
    assert not bool(getattr(P,'use_pti_constraint',False))
    assert np.all(np.asarray(P.owner_ltv_multipliers)==1)
    assert float(P.psi_child)==manifest['economic_contract']['frozen_psi_child']
    assert not bool(getattr(P,'mortgage_origination_only',False))
    assert float(getattr(P,'mortgage_amortization',0.))==0.
    branches=policy_mass_branches(e,P)
    assert len(branches)==2
    total=float(e.g_current.sum())
    assert abs(total-1.)<=MASS_TOL
    # Authenticate existing purchase and origin-specific stayer feasibility too.
    purchase=evaluator.rt['accounting'].audit_purchase_accounting(e,P,shared,bg,model)
    price=float(policy.price[0]);rows=[];checks=[]
    ages=P.age_start+P.da*np.arange(P.J)
    for j,age in enumerate(ages):
        death_possible=(j==P.J-1 or (P.use_age_survival and P.survival_probs[j]<1.))
        for branch_index,tenure,tenures in ((0,'renter',(0,)),
                (0,'buyer',range(1,1+P.n_house)),(1,'owner_stayer',range(1,1+P.n_house))):
            mass_all,saving_all,_=branches[branch_index]
            for ten in tenures:
                mass=np.asarray(mass_all[:,ten,0,j],dtype=float)
                saving=np.asarray(saving_all[:,ten,0,j],dtype=float)
                household_mass=float(mass.sum())
                grid=np.full(mass.shape,float(bg[0]))
                if ten==0:
                    credit=np.broadcast_to(model.renter_borrowing_floor(P,bg,j)[:,None,None,None],mass.shape)
                    death=np.zeros_like(mass)
                    death_enforced_separately=False
                    formula='min(s[j+1]*min(b,0),-D[j+1])'
                    credit_name='renter_unsecured_floor'
                else:
                    house=float(P.H_own[ten-1])
                    collateral=-np.asarray(shared.phi_choice[0,ten],dtype=float)*price*house
                    collateral=collateral[None,None,:,:]
                    if tenure=='buyer':
                        credit=np.broadcast_to(model.owner_borrowing_floor(P,bg[:,None,None,None],collateral,j),mass.shape)
                        formula='-phi[n,m]*p*H (native_purchase_income; next-age multiplier=1)'
                        credit_name='buyer_collateral_floor'
                        death_enforced_separately=False
                    else:
                        credit=np.broadcast_to(model.native_due_owner_floor(bg[:,None,None,None],collateral),mass.shape)
                        formula='min(b,-phi[n,m]*p*H)'
                        credit_name='owner_stayer_principal_floor'
                        death_enforced_separately=death_possible
                    death=np.full(mass.shape,model.native_due_death_floor(P,j,price,house))
                # Only DUE owner-stayer kernel has this explicit death bound.
                # For other branches death solvency is an accounting diagnostic,
                # not an invented extra Bellman constraint.
                effective=np.maximum(credit,grid)
                if death_enforced_separately:effective=np.maximum(effective,death)
                effective_violation=float(mass[saving<effective-VIOLATION_TOL].sum())
                checks.append(dict(age_left=float(age),branch=tenure,owner_rooms=0. if ten==0 else float(P.H_own[ten-1]),
                                   effective_floor_violation_mass=effective_violation))
                if effective_violation>MASS_TOL:raise ValueError('Occupied native saving floor violation: '+str(checks[-1]))
                components=[(credit_name,credit,True),( 'grid_lower_bound',grid,True)]
                if death_possible:components.append(('death_solvency',death,death_enforced_separately))
                for name,floor,enforced in components:
                    slack=saving-floor
                    binding=np.abs(slack)<=BINDING_TOL
                    active=(np.abs(floor-effective)<=BINDING_TOL) if enforced else np.zeros(mass.shape,dtype=bool)
                    violation=float(mass[slack < -VIOLATION_TOL].sum())
                    if violation>MASS_TOL:raise ValueError('Occupied constraint/accounting violation: '+str((age,tenure,name,violation)))
                    rows.append(dict(age_left=float(age),branch=tenure,owner_rooms=0. if ten==0 else float(P.H_own[ten-1]),
                        constraint=name,enforced_as_separate_native_bound=bool(enforced),household_mass=household_mass,
                        at_component_boundary_mass=float(mass[binding].sum()),
                        binding_active_component_mass=float(mass[binding & active].sum()),
                        component_is_effective_floor_mass=float(mass[active].sum()),
                        violation_mass=violation,
                        mean_slack=float(np.sum(mass*slack)/household_mass) if household_mass else None,
                        smallest_occupied_slack=float(slack[mass>1e-12].min()) if np.any(mass>1e-12) else None,
                        fraction_at_boundary=float(mass[binding].sum()/household_mass) if household_mass else None,
                        fraction_binding_active=float(mass[binding & active].sum()/household_mass) if household_mass else None,
                        credit_formula=formula,binding_tolerance=BINDING_TOL))
    table('supplemental_constraint_components.csv',rows)
    # Collapse owner product rows using masses, not an unweighted mean of shares.
    grouped=[]
    for age in ages:
        for tenure in ('renter','buyer','owner_stayer'):
            for kind in ('credit','grid','death'):
                selected=[r for r in rows if r['age_left']==age and r['branch']==tenure and
                    ((kind=='credit' and r['constraint'].endswith('_floor')) or
                     (kind=='grid' and r['constraint']=='grid_lower_bound') or
                     (kind=='death' and r['constraint']=='death_solvency'))]
                if not selected:continue
                mass=sum(r['household_mass'] for r in selected)
                bound=sum(r['at_component_boundary_mass'] for r in selected)
                active=sum(r['binding_active_component_mass'] for r in selected)
                grouped.append(dict(age_left=float(age),branch=tenure,component=kind,
                    household_mass=mass,at_boundary_mass=bound,binding_active_mass=active,
                    fraction_at_boundary=bound/mass if mass else None,
                    fraction_binding_active=active/mass if mass else None))
    table('supplemental_constraints_by_age_tenure.csv',grouped)
    fig,axes=plt.subplots(1,3,figsize=(13,4),sharey=True)
    for ax,tenure in zip(axes,('renter','buyer','owner_stayer')):
        for kind,label in (('credit','Artificial credit limit'),('grid','Grid floor'),('death','Death-solvency boundary')):
            sub=[r for r in grouped if r['branch']==tenure and r['component']==kind]
            field='fraction_at_boundary' if kind=='death' else 'fraction_binding_active'
            ax.plot([r['age_left'] for r in sub],[r[field] for r in sub],label=label)
        ax.set_title(tenure.replace('_',' '));ax.set_xlabel('Age at left of model cell');ax.grid(alpha=.2)
    axes[0].set_ylabel('Share of occupied branch mass');axes[0].legend(fontsize=7)
    fig.suptitle('Supplemental: '+LABEL+'\nBoundary incidence; death boundary is not a separate native constraint for renters/buyers',fontsize=9)
    fig.tight_layout(rect=(0,0,1,.91));fig.savefig(out/'supplemental_constraint_binding.png',dpi=150);plt.close(fig)
    for name,digest in hashes.items():assert sha(SOURCE/name)==digest,name
    write('supplemental_constraints_receipt.json',dict(label=LABEL,status='PASS',model_solves=0,
        source_changes=[],economic_changes=[],elapsed_seconds=time.monotonic()-started,
        slurm_job=os.environ['SLURM_JOB_ID'],script_sha256=sha(__file__),checkpoint_sha256=CHECKPOINT_SHA,
        contract_sha256=CONTRACT_SHA,source_manifest=c['source_manifest'],standard_plots_unchanged=17,
        reference_manifest=str(BASE/'fixed_reference_manifest.json'),common_primary_loss=manifest['common_primary_loss'],
        full_target_table=str(SOURCE/'target_fit.csv'),full_parameter_table=str(SOURCE/'parameters.csv'),
        thresholds=dict(binding_asset_units=BINDING_TOL,violation_asset_units=VIOLATION_TOL,violation_mass=MASS_TOL),
        original_purchase_accounting_replay=purchase,effective_floor_checks=checks,
        definitions=dict(credit='Native component before grid/death max; component equality may overlap other constraints.',
            renter='Date-owned current renter mass and its saving policy; unsecured rollover/taper limit uses next-age indices.',
            buyer='Date-owned current owner mass excluding owner stayers, including purchases by previous owners of another product.',
            stayer='Date-owned owner-stayer mass with bp_pol_stay, not the buyer saving policy.',
            grid='Wealth-grid minimum shown separately; binding means saving at its node and this component attains the effective floor.',
            death='b\u2032 >= -(1-selling_cost)*pH for owners, b\u2032 >= 0 for renters, only where death possible. DUE owner-stayer bound is explicit; others are audited estate solvency.',
            occupancy='All positive mass integrated; smallest slack reported only where cell mass exceeds1e-12.',
            unit='Asset stock in model annual-income units; no annualization of slack or constraint probability.'),
        source_equations=dict(renter='solver.renter_borrowing_floor -> debt_rule_at_age -> unsecured_debt_floor',
            buyer='solver.owner_borrowing_floor with native_purchase_income=True',
            stayer='solver.native_due_owner_floor, then max with native_due_death_floor and grid minimum',
            branch_mass='e5f_overnight_estate_audit.policy_mass_branches'),
        limitations=['No new solution or counterfactual: equality to a bound does not measure a shadow value or causal response.',
            'Components overlap; do not sum binding shares across components.',
            'Purchase entry affordability incidence and excluded alternatives are not implemented. Original exact purchase-feasibility audit is replayed only, including income/R and transaction grid support.',
            'A reported zero purchase-violation mass is not evidence that no household is excluded by affordability.',
            'Only the existing discrete-grid solution is described; off-grid optima and numerical-grid sensitivity remain untested.']))
    print(json.dumps(dict(status='PASS',output=str(out),rows=len(rows),model_solves=0)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
