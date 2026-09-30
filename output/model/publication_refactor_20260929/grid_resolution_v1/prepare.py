"""Zero-solve grid preparation. Proposed numerical projection; never runs a model."""
from __future__ import annotations
import copy
import csv
import json
import math
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'code/model'))
from refactor_lab.inputs import load_inputs, encode, decode, sha256_file
from refactor_lab.engine.utils import make_grid

BUNDLE_SHA = '427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7'
ORIGINAL = ROOT / 'output/model/publication_refactor_20260929/local_export_v1/inputs'

def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')

def rouwenhorst(n, rho, sigma_log):
    p = (1 + rho) / 2
    transition = np.array([[p, 1-p], [1-p, p]])
    for size in range(3, n+1):
        old = transition
        transition = np.zeros((size,size))
        transition[:-1,:-1] += p*old
        transition[:-1,1:] += (1-p)*old
        transition[1:,:-1] += (1-p)*old
        transition[1:,1:] += p*old
        transition[1:-1] *= .5
    transition /= transition.sum(axis=1, keepdims=True)
    weights = np.array([math.comb(n-1,k) for k in range(n)],float) / 2**(n-1)
    levels = np.exp(np.linspace(-sigma_log*math.sqrt(n-1),sigma_log*math.sqrt(n-1),n))
    levels /= weights @ levels
    return levels, weights, transition

def quantile_transport(old_weights, new_weights):
    """CDF-bin overlap, old x new joint mass; no node clipping or mass deletion."""
    for weights in (old_weights,new_weights):
        assert np.isfinite(weights).all() and (weights>0).all()
        np.testing.assert_allclose(weights.sum(),1,rtol=0,atol=3e-16)
    left = np.r_[0.,np.cumsum(old_weights)]
    right = np.r_[0.,np.cumsum(new_weights)]
    result = np.maximum(0.,np.minimum(left[1:,None],right[None,1:])-
                        np.maximum(left[:-1,None],right[None,:-1]))
    np.testing.assert_allclose(result.sum(axis=1),old_weights,rtol=0,atol=3e-16)
    np.testing.assert_allclose(result.sum(axis=0),new_weights,rtol=0,atol=3e-16)
    return result

def wealth_subgrid(grid, conditional, count):
    """Retain all entry atoms, endpoints and zero; evenly retain other indices."""
    required = set(np.flatnonzero(conditional.sum(axis=1)>0).tolist()) | {0,len(grid)-1}
    required.update(np.flatnonzero(grid==0).tolist())
    if len(required)>count: raise ValueError('Too many mandatory wealth atoms')
    # Farthest index insertion covers remaining gaps deterministically. This is
    # a numerical grid choice, never an entry projection.
    selected = set(required)
    while len(selected)<count:
        unused = [i for i in range(len(grid)) if i not in selected]
        selected.add(max(unused,key=lambda i:(min(abs(i-j) for j in selected),-i)))
    return np.array(sorted(selected),dtype=int)

def moments(grid,z,joint):
    wb = joint.sum(axis=1); wz=joint.sum(axis=0)
    eb=float(wb@grid); ez=float(wz@z)
    return dict(total_mass=float(joint.sum()),mean_wealth=eb,mean_income_multiplier=ez,
        negative_wealth_mass=float(wb[grid<0].sum()),
        wealth_second_moment=float(wb@grid**2),income_second_moment=float(wz@z**2),
        wealth_income_covariance=float(np.sum(joint*grid[:,None]*z[None,:])-eb*ez))

def main():
    loaded=load_inputs(ORIGINAL,ROOT,BUNDLE_SHA)
    P=loaded.parameters; grid=loaded.b_grid
    assert len(grid)==P.Nb==160 and len(P.z_grid)==P.Nz==15
    assert P.native_fixed_reference_entry and P.native_explicit_transaction_grid
    assert not P.permanent_income_levels_enabled
    np.testing.assert_array_equal(make_grid(P),grid)
    C=P.fixed_reference_entry_conditional
    np.testing.assert_allclose(C.sum(axis=0),1,rtol=0,atol=2e-12)
    rho=float(P.income_shock_persistence)
    sigma_log=float(np.ptp(np.log(P.z_grid))/(2*math.sqrt(14)))
    z15,w15,t15=rouwenhorst(15,rho,sigma_log)
    np.testing.assert_allclose(z15,P.z_grid,rtol=3e-15,atol=0)
    np.testing.assert_array_equal(w15,P.z_weights)
    np.testing.assert_allclose(t15,P.Pi_z,rtol=3e-15,atol=3e-16)
    z9,w9,t9=rouwenhorst(9,rho,sigma_log)
    np.testing.assert_allclose(t9.sum(axis=1),1,rtol=0,atol=3e-16)
    np.testing.assert_allclose(w9@t9,w9,rtol=0,atol=3e-16)
    assert np.isfinite(t9).all() and (t9>=0).all()
    transport=quantile_transport(w15,w9)
    selected=wealth_subgrid(grid,C,120)
    new_grid=grid[selected]
    old_joint=C*w15[None,:]
    # Tiny probability contraction uses NumPy's ordered non-BLAS path. Some
    # Apple BLAS versions emit invalid floating flags for this finite matmul.
    # Independently compare broadcast summation instead of suppressing warnings.
    new_joint=np.einsum('bi,ij->bj',C[selected],transport,optimize=False)
    direct=np.sum(C[selected,:,None]*transport[None,:,:],axis=1)
    np.testing.assert_allclose(new_joint,direct,rtol=0,atol=3e-16)
    assert np.isfinite(new_joint).all() and (new_joint>=0).all()
    # Exact preservation of all occupied wealth atoms, to floating arithmetic.
    restored=np.zeros(len(grid));restored[selected]=new_joint.sum(axis=1)
    np.testing.assert_allclose(restored,old_joint.sum(axis=1),rtol=0,atol=3e-16)
    newC=new_joint/w9[None,:]
    np.testing.assert_allclose(newC.sum(axis=0),1,rtol=0,atol=2e-12)
    new=copy.deepcopy(P)
    new.Nb=120;new.Nz=9;new.z_grid=z9;new.z_weights=w9;new.Pi_z=t9
    new.earnings_transaction_grid=new_grid.copy()
    new.fixed_reference_entry_grid=new_grid.copy()
    new.fixed_reference_entry_conditional=newC
    # Dormant historical 3x5 metadata is not a meaningful 9-state index. Keep
    # it out of proposed input; feature remains disabled, with explicit receipt.
    del new.permanent_income_group_index
    del new.permanent_income_base_state_index
    np.testing.assert_array_equal(make_grid(new),new_grid)
    assert new.entry_wealth_censor_to_frontier==P.entry_wealth_censor_to_frontier
    # Fixed conditional-entry branch executes before the dormant censor path.
    from refactor_lab.engine.distribution import entry_wealth_grid_weights
    for k,z in enumerate(z9):
        indices,weights=entry_wealth_grid_weights(new_grid,new,z_value=float(z))
        assert (indices>=0).all() and (indices<120).all()
        np.testing.assert_array_equal(weights,newC[indices,k])
    arrays={'b_grid':new_grid,'reference_price':loaded.reference_price.copy()}
    encoded={k:encode(v,k,arrays) for k,v in sorted(vars(new).items())}
    dest=HERE/'proposed_120x9';dest.mkdir(exist_ok=True)
    np.savez(dest/'arrays.npz',**arrays)
    write(dest/'bundle.json',dict(schema='grid_resolution_proposal_v1',status='not_adopted_not_runtime_authorized',
        parameters=encoded,arrays_sha256=sha256_file(dest/'arrays.npz'),
        parent_bundle_sha256=BUNDLE_SHA,parent_identity=loaded.identity,
        changes=['160 to 120 wealth nodes preserving every occupied entry atom',
                 '15 to 9 Rouwenhorst states, same period rho and stationary log variance',
                 'CDF-overlap transport of fixed conditional wealth law, preserving wealth marginal',
                 'remove two inactive historical 3x5 index arrays'],
        credit_choices={'reference_original':None,'corrected_diagnostic':0.14}))
    reference=ROOT/'output/model/fertility_identification_20260928/resume_v1/selected_export/primary'
    for name in ('target_fit.csv','parameters.csv'):
        (HERE/('reference_'+name)).write_bytes((reference/name).read_bytes())
    with (HERE/'reference_target_fit.csv').open() as stream:
        fit_rows=list(csv.DictReader(stream))
    with (HERE/'reference_parameters.csv').open() as stream:
        parameter_rows=list(csv.DictReader(stream))
    assert len(fit_rows)==14 and len(parameter_rows)==31
    assert sum(float(r['weight'] or 0)>0 for r in fit_rows)==10
    proposal_meta=json.loads((dest/'bundle.json').read_text())
    with np.load(dest/'arrays.npz',allow_pickle=False) as archive:
        restored_fields={k:decode(v,archive) for k,v in proposal_meta['parameters'].items()}
        assert restored_fields['Nb']==len(archive['b_grid'])==120
        assert restored_fields['Nz']==len(restored_fields['z_grid'])==9
        assert restored_fields['fixed_reference_entry_conditional'].shape==(120,9)
    outside=(P.z_grid<z9[0]) | (P.z_grid>z9[-1])
    before=moments(grid,P.z_grid,old_joint);after=moments(new_grid,z9,new_joint)
    period=float(P.period_years); annual_rho=rho**(1/period)
    annual_sigma=sigma_log*math.sqrt(1-annual_rho**2)
    source_paths=[HERE/'prepare.py',ROOT/'code/model/refactor_lab/inputs.py',ROOT/'code/model/refactor_lab/engine/utils.py',
        ROOT/'code/model/refactor_lab/engine/distribution.py',
        ROOT/'output/model/fertility_identification_20260928/contract_v1/primary_objective.json',
        ROOT/'output/model/fertility_identification_20260928/two_stream_overnight_v1/config.json',
        ROOT/'output/model/fertility_identification_20260928/two_stream_overnight_v1/search.py',
        ROOT/'code/model/tools/e5f_evening_calibration_runtime.py']
    source_paths+=list((ROOT/'code/model/refactor_lab/engine').glob('*.py'))
    paired=ROOT/'output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed'
    source_paths += [paired/n for n in ('driver.py','single_price.py','phase_a.py','phase_b_ge.py')]
    config=json.loads((ROOT/'output/model/fertility_identification_20260928/two_stream_overnight_v1/config.json').read_text())
    objective=json.loads((ROOT/'output/model/fertility_identification_20260928/contract_v1/primary_objective.json').read_text())
    assert len(objective['target_rows'])==14
    assert sum(float(r['actual_weight'] or 0)>0 for r in objective['target_rows'])==10
    assert {r['parameter'] for r in objective['parameter_restrictions']}==set(config['bounds'])
    write(HERE/'preflight.json',dict(status='zero_solve_proposed_inputs_verified_loop_launch_blocked',
        lifecycle_solves=0,reference_identity=loaded.identity,
        reference_bundle_preserved=sha256_file(ORIGINAL/'bundle.json')==BUNDLE_SHA,
        effective_dimensions={'control':[160,15],'proposal':[len(new_grid),len(z9)]},
        mandatory_entry_wealth_atoms=int(np.count_nonzero(C.sum(axis=1))),
        wealth_atom_max_mass_error=float(np.max(np.abs(restored-old_joint.sum(axis=1)))),
        probability_contraction=dict(method='np.einsum optimize=False',
            independent_check='ordered broadcast sum',
            max_agreement_error=float(np.max(np.abs(new_joint-direct))),
            finite_and_nonnegative=bool(np.isfinite(new_joint).all() and (new_joint>=0).all()),
            cdf_weights_finite_positive_unit_mass=True),
        income_process=dict(period_years=period,period_rho=rho,stationary_log_sd=sigma_log,
            inferred_annual_rho=annual_rho,inferred_annual_innovation_sd=annual_sigma,
            reference_reconstruction_max_z_error=float(np.max(abs(z15-P.z_grid))),
            transition_row_error=float(np.max(abs(t9.sum(axis=1)-1))),
            invariant_error=float(np.max(abs(w9@t9-w9)))),
        entry_moments={'reference':before,'proposal':after,'difference':{k:after[k]-before[k] for k in before}},
        linear_income_projection_failure=dict(outside_support_mass=float(w15[outside].sum()),
            outside_nodes=np.flatnonzero(outside).tolist(),new_support=z9[[0,-1]].tolist(),
            old_support=P.z_grid[[0,-1]].tolist(),
            explanation='Nonnegative linear interpolation cannot retain these old income points inside new support; no clipping was applied.'),
        proposed_bundle_sha256=sha256_file(dest/'bundle.json'),
        source_pins={str(p.relative_to(ROOT)):sha256_file(p) for p in sorted(set(source_paths))},
        calibration_contract=dict(target_weight_fingerprint=config['lanes']['one_birth']['target_fingerprint'],
            free_parameters=config['bounds'],scored_moments=10,normalization_parameter='psi_child',normalization_target=2.1,
            displayed_rows=14,parameter_rows=31,standard_diagnostics=17,
            integration='No existing evaluator directly accepts proposed typed inputs and refactor_lab; isolated factory adapter and full scientific replay are required.'),
        loop_preflight=dict(executable=False,reason='Proposal not lead-reviewed; grid-specific frozen observer adapter and pinned reporting integration absent.',
            proposed_cases=['160x15 full renewal-price/population GE','120x9 full renewal-price/population GE','exact repeat of120x9'],
            held_fixed='Same explicitly selected credit contract in every arm, all reference preferences and primitives, no psi normalization during resolution comparison',
            credit_selection_required=True,maximum_lifecycle_solves=18,total_seconds=2400,case_seconds=300,
            single_cpu=True,checkpoint_after_every_price_trial=True,
            stop=['entry feasibility or estate failure','source/target fingerprint drift','missing14/31/17 outputs','price/population closure failure','time or solve cap']),
        disclosure='CDF transport changes the discretized joint wealth-income association; exact wealth marginal and total mass retained. Numerical approximation only, pending review.'))
    print(json.dumps({'status':'prepared_zero_solves','preflight':str(HERE/'preflight.json'),'entry_moments':{'before':before,'after':after},'outside_support_mass':float(w15[outside].sum())},indent=2))

if __name__=='__main__':main()
