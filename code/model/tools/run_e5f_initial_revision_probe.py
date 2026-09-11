"""Pinned initial-stationary utility/PAYGO probes; no early SMM promotion.

The exact candidate loop, optional 2.1 normalization, independent fiscal and
household checks, checkpoint reload and standard graph packet are exercised.
Legacy stationary moments are saved only as diagnostics: their late-period
target contract is not the proposed early calibration objective.
"""
from __future__ import annotations
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
import copy
import gzip
import json
from pathlib import Path
import pickle
import threading
import time
from types import SimpleNamespace
import numpy as np
import run_e5f_matched_pf_smoke as primitive
import run_e5f_transition_calibration as calibration
from e5f_stationary_paygo import (bind_initial_balanced_pension, certify_initial_pension,
                                rebase_initial_supply, solve_balanced_initial_equilibrium)
from e5f_social_security import fiscal_accounts

ACTIVE_RUN = {}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract',type=Path,required=True)
    parser.add_argument('--contract-sha256',required=True)
    parser.add_argument('--case',choices=('old_old','old_balanced','new_old','new_balanced'),required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    primitive.verify(args.contract,args.contract_sha256)
    c=json.loads(args.contract.read_text())
    if (c.get('schema')!='e5f_parenthood_initial_probe_v1' or c.get('calibrated_smm') is not False
            or c.get('payroll_tax')!=.179 or c.get('housing_supply_elasticity')!=.63
            or c.get('fertility_normalization')!=2.1
            or type(c['normalize']) is not bool or c['repetitions'] not in (1,2)
            or not 1<=c['maximum_stationary_solves_per_repetition']<=23
            or not 1<=c['seconds']<=7200):
        raise ValueError('Invalid explicitly diagnostic initial-probe contract')
    root=Path(__file__).resolve().parents[3]
    required={str(p.relative_to(root)) for p in (root/'code/model').rglob('*.py')}
    if not required or not required.issubset(c['source_sha256']):
        raise ValueError('Initial probe source manifest omits active model Python files')
    for name,pin in c['source_sha256'].items():
        source=(root/name).resolve()
        if not source.is_relative_to(root): raise ValueError('Source path escapes snapshot')
        primitive.verify(source,pin)
    primitive.verify(c['normalized_checkpoint'],c['normalized_checkpoint_sha256'])
    args.output.mkdir(parents=True,exist_ok=False)
    out=args.output
    started=time.monotonic()
    state=dict(phase='loading',case=args.case,completed_repetitions=0,stationary_solves=0)
    finished=threading.Event()
    save_lock=threading.Lock()
    def save(name,value):
        with save_lock:
            primitive.pf.calendar.write_json_atomic(out/name,primitive.pf.calendar.jsonable(value))
    ACTIVE_RUN.update(finished=finished,save=save,state=state,started=started)
    def heartbeat():
        while not finished.wait(min(60.,max(.1,c['seconds']-(time.monotonic()-started)))):
            save('heartbeat.json',dict(state,elapsed_seconds=time.monotonic()-started))
            if time.monotonic()-started>c['seconds']:
                save('timeout.json',dict(state,reason='explicit candidate-loop walltime exhausted'))
                os._exit(124)
    threading.Thread(target=heartbeat,daemon=True).start()
    save('contract.json',dict(c,contract_sha256=args.contract_sha256,case=args.case))
    save('heartbeat.json',state)
    with gzip.open(c['normalized_checkpoint'],'rb') as stream: inherited=pickle.load(stream)
    old=inherited['old']
    if (bool(old.parameters.joint_nested_choice) or float(old.parameters.tau_pay)!=.179
            or not bool(getattr(old.parameters,'exhaustive_saving_control',False))):
        raise ValueError('Expected pinned sequential fixed-payroll-tax seed')
    chain,model=primitive.pf.transition.configure_sequential_model()
    primitive.pf.calendar.apply_fertility=primitive.pf.transition.apply_sequential_fertility
    primitive.pf.calendar.advance_calendar_distribution=primitive.pf.transition.advance_sequential_calendar_distribution
    base,supply_rebase=rebase_initial_supply(old.parameters,
        asset_prices=np.asarray(old.solution.p_eq),elasticity=.63)
    if args.case.startswith('new_'):
        from e5f_parenthood_utility import (initialize_parenthood_utility,validate_parenthood_utility,
                                          validate_parenthood_candidate,bind_parenthood_utility)
        base=initialize_parenthood_utility(base)
        validate_parenthood_utility(base)
    if 'structural_candidate' in c:
        if args.case!='new_balanced' or not c['normalize'] or c.get('observe_early') is not True:
            raise ValueError('Structural probes require normalized balanced new utility and explicit early diagnostics')
        candidate=validate_parenthood_candidate(c['structural_candidate'],require_complete=True)
        base=bind_parenthood_utility(base,candidate)
    if 'initial_psi' in c:
        if not np.isfinite(c['initial_psi']): raise ValueError('Initial preference must be finite')
        base.psi_child=float(c['initial_psi'])
    if 'observe_early' in c and type(c['observe_early']) is not bool:
        raise ValueError('Early diagnostic switch must be an explicit Boolean')
    if c.get('observe_early',False) and args.case!='new_balanced':
        raise ValueError('Early diagnostic panel requires balanced new utility')
    grid=np.asarray(old.b_grid).copy()
    if (len(grid)!=120 or int(base.J)!=17 or int(base.I)!=1
            or not np.array_equal(grid,model.make_grid(base))
            or not 0<float(base.tol_eq)<=2.5e-5):
        raise ValueError('Initial full-grid geometry or strict tolerance changed')
    # Validate the analytic formula against the inherited actual distribution
    # before invoking it for any new utility or parameter point.
    predicted,_=bind_initial_balanced_pension(old.parameters,payroll_tax=.179)
    saved_marginal=certify_initial_pension(old.stationary_g_pre,predicted,
        marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
    save('seed_mapping.json',dict(supply_rebase=supply_rebase,saved_marginal=saved_marginal,
        initial_psi=float(base.psi_child),case=args.case,source_files_verified=len(c['source_sha256'])))
    reference=None
    completed=[]
    for repetition in range(c['repetitions']):
        dest=out/f'repetition_{repetition+1:02d}'
        dest.mkdir()
        records=[]
        def ge(overrides,verbose=False):
            if set(overrides)!={'psi_child'}: raise ValueError('Normalizer changed structural primitives')
            if (len(records)>=c['maximum_stationary_solves_per_repetition']
                    or time.monotonic()-started>=c['seconds']):
                raise TimeoutError('Normalization solve/time budget exhausted')
            P=copy.deepcopy(base)
            P.psi_child=float(overrides['psi_child'])
            record=dict(index=len(records),psi_child=P.psi_child,status='started')
            records.append(record)
            state.update(phase='stationary_equilibrium',repetition=repetition+1,
                stationary_solves=state['stationary_solves']+1,psi_child=P.psi_child)
            save('heartbeat.json',dict(state,elapsed_seconds=time.monotonic()-started))
            tick=time.monotonic()
            if args.case.endswith('_balanced'):
                sol,P,price,fiscal=solve_balanced_initial_equilibrium(model=model,parameters=P,
                    b_grid=grid,initial_prices=old.solution.p_eq,payroll_tax=.179,
                    marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
            else:
                sol,P,price=model.solve_markov_income_equilibrium(old.solution.p_eq,P,grid,verbose=False)
                if not sol.converged or not sol.timings.get('strict_converged',False):
                    raise RuntimeError('Unbalanced diagnostic initial GE failed its unchanged strict gate')
                fiscal=dict(actual_accounts=fiscal_accounts(sol.g,P),fiscal_gate_required=False,
                            production_eligible=False)
            record.update(status='completed',seconds=time.monotonic()-tick,
                price=float(price[0]),market_error=float(sol.timings['best_eq_error']),fiscal=fiscal)
            save(f'repetition_{repetition+1:02d}/stationary_solves.json',records)
            return sol,P,price
        adapter=SimpleNamespace(run_model_cp_dt=ge,extract_moments=chain.extract_moments)
        sol,P,price,_,normalization=calibration.solve_old_steady_state(adapter,{},
            initial_psi=float(base.psi_child),completed_fertility_target=2.1,
            completed_fertility_tolerance=5e-4,normalize=c['normalize'])
        if args.case.startswith('new_'): validate_parenthood_utility(P)
        shared=model.precompute_shared(P,grid)
        P._fert2_probs=np.asarray(sol.fert2_probs).copy()
        policy=primitive.pf.calendar.policy_from_solution(sol,price,P,grid,shared)
        pre,reconstruction=primitive.pf.calendar.reconstruct_stationary_pre_fertility(sol,policy,P,grid,shared)
        operator=primitive.pf.transition.operator_gates(sol,policy,pre,P,grid,shared)
        operator.update(reconstruction)
        for name in ('stationary_post_fertility_nesting_l1','one_step_constant_path_nesting_l1',
                     'mature_flow_abs_error','birth_flow_abs_error','topcode_adjusted_birth_flow_abs_error'):
            if not np.isfinite(operator[name]) or abs(operator[name])>5e-9:
                raise RuntimeError(f'Stationary operator failed: {name}={operator[name]}')
        if (not np.isfinite([operator['zero_entry_mass_accounting_residual'],
                             operator['stationary_feasibility_projection_mass']]).all()
                or abs(operator['zero_entry_mass_accounting_residual'])>2e-8
                or operator['stationary_feasibility_projection_mass']>1e-6):
            raise RuntimeError('Stationary mass/feasibility gate failed')
        # Preserve the estimated initial supply curve. Age reweighting later
        # changes demand; it must not silently reset this supply intercept.
        supply=primitive.pf.calendar.HousingSupplyRule('static-elastic',float(price[0]),
            float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0]),.63)
        evaluation=primitive.pf.calendar.evaluate_period(price,pre,P,grid,shared,
            primitive.pf.calendar.SolveCounter(),supply_rule=supply,supplied_policy=policy)
        budget=primitive.dated_budget(evaluation,P,shared,grid,float(P.user_cost_rate*price[0]))
        if not np.isfinite(evaluation.relative_market_residual) or evaluation.relative_market_residual>2e-4:
            raise RuntimeError('Reconstructed initial market did not retain clearing')
        if args.case.endswith('_balanced'):
            fiscal=certify_initial_pension(evaluation.g_current,P,
                marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
        else: fiscal=dict(actual_accounts=fiscal_accounts(evaluation.g_current,P),production_eligible=False)
        moments=chain.extract_moments(sol,P)
        for name in ('tfr','childless_rate'):
            if not np.isfinite(float(moments[name])):
                raise RuntimeError(f'Required stationary diagnostic is nonfinite: {name}')
        packet=dict(parameters=P,b_grid=grid,evaluation=evaluation,shared=shared,
            supply_rule=supply,solution=sol,stationary_g_pre=pre,
            demographic_seed=inherited.get('demographics'),contract_sha256=args.contract_sha256)
        import run_e5f_independent_numerical_audit as audit
        for name,array in primitive.policy_arrays(policy).items():
            if array is not None and not np.isfinite(np.asarray(array)).all():
                raise RuntimeError(f'Policy array contains nonfinite entries: {name}')
        arrays=audit.policy_array_audit(packet,dest)
        if arrays['occupied_negative_steps']:
            raise RuntimeError('Occupied wealth/value monotonicity gate failed')
        for name,stats in arrays['probabilities'].items():
            if stats['nonfinite'] or stats['minimum'] < 0. or stats['maximum'] > 1.:
                raise RuntimeError(f'Policy probability range gate failed: {name}')
        checkpoint=dest/'initial_state.pkl.gz'
        with gzip.open(checkpoint,'wb',compresslevel=1) as stream: pickle.dump(packet,stream,protocol=5)
        with gzip.open(checkpoint,'rb') as stream: reload=pickle.load(stream)
        np.testing.assert_array_equal(reload['evaluation'].g_current,evaluation.g_current)
        np.testing.assert_array_equal(reload['evaluation'].policy.V,policy.V)
        del reload
        signature=dict(price=np.asarray(price).copy(),moments=moments,
            V=policy.V.copy(),g=evaluation.g_current.copy(),psi=float(P.psi_child))
        if reference is not None:
            for name in ('price','V','g'): np.testing.assert_array_equal(signature[name],reference[name])
            if set(signature['moments'])!=set(reference['moments']) or signature['psi']!=reference['psi']:
                raise RuntimeError('Fresh normalized candidate replay differs')
            for name in moments:
                np.testing.assert_array_equal(np.asarray(moments[name]),
                                              np.asarray(reference['moments'][name]),err_msg=name)
        else: reference=signature
        summary=dict(status='passed_initial_diagnostic',normalization=normalization,
            operator_gates=operator,household_budget=budget,fiscal=fiscal,policy_array_gates=arrays,
            price=float(price[0]),checkpoint_sha256=primitive.digest(checkpoint),
            legacy_stationary_moments= moments,calibrated_smm=False,
            measurement_warning='Legacy stationary observers are diagnostics, not the certified early target system',
            seconds=time.monotonic()-started)
        from e5f_parenthood_utility import PARENTHOOD_SEARCH_DOMAIN
        parameter_rows=[]
        for name,lower,upper,transform in PARENTHOOD_SEARCH_DOMAIN:
            estimate=(float(P.beta)**.25 if name=='beta_annual' else
                float(P.hbar_first_child_jump)+float(P.hbar_child_rooms) if name=='h_P' else
                float(np.asarray(getattr(P,name)).reshape(-1)[0]))
            parameter_rows.append(dict(parameter=name,estimate=estimate,lower=lower,upper=upper,
                transform=transform,near_bound=min(estimate-lower,upper-estimate)<=.01*(upper-lower),
                status='diagnostic candidate; not a certified estimate' if 'structural_candidate' in c else 'inherited diagnostic point; not re-estimated',
                interpretation='mapped first-child requirement' if name=='h_P' else 'structural coordinate'))
        parameter_rows.extend(dict(parameter=name,estimate=float(value),lower=None,upper=None,
            transform='',near_bound=False,status=status,interpretation='') for name,value,status in (
                ('hbar_child_rooms',P.hbar_child_rooms,'zero restriction' if args.case.startswith('new_') else 'old-utility diagnostic'),
                ('psi_child',P.psi_child,'normalized to 2.1' if c['normalize'] else 'fixed inherited intercept'),
                ('payroll_tax',P.tau_pay,'externally fixed'),('pension_period',P.pension,'budget derived' if args.case.endswith('_balanced') else 'old-pension diagnostic'),
                ('housing_supply_elasticity',P.xi_supply[0],'externally fixed'),
                ('tenure_choice_kappa',P.tenure_choice_kappa,'externally fixed'),
                ('alpha_cons',P.alpha_cons,'externally fixed'),('sigma',P.sigma,'externally fixed')))
        primitive.pf.write_csv(dest/'parameters.csv',parameter_rows)
        if c.get('observe_early',False):
            from e5f_initial_fertility_observer import observe_initial_fertility
            from e5f_initial_housing_observer import observe_initial_housing_wealth
            state.update(phase='early_measurement',repetition=repetition+1)
            fertility={projection:observe_initial_fertility(evaluation,P,age_projection=projection)
                for projection in ('uniform_birth_time','constant_post_cell')}
            housing=observe_initial_housing_wealth(evaluation,P,grid,shared,
                diagnostic_enabled=True,age_projection='uniform_within_age_cell',
                diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=True)
            early=dict(fertility=fertility,housing_wealth=housing,calibrated_smm=False,
                empirical_target_contract_activated=False,weights_assigned=False)
            save(f'repetition_{repetition+1:02d}/early_measurement.json',early)
            summary['early_measurement']=early
        save(f'repetition_{repetition+1:02d}/summary.json',summary)
        completed.append(summary)
        save('latest_completed.json',dict(repetition=repetition+1,**summary))
        save('best_so_far.json',dict(selection='same diagnostic seed; no calibrated objective',repetition=repetition+1,**summary))
        state.update(phase='diagnostics',completed_repetitions=repetition+1)
        if repetition==c['repetitions']-1:
            audit.standard_diagnostics(packet,dest,validate_production_young=False)
            if len(list((dest/'standard_diagnostics').glob('*.png')))!=17:
                raise RuntimeError('Standard 17-graph diagnostic packet incomplete')
        del sol,packet,signature
    save('summary.json',dict(status='passed_initial_candidate_loop',case=args.case,
        repetitions=len(completed),stationary_solves=state['stationary_solves'],
        normalized=c['normalize'],calibrated_smm=False,perfect_foresight_solved=False,
        elapsed_seconds=time.monotonic()-started,final=completed[-1]))
    finished.set()


if __name__=='__main__':
    try:
        main()
    except Exception as exc:
        if ACTIVE_RUN:
            ACTIVE_RUN['finished'].set()
            ACTIVE_RUN['save']('failure.json',dict(ACTIVE_RUN['state'],
                status='failed_initial_candidate_loop',error_type=type(exc).__name__,error=str(exc),
                elapsed_seconds=time.monotonic()-ACTIVE_RUN['started']))
        raise
