"""Explicit single-specification stationary objective on a pinned source copy.

The inherited purchase/entry adapters remain authenticated compatibility code.
Preferences, normalization options, observers and the scored target mapping are
ordinary source functions here; no evaluator source text is rewritten.
"""
from __future__ import annotations

import copy
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import pickle
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

COMMON = ('H0','beta_annual','chi','first_birth_fixed_cost','kappa_fert',
          'kappa_fert_continuation','theta0')
FREE = COMMON + ('delta_alpha_jump','child_benefit_curvature')


class NonpositiveNormalizedBenefit(RuntimeError):
    """This parameter proposal cannot retain a positive normalized child benefit."""

    def __init__(self, normalization):
        self.audit=dict(normalization=normalization)
        super().__init__('Normalized child benefit must be positive')


def validate_normalization(normalization, psi_child, solve_count):
    if normalization['stationary_solves']!=solve_count:
        raise RuntimeError('Invalid normalization solve accounting')
    if not math.isfinite(float(psi_child)):
        raise RuntimeError('Nonfinite normalized child benefit')
    if psi_child<=0:
        raise NonpositiveNormalizedBenefit(normalization)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''): h.update(block)
    return h.hexdigest()


def write(path, value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')
    temp.replace(path)


def table(path, rows):
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fields); writer.writeheader(); writer.writerows(rows)


def load_module(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec); sys.modules[name]=result
    spec.loader.exec_module(result)
    return result


def load_seed(c, path):
    """Native Torch loader, or a pinned namespace-only local compatibility loader."""
    if c.get('execution', {}).get('kind', 'slurm') == 'local':
        record = c['files']['local_runtime_adapter']
        if sha(record['path']) != record['sha256']:
            raise RuntimeError('Local compatibility adapter changed')
        adapter = load_module('calibration_local_compatibility', record['path'])
        return adapter.load_checkpoint(path)
    with gzip.open(path, 'rb') as stream:
        return pickle.load(stream)


def score_targets(objective, fertility, housing, recent, completed):
    """Every row is explicit; missing required moments fail before scoring."""
    fertility=fertility[objective['cps_projection']]['moments']
    housing=housing['moments']
    values={
        'initial_normalization':completed,
        'cps_childlessness':fertility['childless_rate_40_44'],
        'cps_exactly_one':fertility['exactly_one_among_mothers_40_44'],
        'nchs_mean_age':fertility['period_mean_age_first_birth'],
        'nchs_share30':fertility['period_share_first_births_age30plus'],
        'early_fertility':fertility['mean_children_ever_born_capped3_age25'],
        'wealth_earnings':housing['aggregate_wealth_to_annual_gross_labor_earnings'],
        'bequest_wealth':housing['annual_bequest_flow_to_aggregate_wealth'],
        'old_dispersion':housing['old_total_wealth_to_annual_income_p90_p50_7684'],
        'mean_rooms':housing['aggregate_mean_occupied_rooms_ahs_uncapped_18_85'],
        'ownership_30_55':housing['own_rate_30_55'],
        'first_birth_rooms':housing['housing_increment_0to1'],
        'family_rooms':housing['prime30_55_model_dependent_3plus_minus_1to2_rooms_capped9'],
        'recent_parent_ownership':recent,
    }
    names=[row['restriction_id'] for row in objective['target_rows']]
    if len(names)!=len(set(names)) or set(names)!=set(values):
        raise ValueError('Target registry must contain each declared moment exactly once')
    rows=[]
    for source in objective['target_rows']:
        name=source['restriction_id']; value=values[name]
        if value is None or not math.isfinite(float(value)):
            raise RuntimeError('Required model moment unavailable: '+name)
        value=float(value); target=float(source['target']); weight=source['actual_weight']
        if not math.isfinite(target) or (weight is not None and
                (not math.isfinite(float(weight)) or float(weight)<=0)):
            raise ValueError('Invalid target/weight: '+name)
        if (weight is None)!=(name=='initial_normalization'):
            raise ValueError('Only the fertility normalization may be unscored')
        gap=value-target
        rows.append(dict(moment=name,target=target,model=value,gap=gap,
            weight='' if weight is None else weight,
            loss_contribution='' if weight is None else float(weight)*gap**2))
    return rows


def setup(c, obj, point, out):
    """Authenticate ancestry, then import and use the new pinned source once."""
    sys.path[:0]=[c['runtime_tools']]
    runner=load_module('calibration_recovery_runner',c['files']['recovery_runner']['path'])
    base=json.loads(Path(c['base_contract']['path']).read_text())
    pair,ancestor,lock,_,_=runner.pair_runtime(base['reference_root'])
    source=Path(c['source_root']).resolve()
    manifest=c['source_manifest']
    if sha(manifest['path'])!=manifest['sha256']:
        raise RuntimeError('New source manifest changed')
    inventory=json.loads(Path(manifest['path']).read_text())
    for relative, expected in inventory['files'].items():
        path=(source/relative).resolve()
        if not path.is_relative_to(source) or sha(path)!=expected:
            raise RuntimeError('New source differs: '+relative)
    sys.path[:0]=[str(source/'code/model/tools'),str(source/'code/model')]
    tax=ancestor.load_tax_driver()
    if sha(pair.PLAN)!=tax.PLAN_SHA or sha(pair.CHECKPOINT)!=tax.CHECKPOINT_SHA:
        raise RuntimeError('Seed ancestry changed')
    plan=json.loads(Path(pair.PLAN).read_text())
    selected=load_seed(c,pair.CHECKPOINT)
    tax.check_checkpoint(selected)
    ancestor.SOURCE=source
    rt=ancestor.setup_runtime(tax,plan,selected,out)
    if Path(rt['model'].__file__).resolve()!=source/'code/model/intergen_eqscale_seq_optimized/solver.py':
        raise RuntimeError('Executing model is not the new source')
    source_view=dict(reference_root=str(source.parent),parent_source_inventory=manifest)
    native,evidence=runner.authenticated_native_failure_type(source_view,rt)
    estate=load_module('calibration_estate_audit',c['files']['estate_audit']['path'])
    inherited_status=runner.failure_status
    runner.failure_status=lambda exc: ('inadmissible_parameter_proposal'
        if type(exc) in (estate.EstateFundingShortfall,NonpositiveNormalizedBenefit)
        else inherited_status(exc))
    runtime=StationaryObjective(c,obj,selected,tax,ancestor,rt,runner.adapter,estate)
    # Zero-solve checks are run before the controller starts an objective.
    P=runtime.bind(point)
    if c['normalization']['warm_price'] and (P.I!=1 or
            str(getattr(P,'markov_equilibrium_method','legacy_damped')).lower()
            not in ('direct','direct_brent','brent')):
        raise ValueError('Warm-price contract requires one-market direct equilibrium')
    shared=rt['model'].precompute_shared(P,selected['b_grid'])
    np.testing.assert_array_equal(shared.type_psi[shared.type_map],shared.psi_flat.ravel())
    rt['estate_funding_shortfall_type']=estate.EstateFundingShortfall
    return runtime,tax,selected,rt,runner,native,evidence,{'runtime':'explicit_native_v1'}


class StationaryObjective:
    write=staticmethod(write)

    def __init__(self,c,obj,selected,tax,ancestor,rt,adapter,estate):
        self.c,self.obj,self.selected,self.tax=c,obj,selected,tax
        self.ancestor,self.rt,self.adapter,self.estate=ancestor,rt,adapter,estate

    def bind(self,point):
        from e5f_parenthood_utility import bind_parenthood_utility
        if set(point)!=set(FREE): raise ValueError('Wrong free coordinates')
        bounds={r['parameter']:(r['lower'],r['upper']) for r in self.obj['parameter_restrictions']}
        if set(bounds)!=set(FREE) or not all(math.isfinite(point[k]) and lo<=point[k]<=hi
                                              for k,(lo,hi) in bounds.items()):
            raise ValueError('Point outside explicit calibration bounds')
        binding={name:point[name] for name in COMMON}
        binding['h_P']=float(self.selected['parameters'].hbar_first_child_jump)
        P=bind_parenthood_utility(self.selected['parameters'],binding)
        P.theta1=self.ancestor.THETA1
        P.delta=1-(1-self.ancestor.ANNUAL_DEP)**int(P.period_years)
        P.tau_H=self.ancestor.ANNUAL_PROPERTY_TAX*int(P.period_years)
        P.user_cost_rate=P.q+P.delta+P.tau_H
        P.adult_entry_clock='split_birth_vintage'
        P.child_room_floor=False
        P.hbar_first_child_jump=P.hbar_child_rooms=0.
        P.delta_alpha=0.; P.delta_alpha_jump=float(point['delta_alpha_jump'])
        P.child_benefit_curvature=float(point['child_benefit_curvature'])
        P.compensated_child_housing_shares=True
        P.utility_reference_rent=float(self.c['fixed']['reference_rent'])
        for name in ('utility_comparison_arm','utility_child_benefit_exponent'):
            if hasattr(P,name): delattr(P,name)
        if (P.sigma!=2 or P.alpha_cons!=.733 or P.tenure_choice_kappa!=.005
                or not np.all(np.asarray(P.phi)==.8) or P.xi_supply[0]!=.63):
            raise ValueError('A retained fixed preference/credit/supply input changed')
        return P

    def normalize(self,point,output,deadline):
        from intergen_eqscale_seq_optimized.adult_entry import require_closed_stationary_renewal
        base=self.bind(point); options=self.c['normalization']
        calibration=self.rt['calibration']; previous=calibration.closure
        ledger=[]; warm={} if options['warm_price'] else None
        payroll_tax,fiscal_rule=self.adapter.pension_tax_from_demographics(base)
        def solve_ge(chain,overrides):
            # The supervising process owns the absolute deadline and process
            # group kill. An inner timeout would race it and look like a fatal
            # scientific error instead of an authenticated censored case.
            if len(ledger)>=options['maximum_stationary_solves']:
                raise RuntimeError('Old-steady-state fertility normalization missed tolerance: 23-solve cap')
            P=copy.deepcopy(base); P.psi_child=float(overrides['psi_child'])
            record=dict(index=len(ledger),psi_child=P.psi_child,status='started',epoch=time.time(),
                        incoming_price_state=copy.deepcopy(warm))
            ledger.append(record); write(output/'stationary_solves.json',ledger)
            start=time.monotonic()
            try:
                sol,P,price,fiscal=self.rt['solve_balanced_initial_equilibrium'](
                    model=self.rt['model'],parameters=P,b_grid=self.selected['b_grid'],
                    initial_prices=self.selected['solution'].p_eq,payroll_tax=payroll_tax,
                    marginal_tolerance=1e-9,fiscal_tolerance=1e-6,warm_price_state=warm)
                self.adapter.verify_pension_ratio(fiscal)
            except Exception as exc:
                record.update(status='failed',error_type=type(exc).__name__,error=str(exc),seconds=time.monotonic()-start)
                write(output/'stationary_solves.json',ledger); raise
            seconds=time.monotonic()-start
            record.update(status='completed',seconds=seconds,price=float(price[0]),
                market_residual=float(sol.timings['best_eq_error']),
                price_evaluations=int(sol.timings['unique_fast_price_evaluations']),
                price_method=sol.timings['equilibrium_method'],
                warm_price_search=sol.timings.get('warm_price_search'))
            write(output/'stationary_solves.json',ledger)
            return sol,P,price,seconds
        calibration.closure=SimpleNamespace(solve_ge=solve_ge)
        try:
            result=calibration.solve_old_steady_state(self.rt['chain'],{},
                initial_psi=options['initial_psi'],initial_step=options['initial_step'],
                completed_fertility_target=2.1,completed_fertility_tolerance=5e-4,normalize=True)
        finally: calibration.closure=previous
        sol,P,price,seconds,norm=result
        validate_normalization(norm,P.psi_child,len(ledger))
        renewal=require_closed_stationary_renewal(sol.entry_rate,sol.adult_entry_adjusted_birth_children,5e-4)
        return result,fiscal_rule,renewal

    def evaluate_point(self,*,tax,objective,selected,runtime,point,output,deadline_epoch,graphs=False,
                       report_population_scale=1.0):
        # Stationary packets have unit household mass. For a closed endpoint,
        # express unchanged absolute supply per household only in this reporter.
        if not math.isfinite(float(report_population_scale)) or report_population_scale <= 0:
            raise ValueError('Report population scale must be finite and positive')
        output.mkdir(parents=True,exist_ok=False)
        (sol,P,price,seconds,normalization),fiscal_rule,renewal=self.normalize(point,output,deadline_epoch)
        rt=self.rt; primitive=rt['primitive']; model=rt['model']; grid=selected['b_grid']
        shared=model.precompute_shared(P,grid); P._fert2_probs=sol.fert2_probs.copy()
        policy=primitive.pf.calendar.policy_from_solution(sol,price,P,grid,shared)
        pre,reconstruction=primitive.pf.calendar.reconstruct_stationary_pre_fertility(sol,policy,P,grid,shared)
        operator=primitive.pf.transition.operator_gates(sol,policy,pre,P,grid,shared)
        operator.update(reconstruction)
        for name in ('stationary_post_fertility_nesting_l1','one_step_constant_path_nesting_l1',
                     'mature_flow_abs_error','birth_flow_abs_error','topcode_adjusted_birth_flow_abs_error'):
            assert abs(operator[name])<=5e-9,name
        assert abs(operator['zero_entry_mass_accounting_residual'])<=2e-8
        assert operator['stationary_feasibility_projection_mass']<=1e-6
        supply=primitive.pf.calendar.HousingSupplyRule('static-elastic',float(price[0]),
            float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0])/report_population_scale,float(P.xi_supply[0]))
        evaluation=primitive.pf.calendar.evaluate_period(price,pre,P,grid,shared,
            primitive.pf.calendar.SolveCounter(),supply_rule=supply,supplied_policy=policy)
        assert evaluation.relative_market_residual<=2e-4
        budget=primitive.dated_budget(evaluation,P,shared,grid,float(P.user_cost_rate*price[0]))
        purchase=rt['accounting'].audit_purchase_accounting(evaluation,P,shared,grid,model)
        fiscal=rt['certify_initial_pension'](evaluation.g_current,P,marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
        estate=self.estate.audit(evaluation,P,grid)
        write(output/'estate_funding.json',estate)
        packet=dict(parameters=P,b_grid=grid,evaluation=evaluation,shared=shared,supply_rule=supply,
            solution=sol,stationary_g_pre=pre,demographic_seed=selected.get('demographic_seed'),
            contract_sha256=self.c['objective']['sha256'],ancestry_contract_sha256=selected.get('contract_sha256'))
        arrays=rt['audit'].policy_array_audit(packet,output)
        assert arrays['occupied_negative_steps']==0
        assert all(not x['nonfinite'] and x['minimum']>=0 and x['maximum']<=1 for x in arrays['probabilities'].values())
        fertility={p:rt['observe_initial_fertility'](evaluation,P,age_projection=p)
                   for p in ('uniform_birth_time','constant_post_cell')}
        housing=rt['observe_initial_housing_wealth'](evaluation,P,grid,shared,
            diagnostic_enabled=True,age_projection='uniform_within_age_cell',
            diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=True)
        checkpoint=output/'initial_state.pkl.gz'
        with gzip.open(checkpoint,'wb',compresslevel=1) as stream: pickle.dump(packet,stream,protocol=5)
        checkpoint_sha=sha(checkpoint)
        recent=rt['observe_recent_parent_flow'](evaluation,P,diagnostic_enabled=True,
            snapshot=rt['SNAPSHOT'],age_projection=rt['AGE_PROJECTION'],diagnostic_allow_residence_proxy=True,
            input_provenance=dict(case_id=output.name,checkpoint_sha256=checkpoint_sha))
        rows=score_targets(objective,fertility,housing,recent['model_value'],normalization['completed_fertility'])
        actual=dict(tax.actual_parameters(P),delta_alpha_jump=P.delta_alpha_jump,
                    child_benefit_curvature=P.child_benefit_curvature)
        params=[]
        for restriction in objective['parameter_restrictions']:
            name=restriction['parameter']; value=actual[name]; lo=restriction['lower']; hi=restriction['upper']
            params.append(dict(parameter=name,estimate=value,lower=lo,upper=hi,
                near_bound=min(value-lo,hi-value)<=.01*(hi-lo),status='free in provisional first fit'))
        for name,value,status in (
            ('psi_child',P.psi_child,'normalized to completed fertility 2.1'),
            ('child_benefit_CRRA_coefficient',(1-P.child_benefit_curvature)*P.psi_child,'derived from normalized one-child benefit'),
            ('theta1',P.theta1,'fixed external restriction'),('sigma',P.sigma,'fixed'),
            ('alpha_cons',P.alpha_cons,'fixed CEX childless expenditure share'),
            ('delta_alpha',P.delta_alpha,'fixed zero later-child loading'),
            ('h_P',0.,'no housing floor'),('utility_reference_rent',P.utility_reference_rent,'fixed substantive utility normalization'),
            ('tenure_choice_kappa',P.tenure_choice_kappa,'inherited provisional; review tomorrow'),
            ('q_annual',(1+P.q)**(1/P.period_years)-1,'inherited provisional; review tomorrow'),
            ('financed_share',P.phi[0],'inherited credit contract'),('housing_supply_elasticity',P.xi_supply[0],'fixed provisional external mapping'),
            ('payroll_tax',P.tau_pay,'derived from adopted pension ratio'),('pension_period',P.pension,'balanced PAYGO'),
            ('annual_depreciation',self.ancestor.ANNUAL_DEP,'adopted'),('period_depreciation',P.delta,'compounded'),
            ('annual_property_tax',self.ancestor.ANNUAL_PROPERTY_TAX,'adopted'),('period_property_tax',P.tau_H,'linear period convention'),
            ('selling_cost',P.psi,'retained'),('rental_cap',P.hR_max,'retained provisional'),
            ('wealth_grid_nodes',len(grid),'retained exact grid'),('income_states',len(P.z_grid),'retained B15')):
            params.append(dict(parameter=name,estimate=float(value),lower='',upper='',near_bound='',status=status))
        table(output/'target_fit.csv',rows); table(output/'parameters.csv',params)
        write(output/'observers.json',dict(fertility=fertility,housing_wealth=housing,recent_parent=recent))
        loss=sum(r['loss_contribution'] for r in rows if r['loss_contribution']!='')
        receipt=dict(status='verified_provisional_calibration_point',loss=loss,point=point,
            normalization=normalization,normalization_inputs=self.c['normalization'],
            target_system_sha256=self.c['objective']['sha256'],source_manifest_sha256=self.c['source_manifest']['sha256'],
            selected_checkpoint_sha256=tax.CHECKPOINT_SHA,case_checkpoint_sha256=checkpoint_sha,
            free_count=len(FREE),weighted_count=sum(r['weight']!='' for r in rows),display_count=len(rows),
            price=float(price[0]),market_residual=evaluation.relative_market_residual,
            fiscal=fiscal,fiscal_rule=fiscal_rule,operator_gates=operator,household_budget=budget,
            purchase_accounting=purchase,policy_array_gates=arrays,estate_funding=estate,
            adult_entry_gate=renewal,chosen_solve_seconds=seconds,
            objective_stationary_solves=normalization['stationary_solves'],
            objective_stationary_solve_seconds=normalization['stationary_solve_seconds'],
            model_observer_warnings=self.c['pending_observer_mismatches'],
            economic_changes=self.c['economic_changes'],
            retained_owner_grid=np.asarray(P.H_own).tolist(),retained_conception_schedule=np.asarray(model.get_fecundity_by_age(P)).tolist())
        if report_population_scale != 1.0:
            receipt['stationary_report_units'] = dict(
                population_scale=float(report_population_scale),
                supply='absolute supply divided by population; economic H0 unchanged',
                distribution='unit household mass; multiply by population for terminal-distance checks')
        receipt=tax.finite_json(primitive.pf.calendar.jsonable(receipt))
        write(output/'receipt.json',receipt)
        if graphs:
            rt['audit'].standard_diagnostics(packet,output,validate_production_young=False)
            assert len(list((output/'standard_diagnostics').glob('*.png')))==17
        return receipt
