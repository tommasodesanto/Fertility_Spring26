"""Native evening stationary calibration: DUE credit and free tenure scale.

New source/target contract required. Ten searched coordinates plus normalized
child benefit; all fourteen moments displayed, three explicit validation rows.
No legacy generated accounting adapters or historical floor binding are used.
"""
from __future__ import annotations
import copy
import csv
import gzip
import importlib
import json
import math
import pickle
import sys
import time
from pathlib import Path
import numpy as np
import e5f_calibration_runtime as base
from e5f_calibration_runtime import write, table, sha
from run_e5f_due_stayer_matched_check import require_no_negative_estates

FREE = base.FREE + ('tenure_choice_kappa',)
VALIDATION = frozenset(('nchs_share30','family_rooms','old_dispersion'))
ROOT = Path(__file__).resolve().parents[3]


def require_current_tool(module,name):
    expected=(ROOT/'code/model/tools'/f'{name}.py').resolve()
    actual=Path(module.__file__).resolve()
    if actual!=expected:
        raise RuntimeError(f'Cached/imported tool namespace is not current: {name}: {actual}')
    return module


require_current_tool(base,'e5f_calibration_runtime')
if Path(require_no_negative_estates.__code__.co_filename).resolve() != (ROOT/'code/model/tools/run_e5f_due_stayer_matched_check.py').resolve():
    raise RuntimeError('Estate gate helper is not from current pinned source')


def require_abs_gate(value,tolerance,label):
    number=float(value)
    if not math.isfinite(number) or abs(number)>tolerance:
        raise RuntimeError(label+' gate failed: '+str(value))


def validate_objective(objective):
    rows=objective['target_rows']; names=[r['restriction_id'] for r in rows]
    if len(names)!=14 or len(set(names))!=14:raise ValueError('All14 unique target rows required')
    scored=0
    for row in rows:
        name=row['restriction_id'];w=row['actual_weight']
        if not math.isfinite(float(row['target'])):raise ValueError('Nonfinite target')
        if name=='initial_normalization':
            if w is not None:raise ValueError('Fertility normalization is a separate equality')
        elif name in VALIDATION:
            if w != 0:raise ValueError('Validation moment must have explicit zero weight: '+name)
        else:
            if w is None or not math.isfinite(float(w)) or float(w)<=0:raise ValueError('Scored moment needs positive weight')
            scored+=1
    restrictions=objective['parameter_restrictions']
    if len(restrictions)!=10 or {r['parameter'] for r in restrictions}!=set(FREE):raise ValueError('Exactly ten searched coordinates required')
    for row in restrictions:
        if not math.isfinite(float(row['lower'])) or not math.isfinite(float(row['upper'])) or row['lower']>=row['upper']:raise ValueError('Invalid parameter bounds')
        if row['parameter']=='tenure_choice_kappa' and row['lower']<=0:raise ValueError('Tenure scale requires positive log-domain bounds')
    if scored!=10:raise ValueError('Ten informative scored rows required')


def score_targets(objective,fertility,housing,recent,completed):
    validate_objective(objective)
    # The mapping and observed values are inherited unchanged. Supply positive
    # temporary weights solely to obtain the complete table, then score the
    # declared lane's weights. No value/moment/sample is changed or omitted.
    observer_view=copy.deepcopy(objective)
    for row in observer_view['target_rows']:
        if row['restriction_id'] in VALIDATION:row['actual_weight']=1.
    rows=base.score_targets(observer_view,fertility,housing,recent,completed)
    weights={r['restriction_id']:r['actual_weight'] for r in objective['target_rows']}
    for row in rows:
        weight=weights[row['moment']]
        row['role']='normalization' if weight is None else 'validation' if row['moment'] in VALIDATION else 'scored'
        row['weight']='' if weight is None else float(weight)
        row['loss_contribution']='' if weight is None else float(weight)*row['gap']**2
    return rows


def verify_sources(contract):
    manifest=contract['source_manifest']
    if sha(manifest['path'])!=manifest['sha256']:raise RuntimeError('Source manifest changed')
    inventory=json.loads(Path(manifest['path']).read_text())['files']
    root=Path(contract['source_root']).resolve()
    if root!=ROOT:raise RuntimeError('Wrong native source root')
    for relative,digest in inventory.items():
        path=(root/relative).resolve()
        if not path.is_relative_to(root) or sha(path)!=digest:raise RuntimeError('Source changed: '+relative)
    return inventory


def setup(contract,objective,output):
    """Fresh worker factory; preparation validates sources and performs no solve."""
    validate_objective(objective);verify_sources(contract)
    objective_pin=contract['objective']
    if sha(objective_pin['path'])!=objective_pin['sha256'] or json.loads(Path(objective_pin['path']).read_text())!=objective:
        raise RuntimeError('Lane objective differs from its complete pinned file')
    import e5f_current_transition_runtime as native
    require_current_tool(native,'e5f_current_transition_runtime')
    ancestry=contract['native_ancestry_contract']
    if sha(ancestry['path'])!=ancestry['sha256']:raise RuntimeError('Native ancestry changed')
    prepared=native.setup(Path(output)/'native_preparation',contract=Path(ancestry['path']),
                          reference=Path(contract['reference_case']),fixed_reference_price=True)
    # Verify target values and original nine bounds against authenticated ancestry.
    old=prepared['objective_definition']
    if objective['cps_projection']!=old['cps_projection']:raise ValueError('Fertility projection changed')
    if not 1 <= int(contract['normalization']['maximum_stationary_solves']) <= 23:raise ValueError('Normalization solve cap changed')
    oldrows={r['restriction_id']:r for r in old['target_rows']}
    if {r['restriction_id'] for r in objective['target_rows']}!=set(oldrows):raise ValueError('Target identity changed')
    for row in objective['target_rows']:
        if row['target']!=oldrows[row['restriction_id']]['target']:raise ValueError('Target value changed')
    oldbounds={r['parameter']:(r['lower'],r['upper']) for r in old['parameter_restrictions']}
    for row in objective['parameter_restrictions']:
        if row['parameter'] in oldbounds and (row['lower'],row['upper'])!=oldbounds[row['parameter']]:raise ValueError('Original parameter bounds changed')
    prepared['parameters'].native_due_stayer_credit=False
    rt=prepared['runtime'];rt['calibration']=require_current_tool(importlib.import_module('run_e5f_transition_calibration'),'run_e5f_transition_calibration')
    purchase=require_current_tool(importlib.import_module('e5f_due_purchase_audit'),'e5f_due_purchase_audit')
    estate=require_current_tool(importlib.import_module('e5f_overnight_estate_audit'),'e5f_overnight_estate_audit')
    if not purchase.SUPPORTS_NATIVE_DUE_STAYER_CREDIT or not estate.SUPPORTS_NATIVE_DUE_STAYER_CREDIT:raise RuntimeError('Origin-specific audit required')
    rt['accounting']=purchase
    result=EveningObjective(contract,objective,prepared,native.EstateAuditContract(estate))
    result.native_failure_type=rt['model'].InfeasibleThetaError
    if rt['model'].DEAD_MASS_TOL!=1e-12 or rt['model'].DEAD_VALUE_CUTOFF!=-1e9:raise RuntimeError('Scientific gates changed')
    verify_sources(contract)
    return result


class EveningObjective:
    def __init__(self,contract,objective,prepared,estate):
        self.c=contract;self.obj=objective;self.prepared=prepared
        self.rt=prepared['runtime'];self.selected=prepared['selected'];self.tax=prepared['tax']
        self.ancestor=prepared['objective'].ancestor;self.adapter=prepared['objective'].adapter;self.estate=estate
        self.due=True;self.fixed_reference=False

    def bind(self,point):
        if set(point)!=set(FREE):raise ValueError('Wrong ten-dimensional point')
        bounds={r['parameter']:(r['lower'],r['upper']) for r in self.obj['parameter_restrictions']}
        if not all(math.isfinite(float(point[k])) and lo<=point[k]<=hi for k,(lo,hi) in bounds.items()):raise ValueError('Point outside pinned bounds')
        P=copy.deepcopy(self.prepared['parameters'])
        for name,value in point.items():
            value=float(value)
            if name=='H0':P.H0=np.full(np.asarray(P.H0).shape,value)
            elif name=='beta_annual':P.beta=value**int(P.period_years);P.rho=1/P.beta-1;P.rho_hat=P.rho
            else:setattr(P,name,value)
        P.eps_fert=P.kappa_fert
        P.native_due_stayer_credit=self.due
        P.native_inherited_distribution_evidence_dir=str(self.case_output/'inherited_state_failures')
        if (P.sigma!=2 or P.alpha_cons!=.733 or P.child_room_floor
                or P.hbar_first_child_jump!=0 or P.hbar_child_rooms!=0 or P.delta_alpha!=0
                or not np.all(np.asarray(P.phi)==.8) or P.xi_supply[0]!=.63
                or getattr(P,'native_solvency_credit',False)
                or P.utility_reference_rent!=self.c['fixed']['reference_rent']):
            raise ValueError('Retained preference/credit/normalization primitives changed')
        return P

    def normalize(self,point,output,deadline):
        if self.fixed_reference:
            if self.due:raise ValueError('Fixed reference replay is DUE-off only')
            # This exact original-parameter replay is a smoke, not a candidate.
            P=self.prepared['parameters']
            original=dict(self.tax.actual_parameters(P),delta_alpha_jump=P.delta_alpha_jump,
                          child_benefit_curvature=P.child_benefit_curvature,tenure_choice_kappa=P.tenure_choice_kappa)
            if set(point)!=set(FREE) or any(abs(float(point[k])-float(original[k]))>1e-13 for k in FREE):
                raise ValueError('Fixed reference smoke must use original de_0093 coordinates')
            result=self.prepared['objective'].normalize(point,output,deadline)
            write(output/'stationary_solves.json',[dict(index=0,status='completed',psi_child=float(P.psi_child),scope='fixed original smoke')])
            return result
        return base.StationaryObjective.normalize(self,point,output,deadline)

    def evaluate(self,point,output,deadline_epoch,graphs=False,due=True,fixed_reference=False):
        self.due=bool(due);self.fixed_reference=bool(fixed_reference);self.case_output=Path(output)
        verify_sources(self.c)
        result=self.evaluate_point(tax=self.tax,objective=self.obj,selected=self.selected,runtime=self.rt,
            point=point,output=self.case_output,deadline_epoch=deadline_epoch,graphs=graphs)
        verify_sources(self.c)
        return result

    def classify_failure(self,exc,context):
        verify_sources(self.c)
        runner=sys.modules['current_reference_recovery']
        evidence=runner.recovery.capture_native_failure(exc,expected_native_type=self.native_failure_type,
            native_gate_tolerance=1e-12,context=runner.recovery.Context(**context))
        status='fatal';classification='unknown_or_integrity_failure'
        if type(exc) in (self.estate.EstateFundingShortfall,base.NonpositiveNormalizedBenefit):
            status='inadmissible';classification='explicit_economic_gate'
        elif runner.failure_status(exc)=='inadmissible_parameter_proposal':
            status='inadmissible';classification='inherited_native_prefix'
        elif evidence.narrow_infeasibility_verified:
            status='inadmissible';classification='authenticated_structured_native_gate'
        import dataclasses
        return dict(status=status,authenticated=(status=='inadmissible'),classification=classification,native=dataclasses.asdict(evidence),
                    audit=getattr(exc,'audit',getattr(exc,'ledger',None)))

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
            require_abs_gate(operator[name],5e-9,name)
        require_abs_gate(operator['zero_entry_mass_accounting_residual'],2e-8,'mass accounting')
        require_abs_gate(operator['stationary_feasibility_projection_mass'],0.,'stationary projection')
        supply=primitive.pf.calendar.HousingSupplyRule('static-elastic',float(price[0]),
            float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0])/report_population_scale,float(P.xi_supply[0]))
        evaluation=primitive.pf.calendar.evaluate_period(price,pre,P,grid,shared,
            primitive.pf.calendar.SolveCounter(),supply_rule=supply,supplied_policy=policy)
        require_abs_gate(evaluation.relative_market_residual,2e-4,'housing market')
        budget=primitive.dated_budget(evaluation,P,shared,grid,float(P.user_cost_rate*price[0]))
        purchase=rt['accounting'].audit_purchase_accounting(evaluation,P,shared,grid,model)
        fiscal=rt['certify_initial_pension'](evaluation.g_current,P,marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
        estate=self.estate.audit(evaluation,P,grid)
        write(output/'estate_funding.json',estate)
        require_no_negative_estates(estate)
        if estate['status'] != 'funded': raise RuntimeError('Estate funding gate failed')
        require_abs_gate(budget['budget_excess_mass'],2e-10,'household budget')
        require_abs_gate(purchase['maximum_occupied_transaction_wealth_error'],1e-9,'transaction wealth')
        for key,value in purchase.items():
            if key.endswith('violation_mass') or key in ('transaction_outside_grid_mass','negative_estate_exposure_mass','saving_outside_grid_mass'):
                require_abs_gate(value,2e-10,'purchase '+key)
        packet=dict(parameters=P,b_grid=grid,evaluation=evaluation,shared=shared,supply_rule=supply,
            solution=sol,stationary_g_pre=pre,demographic_seed=selected.get('demographic_seed'),
            contract_sha256=self.c['objective']['sha256'],ancestry_contract_sha256=selected.get('contract_sha256'))
        if self.fixed_reference:self.last_packet=packet
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
                    child_benefit_curvature=P.child_benefit_curvature,tenure_choice_kappa=P.tenure_choice_kappa)
        params=[]
        for restriction in objective['parameter_restrictions']:
            name=restriction['parameter']; value=actual[name]; lo=restriction['lower']; hi=restriction['upper']
            params.append(dict(parameter=name,estimate=value,lower=lo,upper=hi,
                near_bound=min(value-lo,hi-value)<=.01*(hi-lo),status='free in evening DUE calibration'))
        for name,value,status in (
            ('psi_child',P.psi_child,'normalized to completed fertility 2.1'),
            ('child_benefit_CRRA_coefficient',(1-P.child_benefit_curvature)*P.psi_child,'derived from normalized one-child benefit'),
            ('theta1',P.theta1,'fixed external restriction'),('sigma',P.sigma,'fixed'),
            ('alpha_cons',P.alpha_cons,'fixed CEX childless expenditure share'),
            ('delta_alpha',P.delta_alpha,'fixed zero later-child loading'),
            ('h_P',0.,'no housing floor'),('utility_reference_rent',P.utility_reference_rent,'fixed substantive utility normalization'),
            ('q_annual',(1+P.q)**(1/P.period_years)-1,'author-retained 2% annual real rate'),
            ('financed_share',P.phi[0],'inherited credit contract'),('housing_supply_elasticity',P.xi_supply[0],'fixed provisional external mapping'),
            ('payroll_tax',P.tau_pay,'derived from adopted pension ratio'),('pension_period',P.pension,'balanced PAYGO'),
            ('annual_depreciation',self.ancestor.ANNUAL_DEP,'adopted'),('period_depreciation',P.delta,'compounded'),
            ('annual_property_tax',self.ancestor.ANNUAL_PROPERTY_TAX,'adopted'),('period_property_tax',P.tau_H,'linear period convention'),
            ('selling_cost',P.psi,'retained'),('rental_cap',P.hR_max,'retained provisional'),
            ('wealth_grid_nodes',len(grid),'retained exact grid'),('income_states',len(P.z_grid),'retained B15')):
            params.append(dict(parameter=name,estimate=float(value),lower='',upper='',near_bound='',status=status))
        if len(params)!=31 or len({r['parameter'] for r in params})!=31: raise RuntimeError('Parameter table must retain 31 unique rows')
        table(output/'target_fit.csv',rows); table(output/'parameters.csv',params)
        write(output/'observers.json',dict(fertility=fertility,housing_wealth=housing,recent_parent=recent))
        loss=sum(r['loss_contribution'] for r in rows if r['loss_contribution']!='')
        receipt=dict(status='verified_provisional_calibration_point',loss=loss,point=point,
            normalization=normalization,normalization_inputs=self.c['normalization'],
            target_system_sha256=self.c['objective']['sha256'],source_manifest_sha256=self.c['source_manifest']['sha256'],
            selected_checkpoint_sha256=self.prepared['reference_receipt']['case_checkpoint_sha256'],case_checkpoint_sha256=checkpoint_sha,
            free_count=len(FREE),weighted_count=sum(r['weight']!='' and float(r['weight'])>0 for r in rows),display_count=len(rows),
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
        receipt['credit_rule']='DUE existing-owner with death solvency' if P.native_due_stayer_credit else 'baseline smoke only'
        receipt['free_count']=10; receipt['normalized_count']=1
        receipt['policy_plot_scope']='Standard saving/consumption panels are buyer-conditional; stayer policies retained separately'
        write(output/'receipt.json',receipt)
        if graphs:
            rt['audit'].standard_diagnostics(packet,output,validate_production_young=False)
            assert len(list((output/'standard_diagnostics').glob('*.png')))==17
        return receipt


def baseline_authentication_smoke(contract,objective,output,deadline_epoch,graphs=False):
    """DUE-off original-price/psi replay; compare moments, not changed weights."""
    output=Path(output);evaluator=setup(contract,objective,output)
    def csv_rows(path):
        with Path(path).open(newline='') as stream:return list(csv.DictReader(stream))
    reference=Path(contract['reference_case'])
    old_parameters={r['parameter']:r for r in csv_rows(reference/'parameters.csv')}
    point={name:float(old_parameters[name]['estimate']) for name in FREE}
    receipt=evaluator.evaluate(point,output/'case',deadline_epoch,graphs=graphs,due=False,fixed_reference=True)
    current={r['parameter']:r for r in csv_rows(output/'case/parameters.csv')}
    old_targets={r['moment']:r for r in csv_rows(reference/'target_fit.csv')}
    targets={r['moment']:r for r in csv_rows(output/'case/target_fit.csv')}
    if len(current)!=31 or set(current)!=set(old_parameters):raise RuntimeError('Baseline parameter identities differ')
    if len(targets)!=14 or set(targets)!=set(old_targets):raise RuntimeError('Baseline target identities differ')
    comparison={'parameters':{},'targets':{},'weights_compared':False,
                'reason':'New lane weights intentionally differ; all moments and parameters must nest'}
    for name,row in current.items():
        gap=float(row['estimate'])-float(old_parameters[name]['estimate'])
        comparison['parameters'][name]=dict(original=float(old_parameters[name]['estimate']),current=float(row['estimate']),gap=gap)
    for name,row in targets.items():
        comparison['targets'][name]={key:dict(original=float(old_targets[name][key]),current=float(row[key]),
            gap=float(row[key])-float(old_targets[name][key])) for key in ('target','model','gap')}
    write(output/'baseline_comparison.json',comparison)
    for name,row in comparison['parameters'].items():require_abs_gate(row['gap'],0.,'baseline parameter '+name)
    for name,row in comparison['targets'].items():
        for key,values in row.items():require_abs_gate(values['gap'],0. if key=='target' else 1e-13,'baseline '+name+' '+key)
    import e5f_current_transition_runtime as native
    arrays=native.compare_arrays(evaluator.selected,evaluator.last_packet)
    write(output/'baseline_array_comparison.json',arrays)
    if arrays['array_count']!=107:raise RuntimeError('Baseline array census changed; requires review')
    retained_tail={'evaluation.g_pre','evaluation.g_post_fertility','evaluation.g_current'}
    for name,row in arrays['arrays'].items():
        if name in retained_tail:
            if not row.get('finite'):raise RuntimeError('Nonfinite baseline distribution')
            require_abs_gate(row['l1'],2e-12,'baseline retained-tail '+name)
        elif row.get('exact') is not True:raise RuntimeError('Baseline core array differs: '+name)
    write(output/'baseline_authentication.json',dict(status='passed',due=False,fixed_reference_price=True,
        target_rows=14,parameter_rows=31,array_count=107,core_arrays_exact=True,distribution_tail_tolerance=2e-12,loss_not_compared=True,receipt_sha256=sha(output/'case/receipt.json')))
    return receipt


def main():
    import argparse
    parser=argparse.ArgumentParser(description='Torch-only baseline authentication prerequisite')
    parser.add_argument('--contract',type=Path,required=True);parser.add_argument('--lane',required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--deadline-epoch',type=float,required=True)
    parser.add_argument('--graphs',action='store_true');args=parser.parse_args()
    if time.time()>=args.deadline_epoch:raise TimeoutError('Smoke deadline already expired')
    contract=json.loads(args.contract.read_text());contract=dict(contract,objective=contract['lanes'][args.lane]['objective'])
    objective=json.loads(Path(contract['objective']['path']).read_text())
    baseline_authentication_smoke(contract,objective,args.output,args.deadline_epoch,args.graphs)
if __name__=='__main__':main()
