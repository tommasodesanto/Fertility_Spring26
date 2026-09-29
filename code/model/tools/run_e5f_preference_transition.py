#!/usr/bin/env python3
"""One or four announced preference shocks from the immutable block0506 state.

Preparation, validation and no-shock tests are separate from execution. The
supplied draft plan cannot run a shocked path. No preference normalization or
credit-rule replacement is performed here. Execute imports/numerics on Torch.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import pickle
import shutil
import signal
import sys
import time

ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
LABEL = '2007 stationary reference — block0506, September 28 verified export'
MANIFEST = ROOT/'output/model/fertility_identification_20260928/fixed_reference_manifest.json'
MANIFEST_SHA = '147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4'
SOURCE_NAMES = ('run_e5f_preference_transition.py', 'e5f_four_shock_acceleration.py',
                'e5f_exact_policy_cache.py', 'e5f_social_security_root.py',
                'e5f_matched_pf_path_root.py', 'e5f_ssj_scaled_step_root.py',
                'e5f_ssj_toeplitz_jacobian.py')
GENERATED_FIELDS = {'_bp_pol_stay', '_c_pol_stay', '_entry_censored_mass',
    '_entry_total_mass', '_fert2_probs', '_first_births_by_age',
    '_g_stay_distribution', '_joint_choice', '_second_at_risk_by_age',
    '_second_attempts_by_age', '_second_births_by_age', '_third_at_risk_by_age',
    '_third_attempts_by_age', '_third_births_by_age'}
TERMINAL_KEYS = {'population_relative_gap', 'normalized_distribution_l1',
    'birth_queue_maximum_relative_gap', 'asset_price_relative_gap',
    'renter_price_relative_gap', 'psi_absolute_gap'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def plain(x):
    if hasattr(x, 'tolist'):
        return plain(x.tolist())
    if isinstance(x, dict):
        return {str(k): plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [plain(v) for v in x]
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(plain(value), indent=2, sort_keys=True, allow_nan=False)+'\n')
    temporary.replace(path)


def dump_checkpoint(path, value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temporary=Path(str(path)+'.tmp')
    with gzip.open(temporary,'wb',compresslevel=1) as stream:
        pickle.dump(value,stream,protocol=5)
    temporary.replace(path)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def pinned(record):
    p = Path(record['path'])
    require(sha(p) == record['sha256'], 'Pinned artifact changed: '+str(p))
    return p


def serialized(value):
    if hasattr(value, 'shape') and getattr(value, 'size', 0) > 1024:
        digest = hashlib.sha256()
        for start in range(0, value.size, 65536):
            digest.update(value.flat[start:start+65536].tobytes())
        return dict(serialized_array=True, shape=list(value.shape), dtype=str(value.dtype),
                    size=int(value.size), sha256_c_order_bytes=digest.hexdigest())
    if hasattr(value, 'tolist'):
        return serialized(value.tolist())
    if isinstance(value, dict):
        return {str(k): serialized(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialized(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return dict(nonfinite_float=repr(value))
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError('Unserialized parameter '+str(type(value)))


def shock_path(spec, horizon):
    """Preference levels, all announced at 2007; last level persists forever."""
    require(type(horizon) is int and horizon >= 1, 'Positive integer horizon required')
    require(spec['expectations'] == 'perfect_foresight_announced_at_start',
            'This engine does not implement successive unexpected shocks')
    years = [2007] if spec['kind'] == 'one_permanent' else [2007, 2011, 2015, 2019]
    require(spec['kind'] in ('one_permanent', 'four_announced'), 'One or four shocks required')
    require(spec['years'] == years and spec['period_years'] == 4, 'Shock dates/period units differ')
    levels = spec['levels']
    require(isinstance(levels, list) and len(levels) == len(years), 'Explicit shock levels are missing')
    require(all(type(v) in (int, float) and math.isfinite(v) and v > 0 for v in levels),
            'Shock levels must be finite positive saved-psi units')
    require(horizon >= len(levels), 'Horizon ends before the last shock')
    require(spec.get('interpretation') in ('author_supplied', 'reestimated_current_reference',
                                        'explicit_illustrative_transport'),
            'Shock provenance/interpretation must be explicit')
    return levels + [levels[-1]]*(horizon-len(levels))


def draft_plan(kind='four_announced'):
    return dict(schema='block0506_preference_pf_v1', reference_label=LABEL,
        reference_manifest_sha256=MANIFEST_SHA, execution_enabled=False,
        shocks=dict(kind=kind, years=[2007] if kind == 'one_permanent' else [2007,2011,2015,2019],
                    levels=None, period_years=4, expectations='perfect_foresight_announced_at_start',
                    interpretation=None, provenance=None),
        housing='fixed_stock', credit='saved_reference_unchanged',
        fiscal='fixed_payroll_tax_endogenous_pension', property_rebate=0,
        outside_entry=0, retention=1, birth_to_entry_conversion=1/2.1,
        classifications={
            'preferences':'estimated; frozen except explicit announced psi levels',
            'earnings':'externally estimated; saved B15 unchanged',
            'entry':'empirically normalized conditional entry law; split 16/20-year clock',
            'population':'normalized initial distribution; endogenous subsequent levels',
            'fiscal':'externally normalized payroll tax fixed; dated pension solved',
            'geography':'externally fixed one-market closed population',
            'estate':'provisional net-estate funding of next entrants; settlement outstanding',
            'housing':'author-fixed physical reference stock; elastic_reference is explicit alternative',
            'targets':'frozen reference measurement system'},
        budget=dict(horizon=None, max_evaluations=None, total_seconds=None,
                    case_seconds=None, cache_max_bytes=2*1024**3,
                    observed_mapping_seconds=None, maximum_policy_calls=None),
        numerics=dict(initial_prices=None, initial_pensions=None, price_bounds=None,
            pension_bounds=None, market_tolerance=2e-4, fiscal_tolerance=1e-6,
            final_reproduction_tolerance=1e-10, market_slope=None, fiscal_slope=None,
            max_log_step=None, damping=None, max_condition_number=None, worsening_factor=None,
            terminal_tolerances=None, raw_queue_relative_tolerance=None),
        endpoint=None, readiness_receipt=None, initial_jacobian=None, source_pins={},
        historical_shocks_automatically_imported=False)


def validate_plan(plan, launching=False):
    require(plan['schema'] == 'block0506_preference_pf_v1' and plan['reference_label'] == LABEL,
            'Wrong transition/reference contract')
    require(plan['reference_manifest_sha256'] == MANIFEST_SHA, 'Wrong frozen reference')
    require(plan['credit'] == 'saved_reference_unchanged', 'Credit experiment belongs to the other engine')
    require(plan['housing'] in ('fixed_stock','elastic_reference'), 'Explicit housing closure required')
    require(plan['fiscal'] == 'fixed_payroll_tax_endogenous_pension' and
            plan['property_rebate'] == 0 and plan['outside_entry'] == 0 and plan['retention'] == 1 and
            plan['birth_to_entry_conversion'] == 1/2.1, 'Population/fiscal contract changed')
    require(not plan['historical_shocks_automatically_imported'], 'Old estimates cannot be adopted implicitly')
    required = {'preferences','earnings','entry','population','fiscal','geography','estate','housing','targets'}
    require(set(plan['classifications']) == required and all(plan['classifications'].values()),
            'Every closure object needs its classification')
    missing = []
    for section in ('budget','numerics'):
        missing += [section+'.'+k for k,v in plan[section].items() if v is None]
    missing += [k for k in ('endpoint','readiness_receipt') if not plan[k]]
    missing += ['shocks.'+k for k in ('levels','interpretation','provenance') if not plan['shocks'][k]]
    if not plan['source_pins']:
        missing.append('source_pins')
    if launching:
        require(plan['execution_enabled'] is True, 'Shocked transition execution is disabled')
        require(not missing, 'Unresolved launch inputs: '+', '.join(missing))
        b, n = plan['budget'], plan['numerics']
        shock_path(plan['shocks'], b['horizon'])
        require(type(b['max_evaluations']) is int and b['max_evaluations'] >= 2, 'Finite mapping count required')
        require(all(type(b[k]) in (int,float) and math.isfinite(b[k]) and b[k] > 0
                    for k in ('total_seconds','case_seconds','observed_mapping_seconds')), 'Finite time budgets required')
        require(type(b['cache_max_bytes']) is int and 0 <= b['cache_max_bytes'] <= 2*1024**3,
                'Cache allocation exceeds the 2 GiB contract')
        require(b['maximum_policy_calls'] == 2*b['horizon']*b['max_evaluations'],
                'Conservative solve count must include every backward/forward call and replay')
        require(b['case_seconds'] <= b['total_seconds'], 'Per-mapping budget exceeds total')
        require(0 < n['market_tolerance'] <= 2e-4 and 0 < n['fiscal_tolerance'] <= 1e-6 and
                0 <= n['final_reproduction_tolerance'] <= 1e-10, 'Acceptance gates cannot be loosened')
        require(set(n['terminal_tolerances']) == TERMINAL_KEYS and
                all(math.isfinite(v) and v > 0 for v in n['terminal_tolerances'].values()) and
                math.isfinite(n['raw_queue_relative_tolerance']) and
                n['raw_queue_relative_tolerance'] > 0, 'Complete terminal and raw-queue gates required')
        for key in ('initial_prices','initial_pensions'):
            require(len(n[key]) == b['horizon'] and all(math.isfinite(v) and v > 0 for v in n[key]),
                    'Explicit positive initial paths required')
        for key in ('price_bounds','pension_bounds'):
            require(len(n[key]) == 2 and 0 < n[key][0] < n[key][1] and all(map(math.isfinite,n[key])),
                    'Explicit finite root bounds required')
    return missing


def check_sources(pins):
    require(set(pins) == set(SOURCE_NAMES), 'Complete numerical source pins required')
    for name, expected in pins.items():
        require(sha(Path(__file__).parent/name) == expected, 'Numerical source changed: '+name)


def load_reference(output):
    require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID','').isdigit(), 'Torch Slurm only')
    require(sha(MANIFEST) == MANIFEST_SHA, 'Frozen reference manifest changed')
    m = read(MANIFEST)
    for key in ('contract','objective','source_manifest'):
        pinned(m[key])
    export = Path(m['local_export'])
    require(sha(export/'initial_state.pkl.gz') == m['checkpoint']['sha256'], 'Reference checkpoint changed')
    for name, digest in m['artifact_hashes'].items():
        require(sha(export/name) == digest, 'Reference artifact changed: '+name)
    sys.path[:0] = [str(ROOT/'code/model/tools'), str(ROOT/'tmp/e5f_overnight_local_20260927/portable/tools_v4')]
    import e5f_evening_calibration_runtime as runtime
    evaluator = runtime.setup(dict(read(m['contract']['path']), objective=m['objective']),
                              read(m['objective']['path']), output/'runtime')
    with gzip.open(export/'initial_state.pkl.gz','rb') as stream:
        packet = pickle.load(stream)
    require(serialized(vars(packet['parameters'])) == m['actual_serialized_parameters'],
            'All 260 saved parameter fields, including cached arrays, must match')
    sys.path.insert(0,str(Path(__file__).parent))
    write(output/'reference_identity.json', dict(status='PASS', label=LABEL,
        manifest_sha256=MANIFEST_SHA, checkpoint_sha256=m['checkpoint']['sha256'],
        source_manifest_sha256=m['source_manifest']['sha256'], parameter_fields=260,
        target_rows=len(m['full_target_table']), parameter_rows=len(m['full_parameter_table']),
        unchanged_standard_plots=len(m['standard_diagnostic_names'])))
    return m, packet, evaluator


def supply_rule(packet, pf, housing):
    original = packet['supply_rule']
    if housing == 'elastic_reference':
        return original
    require(housing == 'fixed_stock', 'Unknown housing closure')
    q = float(packet['solution'].p_eq[0])
    return pf.calendar.HousingSupplyRule('fixed-stock', q,
        float(original.quantity([q])[0]), 0.)


def load_endpoint(record, m, packet, evaluator, plan):
    import numpy as np
    from e5f_social_security import fiscal_accounts
    receipt = read(pinned(record))
    require(receipt['schema'] == 'block0506_preference_endpoint_v1' and
            receipt['reference_manifest_sha256'] == MANIFEST_SHA and
            receipt['source_manifest_sha256'] == m['source_manifest']['sha256'], 'Endpoint identity differs')
    require(receipt['housing'] == plan['housing'] and receipt['repeat_verified'] is True and
            receipt['native_one_step_verified'] is True, 'Endpoint closure/replay not verified')
    path = pinned(receipt['checkpoint'])
    with gzip.open(path,'rb') as stream:
        terminal = pickle.load(stream)
    before, after = serialized(vars(packet['parameters'])), serialized(vars(terminal['parameters']))
    allowed = GENERATED_FIELDS | {'psi_child','pension','native_inherited_distribution_evidence_dir'}
    require(set(before) == set(after), 'Endpoint parameter fields differ')
    differences = [k for k in before if k not in allowed and before[k] != after[k]]
    require(not differences, 'Endpoint changed frozen primitives: '+', '.join(differences))
    require(np.array_equal(packet['b_grid'],terminal['b_grid']), 'Endpoint wealth grid differs')
    P = terminal['parameters']
    require(P.psi_child == plan['shocks']['levels'][-1], 'Endpoint preference is not the last shock')
    require(P.psi_child == receipt['psi_child'] and P.pension == receipt['pension'], 'Endpoint levels differ')
    q = float(terminal['solution'].p_eq[0])
    require(q == receipt['price'], 'Endpoint price differs')
    ev = terminal['evaluation']
    pf = evaluator.rt['primitive'].pf
    g = terminal['stationary_g_pre']
    require(np.array_equal(g,ev.g_pre) and abs(float(g.sum())-1) <= 1e-9, 'Endpoint normalized state differs')
    E = float(g[:,:,:,0].sum())
    accounting = pf.transition.calendar_topcode_birth_accounting(ev.g_pre,ev.g_post_fertility,float(ev.births),P)
    renewal = float(accounting['topcode_adjusted_birth_children'])/(2.1*E)-1
    require(abs(renewal) <= 1e-6, 'Endpoint demographic renewal fails; do not renormalize psi')
    fiscal = fiscal_accounts(ev.g_current,P)
    require(abs(fiscal['scaled_pension_budget_residual']) <= 1e-6, 'Endpoint PAYGO fails')
    supply = supply_rule(packet,pf,plan['housing'])
    scale = float(supply.quantity([q])[0])/float(ev.demand_by_loc[0])
    require(math.isfinite(scale) and scale > 0 and
            abs(scale/receipt['population_scale']-1) <= 1e-8, 'Endpoint population scale inconsistent with housing')
    for name in ('target_fit','parameters'):
        with pinned(receipt[name]).open() as stream:
            rows = list(csv.DictReader(stream))
        require(len(rows) == (14 if name == 'target_fit' else 31), 'Complete endpoint tables required')
    return terminal, dict(receipt, population_scale=scale, renewal_relative_gap=renewal)


def mapping(packet, evaluator, terminal, endpoint, prices, pensions, psi_path, housing, output, cache_bytes,
            capture=False):
    import numpy as np
    import e5f_exact_policy_cache as cache_module
    import e5f_overnight_estate_audit as estate
    rt = evaluator.rt
    pf = rt['primitive'].pf
    P = copy.deepcopy(packet['parameters'])
    saved_output = P.native_inherited_distribution_evidence_dir
    P.native_inherited_distribution_evidence_dir = str(output/'inherited_state_evidence')
    write(output/'output_override.json', dict(field='native_inherited_distribution_evidence_dir',
        saved=saved_output, effective=P.native_inherited_distribution_evidence_dir, economic_change=False))
    grid, g = packet['b_grid'], packet['stationary_g_pre']
    state = pf.stationary_initial_state(g,float(g[:,:,:,0].sum()),float(packet['evaluation'].births),P,1/2.1)
    require(np.array_equal(g,state.g_pre), 'Inherited population must not be rescaled')
    rents = pf.rents_from_asset_prices(prices,endpoint['price'],P)
    audit_rows = []
    diagnostic_packets=[]
    def observer(t,ev,parameters,bg,shared,next_cohort):
        date = output/f'date_{t:03d}'
        date.mkdir(parents=True,exist_ok=True)
        budget = rt['primitive'].dated_budget(ev,parameters,shared,bg,float(rents[t]))
        purchase = rt['accounting'].audit_purchase_accounting(ev,parameters,shared,bg,rt['model'])
        funding = estate.audit(ev,parameters,bg,next_entrant_cohort=next_cohort)
        arrays = rt['audit'].policy_array_audit(dict(parameters=parameters,b_grid=bg,evaluation=ev),date)
        gates = dict(budget=budget['budget_excess_mass'] <= 2e-10,
            transaction=abs(purchase['maximum_occupied_transaction_wealth_error']) <= 1e-9,
            funded=funding['status']=='funded',negative_estates=funding['estate']['totals']['net_negative'] <= 1e-10,
            occupied_values=arrays['occupied_negative_steps']==0,
            probabilities=all(not v['nonfinite'] and 0 <= v['minimum'] <= v['maximum'] <= 1
                              for v in arrays['probabilities'].values()))
        for key,value in purchase.items():
            if key.endswith('violation_mass') or key in ('transaction_outside_grid_mass','saving_outside_grid_mass'):
                gates[key] = abs(float(value)) <= 2e-10
        audit_rows.append(dict(period=t,gates=gates,budget=budget,purchase=purchase,estate=funding,arrays=arrays))
        write(output/'dated_audits.json',audit_rows)
        write(output.parent/'heartbeat.json',dict(phase='forward',period=t,epoch=time.time()))
        require(all(gates.values()), 'Dated household/accounting gate failed')
        if capture and t in {0,len(prices)//2,len(prices)-1}:
            path=date/'diagnostic_packet.pkl.gz'
            dump_checkpoint(path,dict(parameters=parameters,b_grid=bg,evaluation=ev,shared=shared,
                                      dated_rent=float(rents[t]),period=t))
            diagnostic_packets.append(dict(path=str(path),sha256=sha(path),period=t,dated_rent=float(rents[t])))
    start=time.monotonic()
    with cache_module.policy_cache(pf,max_bytes=cache_bytes) as cache:
        result=pf.evaluate_path_at_prices(prices=prices,psi_path=psi_path,
            terminal_price=endpoint['price'],terminal_V=terminal['evaluation'].policy.V,
            base_parameters=P,b_grid=grid,initial_state=state,supply_rule=supply_rule(packet,pf,housing),
            birth_to_entry_conversion=1/2.1,transfer_path=np.zeros(len(prices)),
            pension_path=pensions,payroll_tax_path=np.full(len(prices),P.tau_pay),dated_observer=observer)
        cache_stats=cache.snapshot()
    for t,row in enumerate(result.rows):
        row['calendar_year']=2007+4*t
    gates=dict(mass=result.maximum_mass_accounting_error<=2e-8,
        policy_reproduction=result.maximum_policy_reproduction_error<=1e-10,
        projection=result.maximum_feasibility_projection_mass==0,
        dated_audits=len(audit_rows)==len(prices) and all(all(a['gates'].values()) for a in audit_rows))
    record=dict(gates=gates,rows=result.rows,cache=cache_stats,seconds=time.monotonic()-start,
        market_residual=[(r['housing_demand']-r['housing_supply'])/r['housing_supply'] for r in result.rows],
        fiscal_residual=[r['scaled_pension_budget_residual'] for r in result.rows],
        initial_mass=float(g.sum()),final_mass=float(result.terminal_state.g_pre.sum()),
        population_l1=float(np.abs(result.terminal_state.g_pre-g).sum()),
        backward_forward_policy_error=result.maximum_policy_reproduction_error,
        mass_error=result.maximum_mass_accounting_error,projection_mass=result.maximum_feasibility_projection_mass,
        diagnostic_packets=diagnostic_packets)
    write(output/'mapping.json',record)
    return result,record


def terminal_checks(packet,evaluator,terminal,endpoint,result,psi_path,numerics):
    import numpy as np
    pf=evaluator.rt['primitive'].pf
    P=terminal['parameters']; scale=endpoint['population_scale']
    g=terminal['stationary_g_pre']*scale
    E=float(g[:,:,:,0].sum())
    state=pf.stationary_initial_state(g,E,float(terminal['evaluation'].births)*scale,P,1/2.1)
    check=pf.terminal_convergence_diagnostics(evaluation=result,psi_path=psi_path,
        reference_state=state,reference_entry_flow=E,reference_price=endpoint['price'],
        reference_psi=P.psi_child,base_parameters=P,tolerances=numerics['terminal_tolerances'])
    raw_gap=float(np.max(np.abs(pf.birth_queue_values(result.terminal_state.scheduled_raw_entries)-
                               pf.birth_queue_values(state.scheduled_raw_entries))))/max(E,1e-15)
    check['raw_queue_maximum_relative_gap']=raw_gap
    check['raw_queue_pass']=raw_gap<=numerics['raw_queue_relative_tolerance']
    check['all_checks_pass']=check['all_checks_pass'] and check['raw_queue_pass']
    check['status']='passed' if check['all_checks_pass'] else 'not_converged'
    return check


def compare_horizons(short,long,early_periods,tolerance):
    require(long['horizon'] > short['horizon'] >= early_periods >= 1, 'Longer horizon and valid comparison window required')
    require(math.isfinite(tolerance) and tolerance > 0, 'Explicit horizon tolerance required')
    for key in ('reference_manifest_sha256','source_pins','housing','shock_contract'):
        require(short[key] == long[key], 'Horizon comparison changed '+key)
    require(short['root_and_terminal_pass'] and long['root_and_terminal_pass'], 'Both finite paths must first pass')
    gaps={}
    for key in ('asset_price','renter_price','adult_population','birth_children','housing_demand','pension_period_units'):
        a=[row[key] for row in short['rows'][:early_periods]]
        b=[row[key] for row in long['rows'][:early_periods]]
        gaps[key]=max(abs(x-y)/max(abs(x),abs(y),1e-12) for x,y in zip(a,b))
    return dict(passed=all(v<=tolerance for v in gaps.values()),relative_gaps=gaps,
                early_periods=early_periods,tolerance=tolerance)


def render_diagnostics(records,output,audit_module,standard_names):
    receipts=[]
    for record in records:
        with gzip.open(pinned(record),'rb') as stream:
            packet=pickle.load(stream)
        directory=output/f"date_{record['period']:03d}"
        directory.mkdir(parents=True,exist_ok=False)
        audit_module.standard_diagnostics(packet,directory,validate_production_young=False,
                                          dated_rent=record['dated_rent'])
        plots=sorted((directory/'standard_diagnostics').glob('*.png'))
        require({p.name for p in plots}==set(standard_names), 'The standard 17-plot set changed')
        receipts.append(dict(period=record['period'],plots={p.name:sha(p) for p in plots}))
        del packet
    return dict(sampled_dates=receipts,visual_review_pending=True)


def measure_jacobian(evaluate,prices,pensions,perturbed_date,step,output,identity):
    """Five budgeted baseline mappings; no shocked preferences or old derivatives."""
    import numpy as np
    from e5f_four_shock_acceleration import assemble_jacobian
    require(math.isfinite(step) and 0 < step < .01,'Explicit small positive log perturbation required')
    require(len(prices)==len(pensions) and 0 <= perturbed_date < len(prices),'Invalid derivative horizon/date')
    def residual(p,b):
        reply=evaluate(p,b)
        require(reply['mapping_valid'],'A derivative mapping failed its scientific gates')
        return np.r_[reply['market_residual'],reply['fiscal_residual']]
    baseline=residual(np.array(prices,float),np.array(pensions,float))
    columns=[]
    for block in (0,1):
        values=[]
        for sign in (1,-1):
            p,b=np.array(prices,float),np.array(pensions,float)
            (p if block==0 else b)[perturbed_date]*=math.exp(sign*step)
            values.append(residual(p,b))
        columns.append((values[0]-values[1])/(2*step))
    matrix,receipt=assemble_jacobian(columns,len(prices),perturbed_date)
    output.mkdir(parents=True,exist_ok=False)
    np.save(output/'jacobian.npy',matrix,allow_pickle=False)
    write(output/'receipt.json',dict(identity,**receipt,matrix=dict(path=str(output/'jacobian.npy'),
        sha256=sha(output/'jacobian.npy')),residual_units='physical_unscaled',
        coordinate_order=['log_house_price','log_period_pension'],mapping_count=5,
        perturbation_log_step=step,baseline_residual=baseline))
    return matrix


def run(plan,output):
    started=time.monotonic()
    validate_plan(plan,launching=True)
    require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm only')
    import numpy as np
    check_sources(plan['source_pins'])
    readiness=read(pinned(plan['readiness_receipt']))
    require(readiness['schema']=='block0506_pf_preparation_v1' and
            readiness['source_pins']==plan['source_pins'] and readiness['tests_passed'] is True and
            readiness['native_smoke_passed'] is True, 'Matching test/no-shock readiness receipt required')
    output.mkdir(parents=True,exist_ok=False)
    write(output/'plan.json',plan)
    m,packet,evaluator=load_reference(output/'reference')
    terminal,endpoint=load_endpoint(plan['endpoint'],m,packet,evaluator,plan)
    from e5f_four_shock_acceleration import solve_joint_with_acceleration
    from e5f_social_security_root import CandidateDomainError
    b,n=plan['budget'],plan['numerics']; path=shock_path(plan['shocks'],b['horizon'])
    deadline=started+b['total_seconds']; latest={}; count=0
    def evaluate(prices,pensions):
        nonlocal count,latest
        count+=1
        require(count<=b['max_evaluations'], 'Mapping budget exhausted')
        pf=evaluator.rt['primitive'].pf
        try: pf.rents_from_asset_prices(prices,endpoint['price'],packet['parameters'])
        except ValueError as exc: raise CandidateDomainError(str(exc)) from exc
        remaining=min(b['case_seconds'],deadline-time.monotonic())
        require(remaining>0,'Transition deadline reached')
        old=signal.signal(signal.SIGALRM,lambda *_: (_ for _ in ()).throw(TimeoutError('Mapping deadline')))
        signal.setitimer(signal.ITIMER_REAL,remaining)
        try:
            result,record=mapping(packet,evaluator,terminal,endpoint,prices,pensions,path,plan['housing'],
                                  output/f'mapping_{count:04d}',b['cache_max_bytes'],capture=True)
        finally:
            signal.setitimer(signal.ITIMER_REAL,0); signal.signal(signal.SIGALRM,old)
        latest=dict(result=result,record=record,evaluation=count)
        write(output/'latest_completed.json',dict(evaluation=count,**record))
        dump_checkpoint(output/'latest_mapping.pkl.gz',dict(terminal_state=result.terminal_state,
            prices=prices,pensions=pensions,psi_path=path,rows=result.rows,evaluation=count))
        return dict(mapping_valid=all(record['gates'].values()),market_residual=record['market_residual'],
                    fiscal_residual=record['fiscal_residual'],payload=dict(evaluation=count))
    def progress(record):
        write(output/'root_progress.json',record)
        if record.get('new_best'):
            write(output/'best_so_far.json',dict(root=record,mapping=latest['record']))
            shutil.copyfile(output/'latest_mapping.pkl.gz',output/'best_mapping.pkl.gz')
    jacobian=None
    if plan['initial_jacobian']:
        jac=read(pinned(plan['initial_jacobian']))
        require(jac['reference_manifest_sha256']==MANIFEST_SHA and jac['housing']==plan['housing'] and
                jac['source_pins']==plan['source_pins'] and jac['horizon']==b['horizon'], 'Jacobian provenance differs')
        require(jac['residual_units']=='physical_unscaled' and jac['coordinate_order']==['log_house_price','log_period_pension'],
                'Old scaled/three-block Jacobian is not compatible')
        jacobian=np.load(pinned(jac['matrix']),allow_pickle=False)
    root=solve_joint_with_acceleration(closure='fixed_tax',initial_prices=n['initial_prices'],
        initial_fiscal_values=n['initial_pensions'],evaluate=evaluate,
        project_prices=lambda x:np.clip(x,*n['price_bounds']),fiscal_bounds=n['pension_bounds'],
        market_tolerance=n['market_tolerance'],fiscal_tolerance=n['fiscal_tolerance'],
        market_slope=n['market_slope'],fiscal_slope=n['fiscal_slope'],max_log_step=n['max_log_step'],
        damping=n['damping'],max_evaluations=b['max_evaluations'],deadline_monotonic=deadline,
        max_condition_number=n['max_condition_number'],worsening_factor=n['worsening_factor'],
        final_reproduction_tolerance=n['final_reproduction_tolerance'],callback=progress,initial_jacobian=jacobian)
    write(output/'root_result.json',root)
    final=root.get('final'); terminal_gate=None
    if final and latest and final.get('payload',{}).get('evaluation')==latest['evaluation']:
        terminal_gate=terminal_checks(packet,evaluator,terminal,endpoint,latest['result'],path,n)
    diagnostics=None
    if root['converged']:
        require(final and latest and final.get('payload',{}).get('evaluation')==latest['evaluation'],
                'Final replay does not own the retained mapping')
        diagnostics=render_diagnostics(latest['record']['diagnostic_packets'],output/'dated_diagnostics',
                                       evaluator.rt['audit'],m['standard_diagnostic_names'])
        rows=latest['record']['rows']
        with (output/'path.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    for name in ('target_fit','parameters'):
        shutil.copyfile(Path(m['local_export'])/(name+'.csv'),output/('reference_'+name+'.csv'))
        endpoint_record=read(pinned(plan['endpoint']))
        shutil.copyfile(pinned(endpoint_record[name]),output/('endpoint_'+name+'.csv'))
    check_sources(plan['source_pins'])
    write(output/'complete.json',dict(reference_label=LABEL,reference_manifest_sha256=MANIFEST_SHA,
        source_pins=plan['source_pins'],housing=plan['housing'],shock_contract=plan['shocks'],horizon=b['horizon'],
        root_and_terminal_pass=bool(root['converged'] and terminal_gate and terminal_gate['all_checks_pass']),
        terminal=terminal_gate,rows=latest.get('record',{}).get('rows',[]),horizon_extension_verified=False,
        status='finite_path_only_horizon_and_visual_review_still_required',diagnostics=diagnostics,
        reference_fit=m['full_target_table'],reference_parameters=m['full_parameter_table']))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--plan-sha256')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    plan=read(args.plan)
    if not args.execute:
        print(json.dumps(dict(execution_enabled=False,unresolved=validate_plan(plan),
                              note='No model imported and no transition launched')))
        return
    require(args.plan_sha256 and sha(args.plan)==args.plan_sha256,'Exact plan pin required')
    require(args.output is not None,'New output directory required')
    run(plan,args.output)


if __name__=='__main__':
    main()
