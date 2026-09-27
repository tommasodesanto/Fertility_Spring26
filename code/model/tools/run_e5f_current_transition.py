#!/usr/bin/env python3
"""Pinned bounded native fixed-tax price/pension transition root.

No economic defaults: an explicitly approved contract and measured numerical
budget are required. Root convergence is separate from terminal convergence.
"""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'
os.environ['MPLBACKEND'] = 'Agg'
import argparse
import copy
import csv
import gzip
import hashlib
import json
from pathlib import Path
import pickle
import shutil
import signal
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[3]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''): h.update(block)
    return h.hexdigest()


def plain(value):
    if isinstance(value, np.ndarray): return plain(value.tolist())
    if isinstance(value, np.generic): return plain(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    if isinstance(value, dict): return {str(k):plain(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [plain(v) for v in value]
    return value


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(plain(value), indent=2, sort_keys=True, allow_nan=False)+'\n')
    temp.replace(path)


def read(path): return json.loads(Path(path).read_text())


def pinned(record):
    path = Path(record['path'])
    if sha(path) != record['sha256']: raise ValueError('Pinned artifact differs: '+str(path))
    return path


def validate(plan, *, horizon, max_evaluations, deadline, case_seconds, now=None):
    now = time.time() if now is None else now
    if (type(horizon) is not int or horizon < 2 or type(max_evaluations) is not int
            or max_evaluations < 2 or not np.isfinite([deadline,case_seconds]).all()
            or deadline <= now or case_seconds <= 0):
        raise ValueError('Explicit positive horizon/evaluation/absolute time budgets required')
    budget = plan['budget']
    if (horizon != budget['horizon'] or max_evaluations != budget['max_evaluations']
            or deadline != budget['deadline_epoch'] or case_seconds != budget['case_seconds']):
        raise ValueError('CLI budgets must match pinned plan exactly')
    if budget.get('cache_max_bytes') != 2*1024**3:
        raise ValueError('Explicit exact-policy cache budget must be 2 GiB')
    if not budget.get('observed_mapping_seconds') or not budget.get('solve_count_estimate'):
        raise ValueError('Measured mapping cost and total solve-count estimate required')
    if plan['economics'] != {'credit':'experimental_natural_solvency', 'fiscal':'fixed_tax',
            'property_rebates':False, 'outside_entry':0., 'retention':1.,
            'psi':'fixed_reference', 'entry_clock':'split_birth_vintage',
            'estate':'provisional_net_estates_fund_actual_next_cohort_sink_residual'}:
        raise ValueError('Explicit current economic closure differs')
    required = {'preferences','earnings','entry','population','fiscal','geography','estate','credit','targets'}
    if set(plan['classification']) != required:
        raise ValueError('Every closure object requires explicit classification')
    valid = {'estimated','empirically_normalized','externally_fixed','outstanding','experimental','provisional'}
    if any(row['status'] not in valid or not row['description'] for row in plan['classification'].values()):
        raise ValueError('Invalid closure classification')
    n = plan['numerics']
    required_numerics = {'market_tolerance','fiscal_tolerance','final_reproduction_tolerance',
        'initial_prices','initial_pensions','price_bounds','pension_bounds','market_slope',
        'fiscal_slope','max_log_step','damping','max_condition_number','worsening_factor','terminal_tolerances'}
    if set(n) != required_numerics: raise ValueError('Complete explicit numerical configuration required')
    if not n['terminal_tolerances']: raise ValueError('Explicit terminal tolerances required')
    if (not 0 < n['market_tolerance'] <= 2e-4 or not 0 < n['fiscal_tolerance'] <= 1e-6
            or not 0 <= n['final_reproduction_tolerance'] <= 1e-10):
        raise ValueError('Production gates cannot be relaxed')
    for name in ('initial_prices','initial_pensions'):
        a = np.asarray(n[name],float)
        if a.shape != (horizon,) or not np.isfinite(a).all() or np.any(a <= 0):
            raise ValueError('Explicit positive starting paths required')
    for name in ('price_bounds','pension_bounds'):
        bounds = np.asarray(n[name],float)
        if bounds.shape != (2,) or not np.isfinite(bounds).all() or not 0 < bounds[0] < bounds[1]:
            raise ValueError('Explicit numerical bounds required')
    if not plan.get('source_pins'): raise ValueError('Current-source pins required')
    required_sources = [Path(__file__), Path(__file__).with_name('e5f_current_transition_runtime.py'),
        Path(__file__).with_name('e5f_social_security_root.py'),
        Path(__file__).with_name('e5f_matched_pf_path_root.py'),
        Path(__file__).with_name('e5f_exact_policy_cache.py'),
        Path(__file__).with_name('run_e5f_independent_numerical_audit.py')]
    if any(str(path.resolve()) not in plan['source_pins'] for path in required_sources):
        raise ValueError('Root and runtime source pins required')
    for path,expected in plan['source_pins'].items():
        if sha(path) != expected: raise ValueError('Current source changed: '+path)
    approval = read(pinned(plan['approval']))
    if (approval.get('approved') is not True or approval.get('scope') != 'native_current_transition'
            or approval.get('source_pins') != plan['source_pins']
            or approval.get('reference_checkpoint_sha256') != plan['reference_checkpoint_sha256']):
        raise ValueError('Explicit current native/reference approval missing or mismatched')
    replay = read(pinned(approval['reference_replay']))
    if replay.get('status') not in ('native_reference_replay_pass','native_reference_tables_pass_array_review_required'):
        raise ValueError('Native reference replay has not passed')
    if replay.get('status') == 'native_reference_tables_pass_array_review_required' and approval.get('arrays_reviewed') is not True:
        raise ValueError('Explicit lead array review required')
    smoke = read(pinned(approval['smoke_receipt']))
    if not str(smoke.get('status','')).startswith('passed'):
        raise ValueError('Native exact-loop smoke has not passed')
    endpoint = read(pinned(plan['endpoint_complete']))
    if endpoint.get('usable_closed_root') is not True or endpoint.get('repeat_verified') is not True:
        raise ValueError('Authenticated repeated closed endpoint required')
    receipt = read(pinned(plan['endpoint_receipt']))
    checkpoint = pinned(plan['endpoint_checkpoint'])
    if receipt['case_checkpoint_sha256'] != sha(checkpoint):
        raise ValueError('Endpoint checkpoint does not match receipt')
    terminal_approval = read(pinned(plan['terminal_approval']))
    if (terminal_approval.get('approved') is not True
            or terminal_approval.get('status') != 'native_terminal_replay_verified'
            or terminal_approval.get('terminal_checkpoint_sha256') != sha(checkpoint)):
        raise ValueError('Native terminal replay approval missing or mismatched')
    terminal_receipt = Path(terminal_approval['terminal_receipt_path'])
    if terminal_receipt.resolve() != Path(plan['endpoint_receipt']['path']).resolve() or sha(terminal_receipt) != terminal_approval['terminal_receipt_sha256']:
        raise ValueError('Approved native terminal receipt differs')
    terminal_pins = terminal_approval.get('current_source_files',{})
    core = {p:v for p,v in plan['source_pins'].items() if '/intergen_eqscale_seq_optimized/' in p}
    if not core or any(terminal_pins.get(p) != v for p,v in core.items()):
        raise ValueError('Approved native terminal and root model source differ')
    return approval, endpoint, receipt


def dump_checkpoint(path, payload):
    temp = Path(str(path)+'.tmp')
    with gzip.open(temp,'wb',compresslevel=1) as stream: pickle.dump(payload,stream,protocol=5)
    temp.replace(path)


def observed_mapping(prepared, terminal, endpoint, prices, pensions, output, index, *, cache_max_bytes):
    import e5f_exact_policy_cache as exact_cache
    import e5f_overnight_estate_audit as estate
    import e5f_solvency_credit_benchmark as credit_audit
    from run_e5f_current_transition_smoke import audit_gates
    rt = prepared['runtime']; pf = rt['primitive'].pf
    packet = prepared['selected']; P = copy.deepcopy(prepared['parameters'])
    P.native_solvency_credit = True
    if P.property_tax_lump_sum_transfer != 0: raise ValueError('Reference rebate must be zero')
    grid = packet['b_grid']; H = len(prices); terminal_price = float(endpoint['endpoint']['price'])
    rents = pf.rents_from_asset_prices(prices,terminal_price,P)
    initial = pf.stationary_initial_state(packet['stationary_g_pre'],float(packet['stationary_g_pre'][:,:,:,0].sum()),
        packet['evaluation'].births,P,1/2.1)
    np.testing.assert_array_equal(initial.g_pre,packet['stationary_g_pre'])
    audits = []; packet_paths = []
    sampled_dates = sorted({0,H//2,H-1})
    def observer(period,e,parameters,b_grid,shared,next_cohort):
        budget = rt['primitive'].dated_budget(e,parameters,shared,b_grid,float(rents[period]))
        purchase = credit_audit.audit_purchase_accounting(e,parameters,shared,b_grid,rt['model'])
        try: funding = estate.audit(e,parameters,b_grid,next_entrant_cohort=next_cohort)
        except estate.EstateFundingShortfall as exc: funding = exc.audit
        audit = dict(budget=budget,purchase=purchase,estate=funding)
        audit['gates'] = audit_gates(audit)
        for name,values in rt['primitive'].policy_arrays(e.policy).items():
            if 'prob' in name:
                audit['gates']['probability_'+name] = bool(np.isfinite(values).all() and values.min() >= -1e-12 and values.max() <= 1+1e-12)
        audit['gates']['occupied_value_monotonicity'] = not bool(np.any(
            (e.g_pre[:-1] > 1e-12) & (np.diff(e.policy.V,axis=0) < -1e-7)))
        audit['gates']['finite_values'] = bool(np.isfinite(e.policy.V).all())
        audits.append(audit)
        if period in sampled_dates:
            packet_path=output/'candidate_packets'/f'evaluation_{index:04d}'/f'date_{period:04d}.pkl.gz'
            packet_path.parent.mkdir(parents=True,exist_ok=True)
            # Serialize this observed state now; retain only paths, not another population copy.
            dump_checkpoint(packet_path,dict(parameters=parameters,b_grid=b_grid,evaluation=e,
                shared=shared,dated_rent=float(rents[period]),period=period,
                years_from_announcement=period*parameters.period_years))
            packet_paths.append(dict(period=period,path=str(packet_path),sha256=sha(packet_path),
                                     dated_rent=float(rents[period])))
        write(output/'heartbeat.json',dict(epoch=time.time(),phase='mapping',evaluation=index,forward_dates=period+1,horizon=H))
    with exact_cache.policy_cache(pf,max_bytes=cache_max_bytes) as cache:
        result = pf.evaluate_path_at_prices(prices=prices,psi_path=np.full(H,P.psi_child),
            terminal_price=terminal_price,terminal_V=terminal['evaluation'].policy.V,
            base_parameters=P,b_grid=grid,initial_state=initial,supply_rule=packet['supply_rule'],
            birth_to_entry_conversion=1/2.1,transfer_path=np.zeros(H),
            pension_path=pensions,payroll_tax_path=np.full(H,P.tau_pay),dated_observer=observer)
        cache_stats=cache.snapshot()
    gates = dict(mass=result.maximum_mass_accounting_error <= 2e-8,
        policy_reproduction=result.maximum_policy_reproduction_error <= 1e-10,
        feasibility=(result.maximum_feasibility_projection_mass == 0.0
                     if bool(getattr(P,'native_exact_inherited_distribution',False))
                     else result.maximum_feasibility_projection_mass <= 1e-6),
        dated_audits=len(audits)==H and all(all(a['gates'].values()) for a in audits))
    # Economic measurements retain symbolic t=0; native legacy calendar labels do not define this experiment.
    for t,row in enumerate(result.rows):
        row.pop('calendar_year',None); row['period']=t; row['years_from_announcement']=t*P.period_years
    market = np.array([(r['housing_demand']-r['housing_supply'])/r['housing_supply'] for r in result.rows])
    fiscal = np.array([r['scaled_pension_budget_residual'] for r in result.rows])
    return result,audits,gates,market,fiscal,packet_paths,cache_stats


def render_sampled_diagnostics(records, output, audit_module):
    receipts=[]
    for record in records:
        path=pinned(record)
        with gzip.open(path,'rb') as stream: packet=pickle.load(stream)
        destination=Path(output)/f"date_{record['period']:04d}"
        destination.mkdir(parents=True,exist_ok=False)
        audit_module.standard_diagnostics(packet,destination,validate_production_young=False,
                                           dated_rent=float(record['dated_rent']))
        plots=sorted((destination/'standard_diagnostics').glob('*.png'))
        if len(plots)!=17: raise RuntimeError('Standard diagnostic packet must contain exactly 17 PNGs')
        receipts.append(dict(period=record['period'],dated_rent=record['dated_rent'],
            plots={str(p):sha(p) for p in plots},source_packet=record,
            units='actual dated population retained; quantity plot unit recorded by standard reporter'))
        del packet
    return dict(status='generated_17_standard_plots_per_sampled_date',sampled_dates=receipts,
                visual_review_pending=True,horizon_certification='separate terminal and horizon-extension gates')


def terminal_check(pf, result, terminal, endpoint, P, tolerances):
    scale = float(endpoint['endpoint']['population_scale'])
    state = pf.stationary_initial_state(terminal['stationary_g_pre']*scale,
        terminal['solution'].entry_rate*scale,terminal['evaluation'].births*scale,P,1/2.1)
    return pf.terminal_convergence_diagnostics(evaluation=result,
        psi_path=np.full(len(result.prices),P.psi_child),reference_state=state,
        reference_entry_flow=terminal['solution'].entry_rate*scale,
        reference_price=float(endpoint['endpoint']['price']),reference_psi=P.psi_child,
        base_parameters=P,tolerances=tolerances)


def run(args):
    if sha(args.plan) != args.plan_sha256: raise ValueError('Plan SHA differs')
    plan = read(args.plan)
    _,endpoint,endpoint_receipt = validate(plan,horizon=args.horizon,max_evaluations=args.max_evaluations,
        deadline=args.deadline,case_seconds=args.case_seconds)
    args.output.mkdir(parents=True,exist_ok=False); write(args.output/'plan.json',plan)
    import e5f_current_transition_runtime as native
    from e5f_social_security_root import solve_social_security_path, CandidateDomainError
    prepared = native.setup(args.output/'preparation')
    for relative,expected in prepared['preparation']['current_source_files'].items():
        if plan['source_pins'].get(str(ROOT/relative)) != expected:
            raise ValueError('Unapproved loaded current source: '+relative)
    if prepared['preparation']['reference_checkpoint_sha256'] != plan['reference_checkpoint_sha256']:
        raise ValueError('Reference snapshot differs from author contract')
    with gzip.open(plan['endpoint_checkpoint']['path'],'rb') as stream: terminal = pickle.load(stream)
    P = prepared['parameters']; rt = prepared['runtime']; pf = rt['primitive'].pf
    if not np.array_equal(terminal['b_grid'],prepared['selected']['b_grid']): raise ValueError('Endpoint grid differs')
    for name in ('psi_child','tau_pay','tau_H','q','delta','period_years','H0','r_bar','xi_supply',
                 'property_tax_lump_sum_transfer','Pi_z','z_grid','survival_probs','fixed_reference_entry_conditional'):
        if not np.array_equal(np.asarray(getattr(terminal['parameters'],name)),np.asarray(getattr(P,name))):
            raise ValueError('Endpoint primitive differs: '+name)
    tables=[]
    for case in (prepared['reference_path'],Path(plan['endpoint_checkpoint']['path']).parent):
        with (case/'parameters.csv').open() as stream:
            tables.append({r['parameter']:float(r['estimate']) for r in csv.DictReader(stream)})
    if len(tables[0]) != 31 or tables[0] != tables[1]:
        raise ValueError('All 31 endpoint/reference parameter rows must agree')
    n = plan['numerics']; count = 0; latest = {}
    write(args.output/'diagnostic_status.json',dict(standard_17_plots='pending final mapping; sampled first, middle and last dates',
        complete_target_and_parameter_tables='reference and terminal pinned; dated paths exported',
        horizon_certification='pending'))
    def alarm(*_): raise TimeoutError('Mapping or global transition deadline reached')
    def evaluate(prices,pensions):
        nonlocal count,latest
        count += 1
        remaining = min(args.case_seconds,args.deadline-time.time())
        if remaining <= 0: raise TimeoutError('Absolute transition deadline reached')
        try: pf.rents_from_asset_prices(prices,float(endpoint['endpoint']['price']),P)
        except ValueError as exc: raise CandidateDomainError(str(exc)) from exc
        old = signal.signal(signal.SIGALRM,alarm); signal.setitimer(signal.ITIMER_REAL,remaining)
        try:
            result,audits,gates,market,fiscal,packet_paths,cache_stats = observed_mapping(prepared,terminal,endpoint,prices,pensions,args.output,count,cache_max_bytes=plan['budget']['cache_max_bytes'])
        finally:
            signal.setitimer(signal.ITIMER_REAL,0); signal.signal(signal.SIGALRM,old)
            # Restore the outer absolute watchdog without extending its deadline.
            if old is not signal.SIG_DFL:
                signal.setitimer(signal.ITIMER_REAL,max(1e-6,args.deadline-time.time()))
        latest = dict(result=result,audits=audits,gates=gates,prices=prices,pensions=pensions,
                      evaluation=count,plan_sha256=args.plan_sha256,packet_paths=packet_paths,cache=cache_stats)
        dump_checkpoint(args.output/'latest_mapping.pkl.gz',latest)
        write(args.output/f'mapping_{count:04d}.json',dict(evaluation=count,gates=gates,audits=audits,market_residual=market,fiscal_residual=fiscal,packet_paths=packet_paths,cache=cache_stats))
        return dict(mapping_valid=all(gates.values()),market_residual=market,fiscal_residual=fiscal,
                    payload={'evaluation':count})
    def progress(record):
        write(args.output/'latest_completed.json',record)
        if record.get('new_best'):
            write(args.output/'best_so_far.json',record)
            shutil.copyfile(args.output/'latest_mapping.pkl.gz',args.output/'best_mapping.pkl.gz')
        write(args.output/'heartbeat.json',dict(epoch=time.time(),phase='root',evaluation=count,event=record.get('event',record.get('phase'))))
    jacobian = None
    if plan.get('initial_jacobian'):
        jacobian = np.load(pinned(plan['initial_jacobian']),allow_pickle=False)
    started = time.monotonic()
    try:
        result = solve_social_security_path(closure='fixed_tax',initial_prices=n['initial_prices'],
            initial_fiscal_values=n['initial_pensions'],evaluate=evaluate,
            project_prices=lambda p:np.clip(p,*n['price_bounds']),fiscal_bounds=n['pension_bounds'],
            market_tolerance=n['market_tolerance'],fiscal_tolerance=n['fiscal_tolerance'],
            market_slope=n['market_slope'],fiscal_slope=n['fiscal_slope'],max_log_step=n['max_log_step'],
            damping=n['damping'],max_evaluations=args.max_evaluations,
            deadline_monotonic=time.monotonic()+args.deadline-time.time(),
            max_condition_number=n['max_condition_number'],worsening_factor=n['worsening_factor'],
            final_reproduction_tolerance=n['final_reproduction_tolerance'],callback=progress,initial_jacobian=jacobian)
        write(args.output/'root_result.json',result)
        terminal_gate = None; plot_receipt = None
        final = result.get('final')
        if final is not None and latest and final.get('payload',{}).get('evaluation') == latest['evaluation']:
            terminal_gate = terminal_check(pf,latest['result'],terminal,endpoint,P,n['terminal_tolerances'])
            dump_checkpoint(args.output/'final_mapping.pkl.gz',latest)
            rows = latest['result'].rows
            with (args.output/'path.csv').open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        if result['converged'] and latest and final.get('payload',{}).get('evaluation') == latest['evaluation']:
            plot_receipt=render_sampled_diagnostics(latest['packet_paths'],args.output/'dated_diagnostics',rt['audit'])
            write(args.output/'diagnostic_status.json',plot_receipt)
        for path,expected in plan['source_pins'].items():
            if sha(path)!=expected: raise ValueError('Source changed during transition: '+path)
        write(args.output/'complete.json',dict(status='root_and_terminal_pass' if result['converged'] and terminal_gate and terminal_gate['all_checks_pass'] else 'incomplete_transition',
            root_converged=result['converged'],terminal_convergence=terminal_gate,
            evaluations=result['evaluations'],seconds=time.monotonic()-started,
            standard_17_plots_pending=plot_receipt is None,dated_diagnostics=plot_receipt,horizon_extension_verified=False,
            economic_classification=plan['classification'],no_new_calibration=True,
            limitation='Credit is experimental; estate settlement provisional; publication certification requires standard plots and horizon comparison'))
    except Exception as exc:
        write(args.output/'failure.json',dict(error_type=type(exc).__name__,error=str(exc),evaluations=count,epoch=time.time(),classification=getattr(exc,'classification','unexpected_error'),inherited_state_evidence=getattr(exc,'audit',None)))
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True);p.add_argument('--plan-sha256',required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--horizon',type=int,required=True)
    p.add_argument('--max-evaluations',type=int,required=True);p.add_argument('--deadline',type=float,required=True)
    p.add_argument('--case-seconds',type=float,required=True)
    args=p.parse_args()
    if not np.isfinite(args.deadline) or args.deadline <= time.time():
        raise ValueError('Future absolute deadline required')
    def global_timeout(*_): raise TimeoutError('Absolute transition budget exhausted')
    old=signal.signal(signal.SIGALRM,global_timeout)
    signal.setitimer(signal.ITIMER_REAL,args.deadline-time.time())
    try: run(args)
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)

if __name__=='__main__': main()
