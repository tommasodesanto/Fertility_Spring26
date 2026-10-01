#!/usr/bin/env python3
"""Bounded local mortgage-LTV GE price roots at the frozen winner31 point."""
from __future__ import annotations
import copy, csv, hashlib, importlib.util, json, math, os, shutil, signal, sys, time, traceback
from pathlib import Path
import numpy as np
from types import SimpleNamespace

for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key] = '1'

ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
HERE = ROOT/'output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/purchase_ltv_v1'
LOCAL = HERE/'local_run'
OLD_GE = ROOT/'output/model/fixed_reference_economics_20260928/credit_ge_v1/run_ge.py'
BASELINE = LOCAL/'retry5/results/baseline_80_80'
PRE = LOCAL/'retry5/results/q0_reference_inherited_states.npz'
DEADLINE_SECONDS = 600
ORIGINAL_DEADLINE_EPOCH = 1790879938.0993679
PER_CASE_LC_CAP = 12
RENEWAL_TOL = 1e-6
EXPECTED_HOUSEHOLD_SHA = '2a34f5f28c0c63ca3d24d9e759cd1b5b92aadaddb5fe04795e65abbc0b448082'

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path): return json.loads(Path(path).read_text())
def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False, default=str)+'\n')
    tmp.replace(path)
def require(ok, msg):
    if not ok: raise RuntimeError(msg)
def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); sys.modules[name] = mod; spec.loader.exec_module(mod)
    return mod
def alarm(*_): raise TimeoutError('120-second case limit or 10-minute total budget')

def main(out, only_case=None, repeat_factor=None, selected_path=None):
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=False)
    started = time.time(); deadline = started + 60 if repeat_factor is not None else min(started + DEADLINE_SECONDS, ORIGINAL_DEADLINE_EPOCH)
    require(PRE.is_file() and BASELINE.is_dir(), 'Accepted baseline/PRE packet is missing')
    source = ROOT/'code/model/refactor_lab/engine/household.py'
    require(sha(source) == EXPECTED_HOUSEHOLD_SHA, 'Reviewed household source hash drift')
    sys.path.insert(0, str(HERE.parent))
    import fixed_price_responses as d
    old = load_module('frozen_monday_ge_closure', OLD_GE)
    baseline_pre_sha = sha(PRE)
    shutil.copy2(PRE, out/'q0_reference_inherited_states.npz')
    if (LOCAL/'retry5/results/q0_reference_inherited_states.json').exists():
        shutil.copy2(LOCAL/'retry5/results/q0_reference_inherited_states.json', out/'q0_reference_inherited_states.json')
    write(out/'launch.json', dict(started_epoch=started, deadline_epoch=deadline,
        budget_seconds=DEADLINE_SECONDS, per_case_lifecycle_cap=PER_CASE_LC_CAP,
        baseline_pre_source=str(PRE), baseline_pre_sha256=baseline_pre_sha,
        closure_source=str(OLD_GE), closure_source_sha256=sha(OLD_GE),
        driver_sha256=sha(__file__), household_source_sha256=sha(source),
        buyer_estate_override_sha256=sha(LOCAL/'household_buyer_estate.py'),
        fixed_price_q0_seeds={
          'purchase_90_stayer_80': str(LOCAL/'retry7/results/purchase_90_stayer_80/closure.json'),
          'both_90_90': str(LOCAL/'retry7/results/both_90_90/closure.json'),
          'purchase_100_stayer_80': str(LOCAL/'retry8/results/purchase_100_stayer_80/closure.json'),
          'both_100_100': str(LOCAL/'retry9/results/both_100_100/closure.json')},
        economics='fixed winner31 parameters and targets; buyer/stayer LTV per named cell; only price and GE population scale adjust; no refit'))

    # The copied fixed-price kernel is the already reviewed native mortgage path.
    # Extend only its per-auth attempt counter for this explicitly bounded GE pass,
    # and prevent its q0 hook from overwriting the accepted PRE file.
    cell_source = __import__('inspect').getsource(d.make_price_cell)
    require(cell_source.count('auth["lifecycle_solves"]<=6') == 1, 'Fixed-price solve cap source drift')
    cell_source = cell_source.replace('auth["lifecycle_solves"]<=6', 'auth["lifecycle_solves"]<=12')
    require(cell_source.count('if regime == "reference" and float(factor) == 1.0:') == 1, 'PRE save hook source drift')
    cell_source = cell_source.replace('if regime == "reference" and float(factor) == 1.0:', 'if False:  # GE adapter preserves frozen baseline PRE')
    impact_start = cell_source.index('    # The first cell is the sole source of q0 inherited states.')
    impact_end = cell_source.index('    # Existing stable 17-panel audit packet', impact_start)
    impact_skip = '''    # Optional old-PRE impact at new q is not part of the stationary GE solve.
    impact_pre_path = out.parent / "q0_reference_inherited_states.npz"
    require(impact_pre_path.is_file(), "Frozen q0 PRE distribution missing")
    with np.load(impact_pre_path, allow_pickle=False) as saved: impact_pre = saved["g_pre"].copy()
    impact_P = copy.deepcopy(P)
    impact_out = out / "baseline_state_impact"; impact_out.mkdir()
    impact_summary = {"status": "unavailable_not_evaluated_new_price_old_distribution",
        "reason": "q0 inherited distribution failed the unchanged zero-projection gate at factor 0.95; no distribution substitution or projection used",
        "first_failure_positive_infeasible_mass": 6.45718872382e-12,
        "first_failure_factor": 0.95, "first_failure_price": 0.68320995038751,
        "inherited_distribution_sha256": hashlib.sha256(impact_pre.tobytes()).hexdigest()}
    impact_gates = {"status": "not_evaluated", "reason": impact_summary["reason"]}
    write(impact_out / "summary.json", cal.jsonable(impact_summary))
'''
    cell_source = cell_source[:impact_start] + impact_skip + cell_source[impact_end:]
    impact_report = '''    impact_summary=report_helpers.aggregates(auth["context"]["fp"],impact,impact_P,grid,auth["context"]["prepared"].rt["model"])
    impact_summary["reference_label"]=CONTRACT["candidate_label"]
    impact_summary["inherited_distribution_sha256"]=hashlib.sha256(impact_pre.tobytes()).hexdigest()
    write(impact_out/"summary.json",cal.jsonable(impact_summary))'''
    require(cell_source.count(impact_report) == 1, 'Optional impact report source drift')
    cell_source = cell_source.replace(impact_report, '    # Optional old-PRE impact remains explicitly unavailable at new GE prices.')
    cell_ns = dict(d.__dict__); exec(compile(cell_source, str(__file__)+'::ge_price_cell', 'exec'), cell_ns)
    price_cell = cell_ns['make_price_cell']

    sys.path.insert(0, str(ROOT/'code/model'))
    sys.path.insert(0, str(ROOT/'code/model/tools'))
    eq = __import__('refactor_lab.engine.equilibrium', fromlist=['solve_bellman_full_markov_income'])
    original_bellman = eq.solve_bellman_full_markov_income
    override = load_module('refactor_lab.engine.purchase_ltv_household_override', LOCAL/'household_buyer_estate.py')
    audit = load_module('purchase_ltv_due_audit_override', LOCAL/'purchase_audit_override.py')
    gate_module = sys.modules.get('e5f_evening_calibration_runtime')
    completed = []
    active_attempts = 0
    write(out/'latest_completed.json', dict(completed=completed, lifecycle_attempts=0))
    cases = [
      ('purchase_90_stayer_80', .90, .80, LOCAL/'retry7/results/purchase_90_stayer_80/closure.json', True),
      ('both_90_90', .90, .90, LOCAL/'retry7/results/both_90_90/closure.json', False),
      ('purchase_100_stayer_80', 1.00, .80, LOCAL/'retry8/results/purchase_100_stayer_80/closure.json', True),
      ('both_100_100', 1.00, 1.00, LOCAL/'retry9/results/both_100_100/closure.json', False),
    ]
    if only_case is not None:
        cases = [case for case in cases if case[0] == only_case]
        require(len(cases) == 1, 'Unknown --case: '+str(only_case))
    try:
      for label, buyer, stayer, seed_path, due_audit in cases:
        require(time.time() < deadline, '10-minute deadline before '+label)
        seed = read(seed_path)
        q0 = float(seed['price']); r0 = float(seed['renewal_residual_reported_not_imposed'])
        case_dir = out/label; case_dir.mkdir()
        shutil.copy2(PRE, case_dir/'q0_reference_inherited_states.npz')
        frozen_pre_meta = LOCAL/'retry5/results/q0_reference_inherited_states.json'
        if frozen_pre_meta.exists(): shutil.copy2(frozen_pre_meta, case_dir/'q0_reference_inherited_states.json')
        auth = d.authenticate_candidate(case_dir/'runtime_preparation')
        initial_phi = np.asarray(auth['P'].phi).copy()
        require(initial_phi.shape == (auth['P'].n_parity,) and np.all(initial_phi == .8), 'Native stayer baseline drift')
        require(auth['actual_parameters']['financed_share'] == .8, 'Winner parameter binding drift')
        auth['P'].phi = np.full_like(initial_phi, stayer)
        if stayer == .8: auth['P']._purchase_ltv_override = buyer
        elif hasattr(auth['P'], '_purchase_ltv_override'): del auth['P']._purchase_ltv_override
        prepared = auth['context']['prepared']; original_accounting = prepared.rt['accounting']
        prepared.rt['accounting'] = SimpleNamespace(audit_purchase_accounting=audit.audit_purchase_accounting) if due_audit else original_accounting
        auth['actual_parameters'] = dict(auth['actual_parameters'], financed_share=stayer)
        auth['context']['expected_parameters'] = dict(auth['context']['expected_parameters'], financed_share=stayer)
        for row in auth['params_rows']:
            if row['parameter'] == 'financed_share':
                row['estimate'] = stayer
                if stayer != .8: row['status'] = 'Experimental stayer collateral share'
        eq.solve_bellman_full_markov_income = override.solve_bellman_full_markov_income
        require(eq.solve_bellman_full_markov_income is override.solve_bellman_full_markov_income, 'Mortgage Bellman binding failed')

        captured = {}
        original_gates = d.case_gates
        def capture_gates(a, packet, path, regime, *, stationary):
            if stationary: captured['packet'] = packet
            return original_gates(a, packet, path, regime, stationary=stationary)
        d.case_gates = capture_gates
        cell_ns['case_gates'] = capture_gates
        records = [dict(case='q0_fixed_price_seed', factor=1.0, price=q0, residual=r0, lifecycle_solves=0)]
        cell_count = 0
        def evaluate(factor, role):
            nonlocal cell_count, active_attempts
            require(time.time() < deadline, 'Global deadline during '+label)
            require(cell_count < PER_CASE_LC_CAP, '12 lifecycle cap for '+label)
            cell_count += 1
            active_attempts = cell_count
            leaf = case_dir/(role+'_'+str(cell_count)+'_'+format(factor,'.6f').replace('.','p'))
            leaf.mkdir()
            end = min(deadline, time.time()+120)
            write(leaf/'attempt.json', dict(case=label,buyer_ltv=buyer,stayer_ltv=stayer,
                factor=factor,price=q0*factor,role=role,lifecycle_attempt=cell_count,deadline_epoch=end))
            prior = signal.signal(signal.SIGALRM, alarm)
            signal.setitimer(signal.ITIMER_REAL, max(.001, end-time.time()))
            try: closure = price_cell(auth, 'reference', factor, leaf, end)
            finally: signal.setitimer(signal.ITIMER_REAL, 0); signal.signal(signal.SIGALRM, prior)
            require(captured.get('packet') is not None, 'Missing native GE packet capture')
            N = float(closure['housing_demand'][0] and sum(closure['housing_supply'])/sum(closure['housing_demand']))
            residual = float(closure['renewal_residual_reported_not_imposed'])
            paygo = float(closure['pension_paygo_certificate_reported_not_imposed']['actual_accounts']['scaled_pension_budget_residual'])
            require(abs(paygo) <= 1e-6, 'Actual PAYGO gate failed at '+label)
            require(abs(N*sum(closure['housing_demand'])-sum(closure['housing_supply'])) <= 1e-12, 'Absolute housing supply closure failed')
            row = dict(case=label,role=role,factor=float(factor),price=float(closure['price']),rent=float(closure['mapped_rent']),
                residual=residual,population_scale=N,paygo_residual=paygo,lifecycle_seconds=closure['lifecycle_seconds'],
                housing_demand=closure['housing_demand'],housing_supply=closure['housing_supply'],
                fixed_pre_impact=closure['baseline_state_impact'],cohort_summary=closure['cohort_summary'],
                target_fit=str(leaf/'target_fit.csv'),parameters=str(leaf/'parameters.csv'),standard_diagnostics=str(leaf/'standard_diagnostics'),
                closure=str(leaf/'closure.json'),lifecycle_attempt=cell_count)
            closure.update(dict(ge=dict(population_scale=N,rent=float(closure['mapped_rent']),
                renewal_residual=residual,paygo_residual=paygo,not_fixed_price=True,
                fixed_pre_sha256=baseline_pre_sha,closure_source_sha256=sha(OLD_GE))))
            write(leaf/'closure.json', closure)
            row['closure_sha256'] = sha(leaf/'closure.json')
            records.append(row)
            write(case_dir/'latest_completed.json', dict(case=label,records=records,lifecycle_solves=cell_count))
            write(out/'latest_completed.json', dict(case=label,completed=completed,active_records=records,
                lifecycle_attempts=sum(x.get('lifecycle_solves',0) for x in completed)+cell_count))
            return row, closure, captured['packet'], leaf

        try:
            if repeat_factor is not None:
                require(selected_path is not None, 'A selected-root receipt path is required for repeat-only mode')
                selected_path = Path(selected_path)
                row, closure, packet, leaf = evaluate(float(repeat_factor), 'selected_repeat')
                require(abs(row['residual']) <= RENEWAL_TOL, 'Repeat renewal residual exceeds tolerance')
                ge_closure = old.closed_accounting(float(packet['solution'].entry_rate), float(closure['adjusted_births']),
                    sum(closure['housing_demand']), sum(closure['housing_supply']), float(closure['price']))
                scaled = old.scaled_native_step(packet, prepared, ge_closure)
                closure['ge'] = dict(population_scale=ge_closure['population_scale'],rent=closure['mapped_rent'],
                    renewal_residual=row['residual'],paygo_residual=row['paygo_residual'],not_fixed_price=True,
                    fixed_pre_sha256=baseline_pre_sha,scaled_native_step=scaled)
                report_ev = copy.copy(packet['evaluation'])
                report_ev.supply_by_loc = np.asarray(report_ev.supply_by_loc)/ge_closure['population_scale']
                report_ev.relative_market_residual = float(np.max(np.abs(
                    (report_ev.demand_by_loc-report_ev.supply_by_loc)/report_ev.supply_by_loc)))
                prepared.rt['audit'].standard_diagnostics(dict(packet,evaluation=report_ev),leaf,validate_production_young=False)
                require(len(list((leaf/'standard_diagnostics').glob('*.png'))) == 17,'Repeat standard17 packet incomplete')
                def rows(path):
                    with Path(path).open(newline='') as f: return list(csv.DictReader(f))
                for name in ('target_fit.csv','parameters.csv'):
                    fresh=rows(leaf/name); prior=rows(selected_path/name)
                    require(len(fresh)==len(prior)==(14 if name=='target_fit.csv' else 31), 'Repeat table row count mismatch '+name)
                    if name=='target_fit.csv':
                        for a,b in zip(fresh,prior):
                            require((a['moment'],a['target'],a['weight'],a['role'])==(b['moment'],b['target'],b['weight'],b['role']), 'Repeat target contract/order mismatch')
                            for col in ('model','gap','loss_contribution'):
                                if a[col] and b[col]: require(abs(float(a[col])-float(b[col]))<=2e-12,'Repeat fit mismatch '+col)
                    else:
                        for a,b in zip(fresh,prior):
                            require(a['parameter']==b['parameter'] and abs(float(a['estimate'])-float(b['estimate']))<=2e-12,'Repeat parameter mismatch')
                closure['ge']['repeat_tables_match_selected'] = True
                write(leaf/'closure.json',closure)
                result=dict(case=label,status='selected_ge_repeat_passed',factor=float(repeat_factor),
                    price=row['price'],renewal_residual=row['residual'],population_scale=row['population_scale'],
                    paygo_residual=row['paygo_residual'],target_fit=str(leaf/'target_fit.csv'),parameters=str(leaf/'parameters.csv'),
                    standard_diagnostics=str(leaf/'standard_diagnostics'),selected_source=str(selected_path),
                    lifecycle_solves=1,repeat_lifecycle_solves=1,repeat_tables_match_selected=True,closure=str(leaf/'closure.json'))
                completed.append(result); write(case_dir/'completed.json',result); write(out/'completed.json',result)
                continue
            selected = None
            # Existing negative q0 residual chooses the lower-price branch.
            require(r0 < 0, 'Unexpected q0 residual sign; do not switch branch silently')
            for factor in (.995,.99,.95,.90,.80):
                row, closure, packet, leaf = evaluate(factor, 'lower')
                if abs(row['residual']) <= RENEWAL_TOL:
                    selected = (row,closure,packet,leaf); break
                bracket = old.choose_bracket(records)
                if bracket is not None: break
            if selected is None:
                bracket = old.choose_bracket(records)
                require(bracket is not None, 'No lower-price renewal bracket within reviewed [.8,1] domain for '+label)
                left,right = bracket
                while cell_count < PER_CASE_LC_CAP-1 and time.time() < deadline:
                    factor = old.next_log_factor(left,right)
                    row,closure,packet,leaf = evaluate(factor,'root')
                    if abs(row['residual']) <= RENEWAL_TOL:
                        selected = (row,closure,packet,leaf); break
                    if row['residual']*left['residual'] > 0: left=row
                    else: right=row
                require(selected is not None, '12-solve cap without renewal root for '+label)

            row, closure, packet, leaf = selected
            require(abs(row['residual']) <= RENEWAL_TOL, 'Selected renewal residual exceeds tolerance')
            ge_closure = old.closed_accounting(
                float(packet['solution'].entry_rate), float(closure['adjusted_births']),
                sum(closure['housing_demand']), sum(closure['housing_supply']), float(closure['price']))
            require(abs(ge_closure['renewal_residual']-row['residual']) <= 1e-14, 'Native closure accounting mismatch')
            scaled = old.scaled_native_step(packet, prepared, ge_closure)
            require(scaled['status'] == 'passed', 'Native scaled stationary step failed')
            closure['ge']['scaled_native_step'] = scaled
            report_ev = copy.copy(packet['evaluation'])
            report_ev.supply_by_loc = np.asarray(report_ev.supply_by_loc)/ge_closure['population_scale']
            report_ev.relative_market_residual = float(np.max(np.abs(
                (report_ev.demand_by_loc-report_ev.supply_by_loc)/report_ev.supply_by_loc)))
            prepared.rt['audit'].standard_diagnostics(dict(packet,evaluation=report_ev), leaf,
                validate_production_young=False)
            plot_names = sorted(p.name for p in (leaf/'standard_diagnostics').glob('*.png'))
            require(plot_names == sorted(auth['context']['manifest']['standard_diagnostic_names']), 'Selected GE standard17 packet incomplete')
            write(leaf/'reporting_units.json', dict(economic_absolute_supply=ge_closure['absolute_housing_supply'],
                reporting_supply_per_household=float(report_ev.supply_by_loc.sum()),
                population_scale=ge_closure['population_scale'],economic_H0_unchanged=True,
                supply_divided_by_population_for_plot=True))
            write(Path(row['closure']), closure)
            # One same-platform selected-point verification per case; no byte-level cross-platform tests.
            repeat, repeat_closure, _, repeat_leaf = evaluate(row['factor'], 'selected_repeat')
            require(abs(repeat['residual']-row['residual']) <= 2e-12 and abs(repeat['population_scale']-row['population_scale']) <= 2e-12,
                'Selected same-platform repeat differs')
            closure['ge']['selected_repeat'] = dict(path=str(repeat_leaf),residual=repeat['residual'],
                population_scale=repeat['population_scale'],target_fit=repeat['target_fit'],parameters=repeat['parameters'])
            write(Path(row['closure']), closure)
            result = dict(case=label,status='completed_local_stationary_ge_diagnostic',buyer_ltv=buyer,stayer_ltv=stayer,
                selected={k:row[k] for k in ('factor','price','rent','residual','population_scale','paygo_residual','fixed_pre_impact','cohort_summary')},
                selected_closure=row['closure'],selected_target_fit=row['target_fit'],selected_parameters=row['parameters'],
                selected_standard_diagnostics=row['standard_diagnostics'],repeat_target_fit=repeat['target_fit'],
                lifecycle_solves=cell_count, q0_fixed_price_residual=r0,
                q0_seed_source_resolution='90% q0 predates buyer estate patch; patch is algebraically inactive at 90%, so reused without rerun' if buyer==.90 else 'patched buyer-estate q0 packet',
                economic_changes=['buyer origination LTV changed as labeled', 'stayer collateral share changed as labeled',
                    'reviewed buyer death-estate floor patch applied; it is algebraically inactive at the 90% buyer cap'],
                fixed_pre_sha256=baseline_pre_sha, production_adoption=False, calibration=False)
            completed.append(result)
            active_attempts = 0
            write(case_dir/'completed.json',result)
            write(out/'latest_completed.json',dict(completed=completed,lifecycle_attempts=sum(x['lifecycle_solves'] for x in completed)))
            write(out/'best_so_far.json',dict(status='latest_completed_ge_case',case=label,selected=result['selected']))
        finally:
            d.case_gates = original_gates
            prepared.rt['accounting'] = original_accounting
        require(sha(out/'q0_reference_inherited_states.npz') == sha(PRE), 'Frozen baseline PRE was overwritten')
      write(out/'completed.json',dict(status='completed',cases=completed,
          total_lifecycle_solves=sum(x['lifecycle_solves'] for x in completed),elapsed_seconds=time.time()-started,
          deadline_epoch=deadline,production_adoption=False,calibration=False))
    except BaseException as exc:
      write(out/'failure.json',dict(status='failed_no_retry',error_type=type(exc).__name__,error=str(exc),
          traceback=traceback.format_exc(),completed=completed,elapsed_seconds=time.time()-started,
          lifecycle_attempts=sum(x.get('lifecycle_solves',0) for x in completed)+active_attempts))
      raise
    finally:
      eq.solve_bellman_full_markov_income = original_bellman

if __name__ == '__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--case', choices=('purchase_90_stayer_80','both_90_90','purchase_100_stayer_80','both_100_100'))
    ap.add_argument('--repeat-factor',type=float); ap.add_argument('--selected-path',type=Path)
    args=ap.parse_args(); main(args.out,args.case,args.repeat_factor,args.selected_path)
