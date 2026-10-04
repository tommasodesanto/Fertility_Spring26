"""Concrete dated-state adapter over the authenticated, unchanged Estate-A runtime.

Only clock/state bindings and numerical initial guesses differ from NativeAdapter.
The original stationary endpoint, scalar fitter, Jacobian/root algorithms and
physical accounting gates are reused. No model solve occurs at import.
"""
from __future__ import annotations
import copy
import contextlib
import importlib
import importlib.util
import json
import math
from pathlib import Path
from types import SimpleNamespace
import sys
import time
import threading
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'code/model/experiments/transition_readiness'))
import one_shock_floor as retained
from two_shock import require, write, pinned, pin, sha, state_hash
require(Path(retained.__file__).resolve() == ROOT/'code/model/experiments/transition_readiness/one_shock_floor.py',
        'Mutable/cached retained controller import forbidden')


class StageAdapter(retained.NativeAdapter):
    def __init__(self, runtime, plan, inherited_state, start_year, initialization=None):
        super().__init__(runtime,plan)
        self.inherited_state=copy.deepcopy(inherited_state)
        self.start_year=start_year
        self.state_sha256=state_hash(inherited_state,runtime.pf.birth_queue_values)
        self.initial_q,self.initial_b=self.q,self.pension
        self.initialization=initialization
        if initialization is not None:
            self.initial_q=float(initialization['prices'][0]);self.initial_b=float(initialization['pensions'][0])
        # These caches belong to this specific dated state; never share with stage1.
        self.warm={};self.warm_by_psi={};self.endpoints={}

    def identity(self):
        return dict(self.rt.identity(),stage_start_year=self.start_year,inherited_state_sha256=self.state_sha256)

    def _mapping(self, terminal, endpoint, q, b, psi, folder, deadline, initial_state=None):
        call_deadline=min(deadline,time.monotonic()+self.plan['budget']['mapping_seconds'])
        state=self.inherited_state if initial_state is None else initial_state
        remaining=self.plan['budget']['maximum_policy_calls']-self.rt.total_native_calls
        require(remaining>=0,'Shared native call allowance exhausted')
        with retained.watchdog(call_deadline),self.rt.native_budget(call_deadline,remaining):
            reply,record=self.rt.mapping(terminal,endpoint,q,b,psi,folder,
                initial_state=state,start_year=self.start_year)
        self._account(record)
        require(all(math.isfinite(float(v)) for k in ('market_residual','fiscal_residual') for v in record[k]),
                'Finite physical residual vectors required')
        record['two_shock_provenance']=dict(start_year=self.start_year,
            inherited_state_sha256=self.state_sha256,explicit_stationary_state=initial_state is not None)
        write(Path(folder)/'native_record.json',record)
        return reply,record

    def measure_inherited_seed(self,folder,deadline):
        original,_=retained.original_modules()
        h=12;count=0;start=self.calls;init=self.initialization
        require(init is not None and len(init['prices'])==len(init['pensions'])==h,
                'Actual accepted forecast 12-date seed initialization required')
        endpoint=init['endpoint'];terminal=init['terminal_packet'];psi=init['psi']
        def evaluate(q,b):
            nonlocal count
            count+=1
            native,record=self._mapping(terminal,endpoint,q,b,np.full(h,psi),folder/f'map_{count:03d}',deadline)
            valid=bool(record['accounting_valid'] and all(record['gates'].values()))
            if count==1:
                write(folder/'baseline_checks.json',dict(valid=valid,
                    market_residual=record['market_residual'],fiscal_residual=record['fiscal_residual'],
                    stationarity_required=False,reason='Actual 2015 inherited state is nonstationary',identity=self.identity()))
            return dict(mapping_valid=valid,market_residual=record['market_residual'],fiscal_residual=record['fiscal_residual'])
        matrix=original.inner.measure_jacobian(evaluate,np.asarray(init['prices']),np.asarray(init['pensions']),
            5,1e-5,folder/'measured',dict(identity=self.identity(),native_measured=True,
                source_candidate=init['candidate'],initialization_kind=init['kind']))
        receipt=json.loads((folder/'measured/receipt.json').read_text())
        require(count==5,'Original five-map Jacobian required')
        receipt.update(matrix=matrix,policy_calls=self.calls-start,accounting_valid=True,mapping_count=count,
            horizon=h,nonstationary_inherited_state=True,identity=self.identity())
        return receipt

    def evaluate(self, *, psi, start_year, horizon, seed, gates, budget, endpoint_controls, path_controls, deadline, folder):
        import numpy as np
        original, _ = retained.original_modules()
        from e5f_four_shock_acceleration import extend_measured_jacobian, solve_joint_with_acceleration
        folder = Path(folder)
        start = self.calls
        terminal, endpoint = self._endpoint(psi, folder/'endpoint', min(deadline,time.monotonic()+budget['endpoint_seconds']))
        p = path_controls
        psi_path = np.full(horizon,psi)
        q0,b0 = self.initial_q,self.initial_b
        qT,bT = endpoint['price'],terminal['parameters'].pension
        J = extend_measured_jacobian(seed,horizon)
        # Native own-lag slopes are used for the root's emergency reset too.
        slopes = [float(np.median(np.abs(np.diag(J)[i*horizon:(i+1)*horizon]))) for i in range(2)]
        require(all(math.isfinite(s) and s > 0 for s in slopes), 'Native measured own slopes absent; no default fallback')
        warm,warm_kind = self.select_warm(horizon,psi)
        write(folder/'warm_start.json',dict(kind=warm_kind,target_psi_hex=float(psi).hex(),horizon=horizon,
            source_psi_hex=None if warm is None else warm['psi_hex'],identity=self.identity(),
            source_folder=None if warm is None else warm['source_folder'],
            fresh_native_mapping_required=True,residuals_or_fertility_reused=False))
        q = np.linspace(q0,qT,horizon) if warm is None else warm['prices'].copy()
        b = np.linspace(b0,bT,horizon) if warm is None else warm['fiscal_values'].copy()
        latest = {}
        count = 0
        path_deadline = min(deadline,time.monotonic()+budget['path_seconds'])
        def evaluate(q,b):
            nonlocal count
            count += 1
            native,record = self._mapping(terminal,endpoint,q,b,psi_path,folder/f'map_{count:03d}',path_deadline)
            latest.update(native=native,record=record)
            write(folder/'latest_completed.json',record)
            return dict(mapping_valid=all(record['gates'].values()),market_residual=record['market_residual'],fiscal_residual=record['fiscal_residual'])
        root = solve_joint_with_acceleration(closure='fixed_tax',initial_prices=q,initial_fiscal_values=b,evaluate=evaluate,
            project_prices=lambda values:np.clip(values,self.q*p['price_bound_ratios'][0],self.q*p['price_bound_ratios'][1]),
            fiscal_bounds=[self.pension*x for x in p['pension_bound_ratios']],market_tolerance=gates['market_tolerance'],
            fiscal_tolerance=gates['fiscal_tolerance'],market_slope=slopes[0],fiscal_slope=slopes[1],
            max_log_step=p['max_log_step'],damping=p['damping'],max_evaluations=p['max_evaluations'],
            deadline_monotonic=path_deadline,max_condition_number=self.plan['fit']['max_condition_number'],
            worsening_factor=self.plan['fit']['worsening_factor'],final_reproduction_tolerance=1e-10,
            initial_jacobian=J if warm is None else warm['final_jacobian'].copy(),callback=lambda row:write(folder/'root_progress.json',row))
        original.inner.write(folder/'root.json',root)
        record = latest['record']
        terminal_check = self.rt.terminal_checks(terminal,endpoint,latest['native'],psi_path,
            tolerance=1e-3,raw_queue_tolerance=1e-3)
        valid = bool(root['converged'] and terminal_check['all_checks_pass'])
        if valid or (self.plan['mode']=='diagnostic' and root['converged']):
            self.retain_warm(horizon,psi,root,folder)
        return dict(identity=self.identity(),reference_manifest_sha256=self.rt.identity()['reference_sha256'],
            source_pins=self.rt.identity()['source_pins'],housing=self.rt.housing,
            shock_contract=dict(start_year=self.start_year,psi=psi,expectations='permanent_until_next_surprise'),
            psi=psi,horizon=horizon,accounting_valid=bool(record['accounting_valid'] and all(record['gates'].values())),policy_calls=self.calls-start,
            root_and_terminal_pass=valid,stationary_pass=endpoint['stationary_pass'],
            root_pass=bool(root['converged']),terminal_pass=bool(terminal_check['all_checks_pass']),
            replay_pass=bool(root['gates']['market_replay'] and root['gates']['fiscal_replay']) if 'gates' in root else bool(root['converged']),
            stationary_renewal_gap=endpoint['stationary_renewal_gap'],
            market_maximum_residual=max(map(abs,record['market_residual'])),
            fiscal_maximum_residual=max(map(abs,record['fiscal_residual'])),
            replay_maximum_gap=root['final_reproduction_max_abs'] if root['final_reproduction_max_abs'] is not None else float('inf'),
            terminal=terminal_check,rows=record['rows'],fertility=record['fertility'],
            final_mapping_pin=pin(folder/f'map_{count:03d}'/'native_record.json'),path_evaluations=count,
            native_reply=latest['native'],terminal_packet=terminal,endpoint=endpoint,psi_path=psi_path,root=root)



class NativeRuntime:
    def __init__(self,plan,output):
        self.plan=plan;self.output=Path(output);base=json.loads(pinned(plan['base_plan']).read_text())
        overlay=plan.get('legacy_source_overlay')
        if overlay is not None:
            helper=pinned(overlay['helper']);metadata=pinned(overlay['manifest'])
            spec=importlib.util.spec_from_file_location('two_shock_authenticated_overlay',helper)
            helper_module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=helper_module;spec.loader.exec_module(helper_module)
            self.overlay_receipt=helper_module.install_overlay(metadata)
            write(self.output/'legacy_source_overlay.json',self.overlay_receipt)
        sys.path.insert(0,str(ROOT/'code/model'))
        module=importlib.import_module('experiments.birth_count_choice.transition_runtime')
        require(Path(module.__file__).resolve()==pinned(base['source_files']['runtime']).resolve(),
                'Frozen native runtime __file__ differs from pin')
        # All already loaded engine modules must originate in the frozen package.
        for name,loaded in tuple(sys.modules.items()):
            if name.startswith(('experiments.birth_count_choice','experiments.transition_readiness')) and getattr(loaded,'__file__',None):
                require(Path(loaded.__file__).resolve().is_relative_to(ROOT),'Mutable cached model import forbidden: '+name)
        self.rt=module.CurrentEstateARuntime.from_handoff(base['handoff'],self.output/'constructor')
        require(self.rt.identity()==base['identity'],'Actual native input identity differs from pinned base')
        self.adapters={};self.seeds={};self.initialization=None
        self.original,self.fitter=retained.original_modules()
        for key,module in (('original_estimator',self.original),('original_fitter',self.fitter)):
            require(Path(module.__file__).resolve()==pinned(base['source_files'][key]).resolve(),
                    'Frozen original module __file__ differs: '+key)
        contract=self.original.target_contract(pinned(base['target_contract']['blocks']),pinned(base['target_contract']['annual']))
        require(contract==base['target_contract']==plan['target_contract'],'Original annual-builder contract differs')

    def bind_budget(self,deadline,heartbeat):
        self.deadline=deadline;self.heartbeat=heartbeat
        original=self.rt._guard_native_call
        def guarded():
            original();heartbeat('native_call',actual_native_calls=self.rt.total_native_calls)
        self.rt._guard_native_call=guarded

    @contextlib.contextmanager
    def budget_context(self):
        # A long individual native call must not hide progress for five minutes.
        stop=threading.Event()
        def pulse():
            while not stop.wait(240.):
                self.heartbeat('native_in_progress',actual_native_calls=self.rt.total_native_calls)
        thread=threading.Thread(target=pulse,name='two-shock-native-heartbeat',daemon=True)
        thread.start()
        try:
            with self.rt.native_budget(self.deadline,self.plan['budget']['maximum_policy_calls']):yield
        finally:
            stop.set();thread.join(timeout=1.)

    def prepare_reference(self,deadline,folder):
        adapter=retained.NativeAdapter(self.rt,self.plan)
        receipt=adapter.prepare_native_reference(deadline,folder)
        return receipt

    def initial_state(self):
        require(self.rt.reference_verified,'Fresh native reference proof required')
        return copy.deepcopy(self.rt.initial_state)

    def queue_values(self,x): return self.rt.pf.birth_queue_values(x)

    def measure_seed(self,*,stage,start_year,inherited_state,folder,deadline,**unused):
        adapter=StageAdapter(self.rt,self.plan,inherited_state,start_year,self.initialization if stage else None)
        self.adapters[stage]=adapter
        if stage==0:
            receipt=adapter.measure_seed(horizon=12,perturbed_date=5,log_step=1e-5,gates=self.plan['gates'],
                budget=self.plan['budget'],deadline=deadline,folder=folder)
            receipt.update(mapping_count=5,horizon=12,identity=adapter.identity())
        else: receipt=adapter.measure_inherited_seed(Path(folder),deadline)
        receipt['source_evidence']=[pin(path) for path in sorted(Path(folder).glob('map_*/native_record.json'))]
        self.seeds[stage]=copy.deepcopy(receipt);return receipt

    def evaluate_stage(self,*,stage,psi,horizon,start_year,inherited_state,deadline,folder,**unused):
        adapter=self.adapters[stage]
        require(adapter.start_year==start_year and adapter.state_sha256==state_hash(inherited_state,self.queue_values),
                'Stage state/clock changed after measured seed')
        reply=adapter.evaluate(psi=psi,start_year=start_year,horizon=horizon,seed=self.seeds[stage],
            gates=self.plan['gates'],budget=self.plan['budget'],endpoint_controls=self.plan['endpoint'],
            path_controls=self.plan['path'],deadline=deadline,folder=folder)
        native=reply['native_reply'];paths=native.floor_runtime_paths
        reply.update(prices=paths['prices'],pensions=paths['pensions'],values=native.values,dated_states=native.dated_states,
            source_evidence=[pin(Path(folder)/'root.json')]+[pin(path) for path in sorted(Path(folder).rglob('native_record.json'))])
        return reply

    def replay_prefix(self,*,start_year,inherited_state,prices,pensions,values,psi,boundary_price,boundary_value,folder,deadline):
        adapter=self.adapters[0]
        require(start_year==2007 and adapter.state_sha256==state_hash(inherited_state,self.queue_values),'Original prefix state required')
        boundary=dict(evaluation=SimpleNamespace(policy=SimpleNamespace(V=boundary_value)))
        native,record=adapter._mapping(boundary,dict(price=boundary_price),np.asarray(prices),np.asarray(pensions),
            np.asarray(psi),Path(folder)/'implemented',deadline)
        return dict(record,terminal_state=native.terminal_state,values=native.values,native_reply=native)

    def set_stage2_initialization(self,accepted,candidate):
        q=np.asarray(accepted['prices'])[2:14].copy();b=np.asarray(accepted['pensions'])[2:14].copy();kind='accepted_stage1_forecast_slice_2_14'
        if len(q)<12:
            require(self.plan['smoke'] and self.plan.get('smoke_seed_endpoint_padding') is True,
                    'Stage1 forecast insufficient for required stage2 seed')
            q=np.r_[q,np.full(12-len(q),accepted['endpoint']['price'])]
            b=np.r_[b,np.full(12-len(b),accepted['terminal_packet']['parameters'].pension)]
            kind='smoke_only_accepted_slice_then_stationary_endpoint_numerical_padding'
        self.initialization=dict(prices=q,pensions=b,psi=accepted['psi'],endpoint=accepted['endpoint'],
            terminal_packet=accepted['terminal_packet'],candidate=candidate,kind=kind)

    def release_stage1(self):
        self.adapters.pop(0,None);self.seeds.pop(0,None)

    def export_state(self,reply,*,index,year,folder):
        require(year==2023 and index==(year-reply['shock_contract']['start_year'])//4,'Actual dated export index differs')
        return self.rt.export_2023(reply,folder)

    def render_standard(self,reply,folder,deadline):
        with retained.watchdog(deadline): return self.rt.render_standard(reply,folder)


def build_runtime(*,plan,output,smoke=False):
    require(plan['smoke'] is bool(smoke),'Explicit smoke mode differs')
    return NativeRuntime(plan,output)
