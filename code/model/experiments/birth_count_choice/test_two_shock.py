"""Tiny deterministic tests of the actual two-stage controller and native seam."""
import contextlib
import copy
import json
from pathlib import Path
from types import SimpleNamespace as NS
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import two_shock as m
sys.path.insert(0,str(m.ROOT/'code/model/experiments/transition_readiness/pinned_tools'))
import e5f_preference_shock_fit as scalar


def state(year):
    return NS(g_pre=np.array([1.+(year-2007)*.001],dtype=np.float64),
        scheduled_entries=np.array([year/1000.,2.],dtype=np.float64),
        scheduled_raw_entries=np.array([year/900.,3.],dtype=np.float64))


def plan():
    return dict(schema=m.SCHEMA,kind='two_unanticipated_permanent',mode='diagnostic',smoke=False,
      baseline_psi=m.BASELINE_PSI,psi_bound_ratios=[.01,2.],stages=copy.deepcopy(m.STAGES),
      source_pins={},target_contract={},rows=[dict(decision_year=2007+4*i,target=t) for i,t in enumerate(m.TARGETS)],weights=[0,1,0,1],
      horizons=[24,32],gates=m.GATES,seed=dict(horizon=12,perturbed_date=5,log_step=1e-5),
      identity=dict(source_pins={'source':'hash'},reference_sha256='r'),
      stage_starts=[dict(initial=x,bounds=[m.BASELINE_PSI*.01,m.BASELINE_PSI*2]) for x in (.15,.12)],
      budget=dict(total_seconds=21480,maximum_policy_calls=20000,candidate_seconds=100,path_seconds=100,seed_seconds=100,mapping_seconds=100,render_seconds=100,endpoint_seconds=1800),
      fit=dict(max_evaluations=12,fertility_tolerance=.005,log_difference_step=.01,max_log_step=.15,damping=.7,max_condition_number=1e8,worsening_factor=1.5,reproduction_tolerance=1e-8),
      endpoint=dict(max_evaluations=48),path=dict(max_evaluations=12),standard_plot_names=[str(i)+'.png' for i in range(17)],disclosure='test')


class FakeRuntime:
    fitter=scalar
    queue_values=staticmethod(lambda x:x)
    def __init__(self,p):self.p=p;self.calls=[];self.replayed=None;self.selected=None;self.rt=NS(total_native_calls=0)
    def bind_budget(self,deadline,heartbeat):self.deadline=deadline
    def budget_context(self):return contextlib.nullcontext()
    def prepare_reference(self,deadline,folder):return dict(accounting_valid=True,policy_calls=2)
    def initial_state(self):return state(2007)
    def measure_seed(self,**kw):
        self.calls.append(('seed',kw['start_year'],m.state_hash(kw['inherited_state'])))
        if kw['stage']:np.testing.assert_array_equal(kw['inherited_state'].scheduled_raw_entries,state(2015).scheduled_raw_entries)
        return dict(accounting_valid=True,policy_calls=5,mapping_count=5,horizon=12,nonstationary_inherited_state=kw['stage']==1,matrix=np.eye(24))
    def evaluate_stage(self,**kw):
        year=kw['start_year'];h=kw['horizon'];stage=kw['stage'];psi=kw['psi'];self.calls.append(('eval',year,psi,h))
        target=(1.861,1.64575)[stage];center=(.15,.12)[stage]
        fertility=[dict(period_tfr_topcode_adjusted=(2.03 if stage==0 else 1.73)+i*.01) for i in range(h)]
        fertility[1]['period_tfr_topcode_adjusted']=target+.2*np.log(psi/center)
        rows=[dict(calendar_year=year+4*i,asset_price=1.+i*.001,renter_price=.2,adult_population=1.,birth_children=1.,housing_demand=1.,pension_period_units=.1) for i in range(h)]
        paths=dict(prices=np.arange(h,dtype=float)*.001+1.,pensions=np.full(h,.1),psi_path=np.full(h,psi),start_year=year)
        values=[np.array([float(year+4*i)]) for i in range(h+1)]
        states={i:dict(state=state(year+4*i)) for i in range(h)}
        native=NS(dated_states=states,values=values,floor_runtime_paths=paths)
        reply=dict(identity=dict(self.p['identity'],stage_start_year=year,inherited_state_sha256=m.state_hash(kw['inherited_state'])),
          source_pins=self.p['identity']['source_pins'],housing='static-elastic',shock_contract=dict(start_year=year,psi=psi,expectations='permanent_until_next_surprise'),
          psi=psi,horizon=h,accounting_valid=True,policy_calls=10,root_pass=True,replay_pass=True,stationary_pass=True,
          market_maximum_residual=0.,fiscal_maximum_residual=0.,replay_maximum_gap=0.,stationary_renewal_gap=0.,terminal_pass=False,
          rows=rows,fertility=fertility,native_reply=native,prices=paths['prices'],pensions=paths['pensions'],values=values,dated_states=states)
        if 'folder' in kw:
            folder=Path(kw['folder']);folder.mkdir(parents=True);root=folder/'root.json';root.write_text('{}');record=folder/'native_record.json';record.write_text('{}');reply['source_evidence']=[m.pin(root),m.pin(record)]
        self.selected=reply;return reply
    def replay_prefix(self,**kw):
        self.calls.append(('prefix',kw['start_year']));self.replayed=kw
        np.testing.assert_array_equal(kw['boundary_value'],self.selected['values'][2])
        np.testing.assert_array_equal(kw['pensions'],self.selected['pensions'][:2])
        self.assert_original=kw['inherited_state'].scheduled_raw_entries.copy()
        native=NS(floor_runtime_paths=dict(start_year=2007),values=self.selected['values'][:3])
        return dict(accounting_valid=True,policy_calls=2,gates={'accounting':True},market_residual=[0.,0.],fiscal_residual=[0.,0.],
          rows=copy.deepcopy(self.selected['rows'][:2]),fertility=copy.deepcopy(self.selected['fertility'][:2]),
          values=copy.deepcopy(self.selected['values'][:3]),terminal_state=state(2015),native_reply=native)
    def set_stage2_initialization(self,accepted,candidate):self.initialization=(accepted['prices'][2:14].copy(),accepted['psi'])
    def release_stage1(self):pass
    def export_state(self,reply,**kw):
        self.calls.append(('export',reply['shock_contract']['start_year'],kw['index']));folder=Path(kw['folder']);folder.mkdir(parents=True)
        path=folder/'state';path.write_bytes(b'exact-state')
        return dict(m.pin(path),exact_native_state=True,reconstructed_or_rescaled=False,calendar_year=2023,period_index=kw['index'],queue_lags=[16,20],forecast_and_continuation_saved=True)
    def render_standard(self,reply,folder,deadline):
        directory=Path(folder)/'date_000'/'standard_diagnostics';directory.mkdir(parents=True)
        plots={}
        for name in self.p['standard_plot_names']:
            path=directory/name;path.write_bytes(b'plot');plots[name]=m.sha(path)
        return dict(sampled_dates=[dict(period=0,plots=plots)])


class Tests(unittest.TestCase):
    def execute(self,mutation=None):
        with tempfile.TemporaryDirectory() as d:
            p=plan();runtime=FakeRuntime(p)
            if mutation:mutation(runtime)
            controller=m.Controller(p,runtime,Path(d)/'out')
            with patch.object(m,'preflight',return_value={'status':'PASS'}):result=controller.run()
            return result,runtime,controller
    def test_exact_full_loop_clock_lineage_no_final_advance(self):
        result,runtime,c=self.execute();self.assertEqual([x[1] for x in runtime.calls if x[0]=='seed'],[2007,2015])
        self.assertEqual([r['model'] for r in result['fit_table']],[2.03,1.861,1.73,1.64575])
        self.assertEqual([r['target'] for r in result['fit_table']],m.TARGETS)
        self.assertEqual([r['weight'] for r in result['fit_table']],[0,1,0,1])
        self.assertEqual([r['source_stage'] for r in result['fit_table']],[1,1,2,2])
        self.assertEqual([x for x in runtime.calls if x[0]=='prefix'],[('prefix',2007)])
        self.assertEqual(runtime.calls[-1],('export',2015,2));self.assertEqual(c.calls,2+10+2+10*len([x for x in runtime.calls if x[0]=='eval']))
        self.assertEqual(runtime.replayed['boundary_value'][0],2015.);self.assertEqual(runtime.replayed['inherited_state'].g_pre[0],1.)
    def test_prefix_boundary_V_and_both_queue_mismatch(self):
        for key in ('scheduled_entries','scheduled_raw_entries'):
            def mutate(r,key=key):
                original=r.replay_prefix
                def changed(**kw):
                    reply=original(**kw);getattr(reply['terminal_state'],key)[0]+=1.;return reply
                r.replay_prefix=changed
            with self.assertRaisesRegex(ValueError,'handoff'):self.execute(mutate)
    def test_prefix_values_or_pensions_mismatch(self):
        def mutate(r):
            original=r.replay_prefix
            def changed(**kw):
                reply=original(**kw);reply['values'][0][0]+=1.;return reply
            r.replay_prefix=changed
        with self.assertRaisesRegex(ValueError,'values'):self.execute(mutate)
    def test_missing_fertility_length(self):
        p=plan();r=FakeRuntime(p);kw=dict(stage=0,psi=.15,horizon=24,start_year=2007,inherited_state=state(2007))
        a=r.evaluate_stage(**kw);b=r.evaluate_stage(**dict(kw,horizon=32));a['fertility']=a['fertility'][:3]
        with self.assertRaisesRegex(ValueError,'Four'):m.diagnostic_horizon_comparison(a,b)
    def test_housing_identity_and_failed_accounting(self):
        for field,value,match in [('housing','fixed-stock','housing'),('accounting_valid',False,'accounting')]:
            def mutate(r,field=field,value=value):
                original=r.evaluate_stage
                def changed(**kw):reply=original(**kw);reply[field]=value;return reply
                r.evaluate_stage=changed
            with self.assertRaisesRegex(ValueError,match):self.execute(mutate)
    def test_selected_replay_mismatch(self):
        p=plan();c=m.Controller(p,FakeRuntime(p),'/tmp/unused')
        c.last=dict(candidate=2,stage=1,psi=.12,reply={})
        with self.assertRaisesRegex(ValueError,'selected'):
            c.selected(1,dict(converged=True,parameter=dict(estimate=.12),root=dict(final=dict(mapping_valid=True,payload=dict(candidate=1,psi=.12)))))
    def test_shared_call_and_time_cap(self):
        p=plan();c=m.Controller(p,FakeRuntime(p),'/tmp/unused');c.calls=19999
        with self.assertRaisesRegex(ValueError,'cap'):c.account(dict(accounting_valid=True,policy_calls=2))
        c.deadline=0
        with self.assertRaisesRegex(ValueError,'deadline'):c.remaining()
    def test_state_hash_binds_both_queues_and_year_comparator(self):
        a=state(2015);b=copy.deepcopy(a);b.scheduled_raw_entries[0]+=1.
        self.assertNotEqual(m.state_hash(a),m.state_hash(b))
        p=plan();r=FakeRuntime(p);kw=dict(stage=1,psi=.12,horizon=24,start_year=2015,inherited_state=a)
        x=r.evaluate_stage(**kw);y=r.evaluate_stage(**dict(kw,horizon=32));v=m.compare_dated_state(x,y,2)
        self.assertEqual(v['index'],2);self.assertTrue(v['state_value_integrity'])
        y['native_reply'].dated_states[2]['state'].g_pre[0]=-1.
        with self.assertRaisesRegex(ValueError,'nonnegative'):m.compare_dated_state(x,y,2)
    def test_smoke_pin_tamper(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'smoke';path.write_text('{}');pin=m.pin(path);path.write_text('{"changed":true}')
            with self.assertRaisesRegex(ValueError,'changed pin'):m.validate_smoke_pin(plan(),pin)
    def test_stage_schema_targetindex_and_original_absolute_bounds(self):
        p=plan();p['stages'][1]['target_index']=3
        with self.assertRaisesRegex(ValueError,'local target'):m.preflight(p)
        p=plan();p['stage_starts'][1]['bounds'][0]*=2
        with self.assertRaisesRegex(ValueError,'absolute bounds'):m.preflight(p)

class NativeSeamTests(unittest.TestCase):
    @staticmethod
    def measured_seed_receipt():
        horizon=12; date=5
        names=('housing_imbalance<-log_house_price','housing_imbalance<-log_period_pension',
            'pension_imbalance<-log_house_price','pension_imbalance<-log_period_pension')
        profiles={name:[(i+1)*.01*(j+1) for j in range(horizon)] for i,name in enumerate(names)}
        return dict(matrix=np.eye(24),unknown_blocks=['log_house_price','log_period_pension'],
            residual_blocks=['housing_imbalance','pension_imbalance'],residual_units='physical_unscaled',
            coordinate_order=['log_house_price','log_period_pension'],horizon=horizon,perturbed_date=date,
            measured_lags=list(range(-date,horizon-date)),lag_profiles=profiles)

    def test_runtime_preserves_authenticated_receipt_through_real_extension(self):
        import two_shock_runtime as native
        adapter,rt,_=self.adapter(dict(prices=np.full(12,5.),pensions=np.full(12,.7)))
        rt.packet={};rt.reference_price=2.;rt.identity=lambda:plan()['identity']
        runtime=native.NativeRuntime.__new__(native.NativeRuntime)
        runtime.plan=adapter.plan;runtime.rt=rt;runtime.adapters={};runtime.seeds={};runtime.initialization=None
        runtime.queue_values=rt.pf.birth_queue_values
        fixture=self.measured_seed_receipt()
        with patch.object(native.StageAdapter,'measure_seed',return_value=copy.deepcopy(fixture)):
            returned=runtime.measure_seed(stage=0,start_year=2015,inherited_state=state(2015),folder='unused',deadline=m.time.monotonic()+30)
        self.assertIn('lag_profiles',runtime.seeds[0]);self.assertEqual(runtime.seeds[0]['measured_lags'],fixture['measured_lags'])
        import e5f_four_shock_acceleration as acceleration
        matrix=acceleration.extend_measured_jacobian(runtime.seeds[0],24)
        self.assertEqual(matrix.shape,(48,48));self.assertGreater(np.count_nonzero(matrix),0)
        # Adapter evaluation calls the same real extension using runtime-owned seed state.
        runtime.adapters[0]._endpoint=lambda *a:(dict(parameters=NS(pension=.3)),dict(price=2.,stationary_pass=True,stationary_renewal_gap=0.))
        class ReachedRoot(Exception): pass
        observed={}
        def stop_at_root(**kw):
            observed['jacobian']=kw['initial_jacobian'].copy()
            raise ReachedRoot()
        with tempfile.TemporaryDirectory() as d,patch.object(acceleration,'solve_joint_with_acceleration',side_effect=stop_at_root):
            with self.assertRaises(ReachedRoot):
                runtime.evaluate_stage(stage=0,psi=.12,horizon=24,start_year=2015,inherited_state=state(2015),
                    deadline=m.time.monotonic()+30,folder=Path(d))
        np.testing.assert_allclose(observed['jacobian'],matrix)

    def test_real_extension_rejects_missing_seed_provenance(self):
        import e5f_four_shock_acceleration as acceleration
        receipt=self.measured_seed_receipt();receipt.pop('coordinate_order')
        with self.assertRaisesRegex(ValueError,'physical/log convention'):
            acceleration.extend_measured_jacobian(receipt,24)

    def adapter(self,initialization=None):
        import two_shock_runtime as native
        p=plan();p['path'].update(price_bound_ratios=[.05,20.],pension_bound_ratios=[.05,20.],max_log_step=.15,damping=.7)
        p['endpoint'].update(price_bound_ratios=[.05,20.],max_log_step=.15,damping=.7,slope=1.)
        rt=NS(packet={},reference_price=2.,population_scale=1.,total_native_calls=0,
            P=NS(psi_child=m.BASELINE_PSI,pension=.3),pf=NS(birth_queue_values=lambda x:x),housing='static-elastic')
        rt.identity=lambda:p['identity'];rt.native_budget=lambda deadline,remaining:contextlib.nullcontext()
        rt.seen=[]
        def mapping(terminal,endpoint,q,b,psi,folder,initial_state=None,start_year=2007):
            rt.seen.append((initial_state,start_year,q,b,psi));rt.total_native_calls+=1
            record=dict(accounting_valid=True,policy_calls=1,gates=dict(accounting=True,housing=True,fiscal=True),
                market_residual=[0.]*len(q),fiscal_residual=[0.]*len(q),rows=[],fertility=[])
            return NS(),record
        rt.mapping=mapping;rt.terminal_checks=lambda *a,**k:dict(all_checks_pass=False)
        return native.StageAdapter(rt,p,state(2015),2015,initialization),rt,native
    def test_native_mapping_clock_state_and_explicit_endpoint(self):
        adapter,rt,native=self.adapter()
        with tempfile.TemporaryDirectory() as d:
            adapter._mapping({}, {},[2.],[.3],[.12],Path(d)/'dated',m.time.monotonic()+30)
            actual,year,*_=rt.seen[-1];self.assertEqual(year,2015);np.testing.assert_array_equal(actual.g_pre,state(2015).g_pre)
            stationary=state(2035)
            adapter._mapping({}, {},[2.],[.3],[.12],Path(d)/'endpoint',m.time.monotonic()+30,initial_state=stationary)
            self.assertIs(rt.seen[-1][0],stationary);self.assertEqual(rt.seen[-1][1],2015);self.assertEqual(adapter.calls,2)
    def test_state_bound_warm_caches_and_absolute_price_pension_bounds(self):
        init=dict(prices=np.full(12,5.),pensions=np.full(12,.7))
        adapter,rt,native=self.adapter(init);other,_,_=self.adapter()
        self.assertNotEqual(adapter.identity()['inherited_state_sha256'],m.state_hash(state(2007)))
        adapter.warm[24]={'sentinel':True};self.assertFalse(other.warm);adapter.warm.clear()
        import e5f_four_shock_acceleration as acceleration
        observed={}
        def solve(**kw):
            observed.update(q=kw['project_prices'](np.array([0.,100.])),b=kw['fiscal_bounds'],initial_q=kw['initial_prices'][0],initial_b=kw['initial_fiscal_values'][0])
            final=kw['evaluate'](kw['initial_prices'],kw['initial_fiscal_values'])
            self.assertTrue(final['mapping_valid'])
            return dict(converged=True,gates=dict(market_replay=True,fiscal_replay=True),final=dict(prices=kw['initial_prices'],fiscal_values=kw['initial_fiscal_values']),
                final_jacobian=np.eye(48),final_reproduction_max_abs=0.)
        adapter._endpoint=lambda *a:({},dict(price=2.,stationary_pass=True,stationary_renewal_gap=0.))
        # Only the endpoint's pension primitive is needed by the unchanged path algorithm.
        adapter._endpoint=lambda *a:(dict(parameters=NS(pension=.3)),dict(price=2.,stationary_pass=True,stationary_renewal_gap=0.))
        with tempfile.TemporaryDirectory() as d,patch.object(acceleration,'extend_measured_jacobian',return_value=np.eye(48)),patch.object(acceleration,'solve_joint_with_acceleration',side_effect=solve):
            reply=adapter.evaluate(psi=.12,start_year=2015,horizon=24,seed=np.eye(24),gates=m.GATES,budget=adapter.plan['budget'],
                endpoint_controls=adapter.plan['endpoint'],path_controls=adapter.plan['path'],deadline=m.time.monotonic()+30,folder=Path(d))
        np.testing.assert_array_equal(observed['q'],[.1,40.]);np.testing.assert_allclose(observed['b'],[.015,6.])
        self.assertEqual((observed['initial_q'],observed['initial_b']),(5.,.7));self.assertEqual(reply['policy_calls'],1)
        self.assertTrue(reply['root_pass']);self.assertTrue(reply['replay_pass']);self.assertFalse(reply['terminal_pass'])
    def test_export_delegates_actual_2023_local_index(self):
        adapter,rt,native=self.adapter();runtime=native.NativeRuntime.__new__(native.NativeRuntime);runtime.rt=rt
        rt.export_2023=lambda reply,folder:dict(period_index=(2023-reply['shock_contract']['start_year'])//4)
        reply=dict(shock_contract=dict(start_year=2015))
        self.assertEqual(runtime.export_state(reply,index=2,year=2023,folder='/tmp/no-write')['period_index'],2)
        with self.assertRaisesRegex(ValueError,'index'):runtime.export_state(reply,index=4,year=2023,folder='/tmp/no-write')
    def test_shared_guard_blocks_before_native_call(self):
        from floor_runtime import FloorRuntime
        runtime=FloorRuntime();runtime.total_native_calls=0
        with runtime.native_budget(m.time.monotonic()+30,1):
            runtime._guard_native_call()
            with self.assertRaisesRegex(RuntimeError,'allowance'):runtime._guard_native_call()
        self.assertEqual(runtime.total_native_calls,1)
        with runtime.native_budget(0.,1):
            with self.assertRaisesRegex(TimeoutError,'deadline'):runtime._guard_native_call()
        self.assertEqual(runtime.total_native_calls,1)

if __name__=='__main__':unittest.main()
