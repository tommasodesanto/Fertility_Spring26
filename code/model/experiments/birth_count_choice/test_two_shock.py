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
    p=dict(schema=m.SCHEMA,kind='two_unanticipated_permanent',mode='diagnostic',smoke=False,
      baseline_psi=m.BASELINE_PSI,psi_bound_ratios=[.01,2.],stages=copy.deepcopy(m.STAGES),
      source_pins={},target_contract={},rows=[dict(decision_year=2007+4*i,target=t) for i,t in enumerate(m.TARGETS)],weights=[0,1,0,1],
      horizons=[24,32],smoke_seed_endpoint_padding=False,gates=m.GATES,seed=dict(horizon=12,perturbed_date=5,log_step=1e-5),
      identity=dict(source_pins={'source':'hash'},reference_sha256='r'),
      stage_starts=[dict(initial=x,bounds=[m.BASELINE_PSI*.01,m.BASELINE_PSI*2]) for x in (.15,.12)],
      budget=dict(total_seconds=21480,maximum_policy_calls=20000,candidate_seconds=100,path_seconds=100,seed_seconds=100,mapping_seconds=100,render_seconds=100,endpoint_seconds=1800),
      fit=dict(max_evaluations=12,fertility_tolerance=.005,log_difference_step=.01,max_log_step=.15,damping=.7,max_condition_number=1e8,worsening_factor=1.5,reproduction_tolerance=1e-8),
      endpoint=dict(max_evaluations=48,price_bound_ratios=[.05,20.],max_log_step=.15,slope=1.),path=dict(max_evaluations=12,price_bound_ratios=[.05,20.],pension_bound_ratios=[.05,20.],max_log_step=.15),standard_plot_names=[str(i)+'.png' for i in range(17)],disclosure='test')
    p['empirical_controls']=m.numerical_controls(p);p['smoke_protocol']=m.SMOKE_PROTOCOL;p['smoke_targets']=None
    return p


class FakeRuntime:
    fitter=scalar
    queue_values=staticmethod(lambda x:x)
    def __init__(self,p):self.p=p;self.calls=[];self.replayed=None;self.selected=None;self.rt=NS(total_native_calls=0)
    def bind_budget(self,deadline,heartbeat):self.deadline=deadline
    def budget_context(self):return contextlib.nullcontext()
    def prepare_reference(self,deadline,folder):
        self.rt.total_native_calls+=2;folder=Path(folder);folder.mkdir(parents=True)
        packet=folder/'reference';packet.write_bytes(b'fresh-reference')
        receipt=folder/'reference_reconstruction.json';m.write(receipt,dict(status='passed',policy_calls=2,checkpoint_sha256=m.sha(packet)))
        return dict(accounting_valid=True,policy_calls=2,reference_checkpoint=m.pin(packet),reconstruction_receipt=m.pin(receipt))
    def initial_state(self):return state(2007)
    def measure_seed(self,**kw):
        self.calls.append(('seed',kw['start_year'],m.state_hash(kw['inherited_state'])))
        if kw['stage']:np.testing.assert_array_equal(kw['inherited_state'].scheduled_raw_entries,state(2015).scheduled_raw_entries)
        self.rt.total_native_calls+=5;proofs=[]
        for i in range(5):
            path=Path(kw['folder'])/f'map_{i+1:03d}'/'native_record.json'
            m.write(path,dict(accounting_valid=True,gates=dict(accounting=True),
                rows=[dict(calendar_year=kw['start_year']+4*j) for j in range(self.p['seed']['horizon'])],
                two_shock_provenance=dict(start_year=kw['start_year'],inherited_state_sha256=m.state_hash(kw['inherited_state']))))
            proofs.append(m.pin(path))
        return dict(source_evidence=proofs,identity=dict(self.p['identity'],stage_start_year=kw['start_year'],inherited_state_sha256=m.state_hash(kw['inherited_state'])),accounting_valid=True,policy_calls=5,mapping_count=5,horizon=self.p['seed']['horizon'],perturbed_date=self.p['seed']['perturbed_date'],
            perturbation_log_step=self.p['seed']['log_step'],nonstationary_inherited_state=kw['stage']==1,matrix=np.eye(2*self.p['seed']['horizon']))
    def evaluate_stage(self,**kw):
        year=kw['start_year'];h=kw['horizon'];stage=kw['stage'];psi=kw['psi'];self.rt.total_native_calls+=10;self.calls.append(('eval',year,psi,h))
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
          market_maximum_residual=0.,fiscal_maximum_residual=0.,replay_maximum_gap=0.,stationary_renewal_gap=0.,terminal_pass=False,terminal=dict(all_checks_pass=False,raw_queue_pass=False,raw_queue_maximum_relative_gap=.01),
          rows=rows,fertility=fertility,native_reply=native,prices=paths['prices'],pensions=paths['pensions'],values=values,dated_states=states)
        if 'folder' in kw:
            folder=Path(kw['folder']);folder.mkdir(parents=True);root=folder/'root.json'
            m.write(root,dict(converged=True,gates=dict(market_replay=True,fiscal_replay=True),final_reproduction_max_abs=0.,
                final=dict(mapping_valid=True,prices=paths['prices'],fiscal_values=paths['pensions'])))
            record=folder/'native_record.json'
            m.write(record,dict(accounting_valid=True,gates=dict(accounting=True),rows=rows,market_residual=[0.]*h,fiscal_residual=[0.]*h,
                two_shock_provenance=dict(start_year=year,inherited_state_sha256=m.state_hash(kw['inherited_state']))))
            reply.update(source_evidence=[m.pin(root),m.pin(record)],final_mapping_pin=m.pin(record),path_evaluations=2)
        self.selected=reply;return reply
    def replay_prefix(self,**kw):
        self.rt.total_native_calls+=2
        self.calls.append(('prefix',kw['start_year']));self.replayed=kw
        np.testing.assert_array_equal(kw['boundary_value'],self.selected['values'][2])
        np.testing.assert_array_equal(kw['pensions'],self.selected['pensions'][:2])
        self.assert_original=kw['inherited_state'].scheduled_raw_entries.copy()
        native=NS(floor_runtime_paths=dict(start_year=2007),values=self.selected['values'][:3])
        return dict(accounting_valid=True,policy_calls=2,gates={'accounting':True},market_residual=[0.,0.],fiscal_residual=[0.,0.],
          rows=copy.deepcopy(self.selected['rows'][:2]),fertility=copy.deepcopy(self.selected['fertility'][:2]),
          values=copy.deepcopy(self.selected['values'][:3]),terminal_state=state(2015),native_reply=native)
    def set_stage2_initialization(self,accepted,candidate):self.initialization=(accepted['prices'][2:2+self.p['seed']['horizon']].copy(),accepted['psi'])
    def release_stage1(self):pass
    def export_state(self,reply,**kw):
        self.calls.append(('export',reply['shock_contract']['start_year'],kw['index']));folder=Path(kw['folder']);folder.mkdir(parents=True)
        path=folder/'state.pkl.gz';index=kw['index'];actual=reply['dated_states'][index]['state'];h=reply['horizon']
        packet=dict(schema='current_floor_actual_2023_v1',initial_state=actual,calendar_year=2023,period=index,
            continuation_calendar_year=2027,reference_identity=self.p['identity'],initial_population=float(actual.g_pre.sum()),
            scheduled_entries=actual.scheduled_entries.copy(),scheduled_raw_entries=actual.scheduled_raw_entries.copy(),
            current_2023_V=reply['values'][index],continuation_V=reply['values'][index+1],terminal_V=reply['values'][-1],
            forecast_prices=reply['prices'][index:],forecast_pensions=reply['pensions'][index:],forecast_psi=np.full(h-index,reply['psi']),b_grid=np.array([1.]))
        with m.gzip.open(path,'wb') as stream:m.pickle.dump(packet,stream)
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

class ManifestControlTests(unittest.TestCase):
    def prepared_pair(self,directory):
        base=plan();base['target_contract']={'rows':base['rows']}
        path=Path(directory)/'base.json';path.write_text(json.dumps(base))
        return [m.prepare_manifest(m.pin(path),smoke=smoke) for smoke in (True,False)]

    def test_execution_controls_separate_original_empirical_fingerprint(self):
        with tempfile.TemporaryDirectory() as d:smoke,fit=self.prepared_pair(d)
        self.assertEqual(smoke['seed'],dict(horizon=4,perturbed_date=1,log_step=1e-5))
        self.assertEqual(fit['seed'],dict(horizon=12,perturbed_date=5,log_step=1e-5))
        self.assertEqual(smoke['horizons'],[6]);self.assertEqual(smoke['path']['max_evaluations'],3)
        self.assertEqual(smoke['endpoint']['max_evaluations'],8);self.assertEqual(smoke['budget']['total_seconds'],1680)
        self.assertEqual({k:smoke['budget'][k] for k in ('seed_seconds','mapping_seconds','candidate_seconds','path_seconds')},
            dict(seed_seconds=360,mapping_seconds=180,candidate_seconds=360,path_seconds=240))
        self.assertEqual(smoke['budget']['maximum_policy_calls'],400);self.assertIs(smoke['smoke_seed_endpoint_padding'],True)
        self.assertEqual(smoke['gates']['market_tolerance'],.05);self.assertEqual(smoke['gates']['fiscal_tolerance'],.005)
        self.assertEqual(fit['horizons'],[24,32]);self.assertEqual(fit['gates'],m.GATES)
        self.assertEqual(fit['path']['max_evaluations'],12);self.assertEqual(fit['fit']['max_evaluations'],12)
        self.assertEqual(fit['endpoint']['max_evaluations'],48);self.assertEqual(fit['budget']['endpoint_seconds'],1800)
        self.assertEqual(fit['budget']['total_seconds'],21480);self.assertEqual(fit['budget']['maximum_policy_calls'],20000)
        self.assertIs(fit['smoke_seed_endpoint_padding'],False)
        for p in (smoke,fit):m.validate_mode_controls(p)
        self.assertEqual(smoke['empirical_controls'],m.numerical_controls(fit))
        self.assertEqual(m.numerical_controls(fit),plan()['empirical_controls'])
        a,b=m.fingerprints(smoke),m.fingerprints(fit)
        for key in ('source','contract','empirical_controls','smoke_controls'):self.assertEqual(a[key],b[key])
        self.assertNotEqual(a['controls'],b['controls'])

    def test_smoke_cannot_relax_empirical_contract(self):
        with tempfile.TemporaryDirectory() as d:pair=self.prepared_pair(d)
        for mode in pair:
            for block,key,value in [('gates','market_tolerance',.05),('gates','fiscal_tolerance',.005),
                ('budget','total_seconds',1680),('budget','maximum_policy_calls',400),('budget','endpoint_seconds',240),
                ('fit','max_evaluations',3),('path','max_evaluations',3),('endpoint','max_evaluations',8)]:
                p=copy.deepcopy(mode);p['empirical_controls'][block][key]=value
                with self.subTest(smoke=p['smoke'],block=block,key=key),self.assertRaisesRegex(ValueError,'Original empirical'):m.validate_mode_controls(p)
            for key,value in [('seed',dict(horizon=4,perturbed_date=1,log_step=1e-5)),('horizons',[6]),('smoke_seed_endpoint_padding',True)]:
                p=copy.deepcopy(mode);p['empirical_controls'][key]=value
                with self.assertRaisesRegex(ValueError,'Original empirical'):m.validate_mode_controls(p)

    def test_modes_cannot_change_actual_budgets_gates_horizons_padding(self):
        with tempfile.TemporaryDirectory() as d:pair=self.prepared_pair(d)
        for mode in pair:
            for block,key,value in [('gates','market_tolerance',.5),('budget','maximum_policy_calls',401),
                ('budget','total_seconds',1801),('path','max_evaluations',4),('endpoint','max_evaluations',9)]:
                p=copy.deepcopy(mode);p[block][key]=value
                with self.assertRaisesRegex(ValueError,'Exact smoke'):m.validate_mode_controls(p)
            p=copy.deepcopy(mode);p['smoke_protocol']='native_two_stage_execution_only_v1'
            with self.assertRaisesRegex(ValueError,'execution-only'):m.validate_mode_controls(p)
            p=copy.deepcopy(mode);p['smoke_seed_endpoint_padding']=not p['smoke_seed_endpoint_padding']
            with self.assertRaisesRegex(ValueError,'Exact smoke'):m.validate_mode_controls(p)

    def test_old_fitted_smoke_receipt_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'receipt.json';path.write_text(json.dumps(dict(status='matched',schema=m.SCHEMA,smoke=True)))
            with self.assertRaisesRegex(ValueError,'execution-only'):m.validate_smoke_pin(plan(),m.pin(path))


class ExecutionSmokeTests(unittest.TestCase):
    def smoke_plan(self):
        p=plan();p['smoke']=True;p.update(m.smoke_controls(p['empirical_controls']))
        for x in p['stage_starts']:x['initial']=m.BASELINE_PSI
        return p

    def execute(self,directory):
        p=self.smoke_plan();r=FakeRuntime(p);c=m.Controller(p,r,Path(directory)/'out')
        with patch.object(m,'preflight',return_value={'status':'PASS'}),patch.object(r.fitter,'fit_one',side_effect=AssertionError('Scalar fitter forbidden')):
            result=c.run()
        return result,r,c

    def test_fixed_baseline_integration_and_actual_receipt_release(self):
        with tempfile.TemporaryDirectory() as d:
            result,r,c=self.execute(d)
            self.assertEqual([x for x in r.calls if x[0]=='eval'],[('eval',2007,m.BASELINE_PSI,6),('eval',2015,m.BASELINE_PSI,6)])
            self.assertEqual([x for x in r.calls if x[0]=='seed'],[('seed',2007,m.state_hash(state(2007))),('seed',2015,m.state_hash(state(2015)))])
            self.assertEqual(result['seed_controls'],dict(horizon=4,perturbed_date=1,log_step=1e-5))
            self.assertIn('[2:6]',result['initializer']);self.assertIn('no endpoint padding executed',result['initializer'])
            self.assertIn('empirical 12-date Jacobian at perturbed date 5',result['untested'])
            self.assertIn('empirical 12-to-24/32 measured-Jacobian extension',result['untested'])
            self.assertEqual(result['status'],'execution_passed');self.assertIs(result['empirical_fitted'],False)
            self.assertEqual(result['terminal_passes_by_stage'],{'1':False,'2':False});self.assertIs(result['terminal_diagnostics_gating'],False)
            for proof in result['selected_candidate_pins']:
                candidate=json.loads(m.pinned(proof).read_text());self.assertIs(candidate['terminal_pass'],False)
                terminal=json.loads(m.pinned(candidate['terminal_diagnostics_pin']).read_text())
                self.assertIs(terminal['all_checks_pass'],False);self.assertEqual(terminal['raw_queue_maximum_relative_gap'],.01)
            for proof in result['diagnostic_pins']:
                diagnostic=json.loads(m.pinned(proof).read_text())
                self.assertEqual(set(diagnostic['sampled_dates'][0]['plots']),set(plan()['standard_plot_names']))
                self.assertEqual(len(diagnostic['sampled_dates'][0]['plots']),17)
            with m.gzip.open(result['exact_2023_state']['path'],'rb') as stream:packet=m.pickle.load(stream)
            np.testing.assert_array_equal(packet['scheduled_entries'],state(2023).scheduled_entries)
            np.testing.assert_array_equal(packet['scheduled_raw_entries'],state(2023).scheduled_raw_entries)
            self.assertEqual(c.calls,34);self.assertEqual(r.replayed['boundary_price'],1.002)
            self.assertEqual(r.replayed['boundary_value'][0],2015.);np.testing.assert_array_equal(r.replayed['pensions'],[.1,.1])
            self.assertFalse((Path(d)/'out/fertility_fit.csv').exists());self.assertFalse((Path(d)/'out/stage1/fit.json').exists())
            with patch.object(m,'preflight',return_value={'status':'PASS'}):m.validate_smoke_pin(plan(),m.pin(Path(d)/'out/complete.json'))

    def test_controller_rejects_seed_controls_and_missing_actual_maps(self):
        for key,value in [('horizon',12),('perturbed_date',2),('perturbation_log_step',1e-4),('mapping_count',4),('source_evidence',[])]:
            with self.subTest(key=key),tempfile.TemporaryDirectory() as d:
                p=self.smoke_plan();r=FakeRuntime(p);measure=r.measure_seed
                def wrong_seed(**kw):
                    receipt=measure(**kw);receipt[key]=value;return receipt
                r.measure_seed=wrong_seed;c=m.Controller(p,r,Path(d)/'out')
                with patch.object(m,'preflight',return_value={'status':'PASS'}),self.assertRaises(ValueError):c.run()

    def test_receipt_rejects_changed_plots_or_actual_checkpoint(self):
        for kind in ('plot','state','handoff','root','mapping','call_cap','fingerprint','seed_horizon','seed_date','seed_step','seed_maps','seed_proofs'):
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as d:
                result,r,c=self.execute(d);path=Path(d)/'out/complete.json'
                if kind=='plot':
                    pin=result['diagnostic_pins'][-1];plot=Path(pin['path']).parent/'date_000/standard_diagnostics/0.png';plot.write_bytes(b'changed')
                elif kind=='state':
                    exported=result['exact_2023_state']
                    with m.gzip.open(exported['path'],'rb') as stream:packet=m.pickle.load(stream)
                    packet['initial_state'].scheduled_raw_entries[0]+=1.
                    with m.gzip.open(exported['path'],'wb') as stream:m.pickle.dump(packet,stream)
                    # Re-pin all referenced files: content checks must still fail.
                    new=m.pin(exported['path']);exported.update(new)
                    for item in result['evidence_pins']:
                        if item['path']==new['path']:item.update(new)
                    cp=m.pinned(result['selected_candidate_pins'][-1]);candidate=json.loads(cp.read_text());candidate['exact_2023_state'].update(new);m.write(cp,candidate);result['selected_candidate_pins'][-1]=m.pin(cp)
                elif kind=='handoff':
                    cp=m.pinned(result['handoff_pin']);hand=json.loads(cp.read_text());hand['boundary_b']+=1.;m.write(cp,hand)
                    new=m.pin(cp);result['handoff_pin']=new
                    for item in result['evidence_pins']:
                        if item['path']==new['path']:item.update(new)
                elif kind in ('root','mapping'):
                    cp=m.pinned(result['selected_candidate_pins'][-1]);candidate=json.loads(cp.read_text());key='root_pin' if kind=='root' else 'final_mapping_pin'
                    item=m.pinned(candidate[key]);value=json.loads(item.read_text())
                    if kind=='root':value['final_reproduction_max_abs']=1.
                    else:value['accounting_valid']=False
                    m.write(item,value);new=m.pin(item);candidate[key]=new
                    for pins in (candidate['source_evidence'],result['evidence_pins']):
                        for proof in pins:
                            if proof['path']==new['path']:proof.update(new)
                    m.write(cp,candidate);result['selected_candidate_pins'][-1]=m.pin(cp)
                elif kind.startswith('seed_'):
                    cp=m.pinned(result['stage_seed_pins'][-1]);seed=json.loads(cp.read_text())
                    key,value={'seed_horizon':('horizon',12),'seed_date':('perturbed_date',2),'seed_step':('perturbation_log_step',1e-4),
                        'seed_maps':('mapping_count',4),'seed_proofs':('source_evidence',seed['source_evidence'][:4])}[kind]
                    seed[key]=value;m.write(cp,seed);result['stage_seed_pins'][-1]=m.pin(cp)
                elif kind=='call_cap':result['actual_policy_calls']=result['native_actual_policy_calls']=401
                else:result['fingerprints']['empirical_controls']='changed'
                m.write(path,result)
                with patch.object(m,'preflight',return_value={'status':'PASS'}),self.assertRaises(ValueError):m.validate_smoke_pin(plan(),m.pin(path))



class NativeSeamTests(unittest.TestCase):
    def test_real_four_date_jacobian_five_maps_and_six_date_lag_support(self):
        import run_e5f_preference_transition as inner
        import e5f_four_shock_acceleration as acceleration
        profiles=[{-1:.2,0:2.,1:.5,2:.1},{-1:.1,0:.4,1:.3,2:.2},
            {-1:.3,0:.7,1:.2,2:.1},{-1:.2,0:3.,1:.4,2:.3}]
        def matrix(h):
            blocks=[sum((np.diag(np.full(h-abs(lag),value),-lag)
                for lag,value in profile.items()),np.zeros((h,h))) for profile in profiles]
            return np.block([[blocks[0],blocks[1]],[blocks[2],blocks[3]]])
        derivative=matrix(4);calls=[]
        def evaluate(q,b):
            calls.append((q.copy(),b.copy()))
            residual=derivative@np.r_[np.log(q),np.log(b)]
            return dict(mapping_valid=True,market_residual=residual[:4],fiscal_residual=residual[4:])
        with tempfile.TemporaryDirectory() as d:
            measured=inner.measure_jacobian(evaluate,np.ones(4),np.ones(4),1,1e-5,Path(d)/'measured',{})
            receipt=json.loads((Path(d)/'measured/receipt.json').read_text())
        self.assertEqual(len(calls),5);self.assertEqual(receipt['mapping_count'],5)
        self.assertEqual((receipt['horizon'],receipt['perturbed_date']),(4,1))
        self.assertEqual(receipt['measured_lags'],[-1,0,1,2])
        np.testing.assert_allclose(measured,derivative,atol=1e-10)
        extended=acceleration.extend_measured_jacobian(receipt,6)
        np.testing.assert_allclose(extended,matrix(6),atol=1e-10)
        self.assertEqual(extended[5,0],0.) # Lag5 was never measured.
        for index,(q,b) in enumerate(calls[1:]):
            changed=np.flatnonzero(np.r_[q,b]!=1.)
            np.testing.assert_array_equal(changed,[1 if index<2 else 5])

    def test_actual_runtime_routes_pinned_seed_controls_both_stages(self):
        import two_shock_runtime as native
        original,_=native.retained.original_modules()
        for smoke in (False,True):
            with self.subTest(smoke=smoke),tempfile.TemporaryDirectory() as d:
                adapter,rt,_=self.adapter()
                runtime=native.NativeRuntime.__new__(native.NativeRuntime)
                runtime.plan=adapter.plan
                runtime.plan['smoke']=smoke
                if smoke:runtime.plan.update(m.smoke_controls(runtime.plan['empirical_controls']))
                runtime.rt=rt;runtime.adapters={};runtime.seeds={};runtime.initialization=None
                runtime.queue_values=rt.pf.birth_queue_values
                rt.terminal_checks=lambda *a,**k:dict(all_checks_pass=True)
                h=runtime.plan['seed']['horizon'];date=runtime.plan['seed']['perturbed_date']
                accepted=dict(prices=np.arange(6 if smoke else 24.)+1.,
                    pensions=np.full(6 if smoke else 24,.2),psi=m.BASELINE_PSI,
                    endpoint=dict(price=2.),terminal_packet=dict(parameters=NS(pension=.3)))
                with patch.object(original.inner,'measure_jacobian',wraps=original.inner.measure_jacobian) as measure:
                    for stage,year in enumerate((2007,2015)):
                        if stage:runtime.set_stage2_initialization(accepted,1)
                        receipt=runtime.measure_seed(stage=stage,start_year=year,inherited_state=state(year),
                            folder=Path(d)/f'stage{stage}',deadline=m.time.monotonic()+30)
                        self.assertEqual((receipt['horizon'],receipt['perturbed_date']),(h,date))
                        self.assertEqual(receipt['perturbation_log_step'],1e-5)
                        self.assertEqual(receipt['mapping_count'],5);self.assertEqual(len(receipt['source_evidence']),5)
                        self.assertEqual(receipt['measured_lags'],list(range(-date,h-date)))
                        self.assertEqual(measure.call_args.args[3:5],(date,1e-5))
                        self.assertEqual(len(measure.call_args.args[1]),h)
                        if stage:
                            np.testing.assert_array_equal(measure.call_args.args[1],accepted['prices'][2:2+h])
                            np.testing.assert_array_equal(measure.call_args.args[2],accepted['pensions'][2:2+h])
                            self.assertEqual(receipt['initialization_kind'],f'accepted_stage1_forecast_slice_2_{2+h}')
                    self.assertEqual(measure.call_count,2)
                self.assertEqual(rt.total_native_calls,10)
                self.assertEqual([entry[1] for entry in rt.seen],[2007]*5+[2015]*5)


    @staticmethod
    def measured_seed_receipt():
        horizon=12; date=5
        names=('housing_imbalance<-log_house_price','housing_imbalance<-log_period_pension',
            'pension_imbalance<-log_house_price','pension_imbalance<-log_period_pension')
        profiles={name:[(i+1)*.01*(j+1) for j in range(horizon)] for i,name in enumerate(names)}
        return dict(matrix=np.eye(24),mapping_count=5,perturbation_log_step=1e-5,unknown_blocks=['log_house_price','log_period_pension'],
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
    def test_short_initializer_padding_exclusively_smoke(self):
        import two_shock_runtime as native
        accepted=dict(prices=np.arange(6.)+1.,pensions=np.full(6,.2),psi=m.BASELINE_PSI,
            endpoint=dict(price=2.),terminal_packet=dict(parameters=NS(pension=.3)))
        runtime=native.NativeRuntime.__new__(native.NativeRuntime);runtime.plan=plan()
        with self.assertRaisesRegex(ValueError,'forecast insufficient'):runtime.set_stage2_initialization(accepted,1)
        runtime.plan['smoke']=True;runtime.plan.update(m.smoke_controls(runtime.plan['empirical_controls']))
        runtime.set_stage2_initialization(accepted,1)
        np.testing.assert_array_equal(runtime.initialization['prices'],accepted['prices'][2:6])
        np.testing.assert_array_equal(runtime.initialization['pensions'],accepted['pensions'][2:6])
        self.assertEqual(runtime.initialization['kind'],'accepted_stage1_forecast_slice_2_6')
        short=dict(accepted,prices=accepted['prices'][:4],pensions=accepted['pensions'][:4])
        runtime.set_stage2_initialization(short,1)
        np.testing.assert_array_equal(runtime.initialization['prices'],[3.,4.,2.,2.])
        self.assertIn('smoke_only',runtime.initialization['kind'])
        runtime.plan['smoke_seed_endpoint_padding']=False
        with self.assertRaisesRegex(ValueError,'forecast insufficient'):runtime.set_stage2_initialization(short,1)
        runtime.plan=plan()
        accepted.update(prices=np.arange(24.)+1.,pensions=np.arange(24.)*.01+.2)
        runtime.set_stage2_initialization(accepted,1)
        np.testing.assert_array_equal(runtime.initialization['prices'],accepted['prices'][2:14])
        np.testing.assert_array_equal(runtime.initialization['pensions'],accepted['pensions'][2:14])
        self.assertEqual(runtime.initialization['kind'],'accepted_stage1_forecast_slice_2_14')

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
