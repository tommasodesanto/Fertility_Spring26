import ast
import copy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
import numpy as np
import run_e5f_current_transition as driver


class RootContractTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.root=Path(self.temp.name)
        def artifact(name,obj):
            path=self.root/name;path.write_text(json.dumps(obj));return dict(path=str(path),sha256=driver.sha(path))
        self.artifact=artifact
        tools=Path(driver.__file__).parent
        pins={str((tools/name).resolve()):driver.sha(tools/name) for name in
            ('run_e5f_current_transition.py','e5f_current_transition_runtime.py','e5f_social_security_root.py','e5f_matched_pf_path_root.py','e5f_exact_policy_cache.py','run_e5f_independent_numerical_audit.py')}
        pins[str((tools.parent/'intergen_eqscale_seq_optimized/solver.py').resolve())]=driver.sha(tools.parent/'intergen_eqscale_seq_optimized/solver.py')
        replay=artifact('replay.json',{'status':'native_reference_replay_pass'})
        smoke=artifact('smoke.json',{'status':'passed_native_smoke'})
        approval=artifact('approval.json',dict(approved=True,scope='native_current_transition',source_pins=pins,
            reference_checkpoint_sha256='reference',smoke_receipt=smoke,reference_replay=replay))
        checkpoint=artifact('checkpoint.json',{'synthetic':'no pickle needed for pure gate'})
        self.plan=dict(budget=dict(horizon=2,max_evaluations=5,deadline_epoch=100.,case_seconds=10.,
                observed_mapping_seconds=1.,solve_count_estimate=20,cache_max_bytes=2*1024**3),
            economics={'credit':'experimental_natural_solvency','fiscal':'fixed_tax','property_rebates':False,
                'outside_entry':0.,'retention':1.,'psi':'fixed_reference','entry_clock':'split_birth_vintage',
                'estate':'provisional_net_estates_fund_actual_next_cohort_sink_residual'},
            classification={key:dict(status='externally_fixed',description='synthetic') for key in
                ('preferences','earnings','entry','population','fiscal','geography','estate','credit','targets')},
            numerics=dict(market_tolerance=2e-4,fiscal_tolerance=1e-6,final_reproduction_tolerance=1e-10,
                initial_prices=[1.,1.],initial_pensions=[1.,1.],price_bounds=[.5,2.],pension_bounds=[.5,2.],market_slope=1.,fiscal_slope=1.,max_log_step=.1,damping=.5,max_condition_number=1e8,worsening_factor=2.,terminal_tolerances={'population_relative_gap':.01}),
            source_pins=pins,approval=approval,reference_checkpoint_sha256='reference',
            endpoint_complete=artifact('endpoint.json',dict(usable_closed_root=True,repeat_verified=True)),
            endpoint_receipt=artifact('receipt.json',dict(case_checkpoint_sha256=checkpoint['sha256'])),
            endpoint_checkpoint=checkpoint)
        self.plan['terminal_approval']=artifact('terminal_approval.json',dict(approved=True,status='native_terminal_replay_verified',
            terminal_checkpoint_sha256=checkpoint['sha256'],terminal_receipt_path=self.plan['endpoint_receipt']['path'],
            terminal_receipt_sha256=self.plan['endpoint_receipt']['sha256'],current_source_files=pins))

    def tearDown(self):self.temp.cleanup()

    def check(self,plan=None):
        return driver.validate(self.plan if plan is None else plan,horizon=2,max_evaluations=5,
            deadline=100.,case_seconds=10.,now=1.)

    def test_accepts_explicit_pinned_contract(self):self.check()

    def test_rejects_missing_approval_and_changed_source(self):
        p=copy.deepcopy(self.plan);p['approval']=self.artifact('no.json',{'approved':False})
        with self.assertRaises(ValueError):self.check(p)
        p=copy.deepcopy(self.plan);p['source_pins'][str(Path(driver.__file__).resolve())]='wrong'
        with self.assertRaises(ValueError):self.check(p)

    def test_rejects_closure_and_budget_changes(self):
        for section,key,value in [('economics','psi','renormalized'),('economics','property_rebates',True),
                                  ('budget','horizon',3),('numerics','fiscal_tolerance',1e-4)]:
            p=copy.deepcopy(self.plan);p[section][key]=value
            with self.assertRaises(ValueError):self.check(p)

    def test_rejects_unverified_endpoint(self):
        p=copy.deepcopy(self.plan);p['endpoint_complete']=self.artifact('badendpoint.json',dict(usable_closed_root=False,repeat_verified=True))
        with self.assertRaises(ValueError):self.check(p)

    def test_json_arrays_and_nonfinite_records(self):
        path=self.root/'record.json';driver.write(path,dict(x=np.array([1.,np.inf]),v=np.bool_(True)))
        self.assertEqual(json.loads(path.read_text()),dict(x=[1.,None],v=True))

    def test_population_scaling_only_terminal_distance(self):
        class PF:
            def stationary_initial_state(self,g,e,b,P,c):
                self.state=(g,e,b);return 'state'
            def terminal_convergence_diagnostics(self,**kw):self.kw=kw;return {'all_checks_pass':True}
        pf=PF();g=np.array([.4,.6]);V=np.array([2.,3.])
        terminal=dict(stationary_g_pre=g,solution=SimpleNamespace(entry_rate=.1),evaluation=SimpleNamespace(births=.2,policy=SimpleNamespace(V=V)))
        P=SimpleNamespace(psi_child=.4)
        driver.terminal_check(pf,SimpleNamespace(prices=np.ones(2)),terminal,
            {'endpoint':{'population_scale':2.,'price':1.}},P,{})
        np.testing.assert_array_equal(pf.state[0],[.8,1.2]);np.testing.assert_array_equal(g,[.4,.6])
        np.testing.assert_array_equal(V,[2.,3.]);self.assertEqual(pf.kw['reference_entry_flow'],.2)


    def test_render_uses_dated_rent_and_exact_17_packet(self):
        packet=self.root/'packet.pkl.gz'
        driver.dump_checkpoint(packet,{'sentinel':'unchanged'})
        record=dict(path=str(packet),sha256=driver.sha(packet),period=0,dated_rent=.123)
        calls=[]
        class Audit:
            def standard_diagnostics(self,p,out,**kwargs):
                calls.append((p,kwargs));folder=out/'standard_diagnostics';folder.mkdir()
                for i in range(17):(folder/f'{i}.png').write_bytes(b'synthetic')
        receipt=driver.render_sampled_diagnostics([record],self.root/'plots',Audit())
        self.assertEqual(calls[0][1],dict(validate_production_young=False,dated_rent=.123))
        self.assertEqual(len(receipt['sampled_dates'][0]['plots']),17)
        self.assertTrue(receipt['visual_review_pending'])

    def test_dated_rent_validation_before_any_model_work(self):
        source=(Path(driver.__file__).parent/'run_e5f_independent_numerical_audit.py').read_text()
        function=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='standard_diagnostics')
        class Model:
            def compute_markov_statistics(self,*a,**kw):raise RuntimeError('validated')
        namespace={'np':np,'model':Model()}
        exec(compile(ast.Module(body=[function],type_ignores=[]),'<pure-function-test>','exec'),namespace)
        packet={'parameters':None,'b_grid':None,'evaluation':SimpleNamespace(policy=SimpleNamespace(price=np.ones(1)),g_current=None)}
        for bad in (0.,-1.,np.inf,np.nan,[1.,2.],[[1.]]):
            with self.assertRaises(ValueError):namespace['standard_diagnostics'](packet,self.root,dated_rent=bad)
        # A valid rent passes the new validation; remaining model work deliberately stubbed.
        with self.assertRaises(AttributeError):namespace['standard_diagnostics'](packet,self.root,dated_rent=.1)

if __name__=='__main__':unittest.main()
