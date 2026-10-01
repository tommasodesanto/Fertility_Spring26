"""Zero-lifecycle checks of provenance, native binding and population/queue transport."""
import copy, importlib.util, pathlib, sys, tempfile, types, unittest
import numpy as np
HERE=pathlib.Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('floor_runtime',HERE/'floor_runtime.py')
r=importlib.util.module_from_spec(spec);sys.modules[spec.name]=r;spec.loader.exec_module(r)

class RuntimeTests(unittest.TestCase):
    def test_normalized_h0_binding_is_authenticated_once_and_legacy_is_unchanged(self):
        marker=dict(schema='normalized_n0_fixed_h0_v1',N0=1.,H0_source='authenticated_parameter_table',counterfactual_H0_fixed=True)
        handoff=dict(normalized_housing_contract=marker)
        closure=dict(normalized_population=1.,population_scale=1.,H0_derived=6.85157529,H0_bounds=[.2,80.])
        P=types.SimpleNamespace(H0=np.array([7.]),N_target=1.)
        original=P.H0
        self.assertIs(r.bind_normalized_reference_housing(P,{},dict(population_scale=.92),{'H0':7.}),P)
        self.assertIs(P.H0,original)
        r.bind_normalized_reference_housing(P,handoff,closure,{'H0':6.85157529})
        np.testing.assert_array_equal(P.H0,[6.85157529])
        for key,value in (('population_scale',.92),('normalized_population',.92),('H0_derived',7.)):
            bad=dict(closure);bad[key]=value
            with self.assertRaises(RuntimeError):r.bind_normalized_reference_housing(P,handoff,bad,{'H0':6.85157529})
        bad=dict(closure,H0_derived=81.)
        with self.assertRaisesRegex(RuntimeError,'H0 differs'):r.validate_normalized_housing_report(handoff,bad,{'H0':81.})
        P.N_target=.92
        with self.assertRaisesRegex(RuntimeError,'N_target'):r.bind_normalized_reference_housing(P,handoff,closure,{'H0':6.85157529})
        with self.assertRaisesRegex(RuntimeError,'explicit normalized'):r.validate_normalized_housing_report({},closure,{'H0':6.85157529})
        for key,value in (('N0',True),('counterfactual_H0_fixed',False),('H0_source','profile_each_shock')):
            bad=copy.deepcopy(handoff);bad['normalized_housing_contract'][key]=value
            with self.assertRaisesRegex(RuntimeError,'contract differs'):r.normalized_housing_contract(bad)

    def test_normalized_reference_comparison_preserves_all_numeric_values_and_plot_hashes(self):
        import csv
        def table(folder,name,rows):
            with (folder/name).open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        with tempfile.TemporaryDirectory() as temp:
            a,b=[pathlib.Path(temp)/name for name in ('reference','fresh')]
            for folder in (a,b):
                (folder/'standard_diagnostics').mkdir(parents=True)
                for i in range(17):(folder/'standard_diagnostics'/f'{i}.png').write_bytes(str(i).encode())
                table(folder,'target_fit.csv',[dict(moment=str(i),target=1.,model=1.,gap=0.,weight=1.,loss_contribution=0.) for i in range(14)])
                table(folder,'parameters.csv',[dict(parameter='H0' if i==0 else str(i),estimate=6.8,status=('derived calibrated housing supply coefficient at N0=1' if folder==a else 'fixed population scale; not identified by per-household targets') if i==0 else 'fixed',lower='.2',upper='80',near_bound='False') for i in range(31)])
            r.write(a/'closure.json',dict(price=.7,population_scale=1.,normalized_population=1.,H0_derived=6.8,H0_bounds=[.2,80.],housing_supply_coefficient_role='derived calibrated coefficient at N0=1',standard_plot_supply_units='physical housing supply at normalized N0=1'))
            r.write(b/'closure.json',dict(price=.7,population_scale=1.,standard_plot_supply_units='physical supply divided by endogenous household population'))
            r.compare_normalized_reference_reports(a,b)
            with (b/'parameters.csv').open() as stream:self.assertIn('fixed calibrated',list(csv.DictReader(stream))[0]['status'])
            with (b/'target_fit.csv').open() as stream:rows=list(csv.DictReader(stream))
            rows[2]['model']='1.01';table(b,'target_fit.csv',rows)
            with self.assertRaisesRegex(RuntimeError,'numeric value'):r.compare_normalized_reference_reports(a,b)
            rows[2]['model']='1.0';table(b,'target_fit.csv',rows)
            with (b/'parameters.csv').open() as stream:params=list(csv.DictReader(stream))
            params[0]['status']='unexpected H0 metadata';table(b,'parameters.csv',params)
            with self.assertRaisesRegex(RuntimeError,'H0 role'):r.compare_normalized_reference_reports(a,b)
            params[0]['status']='fixed population scale; not identified by per-household targets';table(b,'parameters.csv',params)
            (b/'standard_diagnostics/0.png').write_bytes(b'changed')
            with self.assertRaisesRegex(RuntimeError,'standard plots'):r.compare_normalized_reference_reports(a,b)

    def test_real_handoff_and_complete_source_inventory(self):
        h=r.ROOT/'output/model/transition_readiness_v1/current_floor_handoff/handoff.json'
        handoff,pins,report=r.authenticate_handoff(dict(path=str(h),sha256=r.sha(h)))
        self.assertEqual(len(pins),handoff['checkpoint_and_sources'].get('source_pin_count',111));self.assertEqual(handoff['tables']['parameters_csv']['rows'],31)
    def test_wrong_handoff_hash_rejected(self):
        h=r.ROOT/'output/model/transition_readiness_v1/current_floor_handoff/handoff.json'
        with self.assertRaisesRegex(RuntimeError,'Handoff hash'):r.authenticate_handoff(dict(path=str(h),sha256='0'*64))
    def test_actual_population_and_both_native_queue_inputs(self):
        rt=r.FloorRuntime();seen=[]
        def initial(g,entry,raw,P,conversion):
            seen.append((g.copy(),entry,raw,conversion));return types.SimpleNamespace(g_pre=g,scheduled_entries=('due16',entry,'due20',entry),scheduled_raw_entries=('due16',raw,'due20',raw))
        rt.pf=types.SimpleNamespace(stationary_initial_state=initial)
        g=np.arange(1,9,dtype=float).reshape(1,1,1,2,2,2);g/=g.sum()
        packet=dict(stationary_g_pre=g,evaluation=types.SimpleNamespace(births=.02),parameters=object())
        state=rt.stationary_state(packet,.9222615667)
        self.assertAlmostEqual(state.g_pre.sum(),.9222615667)
        np.testing.assert_array_equal(seen[0][0],g*.9222615667)
        self.assertEqual(seen[0][1],float((g*.9222615667)[:,:,:,0].sum()))
        self.assertEqual(seen[0][2],.02*.9222615667);self.assertEqual(seen[0][3],1/2.1)
        self.assertIn('due20',state.scheduled_entries);self.assertIn('due20',state.scheduled_raw_entries)
    def test_mapping_fails_before_selected_reconstruction(self):
        rt=r.FloorRuntime();rt.reference_verified=False;rt.packet=None
        with self.assertRaisesRegex(RuntimeError,'reconstruction'):rt.mapping(None,None,[1.],[1.],[.1],'.')
    def test_export_copies_actual_clock_state_without_normalization(self):
        rt=r.FloorRuntime();rt.grid=np.array([0.,1.]);rt.identity=lambda:{'test_identity':'pinned'}
        state=types.SimpleNamespace(g_pre=np.array([.3,.4]),scheduled_entries=types.SimpleNamespace(due16=[.1],due20=[.2]),scheduled_raw_entries=types.SimpleNamespace(due16=[.2],due20=[.4]))
        native=types.SimpleNamespace(dated_states={4:dict(state=state,parameters=types.SimpleNamespace(psi_child=.17))},floor_runtime_paths=dict(start_year=2007,prices=np.arange(7)+1.,pensions=np.ones(7),psi_path=np.full(7,.17)),values=[np.ones((2,2))*i for i in range(8)])
        with tempfile.TemporaryDirectory() as folder:
            receipt=rt.export_2023(dict(native_reply=native),folder)
            import gzip,pickle
            with gzip.open(receipt['path'],'rb') as f:packet=pickle.load(f)
            self.assertEqual(receipt['period_index'],4);self.assertEqual(receipt['queue_lags'],[16,20]);self.assertEqual(packet['initial_population'],.7)
            np.testing.assert_array_equal(packet['initial_state'].g_pre,state.g_pre)
            self.assertEqual(packet['initial_state'].scheduled_raw_entries.due20,[.4])
            np.testing.assert_array_equal(packet['current_2023_V'],native.values[4])
            np.testing.assert_array_equal(packet['continuation_V'],native.values[5])
            self.assertEqual(packet['continuation_calendar_year'],2027)

    def test_resume_rejects_identity_before_unpickling(self):
        rt=r.FloorRuntime();rt.identity=lambda:{'reference_sha256':'expected'}
        with tempfile.TemporaryDirectory() as folder:
            path=pathlib.Path(folder)/'receipt.json'
            r.write(path,dict(schema='current_floor_reference_reconstruction_v1',status='passed',identity={'reference_sha256':'other'}))
            with self.assertRaisesRegex(RuntimeError,'numerical identity'):
                rt.load_reconstructed_reference(dict(path=str(path),sha256=r.sha(path)),folder)
    def test_resume_rejects_changed_checkpoint_bytes(self):
        rt=r.FloorRuntime();rt.identity=lambda:{'reference_sha256':'expected'}
        with tempfile.TemporaryDirectory() as folder:
            checkpoint=pathlib.Path(folder)/'state.pkl.gz';checkpoint.write_bytes(b'not a native packet')
            receipt=pathlib.Path(folder)/'receipt.json'
            r.write(receipt,dict(schema='current_floor_reference_reconstruction_v1',status='passed',identity=rt.identity(),lifecycle_calls=2,policy_calls=2,checkpoint=dict(path=str(checkpoint),sha256='0'*64),checkpoint_sha256='0'*64))
            with self.assertRaisesRegex(RuntimeError,'checkpoint hash'):
                rt.load_reconstructed_reference(dict(path=str(receipt),sha256=r.sha(receipt)),folder)

    def test_mapping_passes_exact_cache_and_actual_state_to_native_callback(self):
        import contextlib
        from collections import namedtuple
        State=namedtuple('State','g_pre scheduled_entries scheduled_raw_entries')
        rt=r.FloorRuntime();rt.reference_verified=True;rt.packet={'supply_rule':object()};rt.P=object();rt.grid=np.array([0.,1.]);rt.ctx={'prepared':None}
        rt.initial_state=State(np.array([.3,.5]),[.1,.2],[.2,.4]);seen={}
        rt.native_bindings=lambda:contextlib.nullcontext()
        rt.model=types.SimpleNamespace(solve_bellman_full_markov_income=lambda *a,**k:None)
        rt.pf=types.SimpleNamespace(PFInitialState=State,calendar=types.SimpleNamespace(evaluate_period=lambda *a,**k:types.SimpleNamespace(g_pre=a[1])),transition=types.SimpleNamespace(advance_adult_entry_clock=lambda q,*a,**k:(q[0],q[1:]+[0.])))
        def callback(packet,prepared,terminal,endpoint,prices,pensions,psi_path,housing,output,cache_bytes,**kwargs):
            seen.update(cache_bytes=cache_bytes,state=kwargs['initial_state'])
            rt.pf.calendar.evaluate_period([1.],kwargs['initial_state'].g_pre,rt.P,rt.grid,None,None)
            rt.pf.transition.advance_adult_entry_clock(kwargs['initial_state'].scheduled_entries,0.,1/2.1,None)
            rt.pf.transition.advance_adult_entry_clock(kwargs['initial_state'].scheduled_raw_entries,0.,1/2.1,None)
            return types.SimpleNamespace(),dict(gates={'audits':True})
        rt.scaffold=types.SimpleNamespace(mapping=callback,supply_rule=lambda *a:None)
        with tempfile.TemporaryDirectory() as folder:reply,record=rt.mapping({},dict(price=1.),[1.],[1.],[.17],folder)
        self.assertEqual(seen['cache_bytes'],64*1024**3)
        np.testing.assert_array_equal(seen['state'].g_pre,rt.initial_state.g_pre)
        self.assertEqual(record['policy_calls'],0);self.assertEqual(record['exact_policy_cache_bytes'],64*1024**3)
        self.assertEqual(reply.dated_states[0]['state'].scheduled_raw_entries,[.2,.4])

    def test_old_cache_retry_cannot_start_native_call_after_deadline(self):
        rt=r.FloorRuntime();started=[]
        def native():
            rt._guard_native_call();started.append(True)
        def old_cache_fallback():
            try:return native()
            except Exception:return native()
        with rt.native_budget(r.time.monotonic()-1.,10) as budget:
            with self.assertRaisesRegex(TimeoutError,'before solve'):old_cache_fallback()
        self.assertEqual(started,[]);self.assertEqual(budget['used'],0)
    def test_call_allowance_checks_before_retry_and_nested_limits(self):
        rt=r.FloorRuntime();started=[]
        def native():
            rt._guard_native_call();started.append(True)
        with rt.native_budget(r.time.monotonic()+60.,5) as outer:
            with rt.native_budget(r.time.monotonic()+60.,1) as inner:
                native()
                for _ in range(2):
                    with self.assertRaisesRegex(RuntimeError,'allowance exhausted'):native()
        self.assertEqual(len(started),1);self.assertEqual(inner['used'],1);self.assertEqual(outer['used'],1)

    def test_actual_native_api_and_interp_boundary_values_without_model_solve(self):
        result=r.FloorRuntime.native_api_preflight()
        self.assertEqual(result['policy_calls'],0)
        self.assertEqual(result['callbacks']['interp_indices'],'small_credit_lab.engine.utils')
        from small_credit_lab.engine import solver,utils
        existed=hasattr(solver,'interp_indices');old=getattr(solver,'interp_indices',None)
        rt=r.FloorRuntime();rt.model=solver
        with rt.native_api_bindings():
            self.assertIs(solver.interp_indices,utils.interp_indices)
            indices,weights=solver.interp_indices(np.array([0.,2.,5.]),np.array([-1.,0.,1.,2.,4.,5.,8.]))
            np.testing.assert_array_equal(indices,[0,0,0,1,1,1,1])
            np.testing.assert_allclose(weights,[0.,0.,.5,0.,2/3,1.,1.],rtol=0,atol=1e-15)
        self.assertEqual(hasattr(solver,'interp_indices'),existed)
        if existed:self.assertIs(solver.interp_indices,old)
    def test_native_counter_survives_failure_after_call_entry(self):
        rt=r.FloorRuntime();rt.total_native_calls=0
        def completed_solve_then_observer_failure():
            rt._guard_native_call()
            raise AttributeError('post-solve native report helper missing')
        with self.assertRaises(AttributeError):completed_solve_then_observer_failure()
        self.assertEqual(rt.total_native_calls,1)
        with rt.native_budget(r.time.monotonic()+60.,0):
            with self.assertRaisesRegex(RuntimeError,'allowance'):rt._guard_native_call()
        self.assertEqual(rt.total_native_calls,1)

    def test_actual_indirect_population_helpers_are_native_and_scoped(self):
        r.FloorRuntime.native_api_preflight()
        from small_credit_lab.engine import solver,utils,parameters,household
        rt=r.FloorRuntime();rt.model=solver
        names=('independent_child_maturation_active','get_fecundity_by_age','readiness_settled_state','parent_age_maturation_active')
        old={name:(hasattr(solver,name),getattr(solver,name,None)) for name in names}
        with rt.native_api_bindings() as api:
            for name in names:
                self.assertIs(getattr(solver,name),getattr(parameters,name))
                self.assertEqual(api[name],'small_credit_lab.engine.parameters')
            self.assertIs(solver.birth_destination_child_state,household.birth_destination_child_state)
        for name,(existed,value) in old.items():
            self.assertEqual(hasattr(solver,name),existed)
            if existed:self.assertIs(getattr(solver,name),value)
    def test_completed_native_solve_evidence_saved_without_claiming_verification(self):
        import gzip,pickle
        rt=r.FloorRuntime();rt.total_native_calls=1;rt.identity=lambda:{'native_source':'authenticated'}
        P=types.SimpleNamespace(psi_child=.17);solution=types.SimpleNamespace(V=np.array([1.,2.]))
        with tempfile.TemporaryDirectory() as folder:
            receipt=rt.save_unverified_native_solve(P,np.array([0.,1.]),types.SimpleNamespace(),solution,.72,folder,1)
            self.assertEqual(receipt['policy_calls'],1);self.assertFalse(receipt['reference_verified']);self.assertFalse(receipt['production_ready'])
            self.assertEqual(r.sha(receipt['checkpoint']['path']),receipt['checkpoint']['sha256'])
            with gzip.open(receipt['checkpoint']['path'],'rb') as stream:packet=pickle.load(stream)
            np.testing.assert_array_equal(packet['solution'].V,solution.V)
            self.assertEqual(packet['stage'],'native_solve_completed_reconstruction_pending')
            with self.assertRaisesRegex(RuntimeError,'not verified'):
                rt.load_reconstructed_reference(dict(path=str(pathlib.Path(folder)/'native_solve_unverified.json'),sha256=r.sha(pathlib.Path(folder)/'native_solve_unverified.json')),folder)

    def test_real_population_callback_globals_resolve_native_and_block_reconfigure(self):
        r.FloorRuntime.native_api_preflight()
        sys.path[:0]=[str(r.ROOT/'code/model/tools'),str(r.ROOT/'code/model')]
        import run_e5f_perfect_foresight_transition as pf
        from small_credit_lab.engine import solver
        rt=r.FloorRuntime();rt.model=solver;rt.pf=pf
        old=pf.calendar.model;configure=pf.transition.configure_sequential_model
        rt.rt={'model':old,'primitive':types.SimpleNamespace(model=old),'audit':types.SimpleNamespace(model=old)}
        with rt.native_bindings():
            self.assertIs(pf.transition.apply_sequential_fertility.__globals__['calendar'].model,solver)
            self.assertIs(pf.transition.advance_sequential_calendar_distribution.__globals__['calendar'].model,solver)
            with self.assertRaisesRegex(RuntimeError,'reconfiguration forbidden'):pf.transition.configure_sequential_model()
            self.assertIs(pf.calendar.model,solver)
        self.assertIs(pf.calendar.model,old);self.assertIs(pf.transition.configure_sequential_model,configure)

    def test_observer_adapter_changes_only_dependencies_with_genuine_native_helpers(self):
        r.FloorRuntime.native_api_preflight()
        from small_credit_lab.engine import solver
        rt=r.FloorRuntime();rt.model=solver
        sys.path[:0]=[str(r.ROOT/'code/model/tools'),str(r.ROOT/'code/model')]
        housing=r.load('test_original_housing_observer',r.ROOT/'code/model/tools/e5f_initial_housing_observer.py')
        recent=r.load('test_original_recent_observer',r.ROOT/'code/model/tools/e5f_recent_parent_flow_observer.py')
        def retained_tail(*args,**kwargs):
            kwargs['diagnostic_allow_retained_dead_tail']=True
            return recent.observe_recent_parent_flow(*args,**kwargs)
        rt.rt={'observe_initial_housing_wealth':housing.observe_initial_housing_wealth,'observe_recent_parent_flow':retained_tail}
        with tempfile.TemporaryDirectory() as folder:
            name=solver.__name__;rt.install_observer_adapters(folder)
            receipt=r.read(pathlib.Path(folder)/'transformation.json')
            self.assertEqual(solver.__name__,name)
            self.assertEqual(len(receipt['observers']['housing']['substitutions']),2)
            self.assertEqual(len(receipt['observers']['recent_parent']['substitutions']),1)
            self.assertTrue(all(v['computation_ast_unchanged'] for v in receipt['observers'].values()))
            self.assertEqual(receipt['constants'],dict(DEAD_MASS_TOL=1e-12,DEAD_VALUE_CUTOFF=-1e9))
            function=rt.rt['observe_initial_housing_wealth']
            self.assertIs(function.__globals__['_CURRENT_FLOOR_NATIVE_MODEL'],solver)
            for observer in receipt['observers'].values():
                self.assertEqual(r.sha(observer['adapter']['path']),observer['adapter']['sha256'])
                self.assertEqual(r.sha(observer['original']['path']),observer['original']['sha256'])

if __name__=='__main__':unittest.main()
