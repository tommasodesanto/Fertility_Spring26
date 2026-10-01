"""Zero-lifecycle checks of provenance, native binding and population/queue transport."""
import copy, importlib.util, pathlib, sys, tempfile, types, unittest
import numpy as np
HERE=pathlib.Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('floor_runtime',HERE/'floor_runtime.py')
r=importlib.util.module_from_spec(spec);sys.modules[spec.name]=r;spec.loader.exec_module(r)

class RuntimeTests(unittest.TestCase):
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

if __name__=='__main__':unittest.main()
