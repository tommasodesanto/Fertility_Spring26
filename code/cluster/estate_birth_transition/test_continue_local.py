"""Zero-native contract tests for continue_local.py."""
from __future__ import annotations

import importlib.util
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('continue_local', HERE/'continue_local.py')
mod = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(mod)


class Fixture(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.base = Path(self.tmp.name)
        self.identity = dict(reference_sha256='r'*64, engine_sha256='e'*64,
            entry_sha256='a'*64, grid_sha256='g'*64, effective_parameters_sha256='p'*64,
            source_pins={'source.py':'s'*64})
        self.plan = dict(mode='diagnostic', horizons=[24,32], identity=self.identity,
            gates=dict(final_reproduction_tolerance=1e-10, market_tolerance=2e-4, fiscal_tolerance=2e-5,
                horizon_relative_tolerance=.001),
            budget=dict(total_seconds=21480,seed_seconds=1800,candidate_seconds=7200,endpoint_seconds=1800,
                mapping_seconds=1800,path_seconds=6000,render_seconds=900,maximum_policy_calls=20000),
            fit=dict(max_evaluations=12,fertility_tolerance=.005), path=dict(max_evaluations=12), endpoint=dict(max_evaluations=48),
            initial_psi=.18, psi_bound_ratios=[.01,2.],
            target_contract=dict(rows=[{}, {}, {}, {'target':1.64575}]),
            seed={}, source_files={},handoff={})
        self.plan_path=self.base/'plan.json'; mod.write(self.plan_path,self.plan)
        self.plan_pin=dict(path=str(self.plan_path),sha256=mod.sha(self.plan_path))
        self.panel_path=self.base/'panel.json'; mod.write(self.panel_path,dict(schema='estate_a_transition_panel_v1',
            identity=self.identity,plan=self.plan_pin,guesses=[dict(index=i,psi=.1) for i in range(12)]))
        self.panel_pin=dict(path=str(self.panel_path),sha256=mod.sha(self.panel_path))
        self.failed_path=self.base/'failed.json'; mod.write(self.failed_path,dict(
            schema='current_estate_a_shock_candidate_v1',status='failed',accepted=False,
            contract=mod._contract(self.plan),contract_sha256=mod.fingerprint(mod._contract(self.plan))))
        self.failed_pin=dict(path=str(self.failed_path),sha256=mod.sha(self.failed_path))
        self.candidate=self.base/'candidate_0005'; self.candidate.mkdir()
        mod.write(self.candidate/'complete.json',dict(certified=True,mode='diagnostic',production_ready=False,
            psi=.2,model=1.7,gap=1.7-1.64575,payload=dict(models=[1.7]*4),
            horizon_comparison=dict(passed=True)))
        exp=self.base/'native.gz'; exp.write_bytes(b'checkpoint bytes')
        mod.write(self.candidate/'state_2023_checkpoint/checkpoint_receipt.json',dict(status='complete',
            identity=self.identity,psi=.2,horizon=32,candidate=5,shock_fit_complete=False,
            export=dict(path=str(exp),sha256=mod.sha(exp))))
        for h in (24,32):
            folder=self.candidate/f'horizon_{h:03d}'; folder.mkdir()
            mod.write(folder/'root.json',dict(converged=True,gates=dict(housing=True,mapping=True,
                market_replay=True,fiscal_replay=True,social_security=True),
                final_reproduction_max_abs=0.,market_reproduction_max_abs=0.,fiscal_reproduction_max_abs=0.,
                final=dict(prices=[1.]*h,fiscal_values=[.1]*h,mapping_valid=True,market_gate=True,fiscal_gate=True),
                final_jacobian=[[0.]*(2*h) for _ in range(2*h)]))
            mod.write(folder/'latest_completed.json',dict(accounting_valid=True,
                gates=dict(mass=True,policy_reproduction=True,projection=True,dated_audits=True)))
        self.args=types.SimpleNamespace(plan_pin=self.plan_pin,panel_config_pin=self.panel_pin,
            failed_receipt_pin=self.failed_pin,source_candidate=str(self.candidate),start_psi=.12,index=0,
            output=str(self.base/'continuation'))

    def tearDown(self): self.tmp.cleanup()

    def test_valid_source_authenticates_and_warm_injection_uses_root(self):
        e=mod.authenticate(self.plan_pin,self.panel_pin,self.failed_pin,self.candidate,.12,0)
        class Adapter:
            def __init__(self, runtime, plan): self._identity=plan['identity']; self.warm=[]
            def identity(self): return self._identity
            def retain_warm(self,horizon,psi,root,folder): self.warm.append((horizon,psi,root['converged']))
        class Control: NativeAdapter=Adapter
        roots={h:{'path':r['path'],'sha256':r['sha256']} for h,r in e['roots'].items()}
        base, wrapped=mod._adapter_factory(Control,roots,e['source_psi'])
        obj=wrapped(None,dict(identity=self.identity))
        self.assertIs(base,Adapter); self.assertEqual(obj.warm,[(24,.2,True),(32,.2,True)])

    def test_actual_worker5_candidate5_zero_solve_authentication(self):
        repo=HERE.parents[2]
        base=repo/'output/model/transition_readiness_v1/current_baseline_20261003'
        plan=base/'plans_v7/fit_plan.json'; config=base/'plans_v7/panel_config.json'
        failed=base/'local_workers_v7/worker_5/run/candidate_receipt.json'
        candidate=base/'local_workers_v7/worker_5/run/candidate_0005'
        if not all(p.exists() for p in (plan,config,failed,candidate/'complete.json')):
            self.skipTest('Live v7 source artifacts are not available in this checkout')
        output=self.base/'actual_preflight'
        args=[sys.executable,str(HERE/'continue_local.py'),'--prepare','--output',str(output),
            '--plan-pin',json.dumps(dict(path=str(plan),sha256=mod.sha(plan))),
            '--panel-config-pin',json.dumps(dict(path=str(config),sha256=mod.sha(config))),
            '--failed-receipt-pin',json.dumps(dict(path=str(failed),sha256=mod.sha(failed))),
            '--source-candidate',str(candidate),'--start-psi','.11999694638724082','--index','0']
        completed=subprocess.run(args,check=True,text=True,capture_output=True)
        self.assertIn('prepared',completed.stdout)
        result, rechecked=mod.load_manifest(output/'manifest.json')
        self.assertEqual(result['source_psi'],.12839959665634437)
        self.assertEqual(sorted(rechecked['roots']),[24,32])

    def test_rejects_partial_candidate_and_bad_start_or_failure_receipt(self):
        with self.assertRaises(ValueError): mod.authenticate(self.plan_pin,self.panel_pin,self.failed_pin,
            self.base/'candidate_0006',.12,0)
        with self.assertRaises(ValueError): mod.authenticate(self.plan_pin,self.panel_pin,self.failed_pin,self.candidate,99.,0)
        bad=dict(path=self.failed_pin['path'],sha256='0'*64)
        with self.assertRaises(ValueError): mod.authenticate(self.plan_pin,self.panel_pin,bad,self.candidate,.12,0)

    def test_rejects_source_target_identity_gate_root_and_replay_mismatch(self):
        cp=self.candidate/'state_2023_checkpoint/checkpoint_receipt.json'
        original=json.loads(cp.read_text())
        for key,value in [('identity',dict(wrong=True)),('psi',.3)]:
            d=dict(original);d[key]=value;mod.write(cp,d)
            with self.assertRaises(ValueError): mod.authenticate(self.plan_pin,self.panel_pin,self.failed_pin,self.candidate,.12,0)
        mod.write(cp,original)
        complete=self.candidate/'complete.json'; d=json.loads(complete.read_text());d['gap']=9.;mod.write(complete,d)
        with self.assertRaises(ValueError): mod.authenticate(self.plan_pin,self.panel_pin,self.failed_pin,self.candidate,.12,0)
        d=json.loads(complete.read_text());d['gap']=1.7-1.64575;mod.write(complete,d)
        rootp=self.candidate/'horizon_024/root.json'; original_root=json.loads(rootp.read_text())
        for edit in (dict(converged=False),dict(final_reproduction_max_abs=1.),
                     dict(gates=dict(housing=False)),dict(final_jacobian=[[float('nan')]]),
                     dict(final=dict(prices=[1.],fiscal_values=[.1]*24))):
            d=dict(original_root);d.update(edit)
            rootp.write_text(json.dumps(d,allow_nan=True))
            with self.assertRaises(ValueError): mod.authenticate(self.plan_pin,self.panel_pin,self.failed_pin,self.candidate,.12,0)
        mod.write(rootp,original_root)

    def test_smoke_verifier_requires_both_horizons_and_replay_gates(self):
        out=self.base/'smoke'; candidate=out/'candidate_0001'
        for h in (24,32):
            folder=candidate/f'horizon_{h:03d}'; folder.mkdir(parents=True)
            mod.write(folder/'root.json',dict(converged=True,gates=dict(fiscal_replay=True,housing=True,
                mapping=True,market_replay=True,social_security=True),final_reproduction_max_abs=0.,
                final=dict(mapping_valid=True,market_gate=True,fiscal_gate=True,
                    prices=[1.]*h,fiscal_values=[.1]*h)))
            mod.write(folder/'latest_completed.json',dict(accounting_valid=True,
                gates=dict(mass=True,policy_reproduction=True,projection=True,dated_audits=True)))
        reply=dict(identity=self.identity,psi=.2,horizon=32,root_pass=True,replay_pass=True,
            accounting_valid=True,stationary_pass=True,horizon_comparison=dict(passed=True))
        result=dict(certified=True,payload=dict(models=[1.,2.,3.,4.]),horizon_comparison=dict(passed=True))
        mod.verify_smoke_outputs(out,self.plan,.2,[1.,2.,3.,4.],reply,result)
        result['payload']['models'][0] = 1.0011
        with self.assertRaises(ValueError):
            mod.verify_smoke_outputs(out,self.plan,.2,[1.,2.,3.,4.],reply,result)
        result['payload']['models'][0] = 1.
        bad=mod.read(candidate/'horizon_032/root.json'); bad['gates']['market_replay']=False
        mod.write(candidate/'horizon_032/root.json',bad)
        with self.assertRaises(ValueError):
            mod.verify_smoke_outputs(out,self.plan,.2,[1.,2.,3.,4.],reply,result)

    def test_mock_smoke_leaves_prepared_manifest_immutable(self):
        m=mod.create_manifest(self.args); manifest_path=Path(m['output'])/'manifest.json'
        before=mod.sha(manifest_path); out=self.base/'mock-smoke'; identity=self.identity
        export=self.base/'smoke-export.gz'; export.write_bytes(b'x')
        class Adapter:
            def __init__(self, runtime, plan): self.plan=plan
            def identity(self): return self.plan['identity']
            def retain_warm(self,*args): pass
        class Runner:
            def __init__(self, plan, adapter, folder): self.folder=Path(folder);self.policy_calls=9;self.last=None
            def prepare(self): pass
            def evaluate(self,psi):
                for h in (24,32):
                    f=self.folder/'candidate_0001'/f'horizon_{h:03d}'
                    mod.write(f/'root.json',dict(converged=True,gates=dict(fiscal_replay=True,housing=True,
                        mapping=True,market_replay=True,social_security=True),final_reproduction_max_abs=0.,
                        final=dict(mapping_valid=True,market_gate=True,fiscal_gate=True,
                            prices=[1.]*h,fiscal_values=[.1]*h)))
                    mod.write(f/'latest_completed.json',dict(accounting_valid=True,
                        gates=dict(mass=True,policy_reproduction=True,projection=True,dated_audits=True)))
                mod.write(self.folder/'candidate_0001/complete.json',dict(certified=True,psi=psi,
                    horizon_comparison=dict(passed=True)))
                mod.write(self.folder/'candidate_0001/state_2023_checkpoint/checkpoint_receipt.json',
                    dict(status='complete',candidate=1,horizon=32,psi=psi,identity=identity,
                        export=dict(path=str(export),sha256=mod.sha(export))))
                self.last=dict(identity=identity,psi=psi,horizon=32,root_pass=True,replay_pass=True,
                    accounting_valid=True,stationary_pass=True,horizon_comparison=dict(passed=True))
                return dict(certified=True,payload=dict(models=mod.read(m['complete_pin']['path'])['payload']['models']),
                    horizon_comparison=dict(passed=True))
        fake_driver=types.SimpleNamespace(controller=types.SimpleNamespace(NativeAdapter=Adapter,Controller=Runner),
            preflight=lambda plan:None,
            runtime_module=lambda:types.SimpleNamespace(CurrentEstateARuntime=types.SimpleNamespace(
                from_handoff=lambda handoff,folder:object())))
        import sys
        key='experiments.birth_count_choice.transition';old=sys.modules.get(key);sys.modules[key]=fake_driver
        try: receipt=mod.smoke(manifest_path,out)
        finally:
            if old is None:sys.modules.pop(key,None)
            else:sys.modules[key]=old
        self.assertTrue(receipt['passed']);self.assertEqual(mod.sha(manifest_path),before)
    def test_manifest_is_exclusive_and_detects_changed_root_pin(self):
        m=mod.create_manifest(self.args)
        self.assertEqual(m['new_psi'],.12); self.assertIn('continuation_identity',mod.read(m['panel_config_path']))
        with self.assertRaises(FileExistsError): mod.create_manifest(self.args)
        root=Path(m['roots'][24]['path']); root.write_text('{}')
        with self.assertRaises(ValueError): mod.load_manifest(Path(m['output'])/'manifest.json')

    def test_run_requires_smoke_and_delegates_to_unchanged_panel_once(self):
        m=mod.create_manifest(self.args)
        smcp=self.base/'smoke_checkpoint.json'; mod.write(smcp,dict(status='complete',identity=self.identity,
            psi=.2,candidate=1,horizon=32,export=dict(path=str(self.base/'native.gz'),
            sha256=mod.sha(self.base/'native.gz'))))
        smoke=self.base/'smoke.json'; mod.write(smoke,dict(schema='estate_birth_transition_warm_smoke_v1',
            passed=True,continuation_identity=m['continuation_identity'],wrapper_sha256=m['wrapper']['sha256'],
            source_psi=m['source_psi'],identity=self.identity,horizons=[24,32],
            complete_pin=m['complete_pin'],checkpoint_pin=dict(path=str(smcp),sha256=mod.sha(smcp)),
            fertility_absolute_gaps=[0.,0.,0.,0.],
            root_pins={str(h):dict(path=m['roots'][h]['path'],sha256=m['roots'][h]['sha256']) for h in (24,32)},
            mapped_pins={str(h):dict(path=str(self.candidate/f'horizon_{h:03d}/latest_completed.json'),
                sha256=mod.sha(self.candidate/f'horizon_{h:03d}/latest_completed.json')) for h in (24,32)}))
        manifest_path=Path(m['output'])/'manifest.json'
        manifest_sha=mod.sha(manifest_path)
        calls=[]
        class Adapter: pass
        fake_control=types.SimpleNamespace(NativeAdapter=Adapter)
        fake=types.SimpleNamespace(control=fake_control,evaluate_candidate=lambda *a,**k:(calls.append((a,k)) or {'ok':True}))
        import sys
        old=sys.modules.get('experiments.birth_count_choice.transition_panel')
        sys.modules['experiments.birth_count_choice.transition_panel']=fake
        try:
            result=mod.run(manifest_path,self.base/'fit-output',dict(path=str(smoke),sha256=mod.sha(smoke)))
        finally:
            if old is None: sys.modules.pop('experiments.birth_count_choice.transition_panel',None)
            else: sys.modules['experiments.birth_count_choice.transition_panel']=old
        self.assertEqual(result,{'ok':True}); self.assertEqual(len(calls),1)
        self.assertEqual(mod.sha(manifest_path),manifest_sha)
        self.assertEqual(calls[0][0][1:3],(.12,0))
        with self.assertRaises(ValueError): mod.run(manifest_path,self.base/'bad-fit',dict(path=str(smoke),sha256='0'*64))

    def test_cli_dispatch_accepts_separate_fresh_smoke_and_fit_outputs(self):
        m=mod.create_manifest(self.args); manifest=Path(m['output'])/'manifest.json'
        smoke_out=self.base/'smoke_job'/'run'; fit_out=self.base/'fit_job'/'run'
        smoke_out.parent.mkdir();fit_out.parent.mkdir()
        smoke_file=self.base/'explicit_smoke.json';mod.write(smoke_file,dict(passed=True))
        smoke_pin=dict(path=str(smoke_file),sha256=mod.sha(smoke_file));calls=[]
        old_smoke,old_run=mod.smoke,mod.run
        mod.smoke=lambda *a:(calls.append(('smoke',a)) or dict(passed=True))
        mod.run=lambda *a:(calls.append(('run',a)) or dict(accepted=True))
        import sys
        old_argv=sys.argv[:]
        try:
            sys.argv=['continue_local.py','--manifest',str(manifest),'--smoke','--output',str(smoke_out)]
            with contextlib.redirect_stdout(io.StringIO()):mod.main()
            sys.argv=['continue_local.py','--manifest',str(manifest),'--run','--output',str(fit_out),
                      '--smoke-receipt-pin',json.dumps(smoke_pin)]
            with contextlib.redirect_stdout(io.StringIO()):mod.main()
        finally:
            sys.argv=old_argv;mod.smoke,mod.run=old_smoke,old_run
        self.assertEqual(calls,[('smoke',(str(manifest),str(smoke_out))),
            ('run',(str(manifest),str(fit_out),smoke_pin))])
        self.assertFalse(smoke_out.exists());self.assertFalse(fit_out.exists())


if __name__ == '__main__': unittest.main()
