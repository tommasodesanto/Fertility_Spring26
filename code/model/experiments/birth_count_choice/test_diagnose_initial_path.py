"""Cheap structural tests; no model constructor or solve is imported."""
import contextlib
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import pickle
import sys
import types
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('diagnose_initial_path', HERE / 'diagnose_initial_path.py')
d = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(d)


class FakeArray:
    def __init__(self, values):
        self.values = values if isinstance(values, list) else [values]
        self.dtype = SimpleNamespace(kind='f')
        self.shape = (len(self.values),)
    def tobytes(self): return repr(self.values).encode()
    def __eq__(self, other): return isinstance(other, FakeArray) and self.values == other.values


class ConfigTests(unittest.TestCase):
    def test_analytic_path_has_24_dates(self):
        q = d.expected_paths()
        self.assertEqual(len(q), 24)
        self.assertAlmostEqual(q[0], .7572901438425542)
        self.assertNotEqual(q[-1], d.QT)

    def test_bad_config_rejected_before_import(self):
        with self.assertRaisesRegex(ValueError, 'Config keys differ'):
            d.validate_config({'package_root': '/not-a-package'})

    def test_wrong_path_length_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / 'execution_smoke_v5'
            source = root / 'frozen/source/code/model/experiments/birth_count_choice'
            source.mkdir(parents=True)
            root = root.resolve(); source = source.resolve()
            manifest = root / 'inputs/fit_manifest.json'; manifest.parent.mkdir()
            manifest.write_text('{}')
            receipt = root / 'results/fit_fit_v1/run/reference/reference_reconstruction.json'
            receipt.parent.mkdir(parents=True); receipt.write_text(json.dumps({'checkpoint_sha256': d.REFERENCE_SHA}))
            endpoint = root / 'results/fit_fit_v1/run/stage1/point_013'
            endpoint.mkdir(parents=True)
            pins = {}
            for name, filename in [('checkpoint','native_solve_unverified.pkl.gz'),
                                   ('raw_receipt','native_solve_unverified.json'),
                                   ('stationary','stationary.json'),
                                   ('root','root.json'),
                                   ('latest_completed','latest_completed.json'),
                                   ('one_step_record','one_step/native_record.json')]:
                path = endpoint / filename
                path.parent.mkdir(parents=True, exist_ok=True); path.write_text('x')
                pins[name] = {'path': str(path), 'sha256': d.digest(path)}
            pins['checkpoint']['sha256'] = d.TERMINAL_SHA
            config = dict(package_root=str(root), fit_manifest={'path':str(manifest),'sha256':d.MANIFEST_SHA},
                          driver={'path':str(Path(d.__file__).resolve()),'sha256':d.digest(d.__file__)},
                          reference_receipt={'path':str(receipt),'sha256':d.digest(receipt)},
                          endpoint_pins=pins, q_path=d.expected_paths()[:-1],
                          pension_path=[d.PENSION]*24, psi_path=[d.PSI]*24)
            # Synthetic files deliberately lack real immutable hashes; isolate the path contract.
            with mock.patch.object(d, 'pinned', side_effect=lambda item, root=None: Path(item['path']).resolve()):
                with mock.patch.object(d, 'digest', side_effect=lambda path: d.MANIFEST_SHA if Path(path)==manifest else 'x'):
                    with self.assertRaisesRegex(ValueError, 'Path length differs'):
                        d.validate_config(config)


class OneMapTests(unittest.TestCase):
    def fake_run(self, fail=False):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / 'execution_smoke_v5'
            source = root / 'frozen/source/code/model/experiments/birth_count_choice'
            source.mkdir(parents=True)
            root = root.resolve(); source = source.resolve()
            for name in ('two_shock.py','two_shock_runtime.py'):
                (source/name).write_text('# fake')
            manifest = root / 'manifest.json'; manifest.write_text(json.dumps({'gates':{'market_tolerance':2e-4,'fiscal_tolerance':2e-5}}))
            config = dict(fit_manifest={'path':str(manifest),'sha256':d.digest(manifest)},
                          reference_receipt={'path':'unused','sha256':'unused'},
                          endpoint_pins={name:{'path':'unused','sha256':'unused'} for name in d.PIN_NAMES},
                          q_path=d.expected_paths(), pension_path=[d.PENSION]*24, psi_path=[d.PSI]*24)
            state = SimpleNamespace(g_pre=[1], scheduled_entries=[1], scheduled_raw_entries=[1])
            calls = []
            class RT:
                total_native_calls = 0
                initial_state = state
                grid = [0]
                P = SimpleNamespace(psi_child=0)
                packet = {'parameters':SimpleNamespace(psi_child=0)}
                pf = SimpleNamespace(birth_queue_values=lambda x:x)
                def identity(self): return {'fake':True}
                def load_reconstructed_reference(self, pin, folder): return None
                def _guard_native_call(self): self.total_native_calls += 1
                @contextlib.contextmanager
                def native_budget(self, deadline, remaining):
                    calls.append(('budget', remaining)); yield
                def mapping(self, terminal, endpoint, q, b, psi, folder, **kwargs):
                    calls.append(('map', len(q), len(b), len(psi), kwargs['start_year']))
                    self._guard_native_call()
                    if fail: raise RuntimeError('native failure')
                    record = dict(accounting_valid=True, gates={'dated_audits':True},
                                  policy_calls=1, rows=[{}]*24, market_residual=[.3]*24,
                                  fiscal_residual=[0.]*24)
                    return SimpleNamespace(dated_states={i:None for i in range(24)}),record
            rt = RT()
            fake_driver = SimpleNamespace(__file__=str(source/'two_shock.py'),
                                          preflight=lambda plan:{'native_calls':0},
                                          state_hash=lambda state, values:d.INITIAL_SHA)
            fake_native = SimpleNamespace(__file__=str(source/'two_shock_runtime.py'),
                                          retained=SimpleNamespace(watchdog=lambda deadline:contextlib.nullcontext(),
                                                                   stationary_mapping_valid=lambda record:True),
                                          NativeRuntime=lambda plan, output:SimpleNamespace(rt=rt))
            saved = dict(solution=SimpleNamespace(V=[1]), b_grid=[0],
                         parameters=SimpleNamespace(psi_child=d.PSI))
            fake_array = SimpleNamespace(tobytes=lambda:b'fake')
            patches = [mock.patch.object(d,'INITIAL_SHA',hashlib.sha256(b'fake').hexdigest()),
                       mock.patch.object(d,'TERMINAL_V_SHA',hashlib.sha256(b'fake').hexdigest()),
                       mock.patch.object(d,'validate_config',return_value=(root,source,{})),
                       mock.patch.object(d,'validate_endpoint',return_value=(saved,dict(latest_completed={'population_scale':1},stationary={}))),
                       mock.patch.object(d,'bind_reference_cache_outputs',return_value={}),
                       mock.patch.object(d,'restore_reference_with_bridge_redirect',return_value={}),
                       mock.patch.object(d,'same_public_parameters'),
                       mock.patch.dict(sys.modules,{'two_shock':fake_driver,'two_shock_runtime':fake_native,
                                                    'numpy':SimpleNamespace(array_equal=lambda a,b:a==b,
                                                                            asarray=lambda x:fake_array)})]
            with contextlib.ExitStack() as stack:
                for patch in patches: stack.enter_context(patch)
                out = root / 'diagnostic'
                if fail:
                    with self.assertRaisesRegex(RuntimeError,'native failure'):
                        d.run(config,out)
                    self.assertEqual(json.loads((out/'failure.json').read_text())['phase'],'mapping')
                else:
                    result = d.run(config,out)
                    self.assertEqual(result['status'],'mapping_completed')
                    self.assertFalse(result['empirical_fitted'])
                    self.assertEqual(result['actual_native_calls'],1)
                self.assertEqual(calls,[('budget',64),('map',24,24,24,2007)])

    def test_exactly_one_native_map(self): self.fake_run()
    def test_failure_receipt_retained(self): self.fake_run(fail=True)


class ReferenceBindingTests(unittest.TestCase):
    def test_only_five_outputs_copied_and_primitive_change_rejected(self):
        Array = FakeArray
        fake_numpy = SimpleNamespace(asarray=lambda value:value if isinstance(value,Array) else Array(value),
                                     array_equal=lambda a,b:a==b,
                                     isfinite=lambda value:SimpleNamespace(all=lambda:True))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            saved = SimpleNamespace(psi_child=.17892072066041628, other_primitive=Array([3.]))
            current = SimpleNamespace(psi_child=.17892072066041628, other_primitive=Array([3.]))
            for i, field in enumerate(sorted(d.ENDOGENOUS_FIELDS)):
                setattr(saved, field, Array([float(i+1)]))
                setattr(current, field, Array([0.]))
            checkpoint = root / 'selected_native_packet.pkl.gz'
            with gzip.open(checkpoint,'wb') as stream:
                pickle.dump(dict(b_grid=Array([1.]),parameters=saved),stream)
            actual = d.digest(checkpoint)
            receipt = root / 'reference_reconstruction.json'
            receipt.write_text(json.dumps(dict(schema='current_floor_reference_reconstruction_v1',
                                               status='passed',identity={'reference':1},
                                               checkpoint={'path':str(checkpoint),'sha256':actual},
                                               checkpoint_sha256=actual)))
            pin = {'path':str(receipt),'sha256':d.digest(receipt)}
            runtime = SimpleNamespace(identity=lambda:{'reference':1},grid=Array([1.]),P=current,
                                      total_native_calls=0)
            with mock.patch.dict(sys.modules,{'numpy':fake_numpy}), mock.patch.object(d,'REFERENCE_SHA',actual):
                current.other_primitive = Array([4.])
                with self.assertRaisesRegex(ValueError,'Reference public primitive differs'):
                    d.bind_reference_cache_outputs(runtime,pin,root,root/'failed')
                self.assertTrue(all(getattr(current,field).values == [0.] for field in d.ENDOGENOUS_FIELDS))
                current.other_primitive = Array([3.])
                result = d.bind_reference_cache_outputs(runtime,pin,root,root/'passed')
                self.assertEqual(set(result['fields']),d.ENDOGENOUS_FIELDS)
                self.assertEqual(runtime.total_native_calls,0)
                for field in d.ENDOGENOUS_FIELDS:
                    self.assertEqual(getattr(current,field),getattr(saved,field))
                    self.assertIsNot(getattr(current,field),getattr(saved,field))


class BridgeRedirectTests(unittest.TestCase):
    def test_two_exact_receipts_redirect_and_writer_restored_on_error(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            reports = []
            for i in range(2):
                report = root / f'repeat_{i}/phase_b_ge/selected_root'
                report.mkdir(parents=True)
                original = report / 'saved_platform_bridge.json'
                original.write_text('immutable old content')
                reports.append(report)
            receipt = root / 'reference_reconstruction.json'
            receipt.write_text(json.dumps(dict(schema='current_floor_reference_reconstruction_v1',
                                               status='passed',identity={'id':1},
                                               reports={str(path):{} for path in reports})))
            pin = {'path':str(receipt),'sha256':d.digest(receipt)}
            source = root / 'frozen/source/code/model/experiments/birth_count_choice/transition_runtime.py'
            def original_writer(path,payload):
                path = Path(path); path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(json.dumps(payload,sort_keys=True))
            namespace = {'__file__':str(source),'write':original_writer}
            exec('def compare(self, reference, current):\n'
                 '    payload = dict(schema="current_estate_a_saved_platform_bridge_v1",status="passed",identity={"id":1},value=len(calls))\n'
                 '    calls.append(str(current))\n'
                 '    write(current/"saved_platform_bridge.json", payload)\n'
                 '    if fail[0] and len(calls)==2: raise RuntimeError("original comparison failed")\n',
                 namespace)
            calls = []; fail = [False]
            namespace['calls'] = calls; namespace['fail'] = fail
            runner = SimpleNamespace(compare_repeated=types.MethodType(namespace['compare'],object()))
            class Runtime:
                total_native_calls = 0
                def identity(self): return {'id':1}
                def load_reconstructed_reference(self,pin,folder):
                    for path in reports: runner.compare_repeated(reports[0],path)
                    return {'restored':True}
            rt = Runtime(); rt.runner = runner
            folder = root / 'new_output'
            result = d.restore_reference_with_bridge_redirect(rt,pin,root,folder)
            self.assertEqual(result,{'restored':True})
            self.assertIs(namespace['write'],original_writer)
            self.assertEqual(len(calls),2)
            for i,report in enumerate(reports):
                self.assertEqual((report/'saved_platform_bridge.json').read_text(),'immutable old content')
                self.assertEqual(json.loads((folder/f'repeat_{i}_bridge.json').read_text())['value'],i)
            bridge = json.loads((folder/'bridge_redirect.json').read_text())
            self.assertEqual(len(bridge['writes']),2)
            calls.clear(); fail[0] = True
            with self.assertRaisesRegex(RuntimeError,'original comparison failed'):
                d.restore_reference_with_bridge_redirect(rt,pin,root,root/'failed_output')
            self.assertIs(namespace['write'],original_writer)
            self.assertEqual(len(calls),2)
            self.assertTrue(all((report/'saved_platform_bridge.json').read_text()=='immutable old content'
                                for report in reports))


if __name__ == '__main__':
    unittest.main()
