"""Preset routing checks with no equilibrium solves or reference initialization."""
from pathlib import Path
import importlib.util
import subprocess
import sys
from types import ModuleType, SimpleNamespace

import pytest

MODEL = Path(__file__).resolve().parents[1]


def module_at(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_help_avoids_numeric_and_production_imports():
    for name in ('run_model.py', 'plot_model_policies.py', 'plot_model_aggregates.py',
                 'tools/model_playground.py'):
        script = """import runpy, sys
sys.argv = [sys.argv[1], '--help']
try: runpy.run_path(sys.argv[0], run_name='__main__')
except SystemExit as exc: assert exc.code == 0
assert 'numpy' not in sys.modules
assert 'production.inputs' not in sys.modules
assert 'production.equilibrium' not in sys.modules
"""
        subprocess.run([sys.executable, '-c', script, str(MODEL / name)],
                       check=True, capture_output=True, text=True, timeout=5)


def test_runner_passes_selected_controls_and_identity(monkeypatch, tmp_path):
    runner = module_at('runner_test', MODEL / 'run_model.py')
    config = dict(parameters={'test': 1}, external_inputs={'test': 2},
                  native_overrides={}, price_guess=0.7, budget_seconds=30,
                  closure='fixed_h0', provenance={'role': 'toy'},
                  config_source='/toy.py', config_sha256='abc', config_text='source')
    files = ModuleType('production.parameter_files')
    calls = []
    files.load_parameter_file = lambda path: calls.append(('load', path)) or config
    files.claim_output_root = lambda path: calls.append(('claim', path)) or tmp_path
    flow = ModuleType('production.workflow')
    def stationary(*args, **kwargs):
        calls.append(('run', args, kwargs))
        return SimpleNamespace(price=0.7, closure={}), tmp_path / 'case'
    flow.run_stationary = stationary
    monkeypatch.setitem(sys.modules, 'production.parameter_files', files)
    monkeypatch.setitem(sys.modules, 'production.workflow', flow)
    assert runner.main(['--params', 'toy_params.py']) == tmp_path / 'case'
    assert calls[:2] == [('load', 'toy_params.py'), ('claim', 'toy_params.py')]
    assert calls[2][2]['parameter_file_metadata']['config_sha256'] == 'abc'
    assert calls[2][2]['output_root'] == tmp_path


def test_policy_plotter_selected_cache_and_explicit_directory(monkeypatch, tmp_path):
    plotter = module_at('plotter_test', MODEL / 'plot_model_policies.py')
    runner = ModuleType('run_model')
    runner.PARAMETER_FILE = 'toy_params.py'
    files = ModuleType('production.parameter_files')
    files.output_root_for = lambda path: tmp_path / Path(path).stem
    files.describe_saved_case = lambda directory, preset: str(directory)
    storage = ModuleType('production.storage')
    calls = []
    storage.load_latest = lambda path: calls.append(('latest', path)) or (object(), path)
    storage.load_case = lambda path: calls.append(('case', path)) or (object(), path)
    for name, mod in [('run_model', runner), ('production.parameter_files', files),
                      ('production.storage', storage)]:
        monkeypatch.setitem(sys.modules, name, mod)
    plotter._plot_run = lambda result, path: [path]
    plotter.main()
    plotter.main(parameter_file='best_params.py')
    plotter.RUN_DIRECTORY = tmp_path / 'pinned'
    plotter.main(parameter_file='toy_params.py')
    assert calls == [('latest', tmp_path / 'toy_params'),
                     ('latest', tmp_path / 'best_params'), ('case', tmp_path / 'pinned')]


def test_policy_plotter_missing_toy_never_falls_back(monkeypatch, tmp_path):
    plotter = module_at('missing_plotter_test', MODEL / 'plot_model_policies.py')
    files = ModuleType('production.parameter_files')
    files.output_root_for = lambda path: tmp_path / Path(path).stem
    files.describe_saved_case = lambda directory, preset: str(directory)
    storage = ModuleType('production.storage')
    calls = []
    def missing(path):
        calls.append(path)
        raise FileNotFoundError(path)
    storage.load_latest = missing
    monkeypatch.setitem(sys.modules, 'production.parameter_files', files)
    monkeypatch.setitem(sys.modules, 'production.storage', storage)
    with pytest.raises(FileNotFoundError):
        plotter.main(parameter_file='toy_params.py')
    assert calls == [tmp_path / 'toy_params']


def test_playground_reset_restores_selected_preset(monkeypatch):
    import numpy as np
    playground = module_at('playground_frontend_test', MODEL / 'tools/model_playground.py')
    import production
    import production.parameter_files as files
    inputs = ModuleType('production.inputs')
    calls = []
    def load_inputs(parameters, external, native):
        calls.append((dict(parameters), dict(external), dict(native)))
        return SimpleNamespace(chi=parameters['chi']), np.array([0., 1.])
    inputs.load_inputs = load_inputs
    engine = ModuleType('production.equilibrium')
    monkeypatch.setitem(sys.modules, 'production.inputs', inputs)
    monkeypatch.setitem(sys.modules, 'production.equilibrium', engine)
    monkeypatch.setattr(production, 'equilibrium', engine, raising=False)
    config = dict(parameters={'chi': 1.2}, external_inputs={'sigma': 2.3},
                  native_overrides={'c_min': 0.05}, price_guess=0.6)
    monkeypatch.setattr(files, 'load_parameter_file', lambda path: config)
    model = playground.CanonicalModelPlayground(parameter_file='toy_params.py')
    model.params['chi'] = 5
    model.external_inputs['sigma'] = 5
    model.native_overrides.clear()
    model.reference_price = 5
    model.reset()
    assert model.params is model.parameters
    assert model.params == {'chi': 1.2}
    assert model.external_inputs == {'sigma': 2.3}
    assert model.native_overrides == {'c_min': 0.05}
    assert model.reference_price == 0.6
    assert calls[0] == calls[1]


def test_saved_case_note_discloses_edited_file(monkeypatch, tmp_path):
    import hashlib
    import json
    import production.parameter_files as files
    runner = module_at('note_frontend_test', MODEL / 'run_model.py')
    source = tmp_path / 'toy_params.py'
    source.write_text('original')
    case = tmp_path / 'case'
    case.mkdir()
    metadata = {'config_sha256': hashlib.sha256(source.read_bytes()).hexdigest()}
    (case / 'input_contract.json').write_text(json.dumps({'parameter_file': metadata}))
    monkeypatch.setattr(files, 'resolve_parameter_file', lambda path: source)
    assert 'changed' not in files.describe_saved_case(case, 'toy_params.py')
    source.write_text('edited')
    assert 'last successful saved run' in files.describe_saved_case(case, 'toy_params.py')
    (case / 'input_contract.json').write_text('{}')
    assert 'changed' not in files.describe_saved_case(case, 'toy_params.py')


def test_unsolved_toy_playground_can_initialize_without_solve(monkeypatch, tmp_path, capsys):
    import production.parameter_files as files
    playground = module_at('unsolved_playground_test', MODEL / 'tools/model_playground.py')
    import production.storage as storage
    monkeypatch.setattr(files, 'output_root_for', lambda path: tmp_path / 'toy')
    def missing(path):
        raise FileNotFoundError(path)
    monkeypatch.setattr(storage, 'load_latest', missing)
    calls = []
    class FakePlayground:
        def __init__(self, parameter_file=None):
            calls.append(parameter_file)
            self.P = object()
        def show_parameters(self):
            pass
    monkeypatch.setattr(playground, 'CanonicalModelPlayground', FakePlayground)
    playground.main(['--params', 'toy_params.py'])
    assert playground.sol is None
    assert calls == ['toy_params.py']
    assert 'No saved solution yet' in capsys.readouterr().out


@pytest.mark.parametrize('field,value', [('sigma', 2.01), ('phi', [0.75] * 4)])
def test_direct_playground_edit_replaces_preset_default(monkeypatch, field, value):
    import copy
    import numpy as np
    import production
    import production.parameter_files as files
    playground = module_at('direct_edit_playground_test', MODEL / 'tools/model_playground.py')
    inputs = ModuleType('production.inputs')
    def load_inputs(parameters, external, native):
        for name in set(external) & set(native):
            if not playground._same_value(external[name], native[name]):
                raise ValueError('conflicting controls')
        return SimpleNamespace(**copy.deepcopy({**external, **native})), np.array([0., 1.])
    inputs.load_inputs = load_inputs
    engine = ModuleType('production.equilibrium')
    solves = []
    def solve_at_price(P, grid, price):
        solves.append(P)
        return {'P': P, 'solution': SimpleNamespace()}
    engine.solve_at_price = solve_at_price
    monkeypatch.setitem(sys.modules, 'production.inputs', inputs)
    monkeypatch.setitem(sys.modules, 'production.equilibrium', engine)
    monkeypatch.setattr(production, 'equilibrium', engine, raising=False)
    config = dict(parameters={}, external_inputs={'sigma': 2., 'phi': [0.8] * 4},
                  native_overrides={}, price_guess=0.6)
    monkeypatch.setattr(files, 'load_parameter_file', lambda path: config)
    model = playground.CanonicalModelPlayground(parameter_file='toy_params.py')
    setattr(model.P, field, value)
    result = model.solve()
    assert playground._same_value(getattr(solves[0], field), value)
    assert playground._same_value(getattr(result.P, field), value)
    # The unchanged editable defaults remain available for reset.
    assert model.external_inputs == config['external_inputs']
    model.external_inputs[field] = 2.2 if field == 'sigma' else [0.7] * 4
    with pytest.raises(ValueError, match='Conflicting input'):
        model.solve()
    assert len(solves) == 1
    model.external_inputs[field] = copy.deepcopy(config['external_inputs'][field])
    with pytest.raises(ValueError, match='Conflicting input'):
        model.solve(external_inputs={field: config['external_inputs'][field]})
    with pytest.raises(ValueError, match='Conflicting input'):
        model.solve(native_overrides={field: config['external_inputs'][field]})
    assert len(solves) == 1
