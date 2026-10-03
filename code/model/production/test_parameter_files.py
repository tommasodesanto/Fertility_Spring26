"""Editable-file safety and native calibration export checks; no model solves."""
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from . import parameter_files as files
from .calibration import ANCHOR, check_adapter_without_solves, export_verified_parameters
from .inputs import load_inputs


def _same(left, right):
    assert set(vars(left)) == set(vars(right))
    for key, value in vars(left).items():
        other = getattr(right, key)
        if isinstance(value, np.ndarray):
            assert value.dtype == other.dtype, key
            assert np.array_equal(value, other, equal_nan=True), key
        else:
            assert value == other, key


def test_best_has_exact_default_inputs_and_toy_edits_propagate():
    config = files.load_parameter_file('best_params.py')
    P, grid = load_inputs(config['parameters'], config['external_inputs'], config['native_overrides'])
    base, base_grid = load_inputs()
    _same(P, base)
    assert len(vars(P)) == 245
    assert np.array_equal(grid, base_grid)
    toy = files.load_parameter_file('toy_params.py')
    Q, qgrid = load_inputs(toy['parameters'], toy['external_inputs'], toy['native_overrides'])
    assert Q.beta == (config['parameters']['beta_annual'] - .001) ** P.period_years
    assert Q.beta != P.beta
    toy['external_inputs']['phi'] = [.75] * 4
    Q, _ = load_inputs(toy['parameters'], toy['external_inputs'], toy['native_overrides'])
    assert np.array_equal(Q.phi, [.75] * 4)
    assert np.array_equal(qgrid, grid)
    assert files.load_parameter_file('toy_params.py')['external_inputs']['phi'] == [.8] * 4


def test_fresh_load_relative_path_and_unknown_controls(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    best = files.load_parameter_file('best_params.py')
    assert best['config_source'] == str(files.BEST_FILE.resolve())
    path = tmp_path / 'personal.py'
    text = files.BEST_FILE.read_text()
    path.write_text(text)
    first = files.load_parameter_file(path)
    path.write_text(text.replace('0.9663191380998087', '0.965'))
    second = files.load_parameter_file(path)
    assert second['parameters']['beta_annual'] == .965
    assert first['config_sha256'] != second['config_sha256']
    for addition in ('OUTPUT_ROOT = "/tmp/wrong"', 'EXTERNAL_INPUTZ = {}', 'import os', 'PROVENANCE = dict()'):
        path.write_text(text + '\n' + addition + '\n')
        with pytest.raises(ValueError):
            files.load_parameter_file(path)
    path.write_text(text.replace('"sigma": 2.0', '"unknown_field": 2.0'))
    with pytest.raises(ValueError, match='unknown native field'):
        files.load_parameter_file(path)


def test_output_isolation_and_colliding_sources(tmp_path, monkeypatch):
    monkeypatch.setattr(files, 'PROJECT_ROOT', tmp_path / 'project')
    assert files.output_root_for(files.BEST_FILE) == tmp_path / 'project/output/model/local_solution'
    toy = tmp_path / 'toy.py'
    toy.write_text(files.BEST_FILE.read_text())
    root = files.output_root_for(toy)
    assert not root.exists()
    assert files.claim_output_root(toy) == root
    assert root == tmp_path / 'project/output/model/experiments/toy'
    assert files.output_root_for(toy) == root
    other = tmp_path / 'elsewhere/toy.py'
    other.parent.mkdir()
    other.write_text(toy.read_text())
    with pytest.raises(ValueError, match='already owned'):
        files.output_root_for(other)
    existing = tmp_path / 'existing.py'
    existing.write_text(toy.read_text())
    (tmp_path / 'project/output/model/experiments/existing').mkdir()
    with pytest.raises(ValueError, match='already owned'):
        files.output_root_for(existing)


def _local_receipt(tmp_path):
    """Relocated real historical receipt; deliberately lacks new input identity."""
    receipt = json.loads(ANCHOR.read_text())
    receipt['selected_postcheck']['report'] = str(ANCHOR.parent / 'native_postcheck/selected_postcheck/phase_b_ge/selected_root')
    source = tmp_path / 'completed.json'
    source.write_text(json.dumps(receipt))
    return source


def _synthetic_receipt(tmp_path, P, grid):
    """Explicit mocked solve using real table fixtures; no native verification."""
    import time
    from .calibration import make_evaluator
    receipt = json.loads(_local_receipt(tmp_path).read_text())
    report = Path(receipt['selected_postcheck']['report'])
    def mocked_solver(Q, supplied_grid, **kwargs):
        assert Q.sigma == P.sigma
        assert np.array_equal(Q.phi, P.phi)
        assert np.array_equal(supplied_grid, grid)
        # Identity must be recorded before native output mutation.
        Q.native_inherited_distribution_evidence_dir = '/new/output/directory'
        Q.H0[:] *= 1.01  # Mock a native derived-output update after identity capture.
        return dict(report_directory=str(report), closure=json.loads((report / 'closure.json').read_text()),
                    price=receipt['selected']['price'], lifecycle_solves=2)
    deadline = time.time() + 30
    evaluate = make_evaluator(tmp_path / 'mock', 'floor_s0', P, grid, deadline, solver=mocked_solver)
    receipt['selected_postcheck'] = evaluate('mock_selected', receipt['selected']['parameters'], deadline)
    receipt['test_fixture'] = 'synthetic mocked solver; retained tables/repeat receipt are test fixtures'
    source = tmp_path / 'completed.json'
    source.write_text(json.dumps(receipt))
    return source


def test_synthetic_verified_export_preserves_all_effective_fields(tmp_path):
    original = hashlib.sha256(files.BEST_FILE.read_bytes()).hexdigest()
    P, grid = load_inputs()
    # Caller primitives must survive an export rather than reverting to defaults.
    P.sigma = 2.01
    P.phi[:] = .75
    P.hR_max = 5.9  # Supported advanced primitive, preserved in NATIVE_OVERRIDES.
    receipt = json.loads(ANCHOR.read_text())
    destination = tmp_path / 'run/best_params.py'
    export = export_verified_parameters(destination, _synthetic_receipt(tmp_path, P, grid), P, grid)
    assert export['effective_field_count'] == 245
    assert export['global_default_promoted'] is False
    config = files.load_parameter_file(destination)
    restored, restored_grid = load_inputs(config['parameters'], config['external_inputs'], config['native_overrides'])
    expected = copy.deepcopy(P)
    expected.H0[:] = receipt['selected_postcheck']['H0_derived']
    _same(restored, expected)
    assert np.array_equal(restored_grid, grid)
    assert config['price_guess'] == receipt['selected_postcheck']['price']
    assert config['provenance']['target_fingerprint'] == receipt['target_fingerprint']
    source = Path(config['provenance']['completed_receipt'])
    (destination.parent / 'parameter_file_export.json').write_text(json.dumps(export))
    assert config['provenance']['completed_receipt_sha256'] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert hashlib.sha256(files.BEST_FILE.read_bytes()).hexdigest() == original
    with pytest.raises(ValueError, match='overwrite'):
        export_verified_parameters(destination, _synthetic_receipt(tmp_path, P, grid), P, grid)


@pytest.mark.parametrize('change', ['provisional', 'target', 'historical', 'repeat', 'coordinates'])
def test_export_rejects_unverified_or_other_contracts(tmp_path, change):
    receipt = json.loads(_local_receipt(tmp_path).read_text())
    if change == 'provisional':
        receipt['status'] = 'provisional_search_finished'
    elif change == 'target':
        receipt['target_fingerprint'] = 'other-wealth-target'
    elif change == 'historical':
        receipt['arm'] = 'original'
    elif change == 'repeat':
        receipt['repeat']['status'] = 'unchecked'
    else:
        receipt['selected']['parameters']['beta_annual'] = .965
    source = tmp_path / 'completed.json'
    source.write_text(json.dumps(receipt))
    P, grid = load_inputs()
    with pytest.raises((ValueError, RuntimeError)):
        export_verified_parameters(tmp_path / 'run/best_params.py', source, P, grid)
    assert not (tmp_path / 'run/best_params.py').exists()


def test_export_rejects_changed_structural_grid_and_preserves_adapter_pin(tmp_path):
    P, grid = load_inputs()
    changed_grid = grid.copy()
    changed_grid[10] += .00001
    with pytest.raises(ValueError, match='changed wealth grid'):
        export_verified_parameters(tmp_path / 'run/best_params.py', _synthetic_receipt(tmp_path, P, changed_grid), P, changed_grid)
    assert check_adapter_without_solves()['full_target_pin_rejected']


def test_external_best_exports_have_distinct_output_roots(tmp_path, monkeypatch):
    monkeypatch.setattr(files, 'PROJECT_ROOT', tmp_path / 'project')
    one = tmp_path / 'first_run/best_params.py'
    two = tmp_path / 'second_run/best_params.py'
    for path in (one, two):
        path.parent.mkdir()
        path.write_text(files.BEST_FILE.read_text())
    assert files.output_root_for(one) != files.output_root_for(two)
    assert files.output_root_for(one).name.startswith('best_params__first_run__')
    assert not files.output_root_for(one).exists()


def test_canonical_credit_preparation_and_mock_evaluator_preserve_caller(tmp_path):
    import importlib.util
    import time
    from .calibration import make_evaluator
    driver_path = Path(__file__).resolve().parents[1] / 'experiments/purchase_timing_sandbox/calibrate.py'
    spec = importlib.util.spec_from_file_location('canonical_credit_driver_test', driver_path)
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    P, grid = load_inputs()
    del P.unsecured_credit_limit  # Exact missing field from historical preparation.
    P.sigma = 2.01
    before = copy.deepcopy(vars(P))
    prepared = driver.bind_canonical_credit(P)
    assert prepared is P and prepared.unsecured_credit_limit == 0.0
    assert set(vars(prepared)) - set(before) == {'unsecured_credit_limit'}
    for key, value in before.items():
        if isinstance(value, np.ndarray):
            assert np.array_equal(getattr(prepared, key), value, equal_nan=True), key
        else:
            assert getattr(prepared, key) == value, key
    anchor = json.loads(ANCHOR.read_text())
    report = ANCHOR.parent / 'native_postcheck/selected_postcheck/phase_b_ge/selected_root'
    def mocked_solver(Q, supplied_grid, **kwargs):
        assert Q.unsecured_credit_limit == 0.0
        assert Q.sigma == 2.01
        assert np.array_equal(supplied_grid, grid)
        return dict(report_directory=str(report), closure=json.loads((report / 'closure.json').read_text()),
                    price=anchor['selected']['price'], lifecycle_solves=2)
    deadline = time.time() + 30
    evaluate = make_evaluator(tmp_path / 'mock', 'floor_s0', prepared, grid, deadline, solver=mocked_solver)
    assert evaluate('selected', anchor['selected']['parameters'], deadline)['status'] == 'passed'
    export = export_verified_parameters(tmp_path / 'run/best_params.py', _synthetic_receipt(tmp_path, prepared, grid), prepared, grid)
    assert export['full_input_roundtrip']
    P.unsecured_credit_limit = .1
    with pytest.raises(RuntimeError, match='zero unsecured credit'):
        driver.bind_canonical_credit(P)


def test_real_historical_receipt_cannot_authenticate_changed_fixed_inputs(tmp_path):
    P, grid = load_inputs()
    P.sigma = 2.01
    P.phi[:] = .75
    P.hR_max = 5.9
    with pytest.raises(ValueError, match='lacks authenticated caller-input identity'):
        export_verified_parameters(tmp_path / 'run/best_params.py', _local_receipt(tmp_path), P, grid)
    assert not (tmp_path / 'run/best_params.py').exists()


@pytest.mark.parametrize('change', ['sigma', 'phi', 'hR_max', 'H0', 'grid'])
def test_export_rejects_any_caller_change_after_native_check(tmp_path, change):
    P, grid = load_inputs()
    source = _synthetic_receipt(tmp_path, P, grid)
    if change == 'grid':
        grid[10] += .00001
    elif change in ('phi', 'H0'):
        getattr(P, change)[:] *= .99
    else:
        setattr(P, change, getattr(P, change) + .01)
    with pytest.raises(ValueError, match='differ from the fresh native verified candidate'):
        export_verified_parameters(tmp_path / 'run/best_params.py', source, P, grid)
    assert not (tmp_path / 'run/best_params.py').exists()


def test_input_fingerprint_is_deterministic_for_nonfinite_optional_fields():
    from .calibration import effective_input_fingerprint
    P, grid = load_inputs()
    P.optional_identity_probe = (float('inf'), float('-inf'), float('nan'))
    first = effective_input_fingerprint(P, grid)
    Q = copy.deepcopy(P)
    Q.native_inherited_distribution_evidence_dir = '/another/run/output'
    assert effective_input_fingerprint(Q, grid.copy()) == first
    Q.optional_identity_probe = (float('inf'), float('inf'), float('nan'))
    assert effective_input_fingerprint(Q, grid) != first
