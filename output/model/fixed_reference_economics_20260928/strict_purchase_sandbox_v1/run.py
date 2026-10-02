"""Single fixed-coordinate evaluation of the strict wealth-only purchase experiment."""
from __future__ import annotations
import csv, hashlib, json, os, signal, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SANDBOX = ROOT / 'code/model/experiments/strict_purchase_sandbox/source'
sys.path.insert(0, str(SANDBOX))
from small_credit_lab.engine import solver as sandbox_solver  # noqa: E402
from small_credit_lab.engine import household as sandbox_household, kernels as sandbox_kernels  # noqa: E402
from refactor_lab.engine import household as checked_household, kernels as checked_kernels  # noqa: E402
V2 = HERE.parent / 'normalized_calibration_v2'
sys.path.insert(0, str(V2))
import run_psi as v2  # noqa: E402

LIMIT_SECONDS = 1200
LANE = 'floor_s0'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def verify():
    contract = json.loads((HERE / 'manifest.json').read_text())
    for rel, hashes in contract['source_pairs'].items():
        original, sandbox = rel.split('|')
        assert sha(ROOT / original) == hashes['original'], original
        assert sha(ROOT / sandbox) == hashes['sandbox'], sandbox
    for module in (sandbox_solver, sandbox_household, sandbox_kernels, checked_household, checked_kernels):
        assert Path(module.__file__).resolve().is_relative_to(SANDBOX), module.__file__
    v2.native.verify_sources()
    for rel, digest in json.loads((V2 / 'source_pins.json').read_text()).items():
        assert sha(ROOT / rel) == digest, rel
    incumbent = json.loads((HERE / 'incumbent.json').read_text())
    assert sha(HERE / 'incumbent.json') == contract['incumbent_sha256']
    assert incumbent['chain'] == 2 and incumbent['case'] == '0064_nm'
    assert incumbent['postcheck_status'] == 'selected_numerically_verified'
    assert v2.inputs.canonical(v2.CONFIG['base_target_contract']) == contract['target_fingerprint']
    assert v2.weight_fingerprint({}) == contract['weight_fingerprint']
    assert v2.CONFIG['profiles']['base_control'] == {}
    return incumbent, contract

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    assert not args.out.exists(), 'Refusing existing output directory'
    args.out.mkdir(parents=True)
    t0 = time.time()
    deadline = t0 + LIMIT_SECONDS
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('20-minute cap')))
    signal.setitimer(signal.ITIMER_REAL, LIMIT_SECONDS)
    try:
        incumbent, contract = verify()
        point = incumbent['parameters']
        assert len(point) == 10
        seed, bounds, _ = v2.inputs.seed_and_bounds(LANE)
        bounds = {k: tuple(x) for k, x in bounds.items()}
        bounds['psi_child'] = tuple(v2.CONFIG['psi_bounds'])
        assert set(point) == set(v2.inputs.parameters(LANE)) | {'psi_child'}
        assert all(bounds[k][0] <= point[k] <= bounds[k][1] for k in point)
        v2.inputs.LANES[LANE].update(seed=dict(point), bounds=bounds, free_coordinates=list(point))
        P, grid = v2.inputs.proposal(LANE)
        P, _ = v2.inputs.entry(P, grid, 'nonnegative_mean')
        assert P.R_gross > 1.0 and P.native_purchase_income is True
        assert P.native_due_stayer_credit is True and P.joint_nested_choice is False
        assert P.N_target == 1.0
        v2.write(args.out / 'start.json', dict(status='running', started_epoch=t0,
            deadline_epoch=deadline, selected_price_start=incumbent['selected_price'],
            fixed_parameters=point, experimental_change='Experimental wealth-only purchase origination restriction; current income excluded from eligibility. Budget, interest timing, debt bound and stayer rules unchanged',
            source_manifest_sha256=sha(HERE / 'manifest.json'),
            target_fingerprint=contract['target_fingerprint'], weight_fingerprint=contract['weight_fingerprint']))
        Q = v2.native.utility_checks(P, grid, LANE, args.out)
        evaluate = v2.normalized_objective.make_evaluator(args.out, LANE, Q, grid, deadline,
            incumbent['selected_price'], native_runner=v2.native)
        review = args.out / 'normalization_source_review/receipt.json'
        receipt = json.loads(review.read_text())
        receipt['household_solver_unchanged'] = False
        receipt['isolated_strict_purchase_source_manifest_sha256'] = sha(HERE / 'manifest.json')
        v2.write(review, receipt)
        result = evaluate('strict_purchase', point, deadline)
        assert result['status'] == 'passed', result
        report = Path(result['report'])
        with (report / 'target_fit.csv').open(newline='') as f: targets = list(csv.DictReader(f))
        with (report / 'parameters.csv').open(newline='') as f: parameters = list(csv.DictReader(f))
        assert len(targets) == 14 and len(parameters) == 31
        assert v2.native.target_identity(targets) == v2.CONFIG['base_target_contract']
        assert len(list((report / 'standard_diagnostics').glob('*.png'))) == 17
        repeat = v2.native.compare_repeated(report, report.parent / 'selected_repeat_final')
        result.update(status='fixed_coordinate_strict_purchase_experiment_passed', target_fit=targets,
            parameters=parameters, repeat=repeat,
            loss=float(sum(float(row['loss_contribution'] or 0) for row in targets)),
            source_manifest_sha256=sha(HERE / 'manifest.json'),
            target_fingerprint=contract['target_fingerprint'], weight_fingerprint=contract['weight_fingerprint'],
            all_ten_coordinates_fixed=True, experimental_not_adopted=True,
            elapsed_seconds=time.time() - t0)
        v2.write(args.out / 'completed.json', result)
    except BaseException as exc:
        v2.write(args.out / 'failure.json', dict(type=type(exc).__name__, message=str(exc),
            elapsed_seconds=time.time() - t0))
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)

if __name__ == '__main__': main()
