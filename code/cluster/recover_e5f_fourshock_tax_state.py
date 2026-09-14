"""Recover the announced-four-shock 2023 state, then dispatch its 2% tax arm.

``--prepare`` is login-node-only: it freezes a recovery job behind the original
104-date job.  The compute-node recovery selects a completed native mapping,
uses exactly 5+J backward dates, and submits the already-reviewed tax driver
only after every 2007--2023 row has replayed at native tolerance.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import shlex
import shutil
import subprocess
import sys
import threading
import time
from types import SimpleNamespace as NS
from unittest.mock import patch

for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'

BASE = Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a')
BATCH = BASE / 'batches/announced_original_queue_20260913c'
COMMON = BASE / 'batches/afternoon_original_queue_20260913a/spec.json'
SHARED_TAX2 = BASE / 'batches/inherited_2023_tax_20260913a/tax2/endpoint'
TAX_DRIVER = BASE / 'batches/inherited_2023_tax_20260913b/source/run_e5f_inherited_2023_tax_long.py'
PYTHON = '/share/apps/anaconda3/2025.06/bin/python'
FINAL_PSI = 0.09221854783921073
TOL = 2e-10


def read(path): return json.loads(Path(path).read_text())
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def json_default(value):
    if hasattr(value, 'tolist'):
        return value.tolist()
    raise TypeError(f'Unsupported JSON value: {type(value).__name__}')

def save(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, default=json_default) + '\n'); tmp.replace(path)


def copy_if_different(source, destination):
    """Permit prepare to be invoked either from the repository or batch/source."""
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if source != destination: shutil.copy2(source, destination)


def _coordinates(rows):
    """The only accepted coordinate order is price, pension, equal rebate."""
    if not isinstance(rows, list) or len(rows) != 104:
        raise ValueError('A full announced mapping must contain exactly 104 rows')
    out = []
    for key in ('asset_price', 'pension_period_units', 'equal_transfer_period_units'):
        vals = [float(r[key]) for r in rows]
        if not all(__import__('math').isfinite(x) for x in vals):
            raise ValueError('Nonfinite mapping coordinate: ' + key)
        out.extend(vals)
    return out


def select_mapping(source):
    """Pin the best saved root coordinate if it has a complete native packet.

    Root receipts are not trusted merely by their mapping number.  A candidate
    has to reproduce its advertised 312 coordinates in its own rows.json.
    """
    source = Path(source); out = source / 'output/run'
    best_path = out / 'best_so_far.json'
    advertised = read(best_path).get('prices') if best_path.is_file() else None
    choices = []
    for rows_path in sorted((out / 'mappings').glob('mapping_*/rows.json')):
        try:
            rows = read(rows_path); coords = _coordinates(rows)
        except (OSError, ValueError, KeyError, TypeError):
            continue
        choices.append((rows_path, rows, coords))
    if not choices:
        raise RuntimeError('No complete 104-date native mapping is available after announced job')
    if isinstance(advertised, list) and len(advertised) == 312:
        for rows_path, rows, coords in choices:
            if max(abs(float(a) - float(b)) for a, b in zip(coords, advertised)) <= TOL:
                return rows_path, rows, coords, 'best_so_far_exact_coordinate_match'
    # A complete latest mapping is only a fallback when root best coordinates
    # were not saved alongside an evaluable packet.
    latest = read(out / 'latest_completed_mapping.json') if (out / 'latest_completed_mapping.json').is_file() else {}
    wanted = str(latest.get('path', ''))
    for candidate in choices:
        if str(candidate[0].parent) == wanted:
            return *candidate, 'latest_completed_complete_fallback'
    return *choices[-1], 'lexicographically_latest_complete_fallback'


def prepare(batch):
    batch = Path(batch).resolve()
    if (batch / 'manifest.json').exists(): raise ValueError('Refusing duplicate four-shock tax recovery')
    announced = read(BATCH / 'manifest.json')
    if list(map(float, announced['psi_levels'])) != [0.12891531457859182, 0.11696608375682901, 0.10564290922456478, FINAL_PSI]:
        raise ValueError('Frozen announced four-shock schedule differs')
    one_endpoint = Path(announced['endpoint_pickle']); one_receipt = Path(announced['endpoint_receipt'])
    deps = [BATCH / 'manifest.json', BATCH / 'source/run_e5f_announced_original_queue.py', COMMON, one_endpoint, one_receipt, SHARED_TAX2 / 'terminal.pkl.gz', SHARED_TAX2 / 'root_receipt.json', TAX_DRIVER]
    if not all(p.is_file() for p in deps): raise ValueError('Missing frozen source/input for dispatch')
    if not read(SHARED_TAX2 / 'root_receipt.json').get('verified', False): raise ValueError('Shared 2% endpoint is unverified')
    source = batch / 'source'; source.mkdir(parents=True, exist_ok=True)
    this = source / Path(__file__).name; driver = source / 'run_e5f_inherited_2023_tax_long.py'
    copy_if_different(Path(__file__), this); copy_if_different(TAX_DRIVER, driver)
    m = dict(output=str(batch / 'output'), announced_batch=str(BATCH), announced_manifest=str(BATCH / 'manifest.json'),
        spec=str(COMMON), dependency_job=17711519, recovery_seconds=1800, final_psi=FINAL_PSI,
        psi_levels=announced['psi_levels'], shared_tax_endpoint=str(SHARED_TAX2 / 'terminal.pkl.gz'),
        shared_tax_endpoint_receipt=str(SHARED_TAX2 / 'root_receipt.json'), policy_output=str(batch / 'policy'),
        source_driver=str(driver), source_recovery=str(this),
        file_sha256={str(p): sha(p) for p in deps} | {str(this): sha(this), str(driver): sha(driver)},
        contract=dict(source='announced four-shock job 17711519; no silent substitute', observed_year=2023,
            forward_dates=5, backward_dates=22, root_solves=0, tax=.02, seed_horizon=104,
            information='Four preference levels announced in 2007; 2023 property-tax reform unexpected thereafter'))
    batch.mkdir(parents=True, exist_ok=True); save(batch / 'manifest.json', m)
    folder = batch / 'recovery'; folder.mkdir()
    argv = [PYTHON, str(this), '--run', str(batch / 'manifest.json')]
    script = folder / 'run.sbatch'
    script.write_text('\n'.join(['#!/bin/bash', '#SBATCH --job-name=e5f_fourshock_2023_recover', '#SBATCH --account=torch_pr_570_general', '#SBATCH --cpus-per-task=1', '#SBATCH --mem=32G', '#SBATCH --time=30', f'#SBATCH --dependency=afterany:{m["dependency_job"]}', f'#SBATCH --output={folder}/slurm_%j.log', 'set -euo pipefail', 'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1', shlex.join(argv), '']))
    receipt = dict(manifest_sha256=sha(batch / 'manifest.json'), script_sha256=sha(script), dependency=f'afterany:{m["dependency_job"]}', submitted=False)
    # Preparation intentionally only creates reviewable Slurm material.  The
    # lead submits it after checking the inherited-tax driver interface.
    save(batch / 'dispatch.json', receipt); print(json.dumps(receipt))


def run(manifest):
    import numpy as np
    m = read(manifest)
    for path, digest in m['file_sha256'].items():
        if sha(path) != digest: raise ValueError('Changed pinned input: ' + path)
    recovery = Path(m['output']) / 'recovery'
    if recovery.exists(): raise ValueError('Refusing to overwrite recovery output')
    recovery.mkdir(parents=True); deadline = time.monotonic() + float(m['recovery_seconds']); stop = threading.Event()
    def heartbeat():
        while not stop.wait(60):
            save(recovery / 'heartbeat.json', dict(remaining_seconds=deadline-time.monotonic()))
            if time.monotonic() >= deadline: save(recovery / 'failure.json', dict(error='30-minute recovery deadline')); os._exit(124)
    threading.Thread(target=heartbeat, daemon=True).start()
    try:
        rows_path, fullrows, coords, selection = select_mapping(m['announced_batch'])
        # Record source files/coordinates before scientific imports or work.
        save(recovery / 'seed.json', dict(prices=coords, horizon=104, source_rows=str(rows_path), source_rows_sha256=sha(rows_path), selection=selection))
        save(recovery / 'initial_reference.json', dict(source_rows=str(rows_path), source_rows_sha256=sha(rows_path), coordinate_sha256=sha(recovery / 'seed.json'), selection=selection))
        save(recovery / 'baseline_rows.json', fullrows[4:])
        save(recovery / 'frozen_rows.json', fullrows[:5])
        sys.path.insert(0, str(Path(m['spec']).parent / 'source'))
        sys.path.insert(0, str(Path(m['announced_batch']) / 'source'))
        import run_e5f_original_queue_experiments as runner
        c = runner.load_context(m['spec']); pf = c.joined.pf
        # load_context installs the frozen model runtime before these imports.
        import run_e5f_announced_original_queue as announced
        import run_e5f_transition_calibration as fertility
        from e5f_balanced_terminal import _household_checks
        count, backwards = 5, 5 + int(c.old.parameters.J)
        if backwards != 22: raise ValueError('Expected the 17-age-cell 5+J recovery')
        prices = np.asarray(coords[:104], float); pensions = np.asarray(coords[104:208], float); transfers = np.asarray(coords[208:], float)
        # The original 1% terminal boundary is the one pinned by the announced manifest.
        am = read(m['announced_manifest'])
        with gzip.open(am['endpoint_pickle'], 'rb') as f: endpoint = pickle.load(f)
        psi = np.r_[np.asarray(m['psi_levels'], float), np.full(backwards-len(m['psi_levels']), FINAL_PSI)]
        if psi.shape != (backwards,): raise ValueError('Announced finite recovery psi path must have exactly 22 dates')
        observations, saved = [], {}
        def observe(i, e, P, grid, shared):
            rents = np.asarray(pf.rents_from_asset_prices(prices[:count], float(prices[count]), c.old.parameters), float)
            _, gates = _household_checks(e, P, shared, grid, float(rents[i]), c.primitive, c.audit)
            if not gates or not all(gates.values()): raise RuntimeError('Native household audit failed')
            observations.append(dict(period=i, calendar_year=2007+4*i, **fertility.period_fertility_diagnostics(e, P)))
            if i == 4:
                saved.update(parameters=P, b_grid=grid, evaluation=e, shared=shared, supply_rule=c.old.supply_rule, calendar_year=2023, source='announced_four_shock')
        with c.queue.original_queue_adapter(), c.cache.policy_cache(pf, max_bytes=12*1024**3):
            rents = np.asarray(pf.rents_from_asset_prices(prices[:backwards], float(prices[backwards]), c.old.parameters), float)
            values, _ = pf.backward_value_path(prices=prices[:backwards], rents=rents, psi_path=psi,
                terminal_V=endpoint.policy.V, base_parameters=c.old.parameters, b_grid=c.old.b_grid,
                transfer_path=transfers[:backwards], pension_path=pensions[:backwards], payroll_tax_path=np.full(backwards, .179))
            evaluate, _ = announced.announced_queue_path(c, psi[:count])
            result = evaluate(inherited=c.rebated.InheritedState(2007, c.old.initial_state), old_state=c.old,
                prices=prices[:count], pensions=pensions[:count], transfers=transfers[:count], psi=float(psi[count-1]),
                terminal=NS(parameters=NS(psi_child=FINAL_PSI), policy=NS(V=values[count]), asset_price=float(prices[count])), observer=observe)
        # Strict full-row equivalence catches a wrong announced-news continuation,
        # rather than only checking the three equilibrium coordinates.
        maxgap = 0.
        if len(result.rows) != count: raise RuntimeError('Five-date native replay missing')
        for i, (got, want) in enumerate(zip(result.rows, fullrows[:count])):
            if set(got) != set(want): raise ValueError(f'Row {i} key mismatch')
            for key, value in want.items():
                if isinstance(value, (bool, str)) or value is None:
                    if got[key] != value: raise ValueError(f'Row {i} differs: {key}')
                else:
                    gap = abs(float(got[key])-float(value)); maxgap=max(maxgap,gap)
                    if not np.isfinite(gap) or gap > TOL: raise ValueError(f'Row {i} differs: {key}')
        with gzip.open(recovery / 'native_2023_snapshot.pkl.gz', 'wb', compresslevel=1) as f: pickle.dump(saved, f, protocol=pickle.HIGHEST_PROTOCOL)
        save(recovery / 'rows.json', result.rows); save(recovery / 'fertility.json', observations)
        save(recovery / 'verification.json', dict(status='PASS', observed_year=2023, root_solves=0, forward_dates=count, backward_dates=backwards, numeric_max_abs_gap=maxgap, tolerance=TOL, source_mapping=str(rows_path), source_mapping_sha256=sha(rows_path), model_snapshot_sha256=sha(recovery/'native_2023_snapshot.pkl.gz')))
        dispatch_policy(m, recovery)
    except BaseException as exc:
        save(recovery / 'failure.json', dict(error_type=type(exc).__name__, error=str(exc))); raise
    finally: stop.set()


def dispatch_policy(m, recovery):
    """Write the tax contract only after the native replay is PASS, then sbatch it."""
    proof = read(recovery / 'verification.json')
    if proof.get('status') != 'PASS' or proof.get('root_solves') != 0: raise ValueError('Recovery verification required')
    policy = Path(m['policy_output']); policy.mkdir(parents=True, exist_ok=False)
    announced = read(m['announced_manifest'])
    pins = [Path(m['source_driver']), Path(m['spec']), recovery/'native_2023_snapshot.pkl.gz', recovery/'verification.json', recovery/'rows.json', recovery/'seed.json', recovery/'baseline_rows.json', Path(announced['endpoint_pickle']), Path(announced['endpoint_receipt']), Path(m['shared_tax_endpoint']), Path(m['shared_tax_endpoint_receipt'])]
    pm = dict(spec=m['spec'], recovery=str(recovery), seed=str(recovery/'seed.json'), output=str(policy), endpoint=announced['endpoint_pickle'], endpoint_receipt=announced['endpoint_receipt'], annual_taxes=[.02], start_year=2023, periods=100, psi=FINAL_PSI, seed_horizon=104,
        shared_tax_endpoint=m['shared_tax_endpoint'], shared_tax_endpoint_receipt=m['shared_tax_endpoint_receipt'],
        source_mapping=proof['source_mapping'], source_mapping_sha256=proof['source_mapping_sha256'],
        max_path_evaluations=16, terminal_max_evaluations=24, terminal_seconds=3600, smoke_seconds=1500, path_seconds=8*3600, total_seconds=10*3600,
        expected_mapping_seconds=1900, expected_hours_per_arm=[4,9],
        file_sha256={str(p): sha(p) for p in pins}, closure=dict(preference='Four shocks announced in 2007; final fitted level thereafter', information='Unexpected 2% equally rebated property tax in 2023', initial_distribution='Recovered native announced-four-shock 2023 pre-choice cross-section', population='Original four-vintage birth queue, births/2.1; no immigration or rescaling', fiscal='Equal property-tax rebate to household heads; balanced PAYGO at fixed payroll .179', housing='Fixed calibrated supply curve in asset prices, elasticity .63', production_eligible=False),
        objects=dict(structural_parameters='estimated, retained 2007 calibration', preference_levels='estimated in historical exercise, fixed announced path', housing_supply_scale='empirically normalized at calibration, unchanged', housing_supply_elasticity='externally fixed .63', birth_to_entry_conversion='externally fixed 1/2.1', entry_distribution='retained calibrated entry distribution', migration='externally fixed zero', annual_property_tax='externally fixed policy value .02', payroll_tax='externally fixed .179', pension='endogenous balanced amount', rebate='endogenous equal-per-head balanced amount', inherited_history='outstanding: announced four-shock path convergence is not established', terminal_approach_and_horizon='outstanding until assessed'))
    mp = policy/'manifest.json'; save(mp, pm)
    folder = policy/'tax2'; folder.mkdir(); argv = [PYTHON, m['source_driver'], '--run', str(mp), '--tax', '.02']
    script = folder/'run.sbatch'; script.write_text('\n'.join(['#!/bin/bash', '#SBATCH --job-name=e5f_fourshock_tax2_100', '#SBATCH --account=torch_pr_570_general', '#SBATCH --cpus-per-task=1', '#SBATCH --mem=32G', '#SBATCH --time=610', f'#SBATCH --output={folder}/slurm_%j.log', 'set -euo pipefail', 'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1', shlex.join(argv), '']))
    r = subprocess.run(['sbatch', '--parsable', str(script)], text=True, capture_output=True, check=True)
    receipt = dict(job_id=int(r.stdout.strip().split(';')[0]), tax=.02, manifest_sha256=sha(mp), script_sha256=sha(script), recovery_verification_sha256=sha(recovery/'verification.json'))
    save(policy/'dispatch.json', receipt); save(folder/'submission.json', receipt)


def main():
    p = argparse.ArgumentParser(); g = p.add_mutually_exclusive_group(required=True); g.add_argument('--prepare', type=Path); g.add_argument('--run', type=Path); a = p.parse_args()
    prepare(a.prepare) if a.prepare else run(a.run)


if __name__ == '__main__': main()
