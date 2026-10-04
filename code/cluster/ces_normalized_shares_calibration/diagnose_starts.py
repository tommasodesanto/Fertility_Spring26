"""Four fixed-price single-coordinate diagnostics; no GE search or target scoring.

Run only in the authenticated CES staging container with its dependency overlay.
Each observation uses the production native observer and standard 17-plot packet.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import threading
import time

SEEDS = (.12, .10, .08, .06)
PRICE = .7760569760205563
WRITE_LOCK = threading.Lock()


def write(path, value):
    path = Path(path)
    with WRITE_LOCK:
        temporary = path.with_suffix(path.suffix + '.tmp')
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
        temporary.replace(path)


@contextlib.contextmanager
def deadline_alarm(end):
    def expired(_signum, _frame):
        raise TimeoutError('Native 300-second stage deadline exceeded during diagnostics')
    remaining = end - time.time()
    if remaining <= 0:
        raise TimeoutError('Native stage deadline exhausted before diagnostics')
    prior = signal.signal(signal.SIGALRM, expired)
    timer = signal.setitimer(signal.ITIMER_REAL, remaining)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, *timer)
        signal.signal(signal.SIGALRM, prior)


def standard_packet(context, live, destination):
    """Same native packet construction as observe_price, using normalized live P."""
    import numpy as np
    cal = context['prepared'].rt['primitive'].pf.calendar
    P, grid, sd, sol = (live[k] for k in ('P', 'b_grid', 'sd', 'sol'))
    price = np.asarray(live['price'])
    policy = cal.policy_from_solution(sol, price, P, grid, sd)
    pre, _ = cal.reconstruct_stationary_pre_fertility(sol, policy, P, grid, sd)
    supply = cal.HousingSupplyRule('static-elastic', float(price[0]),
        float(P.H0[0] * (P.user_cost_rate * price[0] / P.r_bar[0]) ** P.xi_supply[0]),
        float(P.xi_supply[0]))
    ev = cal.evaluate_period(price, pre, P, grid, sd, cal.SolveCounter(),
                            supply_rule=supply, supplied_policy=policy)
    packet = dict(parameters=P, b_grid=grid, shared=sd, solution=sol,
                  evaluation=ev, stationary_g_pre=pre, supply_rule=supply,
                  demographic_seed=context['reference'].get('demographic_seed'))
    context['prepared'].rt['audit'].standard_diagnostics(
        packet, destination, validate_production_young=False)
    names = sorted(p.name for p in (destination / 'standard_diagnostics').glob('*.png'))
    if len(names) != 17 or names != sorted(context['manifest']['standard_diagnostic_names']):
        raise RuntimeError('Diagnostic packet must retain exact native 17 plot names')
    return names


def run(args):
    if os.environ.get('CES_NORMALIZED_SHARES_STAGED_CONTEXT') != '1':
        raise RuntimeError('Requires authenticated CES staged container and dependency overlay')
    for key in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        if os.environ.get(key) != '1':
            raise RuntimeError('Single-core environment required: ' + key)
    raw = args.starts_file.read_bytes()
    pin = hashlib.sha256(raw).hexdigest()
    if pin != args.starts_file_sha256:
        raise RuntimeError('Start-plan SHA-256 mismatch')
    plan = json.loads(raw)
    point = dict(plan['starts'][0])
    if len(point) != 11 or set(point) != set(plan['coordinates']):
        raise RuntimeError('Expected exact eleven-coordinate seed0')
    from experiments.ces_normalized_shares import adapter
    from production import equilibrium, inputs, native_price, native_phase_b
    if inputs.DEFAULT_PRICE != PRICE:
        raise RuntimeError('Native default price drift')
    spec = adapter.contract()
    if plan['target_fingerprint'] != spec['target_fingerprint'] or plan['weight_fingerprint'] != spec['weight_fingerprint']:
        raise RuntimeError('CES staged target contract drift')
    coordinate = args.coordinate
    seeds = tuple(float(value) for value in args.seed_values)
    if coordinate not in ('psi_child', 'first_birth_fixed_cost') or coordinate not in point:
        raise ValueError('Unsupported diagnostic coordinate')
    if len(seeds) != 4 or any(not math.isfinite(value) for value in seeds):
        raise ValueError('Exactly four finite seed values required')
    lo, hi = spec['bounds'][coordinate]
    if any(not float(lo) <= value <= float(hi) for value in seeds):
        raise ValueError('Diagnostic seed outside authenticated coordinate bounds')
    started = time.time()
    end = min(float(args.deadline_epoch), started + 2700.)
    if not math.isfinite(end) or end <= started:
        raise ValueError('Positive explicit deadline required')
    args.out.mkdir(parents=True, exist_ok=False)
    completed = []
    stop = threading.Event()
    def heartbeat():
        while not stop.wait(60):
            write(args.out / 'heartbeat.json', dict(time_epoch=time.time(), deadline_epoch=end,
                  completed_cases=len(completed), status='running_diagnostic_only'))
    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    write(args.out / 'plan.json', dict(status='fixed_price_diagnostic_only', no_adoption=True,
        no_GE_search=True, no_optimizer=True, target_scoring=False, price=PRICE, coordinate=coordinate, seed_values=seeds,
        base_seed0=point, starts_file_sha256=pin, deadline_epoch=end,
        total_seconds_cap=2700, per_case_budget_seconds=1200, stage_deadline_seconds=300,
        maximum_lifecycle_solves=4, zero_residual_definition='adjusted births / (2.1 * actual entry) - 1 = 0; birth renewal, not SMM loss',
        economic_changes=f'Only {coordinate} differs from normalized-CES seed0; all other coordinates and fixed inputs retained; diagnostic, not adopted',
        dependency_overlay_inventory_sha256=os.environ.get('CES_DEPENDENCY_OVERLAY_INVENTORY_SHA256')))
    try:
        for index, value in enumerate(seeds):
            if time.time() >= end:
                raise TimeoutError('Global diagnostic deadline exhausted')
            label = f'selected_diagnostic_{index}'  # retains native solution arrays
            case = args.out / label
            parameters = {**point, coordinate: value}
            P, grid = adapter.load_inputs(parameters)
            case_end = min(end, time.time() + 1200.)
            with adapter.install():
                context = equilibrium.build_context(P, grid, case, price_start=PRICE,
                    deadline=case_end, max_lifecycle=2, closure='population_one')
                if (context['selected_d_bar'] != float(P.unsecured_credit_limit)
                    or context['reference_psi'] != float(P.psi_child)
                    or context['closure_mode'] != 'population_one'
                    or not callable(context['write_json'])):
                    raise RuntimeError('Actual reporting context identity drift')
                budget = equilibrium.Budget(case, case_end, 2)
                live = native_price.solve_fixed_price(context, float(P.unsecured_credit_limit),
                    PRICE, budget, label, case / 'stage')
                observed = native_phase_b._observe_with_deadline(context, live, label, final=False)
                write(case / 'renewal_observation.json', observed)
                with deadline_alarm(min(live['case_deadline_epoch'], case_end)):
                    plots = standard_packet(context, live, case / 'diagnostic_packet')
                row = dict(status='fixed_price_diagnostic_only', coordinate=coordinate, value=value,
                    psi_child=float(parameters['psi_child']), parameters=parameters,
                    effective_input_fingerprint=adapter.canonical.effective_input_fingerprint(P, grid),
                    observed=observed, native_summary=live['summary'], plot_names=plots,
                    lifecycle_solves=budget.used_lifecycle,
                    distance_to_zero=abs(float(observed['renewal_residual'])))
                write(case / 'observed.json', row)
                completed.append(row)
                write(args.out / 'latest_completed.json', dict(completed=completed, deadline_epoch=end))
                write(args.out / 'best_so_far.json', min(completed, key=lambda r:r['distance_to_zero']))
                write(args.out / 'heartbeat.json', dict(time_epoch=time.time(), completed_cases=len(completed), status='completed_case'))
        write(args.out / 'completed.json', dict(status='four_fixed_price_diagnostics_completed', cases=completed,
            lifecycle_solves=sum(r['lifecycle_solves'] for r in completed), no_GE_certification=True))
    except Exception as error:
        write(args.out / 'failed.json', dict(status='stopped_no_retry', error_type=type(error).__name__,
            error=str(error), completed_cases=len(completed), time_epoch=time.time()))
        raise
    finally:
        stop.set()
        thread.join(timeout=2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coordinate', choices=('psi_child', 'first_birth_fixed_cost'), default='psi_child')
    parser.add_argument('--seed-values', type=float, nargs=4, default=SEEDS)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--starts-file', type=Path, required=True)
    parser.add_argument('--starts-file-sha256', required=True)
    parser.add_argument('--deadline-epoch', type=float, required=True)
    run(parser.parse_args())


if __name__ == '__main__':
    main()
