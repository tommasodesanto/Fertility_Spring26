#!/usr/bin/env python3
"""Bounded full-objective E5F search in a frozen scientific snapshot.

The population and polishing steps only propose points. Every accepted objective
comes from a separate, hash-pinned full historical evaluation with twelve targets.
"""
from __future__ import annotations
import argparse
import copy
import csv
import json
import math
import os
from pathlib import Path
import random
import shutil
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED

import run_e5f_joint_overnight_case as adapter
import build_e5f_bounded_refinement_plan as planner
planner.adapter = adapter
ROOT = Path(__file__).resolve().parents[3]
N = 11


def transform(u, spec):
    lo, hi, kind = spec['lower'], spec['upper'], spec['transform']
    if kind == 'log': return lo * (hi / lo) ** u
    if kind == 'discount': return lo + (hi-lo) * (1-(1-u)**2)
    if kind == 'softzero': return lo + (hi-lo) * u**2
    if kind == 'asinh': return math.sinh((1-u)*math.asinh(lo)+u*math.asinh(hi))
    raise ValueError(kind)


def unit_to_payload(seed, domain, unit, label):
    if len(unit) != N or any(not math.isfinite(x) or not 0 <= x <= 1 for x in unit):
        raise ValueError('Invalid eleven-dimensional proposal')
    # Preserve the original normalization-coordinate origin, never copy fitted losses.
    out = {'old_psi_child': seed['old_psi_child'], 'best_candidate': {
        'theta': copy.deepcopy(seed['best_candidate']['theta']),
        'old_psi_child': seed['old_psi_child'],
        'new_psi_child': seed['best_candidate']['new_psi_child']}}
    for spec, u in zip(domain, unit):
        name, value = spec['name'], transform(float(u), spec)
        if name == 'beta_annual': out['best_candidate']['theta']['beta'] = value**4
        elif name == 'psi_child_change_2023': out['best_candidate']['new_psi_child'] = out['old_psi_child']+value
        else: out['best_candidate']['theta'][name] = value
    out['proposal'] = {'label': label, 'unit_vector': list(unit), 'status': 'unevaluated full-objective proposal'}
    return out


def repair(parent, mutant):
    return [min(1., max(0., (p+(0. if x < 0 else 1.))/2 if not 0 <= x <= 1 else x))
            for p, x in zip(parent, mutant)]


def de_trial(pop, losses, rng):
    best = min(range(len(pop)), key=lambda i: (losses[i], i))
    trials = []
    for i, parent in enumerate(pop):
        a, b, c = rng.sample([j for j in range(len(pop)) if j != i], 3)
        f = rng.uniform(.5, .9)
        if rng.random() < .8:
            mutant = [p+f*(pop[best][k]-p)+f*(pop[a][k]-pop[b][k]) for k, p in enumerate(parent)]
        else:
            mutant = [pop[a][k]+f*(pop[b][k]-pop[c][k]) for k in range(N)]
        mutant = repair(parent, mutant)
        forced = rng.randrange(N)
        trials.append([mutant[k] if k == forced or rng.random() < .9 else parent[k] for k in range(N)])
    return trials


def select_population(pop, losses, evaluated):
    """Require one actual evaluated candidate for each deferred population slot."""
    ordered = sorted(evaluated, key=lambda item: item[0]['id'])
    if len(ordered) != len(pop) or [r['id'] for r, _ in ordered] != list(range(1, len(pop)+1)):
        raise RuntimeError('Incomplete evaluated generation')
    for i, (row, unit) in enumerate(ordered):
        if not math.isfinite(row['loss']): raise RuntimeError('Nonfinite evaluated loss')
        if row['loss'] < losses[i]: pop[i], losses[i] = list(unit), row['loss']
    return pop, losses


def polish_directions(center, anchor_loss, evaluated):
    direction, gains = [0.]*N, [0.]*N
    for row, unit in evaluated:
        # Each coordinate row identifies its actual dimension in the label.
        j = int(row['label'].split('_')[1])
        gain = anchor_loss-row['loss']
        if gain > gains[j]: gains[j], direction[j] = gain, unit[j]-center[j]
    order = sorted(range(N), key=lambda j: (-gains[j], j))
    return direction, order


def joint_polish_vectors(center, direction, order, seen):
    vectors, labels = [], []
    for name, subset in [('all', set(range(N))), ('top3', set(order[:3])), ('top6', set(order[:6]))]:
        for scale in (.5, 1., 2., 4.):
            u = [min(1., max(0., x+scale*(direction[j] if j in subset else 0.))) for j, x in enumerate(center)]
            key = tuple(u)
            if key in seen: continue
            seen.add(key); vectors.append(u); labels.append(f'{name}_scale_{scale}')
    return vectors, labels


def has_budget(completed, proposed, maximum=360, repeats=2):
    return completed+proposed+repeats <= maximum


def verify_contract(path, sha):
    adapter.verify(path, sha)
    c = adapter.read_json(path)
    if c.get('schema') != 'e5f_joint_overnight_contract_v1': raise RuntimeError('Unknown contract')
    adapter.verify(__file__, c['controller_sha256'])
    adapter.verify(adapter.__file__, c['case_adapter_sha256'])
    adapter.verify(planner.__file__, c['planner_sha256'])
    adapter.verify(c['seed_center'], c['seed_center_sha256'])
    adapter.verify(c['seed_reference'], c['seed_reference_sha256'])
    expected = dict(max_workers=12, case_timeout_seconds=1800, max_histories=360,
                    max_search_seconds=27000, max_total_seconds=30600,
                    max_generations=8, population_size=32, initial_radius=.01,
                    polish_rounds=2, polish_radius=.0025)
    if any(c.get(k) != v for k, v in expected.items()): raise RuntimeError('Changed overnight budget/design')
    if c['search_domain'] != adapter.SEARCH_DOMAIN: raise RuntimeError('Changed parameter restrictions')
    for key in ('source_sha256', 'target_fingerprint', 'code_bundle_sha256'):
        if c[key] != c['base_plan'][key]: raise RuntimeError(f'Mixed {key}')
    return c


def copy_reference(summary_path, dest):
    s = adapter.read_json(summary_path)
    dest.mkdir(parents=True, exist_ok=False)
    for name in ('summary.json', 'target_fit_long.csv', 'parameter_table.csv'):
        shutil.copy2(summary_path.parent/name, dest/name)
    name = s['best_candidate']['candidate']
    (dest/'cases'/name).mkdir(parents=True)
    shutil.copy2(summary_path.parent/'cases'/name/'transition_path.csv', dest/'cases'/name/'transition_path.csv')
    return dest/'summary.json'


def verify_graph_match(reference_case, candidate_case):
    graphs = sorted((reference_case/'standard_diagnostics').glob('*.png'))
    if len(graphs) != 17: raise RuntimeError('Missing reference graph packet')
    for graph in graphs: adapter.verify(candidate_case/'standard_diagnostics'/graph.name, adapter.digest(graph))
    return 17


class Search:
    def __init__(self, contract, mode):
        self.c, self.mode = contract, mode
        self.root = Path(contract['output_root'])/('smoke' if mode == 'smoke' else 'search')
        if self.root.exists(): raise RuntimeError(f'Refusing existing run directory: {self.root}')
        self.root.mkdir(parents=True)
        self.started, self.wall_start = time.monotonic(), time.time()
        self.finish_epoch = min(self.wall_start+contract['max_total_seconds'], contract['absolute_finish_epoch'])
        self.search_epoch = min(self.wall_start+contract['max_search_seconds'], self.finish_epoch-3600)
        self.completed = 0 if mode == 'smoke' else 2
        self.ledger, self.best, self.phase = [], None, 'initializing'
        self.stop = threading.Event()
        self.active, self.lock = {}, threading.Lock()
        self.seed = adapter.read_json(contract['seed_center'])
        self.write_state('running')

    def write_state(self, status, **extra):
        adapter.write_json(self.root/'search_state.json', dict(
            phase=self.phase, status=status, elapsed_seconds=time.monotonic()-self.started,
            epoch=time.time(), completed_histories=self.completed, maximum_histories=self.c['max_histories'],
            best_loss=self.best[0]['loss'] if self.best else None,
            absolute_finish_epoch=self.finish_epoch, production_promoted=False, **extra))

    def kill_owned(self):
        self.stop.set()
        with self.lock: active = list(self.active.values())
        for proc in active:
            if proc.poll() is None:
                try: os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError: pass
        for proc in active:
            try: proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try: os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError: pass
                proc.wait(timeout=5)

    def stage_fits(self, count, repeats=False):
        if not has_budget(self.completed, count, self.c['max_histories'], 0 if repeats else 2): return False
        cutoff = self.finish_epoch if repeats else self.search_epoch
        elapsed = [r['elapsed_seconds'] for r, _ in self.ledger]
        expected = max(600., sorted(elapsed)[len(elapsed)//2]*1.4) if elapsed else 840.
        batches = math.ceil(count/self.c['max_workers'])
        # Include a full timeout after the last launch, and protect verification time.
        return time.time()+max(self.c['case_timeout_seconds'], batches*expected) <= cutoff

    def new_plan(self, stage, vectors=None, labels=None, repeat=None, bridge=False):
        dest = self.root/stage
        dest.mkdir(parents=True, exist_ok=False)
        plan = copy.deepcopy(self.c['base_plan'])
        plan.update(stage=stage, cases=[], input_sha256={},
                    launch_deadline_epoch=self.finish_epoch-1800 if repeat else self.search_epoch,
                    controller_sha256=self.c['controller_sha256'])
        if repeat:
            row, _ = repeat
            source_plan_path = Path(row['plan'])
            source_plan = adapter.load_plan(source_plan_path, row['plan_sha256'])
            original = next(x for x in source_plan['cases'] if x['id'] == row['id'])
            reference = copy_reference(Path(row['summary']), dest/'reference')
            for i in (1, 2):
                center = dest/f'center_{i:03d}.json'
                shutil.copy2(source_plan_path.parent/original['center'], center)
                case = copy.deepcopy(original)
                for k in ('bridge_reference', 'bridge_reference_sha256'): case.pop(k, None)
                case.update(id=i, label=f'exact_repeat_{i}', center=center.name,
                    center_sha256=adapter.digest(center), output=f'task_{i:03d}',
                    reference=str(reference.relative_to(dest)), reference_sha256=adapter.digest(reference))
                plan['cases'].append(case)
            plan['repeat_origin'] = {'plan_sha256': row['plan_sha256'], 'case_id': row['id'],
                                     'center_sha256': original['center_sha256']}
        else:
            for i, (unit, label) in enumerate(zip(vectors, labels, strict=True), 1):
                center = dest/f'center_{i:03d}.json'
                adapter.write_json(center, unit_to_payload(self.seed, self.c['search_domain'], unit, label))
                plan['cases'].append(dict(id=i, label=label, center=center.name,
                    center_sha256=adapter.digest(center), panel_task_id=1, panel_size=1,
                    panel_design='mixed', panel_seed=self.c['seed'], radius=.005, output=f'task_{i:03d}'))
            if bridge:
                reference = copy_reference(Path(self.c['seed_reference']), dest/'bridge_reference')
                for case in plan['cases']:
                    case.update(bridge_reference=str(reference.relative_to(dest)),
                                bridge_reference_sha256=adapter.digest(reference))
        for folder in ('reference', 'bridge_reference'):
            if (dest/folder).exists():
                plan['input_sha256'].update({str(p.relative_to(dest)): adapter.digest(p)
                                            for p in (dest/folder).rglob('*') if p.is_file()})
        path = dest/'plan.json'; adapter.write_json(path, plan)
        return path, adapter.digest(path)

    def run_child(self, plan, sha, case):
        if self.stop.is_set(): raise RuntimeError('Run stopped before child launch')
        deadline = adapter.read_json(plan)['launch_deadline_epoch']
        if time.time() > deadline: raise RuntimeError('Case launch deadline exceeded')
        with (plan.parent/f"case_{case['id']:03d}.log").open('w') as log:
            proc = subprocess.Popen([sys.executable, str(Path(adapter.__file__)), '--plan', str(plan),
                '--plan-sha256', sha, '--case-id', str(case['id'])], stdout=log,
                stderr=subprocess.STDOUT, start_new_session=True)
            with self.lock:
                self.active[proc.pid] = proc
                stopping = self.stop.is_set()
            if stopping:
                try: os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError: pass
            try:
                timeout = min(self.c['case_timeout_seconds'], max(1., self.finish_epoch-time.time()))
                returncode = proc.wait(timeout=timeout)
                if returncode: raise RuntimeError(f"Case {case['id']} exited {returncode}")
                return case['id']
            except subprocess.TimeoutExpired:
                try: os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError: pass
                try: proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL); proc.wait(timeout=5)
                raise RuntimeError(f"Case {case['id']} exceeded its declared time limit")
            finally:
                with self.lock: self.active.pop(proc.pid, None)

    def run_batch(self, path, sha):
        plan = adapter.load_plan(path, sha)
        self.phase = plan['stage']; self.write_state('running')
        cases = iter(plan['cases']); futures = {}; collected = set(); rows = []
        ex = ThreadPoolExecutor(max_workers=self.c['max_workers'])
        try:
            def submit_one():
                case = next(cases, None)
                if case is not None:
                    futures[ex.submit(self.run_child, path, sha, case)] = case
            for _ in range(min(len(plan['cases']), self.c['max_workers'])): submit_one()
            while futures:
                done, _ = wait(futures, timeout=60, return_when=FIRST_COMPLETED)
                if not done:
                    self.write_state('running')
                    if time.time() >= self.finish_epoch: raise RuntimeError('Total wall-time budget exhausted')
                    continue
                for future in done:
                    future.result()  # Failure stops the entire dependent search.
                    del futures[future]
                _, _, rows = planner.collect(path, sha, require_complete=False)
                for row in rows:
                    if row['id'] in collected: continue
                    collected.add(row['id']); self.completed += 1
                    s = adapter.read_json(row['summary'])
                    pair = (row, list(s['panel_design']['unit_vector']))
                    self.ledger.append(pair)
                    if self.best is None or row['loss'] < self.best[0]['loss']: self.best = pair
                self.write_reports()
                self.write_state('running')
                for _ in done: submit_one()
            _, _, rows = planner.collect(path, sha, require_complete=True)
            return [(row, adapter.read_json(row['summary'])['panel_design']['unit_vector']) for row in rows]
        except Exception as error:
            self.kill_owned()
            for future in futures: future.cancel()
            self.write_state('failed', error=str(error))
            adapter.write_json(self.root/'failure.json', {'phase': self.phase, 'error': str(error),
                'completed_histories': self.completed, 'production_promoted': False})
            raise
        finally:
            ex.shutdown(wait=True, cancel_futures=True)

    def evaluate(self, stage, vectors, labels):
        if not has_budget(self.completed, len(vectors), self.c['max_histories']):
            raise RuntimeError('Case budget would consume the two verification slots')
        path, sha = self.new_plan(stage, vectors, labels)
        return self.run_batch(path, sha)

    def write_reports(self):
        if self.best is not None:
            adapter.write_json(self.root/'best_so_far.json', {'loss': self.best[0]['loss'],
                'best': self.best[0], 'production_promoted': False})
        if self.ledger:
            adapter.write_json(self.root/'latest_completed_case.json', {'latest': self.ledger[-1][0], 'phase': self.phase})
            planner.write_csv(self.root/'all_cases.csv', [r for r, _ in self.ledger])

    def finish_report(self, status):
        self.write_reports()
        fits, parameters = [], []
        for row, _ in self.ledger:
            base = Path(row['summary']).parent
            for x in adapter.read_csv(base/'target_fit_long.csv'): fits.append({'stage': base.parent.name, 'case_id': row['id'], **x})
            for x in adapter.read_csv(base/'parameter_table.csv'): parameters.append({'stage': base.parent.name, 'case_id': row['id'], **x})
        planner.write_csv(self.root/'all_target_fits.csv', fits)
        planner.write_csv(self.root/'all_parameters.csv', parameters)
        self.write_state(status)
        lines = ['# Overnight full calibration', '', f'Status: {status}. Production remains retained.', '']
        if self.best:
            row = self.best[0]; base = Path(row['summary']).parent
            lines += [f"Best evaluated loss: {row['loss']:.9f}. First-birth rooms: {row['first_birth_rooms']:.6f} against 0.720246.",
                'All eleven parameters were searched; only the first-child upper bound was enlarged to 2.0.',
                'All twelve empirical targets and weights remain unchanged. No global-optimum claim.', '',
                '| Moment | Target | Model | Gap | Weight | Loss |', '|---|---:|---:|---:|---:|---:|']
            for x in adapter.read_csv(base/'target_fit_long.csv'):
                lines.append('| '+x['moment']+' | '+' | '.join(f"{float(x[k]):.8g}" for k in ['target','model','gap','weight','loss_contribution'])+' |')
            lines += ['', '| Parameter | Estimate | Lower | Upper | Status | Near bound |', '|---|---:|---:|---:|---|---|']
            for x in adapter.read_csv(base/'parameter_table.csv'):
                lines.append('| '+' | '.join(x[k] for k in ['parameter','value','lower_bound','upper_bound','status','near_bound'])+' |')
            lines += ['', f"Selected case: `{base}`.", 'The standard seventeen diagnostic plots and complete ledgers accompany this result.']
        (self.root/'MORNING_SUMMARY.md').write_text('\n'.join(lines)+'\n')

    def smoke(self):
        u = self.seed['panel_design']['unit_vector']
        path, sha = self.new_plan('exact_loop', [u, u], ['bridge_1', 'bridge_2'], bridge=True)
        pairs = self.run_batch(path, sha)
        paths = [Path(r['summary']).parent for r, _ in sorted(pairs, key=lambda x:x[0]['id'])]
        exact = adapter.compare_reference(paths[1], paths[0]/'summary.json')
        exact['exact_standard_pngs'] = verify_graph_match(paths[0], paths[1])
        for base in paths:
            reference = adapter.read_json(base/'case_receipt.json')['reference']
            if not reference.get('bridge_twelve_row_fit') or not reference.get('bridge_parameters'):
                raise RuntimeError('Smoke bridge failed')
        adapter.write_json(self.root/'smoke_verification.json', {'status': 'pass', 'plan_sha256': sha, **exact})
        self.finish_report('smoke_passed')

    def require_smoke(self):
        root = Path(self.c['output_root'])/'smoke'
        proof = adapter.read_json(root/'smoke_verification.json')
        if proof['status'] != 'pass' or proof['exact_standard_pngs'] != 17: raise RuntimeError('Exact loop smoke not passed')
        _, _, rows = planner.collect(root/'exact_loop/plan.json', proof['plan_sha256'])
        if len(rows) != 2: raise RuntimeError('Incomplete smoke')
        paths = [Path(r['summary']).parent for r in sorted(rows, key=lambda x:x['id'])]
        adapter.compare_reference(paths[1], paths[0]/'summary.json')
        verify_graph_match(paths[0], paths[1])
        # The starting incumbent is valid even if later exploratory proposals all fail.
        row = rows[0]
        self.best = (row, adapter.read_json(row['summary'])['panel_design']['unit_vector'])

    def search(self):
        self.require_smoke()
        rng = random.Random(self.c['seed']); base = self.seed['panel_design']['unit_vector']
        hidx = next(i for i,d in enumerate(self.c['search_domain']) if d['name']=='hbar_first_child_jump')
        pop = [base]
        for i in range(31):
            u = [min(1., max(0., x+rng.uniform(-.01,.01))) for x in base]
            h = (.50,.60,.75,.90,1.05,1.20)[i%6]+rng.uniform(-1e-8,1e-8)
            u[hidx] = math.sqrt(h/2.)
            pop.append(u)
        if not self.stage_fits(len(pop)): raise RuntimeError('Insufficient time for initial population')
        pairs = self.evaluate('initial_population', pop, ['seed']+[f'initial_{i}' for i in range(1,32)])
        ordered = sorted(pairs, key=lambda x:x[0]['id'])
        pop, losses = [u for _,u in ordered], [r['loss'] for r,_ in ordered]
        stale = 0
        for generation in range(1,9):
            if not self.stage_fits(32): break
            before = self.best[0]['loss']
            trials = de_trial(pop,losses,rng)
            pairs = self.evaluate(f'generation_{generation:02d}',trials,[f'de_{generation}_{i}' for i in range(32)])
            pop, losses = select_population(pop,losses,pairs)
            spread = max(max(u[j] for u in pop)-min(u[j] for u in pop) for j in range(N))
            stale = stale+1 if (before-self.best[0]['loss'])/max(abs(before),1.) < .001 else 0
            adapter.write_json(self.root/'population.json', {'generation': generation, 'unit_vectors': pop,
                'evaluated_losses': losses, 'range_max': spread, 'consecutive_small_improvements': stale})
            if stale >= 3 and spread < .002: break
        for round_ in range(1,3):
            if not self.stage_fits(22): break
            center, anchor_loss = list(self.best[1]), self.best[0]['loss']
            radius = .0025/2**(round_-1)
            probes, labels = [], []
            for j in range(N):
                for sign, suffix in [(-1,'minus'),(1,'plus')]:
                    probes.append([min(1.,max(0.,x+(sign*radius if k==j else 0.))) for k,x in enumerate(center)])
                    labels.append(f'coordinate_{j}_{suffix}')
            pairs = self.evaluate(f'polish_{round_}_coordinates', probes, labels)
            direction, order = polish_directions(center, anchor_loss, pairs)
            seen = {tuple(center), *(tuple(u) for _,u in self.ledger)}
            vectors, labels = joint_polish_vectors(center,direction,order,seen)
            if vectors and self.stage_fits(len(vectors)):
                self.evaluate(f'polish_{round_}_joint',vectors,labels)
        if self.stage_fits(2,repeats=True):
            selected = copy.deepcopy(self.best)
            path, sha = self.new_plan('final_repeats',repeat=selected)
            pairs = self.run_batch(path,sha)
            for row,_ in pairs:
                proof = adapter.read_json(Path(row['summary']).parent/'case_receipt.json')['reference']
                if not proof or not proof.get('exact_twelve_row_fit'): raise RuntimeError('Missing exact repeat receipt')
                verify_graph_match(Path(selected[0]['summary']).parent,Path(row['summary']).parent)
            adapter.write_json(self.root/'final_verification.json', {'status':'pass','exact_repeats':2,
                'exact_standard_pngs_per_repeat':17,'selected':selected[0],'production_promoted':False})
            self.finish_report('complete_verified')
        else:
            self.finish_report('time_budget_exhausted_without_final_repeats')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--contract',type=Path,required=True); p.add_argument('--contract-sha256',required=True)
    p.add_argument('--mode',choices=('smoke','search'),required=True)
    a=p.parse_args(); c=verify_contract(a.contract,a.contract_sha256)
    run=Search(c,a.mode)
    try: getattr(run,a.mode)()
    except BaseException as error:
        run.kill_owned(); run.phase='stopped'; run.finish_report('failed')
        adapter.write_json(run.root/'failure.json', {'error':str(error),'type':type(error).__name__,
            'completed_histories':run.completed,'production_promoted':False})
        raise

if __name__=='__main__': main()
