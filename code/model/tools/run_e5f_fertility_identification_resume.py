#!/usr/bin/env python3
"""Resume a reviewed, terminated controller without rerunning attempted cases.

The pinned original controller still owns proposals, budgets, solves, gates,
selection and exports. This adapter replays authenticated records through its
finish callback and translates only scheduler stage labels in NEW contexts.
No original artifact is edited. Run only after independent lead review on Torch.
"""
from __future__ import annotations
import argparse
import os
from pathlib import Path
from types import SimpleNamespace
import run_e5f_fertility_identification as original

core = original.core
ALIASES = {'search': 'de', 'jacobian': 'initial'}


def inherit_clock(old, smoke, now):
    assert old == smoke, 'Cannot replace original clock'
    assert old['end'] == 1790630594.4615145, 'Wrong experiment deadline'
    assert old['start'] <= now < old['search_cutoff'], 'No late restart'
    return dict(old)


def authenticate(c, objs, manifest, contract_pin, now):
    """Authenticate every attempted record before replay; fail closed on omissions."""
    assert manifest['status'] == 'lead_reviewed_resume'
    assert manifest['contract_sha256'] == contract_pin
    assert manifest['old_job_terminal'] is True
    pins = manifest['pins']
    assert all(core.sha(p['path']) == p['sha256'] for p in pins)
    pinned = {str(Path(p['path']).resolve()) for p in pins}
    old_path = Path(manifest['old_complete']).resolve()
    assert str(old_path) in pinned
    old = core.read(old_path)
    assert old['status'] == 'incomplete_or_fatal_stop'
    assert old['contract_sha256'] == contract_pin
    approval = core.read(manifest['original_approval'])
    assert str(Path(manifest['original_approval']).resolve()) in pinned
    smoke_pin = approval['smoke_receipt']
    assert core.sha(smoke_pin['path']) == smoke_pin['sha256']
    smoke = core.read(smoke_pin['path'])
    inherit_clock(old['clock'], smoke['clock'], now)
    rows = old['records']; assert len(rows) + 6 <= c['budget']['max_objective_cases']
    indexed = {}; classifications = {}; overlays = manifest.get('rejection_overlays', {})
    assert set(overlays) == {r['case'] for r in rows if r['status'] == 'fatal'}, 'Overlay keys must exactly match fatal cases'
    for row in rows:
        key = row['case']; assert key not in indexed, 'Duplicate completed request'
        assert not key.startswith('repeat_'), 'Only pre-repeat recovery supported'
        req_path = Path(row['request_path']).resolve()
        assert str(req_path) in pinned
        req = core.read(req_path)
        assert req['id'] == key and req['point'] == row['point'] and req['lane'] == row['lane'] and req['design'] == row['design']
        assert req['contract_sha256'] == contract_pin
        assert req['scientific_candidate_id'] == core.identity(c, row['lane'], row['point'])
        status, error, data = row['status'], row['error'], {}
        if status == 'success':
            data = core.validate(Path(row['case_path']).parent, c, objs, req)
            assert all(row[k] == v for k, v in data.items())
        elif status == 'fatal':
            # Review must identify the precise classifier-label failure. An
            # unknown fatal failure must never become a searchable population.
            overlay = overlays[key]
            failure_path = req_path.parent / key / 'failure.json'
            assert str(failure_path.resolve()) in pinned
            failure = core.read(failure_path)
            assert failure['context'] == req['context']
            assert failure['error_type'] == 'NonpositiveNormalizedBenefit'
            assert failure['phase'] == 'objective' and failure['status'] == 'fatal'
            assert failure.get('classifier_error') == 'candidate identity/stage required', 'Not the demonstrated classifier-label failure'
            assert row['error'] == failure, 'Recorded failure differs from pinned evidence'
            assert req['context']['stage'] == 'search'
            assert overlay['status'] == 'inadmissible'
            assert overlay['reason'] == 'verified_nonpositive_benefit_stage_label_only'
            status = 'inadmissible'
            error = dict(original_failure=failure, resume_overlay=overlay,
                         original_record_status=row['status'], original_record_path=str(old_path))
        else:
            assert status in ('inadmissible', 'censored_timeout', 'censored_late_completion'), status
        indexed[key] = dict(row=row, request=req, request_path=str(req_path))
        classifications[key] = (status, data, error)
    attempted = {p.stem.removesuffix('.request') for p in old_path.parent.glob('*.request.json')}
    assert attempted == set(indexed), 'Launched request missing completed record: never rerun it'
    assert sum(k.startswith('jacobian_') for k in indexed) == 40
    assert all(v['row']['status'] == 'success' for k,v in indexed.items() if k.startswith('jacobian_'))
    return indexed, classifications


class ReplayScheduler:
    """Adapter retains the original scheduler for every new child process."""
    def __init__(self, replay, classifications, run_batch):
        self.replay = replay
        self.classifications = classifications
        self.run_batch_original = run_batch
        self.used = set()
        self.new_ids = set()

    def classify(self, folder, c, objs, req, proc, code, fallback):
        key = getattr(proc, '_resume_replay_case', None)
        if key is None:
            return fallback(folder, c, objs, req, proc, code)
        assert key in self.used and req == self.replay[key]['request']
        return self.classifications[key]

    def run_batch(self, requests, **kwargs):
        completed = []; pending = []
        for req in requests:
            key = req['id']
            assert key not in self.used and key not in self.new_ids, 'Duplicate proposal dispatch'
            if key not in self.replay:
                req['context']['stage'] = ALIASES.get(req['context']['stage'], req['context']['stage'])
                pending.append(req); self.new_ids.add(key)
                continue
            saved = self.replay[key]; row = saved['row']; payload = saved['request']
            assert all(req[k] == payload[k] for k in ('id','lane','point','design','context')), 'Deterministic replay diverged'
            req['payload'] = payload; req['request_path'] = saved['request_path']
            self.used.add(key)
            proc = SimpleNamespace(deadline=row['deadline'], _resume_replay_case=key)
            completed.append(kwargs['finish'](req, proc, row['returncode']))
        if pending:
            result = self.run_batch_original(pending, **kwargs)
            result['results'] = completed + result['results']
        else:
            result = dict(complete=True, results=completed)
        return result


def main():
    parser = argparse.ArgumentParser()
    for name in ('contract','output','approval','resume-manifest'):
        parser.add_argument('--'+name, type=Path, required=True)
    for name in ('approval-sha256','resume-manifest-sha256','wrapper-sha256'):
        parser.add_argument('--'+name, required=True)
    a = parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID','').isdigit(), 'Torch Slurm only'
    assert not __import__('sys').flags.optimize
    assert core.sha(__file__) == a.wrapper_sha256
    assert core.sha(a.resume_manifest) == a.resume_manifest_sha256
    manifest = core.read(a.resume_manifest)
    assert Path(manifest['original_approval']).resolve() == a.approval.resolve()
    c, objs = original.verify(a.contract)
    replay, classifications = authenticate(c, objs, manifest, core.sha(a.contract), original.time.time())
    saved_module, saved_classify = core.module, core.classify
    adapters = []
    def load(name, path):
        module = saved_module(name, path)
        if name == 'identification_supervisor':
            assert not adapters
            adapter = ReplayScheduler(replay, classifications, module.run_batch)
            adapters.append(adapter); module.run_batch = adapter.run_batch
        return module
    def classify(folder, c, objs, req, proc, code):
        return adapters[0].classify(folder,c,objs,req,proc,code,saved_classify)
    a.stage = 'run'
    core.module = load; core.classify = classify
    try:
        original.controller(a, c, objs)
        assert adapters[0].used == set(replay), 'Not all inherited attempts replayed'
    finally:
        core.module = saved_module; core.classify = saved_classify
        if a.output.exists():
            core.write(a.output/'resume_provenance.json', dict(
                wrapper_sha256=a.wrapper_sha256, manifest_sha256=a.resume_manifest_sha256,
                original_contract_sha256=core.sha(a.contract), original_complete=manifest['old_complete'],
                replayed_cases=sorted(adapters[0].used) if adapters else [],
                new_request_ids=sorted(adapters[0].new_ids) if adapters else [],
                context_aliases=ALIASES, original_artifacts_modified=False))

if __name__ == '__main__':
    main()
