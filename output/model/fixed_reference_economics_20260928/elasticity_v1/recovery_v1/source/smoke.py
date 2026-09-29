#!/usr/bin/env python3
"""Zero-lifecycle smoke for recovery authentication and orchestration."""
import gzip
import json
import os
from pathlib import Path
import tempfile
import time
from types import SimpleNamespace

import run_recovery as r


def check(ok, message):
    if not ok:
        raise AssertionError(message)


def checkpoint_tests():
    with tempfile.TemporaryDirectory() as d:
        dest = Path(d)/'state.pkl.gz'
        expected = {'full_packet': {'arrays': list(range(100)), 'policy': 'retained'}}
        digest = r.atomic_checkpoint(expected, dest)
        check(r.base.sha(dest) == digest and gzip.open(dest, 'rb').read(), 'Atomic checkpoint missing')
        import pickle
        with gzip.open(dest, 'rb') as f:
            check(pickle.load(f) == expected, 'Atomic checkpoint changed packet')
    class Interrupted:
        def __reduce__(self):
            def fail():
                raise RuntimeError('interrupted pickle write')
            return (fail, ())
    # A source iterator that raises during pickling leaves no published final file.
    class Broken:
        def __reduce_ex__(self, protocol):
            raise RuntimeError('synthetic interruption')
    with tempfile.TemporaryDirectory() as d:
        dest = Path(d)/'state.pkl.gz'
        try:
            r.atomic_checkpoint(Broken(), dest)
        except RuntimeError as exc:
            check('synthetic interruption' in str(exc), 'Wrong interrupted-write branch')
        else:
            raise AssertionError('Interrupted checkpoint unexpectedly passed')
        check(not dest.exists(), 'Incomplete checkpoint was published')


def controller_tests(plan):
    original_verify, original_popen, original_comparison = r.verify_plan, r.subprocess.Popen, r.make_comparison
    original_kill = r.os.killpg
    original_time, original_sleep = r.time.time, r.time.sleep
    p = dict(plan)
    p['total_seconds'] = 2100
    p['case_seconds'] = 900
    calls = []
    r.verify_plan = lambda path, full_inputs=False: p
    r.make_comparison = lambda out, records, p: dict(comparison_sha256='test-only',
        elasticities_sha256='test-only', outcome_rows=84, elasticity_rows=28)
    class FakeProcess:
        pid = 99999999
        def __init__(self, command, stdout, stderr, env, start_new_session):
            case = command[command.index('--case')+1]
            out = Path(command[command.index('--output')+1])
            calls.append(case)
            if mode == 'passed':
                out.mkdir()
                (out/'conditional_cohort_state.pkl.gz').write_bytes(b'synthetic-only')
                r.base.write(out/'receipt.json', dict(status='passed', lifecycle_solves=1,
                    plan_sha256=r.base.sha(plan_path), standard_plot_count=17,
                    checkpoint=dict(sha256=r.base.sha(out/'conditional_cohort_state.pkl.gz')),
                    cohort_gates=dict(relative_market_residual=0.01)))
                self.returncode = 0
            elif mode == 'failed_child':
                self.returncode = 1
            else:
                self.returncode = None
        def poll(self): return self.returncode
        def wait(self, timeout=None):
            self.returncode = -15
            return self.returncode
    r.subprocess.Popen = FakeProcess
    r.os.killpg = lambda pid, sign: None
    try:
        with tempfile.TemporaryDirectory() as d:
            global plan_path, mode
            plan_path = Path(d)/'plan.json'; plan_path.write_text('{}')
            os.environ['RECOVERY_STARTED_EPOCH'] = str(time.time())
            os.environ['SLURM_JOB_ID'] = '999999'
            mode = 'passed'
            r.controller(SimpleNamespace(output=Path(d)/'passed', plan=plan_path), p)
            done = r.base.read(Path(d)/'passed/completed.json')
            check(calls == ['reference_1010','credit_1010'] and done['complete_three_price'] and
                  not done['complete_five_price'] and done['new_lifecycle_solves'] == 2,
                  'Two-case controller did not complete exactly')
            calls.clear(); mode = 'failed_child'
            try: r.controller(SimpleNamespace(output=Path(d)/'failed', plan=plan_path), p)
            except RuntimeError as exc: check('no retry' in str(exc), 'Wrong failed-child branch')
            else: raise AssertionError('Failed child was accepted')
            check(calls == ['reference_1010'], 'Failed child retried or next case started')
            calls.clear(); mode = 'timeout'
            os.environ['RECOVERY_STARTED_EPOCH'] = str(time.time()-2040)
            try: r.controller(SimpleNamespace(output=Path(d)/'budget', plan=plan_path), p)
            except RuntimeError as exc: check('budget' in str(exc), 'Wrong budget branch')
            else: raise AssertionError('Expired budget was accepted')
            check(not calls, 'Expired budget started a child')
            calls.clear(); mode = 'timeout'
            clock = [original_time()]
            os.environ['RECOVERY_STARTED_EPOCH'] = str(clock[0])
            r.time.time = lambda: clock[0]
            r.time.sleep = lambda _: clock.__setitem__(0, clock[0]+901)
            try: r.controller(SimpleNamespace(output=Path(d)/'timeout', plan=plan_path), p)
            except TimeoutError as exc: check('deadline' in str(exc), 'Wrong timed-out child branch')
            else: raise AssertionError('Timed-out child was accepted')
            check(calls == ['reference_1010'], 'Timed-out child retried or next case started')
    finally:
        r.verify_plan, r.subprocess.Popen, r.make_comparison = original_verify, original_popen, original_comparison
        r.os.killpg = original_kill
        r.time.time, r.time.sleep = original_time, original_sleep


if __name__ == '__main__':
    plan_path = Path('/work/recovery_results/plan_recovery.json')
    p = r.verify_plan(plan_path, full_inputs=True)
    real_plan_digest = r.base.sha(plan_path)
    checkpoint_tests()
    controller_tests(p)
    print(json.dumps(dict(status='passed', zero_lifecycle_solves=True,
        six_prior_inputs_authenticated=True, full_plan_sha256=real_plan_digest,
        test_cases=['two_case_flow','failed_child_no_retry','budget_deadline','atomic_checkpoint','interrupted_write'])))
