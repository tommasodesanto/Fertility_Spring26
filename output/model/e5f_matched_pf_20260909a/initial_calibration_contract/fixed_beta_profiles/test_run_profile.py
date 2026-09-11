import contextlib
import copy
import io
import json
import math
import tempfile
from pathlib import Path
import unittest
import numpy as np
import run_profile as p


def restrictions():
    return {n: dict(lower=.01, upper=100.) for n in p.ALL_NAMES}


def score(item, beta):
    # Smooth 12-row objective with eight independently informative directions.
    x = np.log([item['parameters'][n] for n in p.FREE_NAMES])
    r = np.r_[x + np.linspace(.01, .08, 8), np.ones(4) * .1]
    fits = [dict(restriction_id='initial_normalization', target=2.1, model=2.1, gap=0., scored=False)]
    fits += [dict(restriction_id=f'target_{i}', target=0., model=float(v), gap=float(v),
                  actual_weight=1., loss_contribution=float(v*v), scored=True) for i, v in enumerate(r)]
    rows = [dict(parameter=n, estimate=item['parameters'][n], lower=.01, upper=100., status='raw original') for n in p.ALL_NAMES]
    # Deliberately different final and input psi: repetition must retain the INPUT.
    rows.append(dict(parameter='psi_child', estimate=item['initial_psi'] + .123, status='normalization'))
    return dict(contract_sha256=p.APPROVED_OBJECTIVE, free_parameter_count=9,
                loss=math.fsum(float(v*v) for v in r), target_fit=fits, parameters=rows, normalization={'completed_fertility': 2.1})


def plan(beta=.98):
    return dict(fixed_beta=beta, workers=8, rounds=2, maximum_search_cases=56,
                search_seconds=100000, resume_proposal=dict(parameters={n: 1. for n in p.ALL_NAMES}, initial_psi=.2))


class ProfileTests(unittest.TestCase):
    def execute(self, plan_value=None, reject=None, seed_only=False):
        config = plan() if plan_value is None else plan_value
        inherited = score(dict(parameters={n: 1. for n in p.ALL_NAMES}, initial_psi=.2), .999)
        inherited['loss'] = -100.  # Impossible best: proves inherited loss is not eligible.
        calls = []
        def worker(item):
            calls.append(copy.deepcopy(item))
            status = 'rejected_equilibrium' if reject and reject(item) else 'verified'
            s = score(item, config['fixed_beta'])
            return dict(case_id=item['case_id'], proposal=item, status=status, score=s,
                        loss=s['loss'], output='/nonexistent/synthetic')
        with tempfile.TemporaryDirectory() as temp, contextlib.redirect_stdout(io.StringIO()):
            out = Path(temp)
            summary, best = p.search(config, inherited, restrictions(), worker, out, seed_only=seed_only)
            jacobians = [json.loads(f.read_text()) for f in out.glob('jacobian_round_*.json')]
            p.save_profile_tables(out, best, config['fixed_beta'])
            overlay = p.read(out/'selected_profile_score.json') if best else None
        return summary, best, calls, jacobians, overlay

    def test_exact_two_round_production_loop_both_betas(self):
        for beta in (.98, .99):
            summary, best, calls, js, overlay = self.execute(plan(beta))
            self.assertEqual(summary['status'], 'completed')
            self.assertEqual(summary['search_attempted_cases'], 56)
            self.assertEqual(len(calls), 58)  # fixed seed + 56 candidates + 2-repetition call
            self.assertTrue(summary['selected_exact_repetitions_verified'])
            self.assertEqual(calls[0]['case_id'], 'fixed_beta_seed')
            self.assertEqual(calls[-1]['repetitions'], 2)
            for item in calls:
                self.assertEqual(set(item['parameters']), set(p.ALL_NAMES))
                self.assertEqual(item['parameters']['beta_annual'], beta)
            self.assertEqual(len(js), 2)
            for j in js:
                self.assertEqual(np.shape(j['weighted_jacobian']), (12, 8))
                self.assertEqual(j['names'], list(p.FREE_NAMES))
            self.assertGreater(best['loss'], -100.)
            self.assertEqual(overlay['free_parameter_count'], 8)
            self.assertEqual(overlay['raw_scorer_free_parameter_count'], 9)
            self.assertEqual(best['score']['free_parameter_count'], 9)
            b = next(r for r in overlay['parameters'] if r['parameter']=='beta_annual')
            self.assertFalse(b['profile_free_parameter']); self.assertEqual(b['profile_fixed_value'], beta)
            # Find selected original proposal by exact structural values and ensure input psi matches.
            original = next(c for c in calls[:-1] if c['parameters'] == calls[-1]['parameters'])
            self.assertEqual(calls[-1]['initial_psi'], original['initial_psi'])

    def test_one_rejected_derivative_continues_with_one_sided_column(self):
        summary, _, _, js, _ = self.execute(reject=lambda i: i['case_id']=='r0_d0_-1')
        self.assertEqual(summary['status'], 'completed'); self.assertEqual(len(js), 2)

    def test_rejected_primary_seed_uses_two_bounded_fixed_seeds(self):
        summary, _, calls, _, _ = self.execute(reject=lambda i: i['case_id']=='fixed_beta_seed')
        self.assertEqual(summary['status'], 'completed')
        self.assertEqual(len([c for c in calls if c['case_id'].startswith('fixed_beta_seed')]), 3)
        self.assertEqual(summary['attempted_cases'], 60)

    def test_all_seed_rejections_never_select_unrestricted(self):
        summary, best, calls, _, _ = self.execute(reject=lambda i: i['case_id'].startswith('fixed_beta_seed'))
        self.assertEqual(summary['status'], 'stopped_for_review'); self.assertIsNone(best)
        self.assertEqual(len(calls), 3); self.assertIsNone(summary['best_loss'])

    def test_no_time_for_stage_still_repeats_fixed_seed(self):
        config = plan(); config['search_seconds'] = 1
        summary, _, calls, _, _ = self.execute(config)
        self.assertEqual(len(calls), 2); self.assertTrue(summary['selected_exact_repetitions_verified'])
        self.assertEqual(summary['stop_reason'], 'stage_deadline_or_case_budget')

    def test_seed_only_reuses_exact_first_stage(self):
        summary, _, calls, _, _ = self.execute(seed_only=True)
        self.assertEqual(summary['status'], 'seed_smoke_passed'); self.assertEqual(len(calls), 1)
        self.assertFalse(summary['selected_exact_repetitions_verified'])

    def test_wrong_beta_raw_metadata_or_normalization_is_rejected(self):
        item = dict(parameters=p.full_parameters({n:1. for n in p.ALL_NAMES}, .98), initial_psi=.2)
        raw = score(item, .98); p.validate_score(raw, item, .98)
        for change in ('beta', 'metadata', 'normalization'):
            modified = copy.deepcopy(raw)
            if change == 'beta': modified['parameters'][0]['estimate'] = .99
            if change == 'metadata': modified['free_parameter_count'] = 8
            if change == 'normalization': modified['target_fit'][0]['gap'] = .001
            with self.assertRaises(ValueError): p.validate_score(modified, item, .98)

    def test_feasible_bounds_and_no_fabricated_column(self):
        self.assertEqual(p.feasible_steps(-.02, 0), [(-1, -.02)])
        center=np.arange(12.); slope=np.ones(12)
        np.testing.assert_allclose(p.derivative_column([(.02,center+.02*slope)],center),slope)
        with self.assertRaises(RuntimeError): p.derivative_column([],center)

    def test_repetition_signature_ignores_timing_but_checks_economics(self):
        item=dict(parameters=p.full_parameters({n:1. for n in p.ALL_NAMES},.98),initial_psi=.2)
        first=score(item,.98); second=copy.deepcopy(first)
        first['normalization']['stationary_solve_seconds']=[10.,20.]
        second['normalization']['stationary_solve_seconds']=[12.,24.]
        self.assertEqual(p.numeric_signature(first),p.numeric_signature(second))
        second['normalization']['completed_fertility'] += 1e-8
        self.assertNotEqual(p.numeric_signature(first),p.numeric_signature(second))

    def test_two_wave_limit(self):
        with self.assertRaises(ValueError): p.batch(list(range(17)),lambda x:x,lambda x:None,8)

    def test_only_exact_preflighted_housing_error_is_recoverable(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t); (d/'raw').mkdir()
            f=dict(error_type='RuntimeError',phase='stationary_equilibrium',error='Initial housing equilibrium failed its unchanged strict gate')
            p.write(d/'raw/failure.json',f)
            self.assertEqual(p.failure_status(d),'failed')
            p.write(d/'preflight.json',{})
            self.assertEqual(p.failure_status(d),'rejected_equilibrium')
            f['error']='Source mismatch';p.write(d/'raw/failure.json',f)
            self.assertEqual(p.failure_status(d),'failed')


if __name__ == '__main__': unittest.main()
