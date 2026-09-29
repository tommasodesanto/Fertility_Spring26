"""Torch-only synthetic checks for the isolated solver transformations.

No Bellman/grid imports: extract the actual transformed choice block and pure
helpers, then check those exact bodies with small numerical arrays. A complete
flag-off fixed-price solve/KFE replay is an additional launch gate owned by the
lead, not something these synthetic checks certify.
"""
import ast
import os
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

assert sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit()
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import patch_solver

ROOT = Path(os.environ.get('TWO_BIRTH_SOURCE_ROOT', '/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26'))
ORIGINAL = (ROOT / patch_solver.SOLVER).read_text()
PATCHED = patch_solver.patch_text(ORIGINAL)


def functions(source, names):
    nodes = [node for node in ast.parse(source).body
             if isinstance(node, ast.FunctionDef) and node.name in names]
    if {node.name for node in nodes} != set(names):
        raise AssertionError('Requested helper definitions missing')
    return ast.Module(body=nodes, type_ignores=[])


def make_namespace():
    namespace = dict(np=np, SimpleNamespace=SimpleNamespace, DEAD_VALUE_CUTOFF=-1e9,
        independent_child_maturation_active=lambda P: P.child_state_mode == 'independent_count',
        parent_age_maturation_active=lambda P: bool(getattr(P, 'parent_age_mode', False)),
        readiness_gate_active=lambda P: bool(getattr(P, 'readiness_mode', False)),
        readiness_settled_state=lambda P: 0,
        birth_destination_child_state=lambda P, cs: cs + 1)
    utils = (ROOT / 'code/model/intergen_eqscale_seq_optimized/utils.py').read_text()
    exec(compile(functions(utils, ['logsumexp']), '<native logsumexp>', 'exec'), namespace)
    exec(compile(functions(PATCHED, ['two_births_active', 'two_birth_last_opportunity',
         'two_birth_fertility_step', 'snapshot_extra_birth_probabilities']), '<transformed helpers>', 'exec'), namespace)
    return namespace


NS = make_namespace()


def parameters(**overrides):
    values = dict(sequential_births=True, joint_nested_choice=False,
        child_state_mode='independent_count', n_parity=4, n_child_states=4,
        kappa_fert=.175, kappa_fert_continuation=.332, first_birth_fixed_cost=.621)
    values.update(overrides)
    return SimpleNamespace(**values)


def choice_body(source):
    bellman = next(node for node in ast.parse(source).body
                   if isinstance(node, ast.FunctionDef) and node.name == 'solve_bellman_full_markov_income')
    found = [node for node in ast.walk(bellman) if isinstance(node, ast.If)
             and ast.unparse(node.test) == "bool(getattr(P, 'sequential_births', False))"]
    if len(found) != 1:
        raise AssertionError(f'Expected one sequential choice block, found {len(found)}')
    return ast.Module(body=found[0].body, type_ignores=[])


def execute_choice(source, two=False, pi=.83, P=None):
    namespace = dict(NS)
    shape = (2, 1, 1, 1, 1, 4, 4)
    values = np.arange(32, dtype=float).reshape((2, 1, 1, 4, 4)) / 19. - 2.
    namespace.update(P=P or parameters(two_births_per_period=two), Nb=2, nt=1, I=1,
        J=1, Nz=1, npar=4, ncs=4, j=0, zz=0, pi_j=pi, two_births=two,
        VI=values, VI_ex=None, V=np.zeros(shape),
        fert_probs=np.zeros(shape[:5] + (4,)), fert_value=np.zeros(shape[:5]),
        fert2_probs=np.zeros(shape[:5] + (2, 2, 4)),
        fert_extra_probs=np.zeros(shape[:5] + (2, 2, 4)) if two else None)
    exec(compile(choice_body(source), '<actual sequential Bellman block>', 'exec'), namespace)
    return namespace


class SolverPatchTests(unittest.TestCase):
    def test_exact_transform_and_syntax(self):
        compile(PATCHED, '<patched solver>', 'exec')
        self.assertEqual(patch_solver.patch_sources({patch_solver.SOLVER: ORIGINAL})[patch_solver.SOLVER], PATCHED)
        with self.assertRaises(ValueError):
            patch_solver.patch_text(PATCHED)
        with self.assertRaises(ValueError):
            patch_solver.patch_text(ORIGINAL.replace('    P._joint_choice = joint\n', ''))

    def test_flag_off_actual_choice_block_exact(self):
        for pi in (0., .5, 1.):
            old = execute_choice(ORIGINAL, pi=pi)
            new = execute_choice(PATCHED, pi=pi)
            for name in ('V', 'fert_probs', 'fert_value', 'fert2_probs'):
                np.testing.assert_array_equal(old[name], new[name], err_msg=name)

    def test_bellman_two_levels_against_direct_formula(self):
        state = execute_choice(PATCHED, two=True)
        VI, P, pi = state['VI'], state['P'], state['pi_j']
        for n, m in ((0, 0), (1, 0), (1, 1)):
            W = VI[..., n + 1, m + 1]
            T = pi * VI[..., n + 2, m + 2] + (1-pi) * W
            k = P.kappa_fert_continuation
            inner = np.logaddexp(W/k, T/k) * k
            inner_try = np.exp(T/k - np.logaddexp(W/k, T/k))
            np.testing.assert_allclose(state['fert_extra_probs'][:, :, :, 0, 0, 1, n, m], inner_try, atol=2e-15)
            outer_try = pi * inner + (1-pi) * VI[..., n, m]
            if n == 0:
                outer_try -= pi * P.first_birth_fixed_cost
            scale = P.kappa_fert if n == 0 else k
            expected = scale * np.logaddexp(VI[..., n, m]/scale, outer_try/scale)
            np.testing.assert_allclose(state['V'][:, :, :, 0, 0, n, m], expected, atol=2e-15)
        # A household starting with two children gets one remaining birth only.
        old = execute_choice(ORIGINAL)
        np.testing.assert_array_equal(state['V'][:, :, :, 0, 0, 2, :3], old['V'][:, :, :, 0, 0, 2, :3])

    def test_inner_pi_edges_and_deterministic_limit(self):
        last = NS['two_birth_last_opportunity']
        wait = np.array([0., 2.]); success = np.array([4., -1.]); k = .332
        value, prob = last(wait, success, 0., k)
        np.testing.assert_allclose(value, wait + k * np.log(2.), atol=1e-15)
        np.testing.assert_allclose(prob, .5)
        value, prob = last(wait, success, 1., 1e-9)
        np.testing.assert_allclose(value, np.maximum(wait, success), atol=1e-8)
        np.testing.assert_array_equal(np.argmax(prob, axis=-1), [1, 0])
        with self.assertRaises(ValueError):
            last(wait, success, 1., 0.)
        # The additional inclusive value is multiplied by zero when no first
        # success is possible, leaving the outer choice exactly unchanged.
        old = execute_choice(ORIGINAL, pi=0.)
        new = execute_choice(PATCHED, two=True, pi=0.)
        np.testing.assert_array_equal(old['V'], new['V'])

    def test_large_live_values_normalize_without_changing_inclusive_value(self):
        wait = np.array([-8e8, -1e8, -1e6, -10., -2e9])
        success = wait + np.array([0., 5., -7., 2., 0.])
        pi, scale = .83, .332
        menu = np.stack((wait, pi * success + (1-pi) * wait), axis=-1)
        old_value, old_prob = NS['logsumexp'](menu / scale, axis=-1)
        value, prob = NS['two_birth_last_opportunity'](wait, success, pi, scale)
        live = np.max(menu, axis=-1) > -1e9
        self.assertGreater(float(np.max(abs(old_prob[live].sum(axis=-1)-1))), 1e-12)
        self.assertLess(float(np.max(abs(prob[live].sum(axis=-1)-1))), 1e-12)
        np.testing.assert_array_equal(value, scale * old_value)
        np.testing.assert_array_equal(prob[~live], 0.)
        self.assertTrue(np.isfinite(prob).all() and prob.min() >= 0 and prob.max() <= 1)

    def test_gates_and_cache_snapshot(self):
        active = NS['two_births_active']
        self.assertFalse(active(SimpleNamespace()))
        self.assertTrue(active(parameters(two_births_per_period=True)))
        for bad in ({'joint_nested_choice': True}, {'readiness_mode': True},
                    {'parent_age_mode': True}, {'child_state_mode': 'historical'}, {'n_parity': 3}):
            with self.assertRaises(ValueError):
                active(parameters(two_births_per_period=True, **bad))
        P = parameters(two_births_per_period=True)
        with self.assertRaises(RuntimeError):
            NS['snapshot_extra_birth_probabilities'](P)
        P._fert_extra_probs = np.full((2, 1, 1, 1, 1, 2, 2, 4), .25)
        snapshot = NS['snapshot_extra_birth_probabilities'](P)
        P._fert_extra_probs.fill(.75)
        np.testing.assert_array_equal(snapshot, .25)
        self.assertFalse(np.shares_memory(snapshot, P._fert_extra_probs))
        P.two_births_per_period = False
        self.assertIsNone(NS['snapshot_extra_birth_probabilities'](P))

    def test_actual_accepted_price_payload_restores_owned_cache(self):
        # Run the actual retained-payload statements and upgrade function with
        # a no-computation forward stub. A rejected later price mutates P.
        tree = ast.parse(PATCHED)
        fixed = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                     and n.name == 'solve_markov_income_at_prices')
        retain = next(n for n in ast.walk(fixed) if isinstance(n, ast.If)
                      and ast.unparse(n.test) == 'fast_stats and retain_payload')
        ns = dict(NS)
        P = parameters(two_births_per_period=True)
        P._fert2_probs = np.full((1, 2), .125)
        P._fert_extra_probs = np.full((1, 2, 2, 4), .25)
        ns.update({name: np.zeros((1,)) for name in
                   ('V', 'c_pol', 'hR_pol', 'bp_pol', 'tc', 'tp', 'lp_j', 'fp', 'fv', 'r', 'p')})
        ns.update(P=P, sol=SimpleNamespace())
        exec(compile(ast.Module(body=retain.body, type_ignores=[]), '<actual payload>', 'exec'), ns)
        retained = ns['sol']
        P._fert2_probs.fill(.875)
        P._fert_extra_probs.fill(.75)
        np.testing.assert_array_equal(retained._model_payload[-2], .125)
        np.testing.assert_array_equal(retained._model_payload[-1], .25)
        captured = {}
        def forward(*args, **kwargs):
            captured['extra'] = P._fert_extra_probs.copy()
            captured['ordinary'] = P._fert2_probs.copy()
            return None, SimpleNamespace()
        def pack(*args):
            return SimpleNamespace(fert_extra_probs=NS['snapshot_extra_birth_probabilities'](P))
        ns.update(forward_distribution_markov_income=forward,
                  pack_solution_markov_income=pack,
                  time=SimpleNamespace(perf_counter=lambda: 0.))
        P.w_hat = np.ones(1)
        exec(compile(functions(PATCHED, ['upgrade_fast_markov_solution']), '<actual upgrade>', 'exec'), ns)
        upgraded = ns['upgrade_fast_markov_solution'](retained, P, np.ones(1), SimpleNamespace())
        np.testing.assert_array_equal(captured['extra'], .25)
        np.testing.assert_array_equal(captured['ordinary'], .125)
        np.testing.assert_array_equal(upgraded.fert_extra_probs, .25)
        self.assertFalse(np.shares_memory(P._fert_extra_probs, retained._model_payload[-1]))

    def test_flow_extremes_and_no_three_birth_chain(self):
        step = NS['two_birth_fertility_step']
        leading = (2, 1, 1, 3)
        pre = np.zeros(leading + (4, 4))
        first = np.ones(leading)
        continuation = np.ones(leading + (2, 2, 4))
        extra = np.ones_like(continuation)
        for n in range(4):
            pre.fill(0.); pre[..., n, n] = 1.
            post, births, attempts, risk, tagged = step(pre, first, continuation, extra, 1.)
            dest = min(n + 2, 3)
            np.testing.assert_array_equal(post[..., dest, dest], 1.)
            np.testing.assert_array_equal(post.sum(axis=(-2, -1)), 1.)
            np.testing.assert_array_equal(births.sum(axis=-1), dest-n)
            np.testing.assert_array_equal(births, attempts)
            self.assertTrue(np.all(attempts <= risk))
            self.assertEqual(float(tagged.sum()), float(pre.sum()) if n == 0 else 0.)
        pre.fill(0.); pre[..., 0, 0] = 1.
        post, births, *_ = step(pre, first, continuation, extra, .5)
        np.testing.assert_array_equal(post[..., 0, 0], .5)
        np.testing.assert_array_equal(post[..., 1, 1], .25)
        np.testing.assert_array_equal(post[..., 2, 2], .25)
        for pi, attempt in ((0., 1.), (1., 0.)):
            post, births, *_ = step(pre, np.full(leading, attempt), continuation, extra, pi)
            np.testing.assert_array_equal(post, pre)
            self.assertEqual(float(births.sum()), 0.)

    def test_random_flow_stock_identity_and_zero_extra(self):
        rng = np.random.default_rng(612)
        leading = (3, 2, 1, 2)
        pre = rng.uniform(0., 1., leading + (4, 4))
        for n in range(4):
            pre[..., n, n+1:] = 0.
        first = rng.uniform(size=leading)
        continuation = rng.uniform(size=leading+(2, 2, 4))
        extra = rng.uniform(size=continuation.shape)
        step = NS['two_birth_fertility_step']
        post, births, attempts, risk, tagged = step(pre, first, continuation, extra, .79)
        self.assertGreaterEqual(float(post.min()), 0.)
        np.testing.assert_allclose(post.sum(axis=(-2,-1)), pre.sum(axis=(-2,-1)), atol=1e-14)
        before = pre.sum(axis=-1); after = post.sum(axis=-1)
        np.testing.assert_allclose(np.cumsum(before-after, axis=-1)[..., :3], births, atol=2e-15)
        np.testing.assert_allclose(births, .79*attempts, atol=2e-15)
        self.assertTrue(np.all(attempts <= risk + 1e-14))
        np.testing.assert_allclose(tagged.sum(axis=(-2,-1)), births[..., 0], atol=2e-15)
        extra.fill(0.)
        post, births, *_ = step(pre, first, continuation, extra, .79)
        expected = pre.copy()
        for n in range(3):
            for m in range(n+1):
                p = first if n == 0 else continuation[..., 1, n-1, m]
                success = pre[..., n, m] * p * .79
                expected[..., n, m] -= success
                expected[..., n+1, m+1] += success
        np.testing.assert_array_equal(post, expected)


if __name__ == '__main__':
    unittest.main()
