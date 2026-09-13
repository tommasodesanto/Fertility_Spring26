from __future__ import annotations

from types import SimpleNamespace
import unittest

import numpy as np

import e5f_exact_policy_cache as c


def policy_function(*, price, rent, P, b_grid, shared, continuation_V):
    policy_function.calls += 1
    if getattr(P, "fail", False): raise RuntimeError("native failure")
    return SimpleNamespace(value=np.array([price, rent, P.psi, shared.scale,
        b_grid.sum(), continuation_V.sum()], dtype=float))


policy_function.calls = 0


def arguments():
    return dict(price=1.0, rent=.1, P=SimpleNamespace(psi=.2, nested=SimpleNamespace(x=1)),
        b_grid=np.arange(6.0), shared=SimpleNamespace(scale=2.0, array=np.arange(3.0)),
        continuation_V=np.arange(8.0).reshape(2, 4))


class ExactPolicyCacheTests(unittest.TestCase):
    def setUp(self):
        policy_function.calls = 0
        self.module = SimpleNamespace(solve_date_policy=policy_function)

    def test_key_includes_every_required_argument(self):
        base = arguments(); key = c.exact_call_key(policy_function, (), base)
        mutations = []
        for name, value in (("price", 1.1), ("rent", .11)):
            item = arguments(); item[name] = value; mutations.append(item)
        item = arguments(); item["P"].nested.x = 2; mutations.append(item)
        item = arguments(); item["b_grid"][0] = 9; mutations.append(item)
        item = arguments(); item["shared"].array[0] = 9; mutations.append(item)
        item = arguments(); item["continuation_V"][0, 0] = 9; mutations.append(item)
        for changed in mutations:
            self.assertNotEqual(key, c.exact_call_key(policy_function, (), changed))

    def test_noncontiguous_and_contiguous_arrays_hash_by_exact_values(self):
        left = arguments(); left["continuation_V"] = np.arange(12.0).reshape(3, 4)[:, ::2]
        right = arguments(); right["continuation_V"] = np.ascontiguousarray(left["continuation_V"])
        self.assertFalse(left["continuation_V"].flags.c_contiguous)
        self.assertEqual(c.exact_call_key(policy_function, (), left),
                         c.exact_call_key(policy_function, (), right))

    def test_scalar_types_are_distinct_and_hidden_object_state_bypasses(self):
        left=arguments();right=arguments();right['price']=np.float64(left['price'])
        self.assertNotEqual(c.exact_call_key(policy_function,(),left),
                            c.exact_call_key(policy_function,(),right))
        args=arguments();args['P'].callback=lambda x:x
        with c.policy_cache(self.module,max_bytes=10000) as stats:
            self.module.solve_date_policy(**args)
            self.module.solve_date_policy(**args)
            self.assertEqual(stats.snapshot()['hits'],0)
            self.assertEqual(stats.snapshot()['serialization_bypasses'],2)

    def test_hit_returns_fresh_unaliased_pickle(self):
        original = self.module.solve_date_policy
        with c.policy_cache(self.module, max_bytes=10000) as stats:
            first = self.module.solve_date_policy(**arguments())
            first.value[:] = -99
            second = self.module.solve_date_policy(**arguments())
            second.value[:] = -88
            third = self.module.solve_date_policy(**arguments())
            self.assertNotEqual(float(third.value[0]), -99)
            self.assertNotEqual(float(third.value[0]), -88)
            self.assertEqual(stats.snapshot()["hits"], 2)
            self.assertEqual(stats.snapshot()["actual_solves"], 1)
        self.assertIs(self.module.solve_date_policy, original)

    def test_lru_eviction_respects_byte_bound(self):
        with c.policy_cache(self.module, max_bytes=300) as stats:
            for price in (1.0, 1.1, 1.2, 1.0):
                args = arguments(); args["price"] = price
                self.module.solve_date_policy(**args)
            snapshot = stats.snapshot()
            self.assertLessEqual(snapshot["bytes"], 300)
            self.assertGreater(snapshot["evictions"], 0)

    def test_native_exceptions_are_not_cached_and_original_restored(self):
        original = self.module.solve_date_policy
        args = arguments(); args["P"].fail = True
        with self.assertRaises(RuntimeError):
            with c.policy_cache(self.module, max_bytes=10000) as stats:
                for _ in range(2):
                    with self.assertRaisesRegex(RuntimeError, "native failure"):
                        self.module.solve_date_policy(**args)
                self.assertEqual(stats.snapshot()["actual_solves"], 2)
                self.assertEqual(stats.snapshot()["exceptions"], 2)
                raise RuntimeError("leave context")
        self.assertIs(self.module.solve_date_policy, original)

    def test_unserializable_key_bypasses_without_changing_result(self):
        args = arguments(); args["P"].cycle = args["P"]
        with c.policy_cache(self.module, max_bytes=10000) as stats:
            result = self.module.solve_date_policy(**args)
            self.assertEqual(result.value[0], 1.0)
            self.assertEqual(stats.snapshot()["serialization_bypasses"], 1)
            self.assertEqual(stats.snapshot()["actual_solves"], 1)


if __name__ == "__main__": unittest.main()
