"""Exact-loop mock tests: no model evaluations."""
from types import SimpleNamespace as NS
import math
import time
import unittest
from unittest import mock

import e5f_matched_pf_price_root as root


def endpoint(residual, *, valid=True, accepted=None, regime="fixed_transfer", price=1.):
    return NS(contract={"fiscal_regime": regime, "asset_price": price}, root_coordinates=("log_asset_price",),
              residuals={"housing_relative": residual}, root_residuals=[residual],
              mapping_valid=valid, accepted=abs(residual) <= 2e-4 if accepted is None else accepted,
              gates={"inner": valid, "market": abs(residual) <= 2e-4})


class PriceRootTest(unittest.TestCase):
    def args(self, evaluate, **kw):
        return dict(fresh_evaluate=evaluate, start_price=1., bound_ratios=(.5, 2.),
                    maximum_evaluations=20, deadline_monotonic=time.monotonic() + 60,
                    market_tolerance=2e-4, replay_tolerance=2e-10, **kw)

    def test_direction_bracket_root_and_uncached_final_replay(self):
        for desired in (.8, 1.2):
            calls, records = [], []
            def evaluate(price):
                calls.append(price)
                return endpoint(math.log(desired / price), price=price)
            result = root.solve_price_root(**self.args(evaluate, progress_callback=lambda latest, best: records.append((latest, best))))
            self.assertTrue(result.converged)
            self.assertEqual(result.evaluations, len(calls))
            self.assertEqual(calls[-1], calls[-2])
            self.assertTrue(calls[1] < calls[0] if desired < 1 else calls[1] > calls[0])
            self.assertAlmostEqual(calls[-1], desired)
            self.assertEqual(records[-1][0]["phase"], "final_replay")
            self.assertTrue(all(best is not None for _, best in records))

    def test_starting_root_still_requires_fresh_replay(self):
        evaluate = mock.Mock(side_effect=[endpoint(0.), endpoint(0.)])
        result = root.solve_price_root(**self.args(evaluate))
        self.assertTrue(result.converged)
        self.assertEqual(evaluate.call_count, 2)
        self.assertEqual(result.evaluations, 2)

    def test_final_replay_can_fail_despite_cached_search_success(self):
        evaluate = mock.Mock(side_effect=[endpoint(0.), endpoint(1e-3)])
        result = root.solve_price_root(**self.args(evaluate))
        self.assertEqual(result.status, "incomplete_replay_failure")
        self.assertFalse(result.converged)

    def test_replay_requires_reproduction_even_inside_market_gate(self):
        evaluate = mock.Mock(side_effect=[endpoint(0.), endpoint(1e-5)])
        result = root.solve_price_root(**self.args(evaluate))
        self.assertEqual(result.status, "incomplete_replay_failure")
        self.assertEqual(result.replay_residual_gap, 1e-5)

    def test_discontinuous_nonzero_residual_never_passes_on_bracket_width(self):
        evaluate = lambda price: endpoint(-.01 if price > .812345 else .01, price=price)
        args = self.args(evaluate)
        args["maximum_evaluations"] = 90
        result = root.solve_price_root(**args)
        self.assertFalse(result.converged)
        self.assertIn(result.status, {"incomplete_price_resolution", "incomplete_budget"})
        self.assertEqual(abs(result.best_record["housing_residual"]), .01)

    def test_budget_reserves_replay_slot_and_retains_best(self):
        calls = mock.Mock(side_effect=lambda price: endpoint(math.log(.8 / price), price=price))
        args = self.args(calls)
        args["maximum_evaluations"] = 2
        result = root.solve_price_root(**args)
        self.assertEqual(result.status, "incomplete_budget")
        self.assertEqual(calls.call_count, 1)
        self.assertIsNotNone(result.best_endpoint)

    def test_invalid_mapping_stops_without_trying_other_prices(self):
        evaluate = mock.Mock(return_value=endpoint(.1, valid=False))
        result = root.solve_price_root(**self.args(evaluate))
        self.assertEqual(result.status, "incomplete_invalid_mapping")
        self.assertEqual(evaluate.call_count, 1)

    def test_economic_nonexistence_records_failure_without_retry(self):
        evaluate = mock.Mock(side_effect=RuntimeError("no finite demographic root"))
        result = root.solve_price_root(**self.args(evaluate))
        self.assertEqual(result.status, "incomplete_evaluation_failure")
        self.assertIn("no finite demographic root", result.records[-1]["error"])
        self.assertEqual(result.evaluations, 1)

    def test_no_bracket_is_not_success(self):
        result = root.solve_price_root(**self.args(lambda price: endpoint(-.1, price=price)))
        self.assertEqual(result.status, "incomplete_no_bracket")
        self.assertEqual(result.evaluations, 3)

    def test_scalar_cache_uses_executable_price_not_log_roundoff(self):
        evaluate = mock.Mock(side_effect=lambda price: endpoint(-.1, price=price))
        args = self.args(evaluate)
        args["bound_ratios"] = (math.nextafter(1., 0.), 2.)
        result = root.solve_price_root(**args)
        self.assertEqual(result.status, "incomplete_no_bracket")
        self.assertEqual(result.cache_hits, 1)
        self.assertEqual(evaluate.call_count, 2)

    def test_wrong_evaluation_price_is_not_a_valid_root_trial(self):
        result = root.solve_price_root(**self.args(lambda price: endpoint(0., price=2.)))
        self.assertEqual(result.status, "incomplete_evaluation_failure")

    def test_market_gate_cannot_be_relaxed(self):
        args = self.args(lambda price: endpoint(0.))
        args["market_tolerance"] = 1e-3
        with self.assertRaises(ValueError):
            root.solve_price_root(**args)

    def test_expired_deadline_never_starts_evaluation(self):
        evaluate = mock.Mock()
        args = self.args(evaluate)
        args["deadline_monotonic"] = 0.
        result = root.solve_price_root(**args)
        self.assertEqual(result.status, "incomplete_deadline")
        evaluate.assert_not_called()

    def test_fiscal_regime_cannot_change(self):
        result = root.solve_price_root(**self.args(lambda price: endpoint(0., regime="equal_rebate")))
        self.assertEqual(result.status, "incomplete_evaluation_failure")

    def test_endpoint_acceptance_is_not_replaced_by_market_residual(self):
        result = root.solve_price_root(**self.args(lambda price: endpoint(0., accepted=False)))
        self.assertEqual(result.status, "incomplete_endpoint_acceptance")
        self.assertEqual(result.evaluations, 1)


if __name__ == "__main__":
    unittest.main()
