import unittest
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np

import run_e5f_successive_surprises_overnight as runner


class FiniteCarryTests(unittest.TestCase):
    def _result(self, *, next_state=None, converged=True):
        return SimpleNamespace(
            path=SimpleNamespace(rows=[{"renter_price": 1.25, "period": 7}]),
            root_receipt={
                "finite_horizon_market_fiscal_converged": converged,
                "terminal_distance_passed": False,
                "production_eligible": True,
                "horizon_verified": True,
                "final": {"prices": [2.0, 2.5], "fiscal_values": [3.0, 3.5]},
            },
            next_state=next_state,
            realized_row=None,
        )

    def test_disabled_keeps_next_state_none(self):
        result = self._result()
        module = SimpleNamespace(first_period_state=Mock())
        self.assertIs(
            runner.carry_finite_diagnostic(
                result, enabled=False, inherited=SimpleNamespace(year=2007),
                old="old", demographics="demo", psi=.4, module=module,
            ), result
        )
        module.first_period_state.assert_not_called()

    def test_failed_finite_horizon_cannot_carry(self):
        result = self._result(converged=False)
        module = SimpleNamespace(first_period_state=Mock())
        self.assertIs(
            runner.carry_finite_diagnostic(
                result, enabled=True, inherited="inherited", old="old",
                demographics="demo", psi=.4, module=module,
            ), result
        )
        module.first_period_state.assert_not_called()

    def test_finite_success_replays_exact_first_period_and_flags_diagnostic(self):
        result = self._result()
        replayed = object()
        module = SimpleNamespace(first_period_state=Mock(return_value=replayed), SurpriseResult=None)
        # Supply the constructor used by the runner without importing the full
        # surprise module; this keeps the test independent of model runtime.
        module.SurpriseResult = lambda path, receipt, next_state, realized: SimpleNamespace(
            path=path, root_receipt=receipt, next_state=next_state, realized_row=realized
        )
        inherited = SimpleNamespace(year=2007)
        carried = runner.carry_finite_diagnostic(
            result, enabled=True, inherited=inherited, old="old",
            demographics="demo", psi=.4, module=module,
        )
        module.first_period_state.assert_called_once_with(
            inherited=inherited, old_state="old", demographics="demo",
            path=result.path, prices=[2.0, 2.5], pensions=[3.0, 3.5], psi=.4,
        )
        self.assertIs(carried.next_state, replayed)
        self.assertFalse(carried.root_receipt["terminal_distance_passed"])
        self.assertTrue(carried.root_receipt["diagnostic_finite_horizon_state_carry"])
        self.assertFalse(carried.root_receipt["production_eligible"])
        self.assertFalse(carried.root_receipt["horizon_verified"])
        self.assertEqual(carried.realized_row["forecast_vintage_year"], 2007)

    def test_replay_exception_propagates(self):
        result = self._result()
        error = RuntimeError("replay failed")
        module = SimpleNamespace(first_period_state=Mock(side_effect=error))
        with self.assertRaises(RuntimeError) as raised:
            runner.carry_finite_diagnostic(
                result, enabled=True, inherited="inherited", old="old",
                demographics="demo", psi=.4, module=module,
            )
        self.assertIs(raised.exception, error)

    def test_existing_next_state_is_preserved(self):
        existing = object()
        result = self._result(next_state=existing)
        module = SimpleNamespace(first_period_state=Mock())
        carried = runner.carry_finite_diagnostic(
            result, enabled=True, inherited="inherited", old="old",
            demographics="demo", psi=.4, module=module,
        )
        self.assertIs(carried, result)
        self.assertIs(carried.next_state, existing)
        module.first_period_state.assert_not_called()


class LargeOwnerObservationTests(unittest.TestCase):
    def test_extracts_tenure_and_child_axes_and_excludes_renter_small_owner(self):
        nb, tenure, locations, ages, zones, wealth, children = (1, 5, 1, 2, 1, 1, 3)
        g = np.zeros((nb, tenure, locations, ages, zones, wealth, children))
        # For each age, put mass into renter, small-owner, and two qualifying
        # owner rungs. Only the latter two should enter the sufficient stats.
        g[0, 0, 0, 0, 0, 0, :] = [100, 100, 100]
        g[0, 1, 0, 0, 0, 0, :] = [200, 200, 200]
        g[0, 3, 0, 0, 0, 0, :] = [1, 2, 3]
        g[0, 4, 0, 0, 0, 0, :] = [4, 5, 6]
        g[0, 3, 0, 1, 0, 0, :] = [10, 20, 30]
        g[0, 4, 0, 1, 0, 0, :] = [40, 50, 60]
        e = SimpleNamespace(g_current=g)
        p = SimpleNamespace(
            child_state_mode="independent_count", H_own=np.array([2.0, 5.0, 6.0, 8.0]),
            J=ages, age_start=25, da=4,
        )
        self.assertEqual(
            runner.large_owner_observation(e, p),
            [
                {"age": 25.0, "age_width": 4.0, "without_children": 5.0, "with_children": 16.0},
                {"age": 29.0, "age_width": 4.0, "without_children": 50.0, "with_children": 160.0},
            ],
        )


if __name__ == "__main__":
    unittest.main()
