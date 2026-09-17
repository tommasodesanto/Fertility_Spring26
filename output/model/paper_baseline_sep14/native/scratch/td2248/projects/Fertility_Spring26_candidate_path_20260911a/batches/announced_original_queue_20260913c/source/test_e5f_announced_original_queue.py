"""Pure routing tests for the announced original-queue adapter; no model solve."""
from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace as NS
import unittest
from unittest.mock import patch

import numpy as np


MODULE = Path(__file__).with_name("run_e5f_announced_original_queue.py")
spec = importlib.util.spec_from_file_location("announced_driver", MODULE)
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)


class FakeCalendar:
    def evaluate_period(self, *args, **kwargs):
        return NS()


class FakePF:
    def __init__(self, *, bad_gate=False, missing_psi=False):
        self.calendar = FakeCalendar()
        self.calls = []
        self.bad_gate, self.missing_psi = bad_gate, missing_psi

    def evaluate_path_at_prices(self, **kw):
        self.calls.append(kw)
        psi = np.asarray(kw["psi_path"], dtype=float)
        # The adapter's observer wraps this exact native calendar entry point.
        for value in psi:
            self.calendar.evaluate_period(None, None, NS(psi_child=float(value)), None, None)
        rows = [({} if self.missing_psi else {"psi_child": float(x)}) for x in psi]
        terminal = NS(g_pre=np.array([1.0]), scheduled_entries=[1.0] * 4,
                      scheduled_raw_entries=[1.0] * 4)
        return NS(rows=rows, values=[np.asarray(kw["terminal_V"])] * (len(psi) + 1),
            terminal_state=terminal, bellman_solves=0, maximum_market_residual=0.0,
            maximum_mass_accounting_error=3e-8 if self.bad_gate else 0.0,
            maximum_policy_reproduction_error=0.0, maximum_feasibility_projection_mass=0.0,
            elapsed_seconds=0.0)


def fake_context(pf):
    def paths(prices, pensions, transfers):
        return tuple(np.asarray(x, dtype=float) for x in (prices, pensions, transfers))
    queue = NS(_validated_paths=paths, BIRTH_TO_ENTRY_CONVERSION=1 / 2.1,
               PAYROLL_TAX=.179, annotate_original_queue_metadata=lambda: {})
    rebated = NS(_terminal_parts=lambda t: (t.parameters, t.policy, t.asset_price, None),
                 InheritedState=lambda year, households: NS(year=year, households=households))
    return NS(queue=queue, rebated=rebated, joined=NS(pf=pf))


def inputs(n, final):
    old = NS(parameters=NS(), b_grid=np.array([0.]), supply_rule=NS())
    inherited = NS(year=2007, households=NS())
    terminal = NS(parameters=NS(psi_child=float(final)), policy=NS(V=np.array([7.])), asset_price=2.)
    return dict(inherited=inherited, old_state=old, prices=np.arange(1, n + 1.),
                pensions=np.ones(n), transfers=np.zeros(n), psi=float(final), terminal=terminal)


class AnnouncedQueueRoutingTests(unittest.TestCase):
    def test_full_vector_reaches_native_without_flattening(self):
        psi = np.array([.128, .117, .106, .092, .092, .092])
        pf = FakePF(); evaluate, _ = driver.announced_queue_path(fake_context(pf), psi)
        out = evaluate(**inputs(6, psi[-1]))
        np.testing.assert_array_equal(pf.calls[0]["psi_path"], psi)
        self.assertEqual([r["announced_psi"] for r in out.rows], psi.tolist())
        self.assertEqual(len(out.values), 7)

    def test_first_replay_uses_first_psi_and_next_price_value(self):
        psi = np.array([.128, .117, .106, .092, .092, .092])
        pf = FakePF(); _, first = driver.announced_queue_path(fake_context(pf), psi)
        path = NS(values=[np.array([9.]), np.array([9.])],
                  rows=[dict(psi_child=float(psi[0]),pension_period_units=1.)])
        state = first(path=path, prices=[1., 1.5], pensions=[1., 1.], transfers=[0., 0.],
                      inherited=NS(year=2007, households=NS()), old_state=NS(parameters=NS(), b_grid=np.array([0.]), supply_rule=NS()),
                      demographics=None, psi=psi[-1])
        self.assertEqual(state.year, 2011)
        self.assertEqual(len(pf.calls), 1)
        np.testing.assert_array_equal(pf.calls[0]["psi_path"], [psi[0]])
        self.assertEqual(float(pf.calls[0]["terminal_price"]), 1.5)
        np.testing.assert_array_equal(pf.calls[0]["terminal_V"], [9.])

    def test_bad_native_gate_rejects(self):
        psi = np.full(6, .092)
        evaluate, _ = driver.announced_queue_path(fake_context(FakePF(bad_gate=True)), psi)
        with self.assertRaisesRegex(RuntimeError, "numerical gates"):
            evaluate(**inputs(6, psi[-1]))

    def test_missing_native_psi_rejects(self):
        psi = np.full(6, .092)
        evaluate, _ = driver.announced_queue_path(fake_context(FakePF(missing_psi=True)), psi)
        with self.assertRaisesRegex(RuntimeError, "Native row psi"):
            evaluate(**inputs(6, psi[-1]))

    def test_nonfinite_native_gate_rejects(self):
        pf=FakePF(); original=pf.evaluate_path_at_prices
        def invalid(**kwargs):
            result=original(**kwargs); result.maximum_mass_accounting_error=float("nan")
            return result
        pf.evaluate_path_at_prices=invalid
        evaluate,_=driver.announced_queue_path(fake_context(pf),np.full(6,.092))
        with self.assertRaisesRegex(RuntimeError,"numerical gates"):
            evaluate(**inputs(6,.092))

    def test_nonfinite_or_mismatched_reproduction_objects_reject(self):
        for a,b in (([float("nan")],[0.]),([1.,2.],[1.])):
            with self.assertRaises(ValueError): driver.exact_gap(a,b)

    def test_initial_guess_extends_saved_100_dates_toward_endpoint(self):
        seed = np.arange(300., dtype=float).reshape(3, 100)
        endpoint = NS(coordinates=np.array([400., 500., 600.]))
        c = NS()
        with patch.object(driver, "read", return_value={"prices": seed.reshape(-1).tolist()}):
            guess = driver.initial_guess(c, endpoint, {"initial_seed_json": "unused"}, 104)
        np.testing.assert_array_equal(guess[:, :100], seed)
        for i, weight in enumerate((.2, .4, .6, .8), 100):
            np.testing.assert_allclose(guess[:, i], (1 - weight) * seed[:, -1] + weight * endpoint.coordinates)


if __name__ == "__main__":
    unittest.main()
