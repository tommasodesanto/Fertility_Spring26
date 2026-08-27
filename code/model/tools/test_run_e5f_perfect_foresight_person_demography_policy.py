from __future__ import annotations

import unittest
from types import SimpleNamespace

import numpy as np

import run_e5f_perfect_foresight_person_demography_policy as policy


def state(*, persons: np.ndarray, heads: np.ndarray, g_pre: np.ndarray, year: int):
    return SimpleNamespace(
        persons=SimpleNamespace(persons=persons, heads=heads, year=year),
        g_pre=g_pre,
    )


class TerminalConvergenceDiagnosticsTest(unittest.TestCase):
    def test_normalized_l1_is_invariant_to_total_mass(self) -> None:
        reference = np.array([1.0, 2.0, 3.0])
        self.assertAlmostEqual(policy.normalized_l1(7.0 * reference, reference), 0.0)

    def test_all_declared_terminal_checks_are_jointly_enforced(self) -> None:
        reference = state(
            persons=np.array([[2.0, 3.0], [4.0, 1.0]]),
            heads=np.array([[0.5, 1.5], [2.0, 0.5]]),
            g_pre=np.array([[[1.0, 2.0], [3.0, 4.0]]]),
            year=2100,
        )
        terminal = SimpleNamespace(
            state=reference,
            asset_price=2.0,
            renter_price=0.5,
            equal_transfer=0.25,
            psi_child=-0.07,
        )
        passing = SimpleNamespace(
            terminal_state=state(
                persons=reference.persons.persons * 1.005,
                heads=reference.persons.heads * 1.005,
                g_pre=reference.g_pre * 1.005,
                year=2535,
            ),
            prices=np.array([2.0 * 1.005]),
            rents=np.array([0.5 * 1.005]),
            rows=[{"equal_transfer_period_units": 0.25 * 1.005}],
        )
        result = policy.terminal_convergence_diagnostics(
            passing, terminal=terminal, psi_path=[-0.0695]
        )
        self.assertTrue(result["all_checks_pass"])
        self.assertTrue(all(result["checks"].values()))

        failing = SimpleNamespace(
            terminal_state=passing.terminal_state,
            prices=np.array([2.0 * 1.011]),
            rents=passing.rents,
            rows=passing.rows,
        )
        result = policy.terminal_convergence_diagnostics(
            failing, terminal=terminal, psi_path=[-0.0695]
        )
        self.assertFalse(result["all_checks_pass"])
        self.assertFalse(result["checks"]["asset_price_relative_gap"])
        self.assertTrue(
            all(
                passed
                for name, passed in result["checks"].items()
                if name != "asset_price_relative_gap"
            )
        )


if __name__ == "__main__":
    unittest.main()
