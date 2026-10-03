"""Focused checks for saved-solution aggregates and raw policy plots."""
from __future__ import annotations

import unittest
import sys
from types import SimpleNamespace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from model_policy_tools import aggregate_solution, plot_policy


def synthetic_solution() -> SimpleNamespace:
    shape = (2, 3, 1, 1, 1, 1, 1)
    g = np.zeros(shape)
    g[0, 0, 0, 0, 0, 0, 0] = 0.2  # renter
    g[1, 1, 0, 0, 0, 0, 0] = 0.3  # owner 1, including a stayer
    g[0, 2, 0, 0, 0, 0, 0] = 0.5  # owner 2, including a stayer
    stay = np.zeros(shape)
    stay[1, 1, 0, 0, 0, 0, 0] = 0.1
    stay[0, 2, 0, 0, 0, 0, 0] = 0.2
    beginning = np.zeros(shape)
    beginning[1, 0, 0, 0, 0, 0, 0] = 0.2
    beginning[0, 1, 0, 0, 0, 0, 0] = 0.3
    beginning[1, 2, 0, 0, 0, 0, 0] = 0.5
    c, bp, cs, bps, h = [np.zeros(shape) for _ in range(5)]
    c[0, 0, 0, 0, 0, 0, 0] = 2.0
    c[1, 1, 0, 0, 0, 0, 0] = 3.0
    cs[1, 1, 0, 0, 0, 0, 0] = 4.0
    c[0, 2, 0, 0, 0, 0, 0] = 5.0
    cs[0, 2, 0, 0, 0, 0, 0] = 6.0
    bp[0, 0, 0, 0, 0, 0, 0] = 0.1
    bp[1, 1, 0, 0, 0, 0, 0] = 0.4
    bps[1, 1, 0, 0, 0, 0, 0] = 0.5
    bp[0, 2, 0, 0, 0, 0, 0] = 0.6
    bps[0, 2, 0, 0, 0, 0, 0] = 0.7
    h[0, 0, 0, 0, 0, 0, 0] = 1.5
    return SimpleNamespace(
        b_grid=np.array([0.0, 2.0]), g=g, g_stay_distribution=stay,
        g_beginning_distribution=beginning, c_pol=c, bp_pol=bp,
        c_pol_stay=cs, bp_pol_stay=bps, hR_pol=h,
        V=np.ones(shape), houses=np.array([2.0, 4.0]),
    )


class ModelPolicyToolsTests(unittest.TestCase):
    def test_aggregate_solution_hand_calculation(self):
        out = aggregate_solution(synthetic_solution(), age_start=18, period_years=4)
        self.assertAlmostEqual(out["population_mass"], 1.0)
        overall = out["overall"]
        self.assertAlmostEqual(overall["mean_consumption"], 4.1)
        self.assertAlmostEqual(overall["mean_inherited_assets"], 1.4)
        self.assertAlmostEqual(overall["mean_next_assets"], 0.47)
        self.assertAlmostEqual(overall["mean_asset_change"], -0.93)
        self.assertAlmostEqual(overall["mean_rooms"], 2.9)
        self.assertAlmostEqual(overall["ownership_rate"], 0.8)
        self.assertEqual(out["by_age"][0]["age"], 18)
        np.testing.assert_allclose(out["inherited_asset_distribution"]["mass"], [0.3, 0.7])

    def test_distribution_gates_and_occupied_policy_gate(self):
        sol = synthetic_solution()
        sol.g_stay_distribution[0, 0, 0, 0, 0, 0, 0] = 0.01
        with self.assertRaisesRegex(ValueError, "renter stayer"):
            aggregate_solution(sol)
        sol = synthetic_solution()
        sol.c_pol[0, 0, 0, 0, 0, 0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            aggregate_solution(sol)
        sol = synthetic_solution()
        sol.V[0, 0, 0, 0, 0, 0, 0] = -1.0e12
        with self.assertRaisesRegex(ValueError, "infeasible"):
            aggregate_solution(sol)

    def test_raw_policy_plot_returns_axes_without_showing(self):
        import matplotlib
        matplotlib.use("Agg")
        ax = plot_policy(synthetic_solution(), variable="next_assets", age=18,
                         income=0, branch=1, owner_policy="staying")
        np.testing.assert_allclose(ax.lines[0].get_ydata(), [0.0, 0.5])
        np.testing.assert_allclose(ax.lines[0].get_xdata(), [0.0, 2.0])
        self.assertEqual(ax.get_lines()[0].get_marker(), ".")
        self.assertEqual(ax.get_xlabel(), "Assets b (model units)")
        self.assertIn("mean annual gross-earnings units", ax.get_title())
        ax_zoom = plot_policy(synthetic_solution(), variable="next_assets", age=18,
                              income=0, branch=1, owner_policy="staying", xlim=(0.0, 0.1))
        self.assertLess(max(ax_zoom.get_ylim()), 0.5)
        sol = synthetic_solution()
        del sol.houses
        ax_node = plot_policy(sol, variable="housing", age=18,
                              income=0, branch=2, axis="node_index",
                              houses=[2.5, 4.5])
        np.testing.assert_allclose(ax_node.lines[0].get_xdata(), [1, 2])
        np.testing.assert_allclose(ax_node.lines[0].get_ydata(), [4.5, 4.5])
        self.assertEqual(ax_node.get_lines()[0].get_marker(), ".")

    def test_aggregate_plot_defaults_to_central_wealth_nodes(self):
        from model_policy_tools import plot_aggregates
        shape = (5, 3, 1, 1, 1, 1, 1)
        g = np.zeros(shape)
        g[0, 1, 0, 0, 0, 0, 0] = 0.00001
        g[2, 0, 0, 0, 0, 0, 0] = 0.99998
        g[4, 2, 0, 0, 0, 0, 0] = 0.00001
        zero = np.zeros(shape)
        sol = SimpleNamespace(
            b_grid=np.arange(5.0), g=g, g_stay_distribution=zero.copy(),
            g_beginning_distribution=g.copy(), c_pol=zero.copy(), bp_pol=zero.copy(),
            c_pol_stay=zero.copy(), bp_pol_stay=zero.copy(), hR_pol=zero.copy(),
            V=np.ones(shape), houses=[2.0, 4.0],
        )
        _, axes = plot_aggregates(sol)
        self.assertEqual(tuple(axes[1, 1].get_xlim()), (1.0, 3.0))
        np.testing.assert_allclose(axes[1, 1].lines[0].get_xdata(), [1.0, 2.0, 3.0])


if __name__ == "__main__":
    unittest.main()
