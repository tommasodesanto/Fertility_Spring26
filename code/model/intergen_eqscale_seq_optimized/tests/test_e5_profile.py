from __future__ import annotations

import os
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from intergen_eqscale_seq_optimized import e5_profile
from intergen_eqscale_seq_optimized.solver import add_aggregate_wealth_bequest_flow_moments


class E5ProfileTests(unittest.TestCase):
    def test_pending_timing_targets_refuse_to_build(self) -> None:
        with self.assertRaisesRegex(ValueError, "pending"):
            e5_profile.e5_target_system()

    def test_dummy_timing_targets_complete_twelve_row_contract(self) -> None:
        with patch.dict(e5_profile.E5_TARGETS, {
            "mean_age_first_birth": 26.5,
            "share_first_births_age30plus": 0.35,
        }):
            system = e5_profile.e5_target_system()
        self.assertEqual(system.count, 12)
        self.assertEqual(set(system.targets_dict()), set(system.weights_dict()))
        self.assertEqual(len(e5_profile.E5_DOMAIN), 10)

    def test_wealth_flow_uses_gross_labor_earnings_without_payroll_wedge(self) -> None:
        p = SimpleNamespace(
            J=2, I=1, J_R=1, period_years=4.0, da=4.0,
            use_age_survival=True, survival_probs=np.array([0.5]),
            income=np.array([[4.0, 4.0]]), tau_pay=0.2,
            H_own=np.array([2.0]), z_grid=np.array([1.0]),
        )
        wealth_g = np.zeros((2, 2, 1, 2, 1, 1, 1))
        wealth_g[0, 0, 0, 0, 0, 0, 0] = 1.0
        wealth_g[1, 1, 0, 1, 0, 0, 0] = 1.0
        death_g = wealth_g.copy()
        bp = np.broadcast_to(np.array([-1.0, 3.0]).reshape(2, 1, 1, 1, 1, 1, 1), wealth_g.shape).copy()
        stats = SimpleNamespace()
        add_aggregate_wealth_bequest_flow_moments(stats, wealth_g, death_g, bp, p, np.array([-1.0, 3.0]), np.array([1.0]))
        # Wealth is (-1) + (3 + 2) = 4.  P.income is stored AFTER the payroll
        # wedge, so the gross denominator is 4 / 4 / (1 - 0.2) = 1.25 annually
        # and the ratio is 4 / 1.25 = 3.2.
        self.assertAlmostEqual(stats.aggregate_wealth_to_annual_gross_labor_earnings, 3.2)
        self.assertAlmostEqual(stats.annual_bequest_flow_to_aggregate_wealth, 0.3125)

    def test_e5_smoke_metadata_with_injected_timing_targets(self) -> None:
        # Import after setting E5 so the runner selects its gated domain.
        with patch.dict(os.environ, {"E5": "1"}, clear=False):
            import importlib
            import intergen_eqscale_seq_optimized.run_e1_chain as chain
            chain = importlib.reload(chain)
            with patch.dict(e5_profile.E5_TARGETS, {
                "mean_age_first_birth": 26.5,
                "share_first_births_age30plus": 0.35,
            }), tempfile.TemporaryDirectory() as tmp:
                # This is a real --smoke runner invocation, with the model
                # evaluator replaced only to keep the metadata gate unit-fast.
                def fake_runtime() -> None:
                    import numpy as local_np
                    chain.np = local_np
                    chain.base_overrides = lambda **_kw: {}
                    chain.diagnostic_loss = lambda *_a, **_kw: 0.0
                    chain.external_entry_wealth_overrides_1824 = lambda: {}
                    chain.extract_moments = lambda *_a: {}
                    chain.get_target_set = lambda *_a: ({}, {})
                    chain.income_process_overrides = lambda *_a: {}
                    from intergen_eqscale_seq_optimized.local_panel import jsonable
                    chain.jsonable = jsonable
                    chain.production_profile_overrides = lambda: {}
                    chain.flhsv_income_overrides = lambda: {}
                    chain.InfeasibleThetaError = type("FakeInfeasibleThetaError", (Exception,), {})
                    chain.run_model_cp_dt = lambda *_a, **_kw: (_ for _ in ()).throw(ValueError("stop after metadata"))
                with patch.object(chain, "load_runtime", fake_runtime), patch.object(chain, "build_seed_theta", return_value={
                    "beta": 0.96 ** 4, "delta_alpha": 0.05, "delta_alpha_jump": 0.10,
                    "psi_child": 0.0, "kappa_fert": 1.0, "kappa_fert_continuation": 0.3,
                    "chi": 1.0, "H0": 4.0, "theta0": 0.1, "theta1": 0.2,
                    **e5_profile.E5_FIXED,
                }), patch("sys.argv", ["run_e1_chain.py", "--outdir", tmp, "--smoke", "--minutes", "1", "--max-evals", "1"]):
                    chain.main()
                metadata = __import__("json").loads((Path(tmp) / "metadata.json").read_text())
        self.assertEqual(metadata["arm"], "E5")
        self.assertEqual(metadata["free_parameter_count"], 10)
        self.assertEqual(metadata["target_count"], 12)

    def test_collector_accepts_e5_contract(self) -> None:
        from intergen_eqscale_seq_optimized import collect_e1
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "chains"
            chain = root / "chain_001"
            chain.mkdir(parents=True)
            tight = {"strict_converged": True, "rank_loss": 1.0, "target_fit": []}
            summary = {
                "metadata": {"arm": "E5", "free_parameter_count": 10, "target_count": 12, "seed": 1},
                "best_tight": tight,
                "tight_repeat_check": {"both_strict": True, "loss_abs_difference": 0.0, "max_abs_moment_difference": 0.0},
                "n_cases_completed": 1, "n_strict": 1,
            }
            (chain / "summary.json").write_text(json.dumps(summary))
            outdir = Path(tmp) / "report"
            with patch("sys.argv", ["collect_e1.py", "--results-root", str(root), "--outdir", str(outdir)]):
                collect_e1.main()
            self.assertTrue((outdir / "results.json").exists())


if __name__ == "__main__":
    unittest.main()
