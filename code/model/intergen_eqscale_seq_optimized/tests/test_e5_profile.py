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
from intergen_eqscale_seq_optimized.diagnostics import completed_fertility_measure
from intergen_eqscale_seq_optimized.solver import add_aggregate_wealth_bequest_flow_moments


class E5ProfileTests(unittest.TestCase):
    def test_pending_timing_targets_refuse_to_build(self) -> None:
        # The shipped profile carries the NCHS values; the refusal mechanism
        # must still fire if any target is unset.
        with patch.dict(e5_profile.E5_TARGETS, {"mean_age_first_birth": None}):
            with self.assertRaisesRegex(ValueError, "pending"):
                e5_profile.e5_target_system()

    def test_shipped_timing_targets_are_the_nchs_build_values(self) -> None:
        self.assertAlmostEqual(e5_profile.E5_TARGETS["mean_age_first_birth"], 25.310560799362)
        self.assertAlmostEqual(e5_profile.E5_TARGETS["share_first_births_age30plus"], 0.270062376851342)
        self.assertEqual(len(e5_profile.e5_target_system().targets_dict()), 12)

    def test_rooms_target_uses_idfe_estimate_and_measured_se(self) -> None:
        name = "housing_increment_0to1"
        self.assertEqual(e5_profile.E5_TARGET_SET, "e5_idfe_review_20260809")
        self.assertAlmostEqual(e5_profile.E5_TARGETS[name], 0.80494368)
        self.assertAlmostEqual(e5_profile.E5_MEASURED_SES[name], 0.16728361)
        self.assertAlmostEqual(
            e5_profile.E5_WEIGHTS[name],
            1.0 / e5_profile.E5_MEASURED_SES[name] ** 2,
        )

    def test_standard_diagnostic_uses_literal_top_group_measurement(self) -> None:
        sol = SimpleNamespace(
            mean_completed_fertility=1.5,
            parity_dist=np.array([0.2, 0.3, 0.3, 0.2]),
        )
        parameters = SimpleNamespace(
            fertility_units="literal_topcode",
            tfr_top_bin_weight=3.6,
        )
        self.assertAlmostEqual(
            completed_fertility_measure(sol, parameters),
            0.3 + 2.0 * 0.3 + 3.6 * 0.2,
        )

    def test_dummy_timing_targets_complete_twelve_row_contract(self) -> None:
        with patch.dict(e5_profile.E5_TARGETS, {
            "mean_age_first_birth": 26.5,
            "share_first_births_age30plus": 0.35,
        }):
            system = e5_profile.e5_target_system()
        self.assertEqual(system.count, 12)
        self.assertEqual(set(system.targets_dict()), set(system.weights_dict()))
        self.assertEqual(len(e5_profile.E5_DOMAIN), 10)
        self.assertEqual(e5_profile.E5_DOMAIN[0][:3], ("beta_annual", 0.94, 0.9995))

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
            theta = {
                ("beta" if name == "beta_annual" else name): (
                    0.96**4 if name == "beta_annual" else (lower + upper) / 2.0
                )
                for name, lower, upper, _transform in e5_profile.E5_DOMAIN
            }
            tight = {
                "strict_converged": True,
                "rank_loss": 1.0,
                "target_fit": [],
                "theta": theta,
            }
            summary = {
                "metadata": {
                    "arm": "E5",
                    "free_parameter_count": 10,
                    "target_count": 12,
                    "seed": 1,
                    "active_domain": [
                        {
                            "name": name,
                            "lower": lower,
                            "upper": upper,
                            "transform": transform,
                        }
                        for name, lower, upper, transform in e5_profile.E5_DOMAIN
                    ],
                },
                "best_tight": tight,
                "tight_repeat_check": {"both_strict": True, "loss_abs_difference": 0.0, "max_abs_moment_difference": 0.0},
                "n_cases_completed": 1, "n_strict": 1,
            }
            (chain / "summary.json").write_text(json.dumps(summary))
            outdir = Path(tmp) / "report"
            with patch("sys.argv", ["collect_e1.py", "--results-root", str(root), "--outdir", str(outdir)]):
                collect_e1.main()
            self.assertTrue((outdir / "results.json").exists())

    def test_e6a_external_has_no_added_free_parameter(self) -> None:
        from intergen_eqscale_seq_optimized.e6a_profile import (
            e6a_metadata,
            e6a_overrides,
        )

        overrides = e6a_overrides()
        metadata = e6a_metadata()
        self.assertEqual(len(e5_profile.E5_DOMAIN), 10)
        self.assertAlmostEqual(overrides["fecundity_omega1"], 0.01331)
        self.assertAlmostEqual(overrides["fecundity_omega2"], 0.14960)
        self.assertAlmostEqual(
            metadata["terminal_success_probability_before_hard_close"], 0.03
        )
        self.assertGreater(overrides["fecundity_terminal_decay"], 0.0)


if __name__ == "__main__":
    unittest.main()
