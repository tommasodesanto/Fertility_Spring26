"""Regression tests for the empirical outside-entry target contract."""

from pathlib import Path
from types import SimpleNamespace
import sys
import unittest


TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from entry_target_contract import (
    entry_target_from_outside_origin_share,
    quota_closure_from_outside_origin_share,
    quota_population_scale,
)


class EntryTargetContractTests(unittest.TestCase):
    def test_outside_origin_share_maps_to_candidate_specific_entry_probability(self) -> None:
        solution = SimpleNamespace(
            entry_rate=0.06173345618074889,
            entrants_mature_total=0.05292942076730084,
        )
        parameters = SimpleNamespace(local_birth_entry_weight=1.0)

        contract = entry_target_from_outside_origin_share(
            solution,
            parameters,
            outside_origin_share=0.169,
        )

        self.assertAlmostEqual(contract["target_qbar"], 0.9692247023019598)
        self.assertAlmostEqual(contract["outside_origin_entrant_share_identity_check"], 0.169)

    def test_infeasible_empirical_share_is_rejected(self) -> None:
        solution = SimpleNamespace(entry_rate=1.0, entrants_mature_total=0.2)
        parameters = SimpleNamespace(local_birth_entry_weight=1.0)

        with self.assertRaisesRegex(RuntimeError, "infeasible model entry probability"):
            entry_target_from_outside_origin_share(
                solution,
                parameters,
                outside_origin_share=0.169,
            )


class QuotaClosureContractTests(unittest.TestCase):
    def test_floor_and_tilt_baselines_reproduce_scale_and_outside_share(self) -> None:
        candidates = (
            (0.061733456180748894, 0.05292942076730084, 0.9692247023019599),
            (0.06173345618074878, 0.05289560240559147, 0.9698443680221587),
        )
        parameters = SimpleNamespace(local_birth_entry_weight=1.0)
        for entry_flow, mature_flow, expected_retention in candidates:
            with self.subTest(expected_retention=expected_retention):
                solution = SimpleNamespace(
                    entry_rate=entry_flow,
                    entrants_mature_total=mature_flow,
                    total_mass=1.0,
                )
                quota = quota_closure_from_outside_origin_share(
                    solution,
                    parameters,
                    outside_origin_share=0.169,
                )
                scale, nodes = quota_population_scale(solution, parameters, quota)

                self.assertAlmostEqual(quota["retention_rate"], expected_retention)
                self.assertAlmostEqual(scale.scale_factor, 1.0)
                self.assertAlmostEqual(scale.outside_origin_entrant_share, 0.169)
                self.assertAlmostEqual(scale.entry_residual, 0.0)
                self.assertEqual(nodes, [])

    def test_policy_scale_uses_fixed_retention_and_inflow_identity(self) -> None:
        baseline = SimpleNamespace(
            entry_rate=1.0,
            entrants_mature_total=1.0,
            total_mass=1.0,
        )
        parameters = SimpleNamespace(local_birth_entry_weight=1.0)
        quota = quota_closure_from_outside_origin_share(
            baseline,
            parameters,
            outside_origin_share=0.2,
        )
        policy = SimpleNamespace(entrants_mature_total=1.01, total_mass=1.0)

        scale, _ = quota_population_scale(policy, parameters, quota)

        self.assertAlmostEqual(quota["retention_rate"], 0.8)
        self.assertAlmostEqual(quota["outside_inflow"], 0.2)
        self.assertAlmostEqual(scale.scale_factor, 0.2 / (1.0 - 0.8 * 1.01))
        self.assertAlmostEqual(scale.entry_residual, 0.0)

    def test_quota_mode_does_not_read_values_or_logit_objects(self) -> None:
        class PoisonSolution:
            entry_rate = 1.0
            entrants_mature_total = 1.0
            total_mass = 1.0

            @property
            def V(self):  # pragma: no cover - access is the test failure
                raise AssertionError("quota closure must not read the value function")

        solution = PoisonSolution()
        parameters = SimpleNamespace(local_birth_entry_weight=1.0)
        quota = quota_closure_from_outside_origin_share(
            solution,
            parameters,
            outside_origin_share=0.2,
        )

        self.assertNotIn("outside_value", quota)
        self.assertNotIn("entry_taste_scale", quota)
        scale, _ = quota_population_scale(solution, parameters, quota)
        self.assertTrue(scale.finite)

    def test_retention_above_one_fails_loudly(self) -> None:
        solution = SimpleNamespace(entry_rate=1.0, entrants_mature_total=0.2)
        parameters = SimpleNamespace(local_birth_entry_weight=1.0)

        with self.assertRaisesRegex(RuntimeError, "Rbar exceeds one"):
            quota_closure_from_outside_origin_share(
                solution,
                parameters,
                outside_origin_share=0.169,
            )


if __name__ == "__main__":
    unittest.main()
