#!/usr/bin/env python3
"""Focused tests for the fixed-share conditional transition wrapper."""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np

import run_e5f_person_policy_tenure_conditional_transition as conditional


class ConditionalHelpersTests(unittest.TestCase):
    def test_extract_atom_options_and_choices(self) -> None:
        atom = conditional.AtomSpec(
            calendar_year=2023,
            wealth_index=0,
            origin_tenure=0,
            market=0,
            age_index=1,
            child_state=0,
            renter_tenure=0,
            owner_tenure=2,
        )
        choices = np.zeros((2, 3, 1, 2, 2, 1, 1), dtype=np.int16)
        choices[0, 0, 0, 1, 1, 0, 0] = 2
        policy = SimpleNamespace(tenure_choice=choices)
        slices = {
            z: np.full((2, 3, 1, 1, 1, 3), float(z)) for z in range(2)
        }
        options, selected = conditional.extract_atom_options_and_choices(
            policy,
            slices,
            number_of_income_states=2,
            atom=atom,
        )
        self.assertEqual(options.shape, (2, 1, 3))
        np.testing.assert_array_equal(options[1], 1.0)
        np.testing.assert_array_equal(selected[:, 0], [0, 2])

    def test_path_array_hashes_use_only_prices_and_transfers(self) -> None:
        rows = [
            {"asset_price": 1.0, "equal_transfer_period_units": 0.1, "other": 7},
            {"asset_price": 2.0, "equal_transfer_period_units": 0.2, "other": 8},
        ]
        first = conditional.path_array_hashes(rows)
        rows[0]["other"] = 99
        self.assertEqual(first, conditional.path_array_hashes(rows))
        rows[0]["asset_price"] = 1.1
        self.assertNotEqual(first, conditional.path_array_hashes(rows))

    def test_atom_defaults_select_distinct_tenures(self) -> None:
        atom = conditional.AtomSpec()
        self.assertEqual(atom.owner_tenure, 2)
        self.assertNotEqual(atom.owner_tenure, atom.renter_tenure)

    def test_atomic_writers_leave_complete_files(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            json_path = root / "checkpoint.json"
            csv_path = root / "checkpoint.csv"
            conditional.write_json_atomic(json_path, {"value": 3})
            conditional.write_csv_atomic(csv_path, [{"a": 1, "b": 2}])
            self.assertEqual(json_path.read_text(encoding="utf-8"), '{\n  "value": 3\n}\n')
            self.assertEqual(csv_path.read_text(encoding="utf-8"), "a,b\n1,2\n")
            self.assertFalse((root / "checkpoint.json.tmp").exists())
            self.assertFalse((root / "checkpoint.csv.tmp").exists())

    def test_production_packet_audit_passes_and_rejects_wrong_flag(self) -> None:
        rows = [
            {
                "period": period,
                "asset_price": 1.0 + period,
                "equal_transfer_period_units": 0.1,
                "relative_market_residual": 1.0e-4,
                "government_budget_residual": 1.0e-5,
                "feasibility_frontier_projection_mass": 0.0,
            }
            for period in range(2)
        ]
        final_record = {
            "maximum_market_residual": 1.0e-4,
            "maximum_fiscal_residual": 1.0e-5,
            "maximum_feasibility_projection_mass": 0.0,
            **conditional.path_array_hashes(rows),
        }
        summary = {
            "maximum_market_residual": 1.0e-4,
            "maximum_fiscal_residual": 1.0e-5,
            "path_converged": True,
            "history_reproduction_status": "passed",
            "accounting_gates": {"all": True},
            "terminal_convergence": {"all_checks_pass": True},
        }
        audit = conditional.audit_production_packet(
            summary=summary,
            rows=rows,
            final_record=final_record,
            expected_horizon=2,
        )
        self.assertTrue(audit["path_converged"])
        summary["path_converged"] = False
        with self.assertRaisesRegex(RuntimeError, "convergence flag"):
            conditional.audit_production_packet(
                summary=summary,
                rows=rows,
                final_record=final_record,
                expected_horizon=2,
            )

    def test_default_off_wrapper_delegates_without_changing_result(self) -> None:
        original_solve = conditional.pf.solve_date_policy
        original_evaluate = (
            conditional.person_pf.evaluate_path_at_prices_person_demography
        )
        sentinel_policy = object()
        sentinel_evaluation = object()

        def fake_solve(**_kwargs: object) -> object:
            return sentinel_policy

        def fake_evaluate(*_args: object, **_kwargs: object) -> object:
            return sentinel_evaluation

        conditional.pf.solve_date_policy = fake_solve
        conditional.person_pf.evaluate_path_at_prices_person_demography = fake_evaluate
        try:
            with TemporaryDirectory() as temporary:
                context = conditional.ConditionalCapture(
                    output_dir=Path(temporary),
                    atom=conditional.AtomSpec(),
                    owner_share=0.5,
                )
                restore = conditional.install_conditional_share(context)
                self.assertIs(conditional.pf.solve_date_policy(test=True), sentinel_policy)
                restore()
                self.assertIs(conditional.pf.solve_date_policy, fake_solve)
                self.assertIs(
                    conditional.person_pf.evaluate_path_at_prices_person_demography,
                    fake_evaluate,
                )
        finally:
            conditional.pf.solve_date_policy = original_solve
            conditional.person_pf.evaluate_path_at_prices_person_demography = (
                original_evaluate
            )

    def test_source_hashes_pin_active_sequential_solver(self) -> None:
        expected = conditional.file_sha256(
            conditional.MODEL_ROOT / "intergen_eqscale_seq_optimized/solver.py"
        )
        self.assertEqual(conditional.source_hashes()["solver"], expected)


if __name__ == "__main__":
    unittest.main()
