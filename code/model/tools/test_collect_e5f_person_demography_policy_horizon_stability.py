"""Unit tests for the person-demography policy horizon audit."""

from __future__ import annotations

import copy
import unittest
from pathlib import Path
from typing import Any

import collect_e5f_person_demography_policy_horizon_stability as audit


def _row(period: int, level: float, owner_rate: float = 0.55) -> dict[str, str]:
    return {
        "period": str(period),
        "calendar_year": str(2023 + 4 * period),
        "asset_price": str(level),
        "renter_price": str(0.2 * level),
        "equal_transfer_period_units": str(0.3 * level),
        "birth_children_topcode_adjusted": str(0.08 * level),
        "resident_persons_actual_scale": str(300_000_000.0 * level),
        "resident_persons": str(3.0 * level),
        "household_heads": str(level),
        "owner_rate": str(owner_rate),
        "relative_market_residual": "0.0001",
        "government_budget_residual": "0.00001",
    }


def _run(
    case: str, horizon: int, reform_scale: float = 1.0, owner_rate: float = 0.55
) -> dict[str, Any]:
    rows = [
        _row(period, reform_scale * (1.0 + 0.01 * period), owner_rate)
        for period in range(horizon)
    ]
    return {
        "case": case,
        "horizon": horizon,
        "rows": rows,
        "summary": {
            "maximum_market_residual": 0.0001,
            "maximum_fiscal_residual": 0.00001,
            "path_converged": True,
            "iterations": 4,
            "elapsed_seconds": 10.0,
            "terminal_convergence": {"all_checks_pass": True},
            "source_hashes": {"driver": "same-source"},
            "terminal_root": {"case": case, "root": "same-root"},
        },
        "summary_path": Path("summary.json"),
        "summary_sha256": "summary",
        "path_file": Path("transition_path.csv"),
        "path_sha256": "path",
    }


class PersonDemographyHorizonAuditTests(unittest.TestCase):
    def test_path_gate_is_recomputed_from_saved_rows(self) -> None:
        run = _run("rebated-tax1-baseline", 2)
        self.assertTrue(audit.path_gate_audit(run)["path_checks_pass"])
        inconsistent = copy.deepcopy(run)
        inconsistent["summary"]["maximum_market_residual"] = 0.00015
        result = audit.path_gate_audit(inconsistent)
        self.assertFalse(result["summary_matches_saved_path"])
        self.assertFalse(result["path_checks_pass"])

    def test_level_stability_uses_declared_relative_gate(self) -> None:
        short = _run("rebated-tax1-baseline", 2)
        reference = _run("rebated-tax1-baseline", 3)
        self.assertEqual(
            audit.level_stability([short, reference], 2, 1e-3)["status"],
            "passed",
        )
        reference["rows"][0]["asset_price"] = "1.002"
        self.assertEqual(
            audit.level_stability([short, reference], 2, 1e-3)["status"],
            "failed",
        )

    def test_policy_effect_stability_catches_ownership_gap(self) -> None:
        baseline = [
            _run("rebated-tax1-baseline", 2),
            _run("rebated-tax1-baseline", 3),
        ]
        reform = [
            _run("rebated-tax2-reform", 2, reform_scale=1.01, owner_rate=0.56),
            _run("rebated-tax2-reform", 3, reform_scale=1.01, owner_rate=0.56),
        ]
        self.assertEqual(
            audit.effect_stability(baseline, reform, 2, 0.02)["status"],
            "passed",
        )
        reform[1]["rows"][0]["owner_rate"] = "0.5605"
        result = audit.effect_stability(baseline, reform, 2, 0.02)
        self.assertEqual(result["status"], "failed")
        self.assertGreater(result["maximum_absolute_gap"], 0.02)

    def test_common_contract_requires_matching_sources_roots_and_initial_state(self) -> None:
        runs = [
            _run("rebated-tax1-baseline", 2),
            _run("rebated-tax1-baseline", 3),
            _run("rebated-tax2-reform", 2),
            _run("rebated-tax2-reform", 3),
        ]
        self.assertEqual(audit.contract_gate(runs)["status"], "passed")

        wrong_source = copy.deepcopy(runs)
        wrong_source[-1]["summary"]["source_hashes"] = {"driver": "changed"}
        self.assertEqual(audit.contract_gate(wrong_source)["status"], "failed")

        wrong_root = copy.deepcopy(runs)
        wrong_root[1]["summary"]["terminal_root"] = {
            "case": "rebated-tax1-baseline",
            "root": "changed",
        }
        self.assertEqual(audit.contract_gate(wrong_root)["status"], "failed")

        wrong_initial_state = copy.deepcopy(runs)
        wrong_initial_state[2]["rows"][0]["resident_persons"] = "3.1"
        self.assertEqual(
            audit.contract_gate(wrong_initial_state)["status"], "failed"
        )


if __name__ == "__main__":
    unittest.main()
