from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest


MODULE_PATH = Path(__file__).with_name("analyze_e5f_transition_calibration_panel.py")
SPEC = importlib.util.spec_from_file_location("calibration_panel_diagnostic", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def task(task_id: int, unit: list[float], residual: list[float]) -> object:
    return MODULE.Task(
        task_id=task_id,
        summary_path=Path("summary.json"),
        fit_path=Path("fit.csv"),
        summary={},
        panel={"domain": [{"name": "x"}]},
        candidate={},
        unit=np.asarray(unit, dtype=float),
        residual=np.asarray(residual, dtype=float),
        model=np.zeros(len(residual)),
        loss=float(np.asarray(residual) @ np.asarray(residual)),
    )


def test_central_jacobian_recovers_linear_derivative() -> None:
    anchor = task(1, [0.5], [1.0, -2.0])
    minus = task(2, [0.4], [0.7, -1.5])
    plus = task(3, [0.6], [1.3, -2.5])
    jacobian, sides, backward, forward = MODULE.central_jacobian(anchor, [(minus, plus)])
    np.testing.assert_allclose(jacobian[:, 0], [3.0, -5.0])
    np.testing.assert_allclose(backward[:, 0], [3.0, -5.0])
    np.testing.assert_allclose(forward[:, 0], [3.0, -5.0])
    assert sides[0]["side_disagreement_ratio"] < 1e-12
    assert sides[0]["forward_backward_cosine"] == pytest.approx(1.0)


def test_rank_report_flags_weak_direction() -> None:
    report = MODULE.rank_report(np.diag([1.0, 1e-5]))
    assert report["relative_ranks"]["relative_0.001"] == 1
    assert report["numerical_rank"] == 2


def test_target_parser_rejects_inconsistent_gap() -> None:
    rows = []
    for index in range(12):
        rows.append(
            {
                "candidate": "task_001",
                "moment": f"m{index}",
                "target": "0",
                "model": "1",
                "gap": "1",
                "weight": "4",
                "loss_contribution": "4",
                "standardized_gap": "2",
            }
        )
    rows[4]["gap"] = "0"
    with pytest.raises(MODULE.ContractError, match="gap is inconsistent"):
        MODULE.parse_target_rows(rows, "task_001")


def test_target_parser_accepts_overidentified_row_count() -> None:
    rows = []
    for index in range(13):
        rows.append(
            {
                "candidate": "task_001",
                "moment": f"m{index}",
                "target": "0",
                "model": "1",
                "gap": "1",
                "weight": "4",
                "loss_contribution": "4",
                "standardized_gap": "2",
            }
        )
    moments, _, _, _, residual, loss = MODULE.parse_target_rows(rows, "task_001")
    assert len(moments) == 13
    assert residual.shape == (13,)
    assert loss == pytest.approx(52.0)
