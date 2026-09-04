#!/usr/bin/env python3
"""Audit an E5F coordinate panel and diagnose calibration gaps/identification.

This is a read-only post-processor.  It never proposes or evaluates a new
parameter vector.  The local Jacobian is with respect to the driver's bounded
unit coordinates, which makes differently scaled structural parameters
comparable but remains a local numerical diagnostic rather than an economic
elasticity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


class ContractError(RuntimeError):
    """Raised when an input artifact does not satisfy the audit contract."""


@dataclass(frozen=True)
class Task:
    task_id: int
    summary_path: Path
    fit_path: Path
    summary: dict[str, Any]
    panel: dict[str, Any]
    candidate: dict[str, Any]
    unit: np.ndarray
    residual: np.ndarray
    model: np.ndarray
    loss: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coordinate-results-dir", type=Path, required=True)
    parser.add_argument("--selected-report-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ContractError(f"Cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise ContractError(f"Expected a JSON object in {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))
    except OSError as error:
        raise ContractError(f"Cannot read CSV {path}: {error}") from error


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise ContractError(f"Refusing to write an empty table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ContractError(f"{label} is not numeric: {value!r}") from error
    if not math.isfinite(result):
        raise ContractError(f"{label} is not finite: {value!r}")
    return result


def parse_target_rows(
    rows: list[dict[str, str]], expected_candidate: str
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    if not rows:
        raise ContractError("Expected at least one target row")
    moments: list[str] = []
    targets: list[float] = []
    models: list[float] = []
    weights: list[float] = []
    residuals: list[float] = []
    seen: set[str] = set()
    for row in rows:
        if row.get("candidate") != expected_candidate:
            raise ContractError(
                f"Target row candidate={row.get('candidate')!r}; "
                f"expected {expected_candidate!r}"
            )
        moment = str(row.get("moment") or "")
        if not moment or moment in seen:
            raise ContractError(f"Duplicate or empty moment {moment!r}")
        seen.add(moment)
        target = finite(row.get("target"), f"{moment} target")
        model = finite(row.get("model"), f"{moment} model")
        weight = finite(row.get("weight"), f"{moment} weight")
        gap = finite(row.get("gap"), f"{moment} gap")
        standardized = finite(row.get("standardized_gap"), f"{moment} standardized gap")
        contribution = finite(row.get("loss_contribution"), f"{moment} contribution")
        if weight <= 0.0:
            raise ContractError(f"{moment} has nonpositive weight")
        expected_gap = model - target
        expected_standardized = math.sqrt(weight) * expected_gap
        if not math.isclose(gap, expected_gap, rel_tol=1e-11, abs_tol=1e-12):
            raise ContractError(f"{moment} gap is inconsistent")
        if not math.isclose(
            standardized, expected_standardized, rel_tol=1e-10, abs_tol=1e-11
        ):
            raise ContractError(f"{moment} standardized gap is inconsistent")
        if not math.isclose(
            contribution, expected_standardized**2, rel_tol=1e-9, abs_tol=1e-10
        ):
            raise ContractError(f"{moment} loss contribution is inconsistent")
        moments.append(moment)
        targets.append(target)
        models.append(model)
        weights.append(weight)
        residuals.append(expected_standardized)
    residual = np.asarray(residuals, dtype=float)
    return (
        moments,
        np.asarray(targets, dtype=float),
        np.asarray(models, dtype=float),
        np.asarray(weights, dtype=float),
        residual,
        float(residual @ residual),
    )


def rank_report(jacobian: np.ndarray) -> dict[str, Any]:
    singular = np.linalg.svd(jacobian, compute_uv=False)
    largest = float(singular[0]) if singular.size else 0.0
    relative = singular / largest if largest > 0.0 else np.zeros_like(singular)
    thresholds = (1e-2, 1e-3, 1e-4)
    positive = singular[singular > np.finfo(float).eps * max(jacobian.shape) * largest]
    condition = float(largest / positive[-1]) if positive.size else math.inf
    return {
        "shape": list(jacobian.shape),
        "singular_values": [float(value) for value in singular],
        "relative_singular_values": [float(value) for value in relative],
        "relative_ranks": {
            f"relative_{threshold:g}": int(np.count_nonzero(relative >= threshold))
            for threshold in thresholds
        },
        "numerical_rank": int(np.linalg.matrix_rank(jacobian)),
        "condition_number_numerical_rank": condition,
    }


def central_jacobian(
    anchor: Task, pairs: list[tuple[Task, Task]]
) -> tuple[np.ndarray, list[dict[str, Any]], np.ndarray, np.ndarray]:
    columns: list[np.ndarray] = []
    forward_columns: list[np.ndarray] = []
    backward_columns: list[np.ndarray] = []
    side_rows: list[dict[str, Any]] = []
    domain = list(anchor.panel["domain"])
    for dimension, (minus, plus) in enumerate(pairs):
        x0 = float(anchor.unit[dimension])
        xm = float(minus.unit[dimension])
        xp = float(plus.unit[dimension])
        if not xm < x0 < xp:
            raise ContractError(f"Coordinate {dimension} is not an ordered central pair")
        expected_other = np.delete(anchor.unit, dimension)
        if not np.allclose(np.delete(minus.unit, dimension), expected_other, atol=2e-12, rtol=0.0):
            raise ContractError(f"Coordinate {dimension} minus task changes another dimension")
        if not np.allclose(np.delete(plus.unit, dimension), expected_other, atol=2e-12, rtol=0.0):
            raise ContractError(f"Coordinate {dimension} plus task changes another dimension")
        central = (plus.residual - minus.residual) / (xp - xm)
        backward = (anchor.residual - minus.residual) / (x0 - xm)
        forward = (plus.residual - anchor.residual) / (xp - x0)
        columns.append(central)
        backward_columns.append(backward)
        forward_columns.append(forward)
        denominator = max(float(np.linalg.norm(backward) + np.linalg.norm(forward)), 1e-15)
        side_disagreement = float(np.linalg.norm(forward - backward) / denominator)
        cosine_denominator = float(np.linalg.norm(forward) * np.linalg.norm(backward))
        cosine = (
            float(forward @ backward / cosine_denominator)
            if cosine_denominator > 0.0
            else 0.0
        )
        side_rows.append(
            {
                "dimension": dimension,
                "parameter": str(domain[dimension]["name"]),
                "unit_minus": xm,
                "unit_anchor": x0,
                "unit_plus": xp,
                "loss_minus": minus.loss,
                "loss_anchor": anchor.loss,
                "loss_plus": plus.loss,
                "side_disagreement_ratio": side_disagreement,
                "forward_backward_cosine": cosine,
            }
        )
    return (
        np.column_stack(columns),
        side_rows,
        np.column_stack(backward_columns),
        np.column_stack(forward_columns),
    )


def load_panel(root: Path) -> tuple[list[Task], list[str], np.ndarray, np.ndarray, dict[str, Any]]:
    root = root.resolve()
    task_dirs = sorted(path for path in root.glob("task_*" ) if path.is_dir())
    if [path.name for path in task_dirs] != [f"task_{index:03d}" for index in range(1, 24)]:
        raise ContractError("Coordinate panel must contain exactly task_001 through task_023")
    tasks: list[Task] = []
    common: dict[str, Any] | None = None
    common_moments: list[str] | None = None
    common_targets: np.ndarray | None = None
    common_weights: np.ndarray | None = None
    source_hashes: dict[str, Any] = {}
    for task_id, task_dir in enumerate(task_dirs, start=1):
        summary_path = task_dir / "summary.json"
        fit_path = task_dir / "target_fit_long.csv"
        summary = read_json(summary_path)
        if summary.get("status") != "complete_transition_calibration_panel_task":
            raise ContractError(f"task {task_id} is not complete")
        panel = dict(summary.get("panel_design") or {})
        candidate = dict(summary.get("best_candidate") or {})
        if int(panel.get("task_id", -1)) != task_id or int(panel.get("panel_size", -1)) != 23:
            raise ContractError(f"task {task_id} panel identity mismatch")
        if panel.get("panel_design") != "coordinate" or len(panel.get("domain") or []) != 11:
            raise ContractError(f"task {task_id} is not the declared 11-dimensional coordinate panel")
        expected_design = "anchor" if task_id == 1 else (
            f"coordinate_{(task_id - 2) // 2}_{'minus' if (task_id - 2) % 2 == 0 else 'plus'}"
        )
        if panel.get("design") != expected_design:
            raise ContractError(f"task {task_id} design mismatch")
        unit = np.asarray(panel.get("unit_vector"), dtype=float)
        if unit.shape != (11,) or not np.all(np.isfinite(unit)):
            raise ContractError(f"task {task_id} has an invalid unit vector")
        contract = {
            "source_sha256": summary.get("source_sha256"),
            "target_set": summary.get("target_set"),
            "target_fingerprint": summary.get("target_fingerprint"),
            "code_bundle_sha256": dict(summary.get("code_fingerprints") or {}).get("bundle_sha256"),
            "model_profile": dict(summary.get("model_profile") or {}).get("name"),
            "housing_supply_elasticity": dict(summary.get("external_closure_contract") or {}).get("housing_supply_elasticity"),
            "tenure_choice_kappa": dict(summary.get("external_closure_contract") or {}).get("tenure_choice_kappa"),
            "panel_seed": panel.get("panel_seed"),
            "center_sha256": panel.get("center_sha256"),
            "local_radius": panel.get("local_radius"),
            "domain": panel.get("domain"),
        }
        if common is None:
            common = contract
        elif json.dumps(contract, sort_keys=True) != json.dumps(common, sort_keys=True):
            raise ContractError(f"task {task_id} scientific/panel contract differs")
        candidate_label = str(candidate.get("candidate") or "")
        rows = read_csv(fit_path)
        moments, targets, models, weights, residual, loss = parse_target_rows(rows, candidate_label)
        reported_loss = finite(candidate.get("transition_loss"), f"task {task_id} loss")
        if not math.isclose(loss, reported_loss, rel_tol=1e-10, abs_tol=1e-9):
            raise ContractError(f"task {task_id} loss does not reconstruct")
        if common_moments is None:
            common_moments = moments
            common_targets = targets
            common_weights = weights
        elif moments != common_moments or not np.array_equal(targets, common_targets) or not np.array_equal(weights, common_weights):
            raise ContractError(f"task {task_id} target contract differs")
        tasks.append(Task(task_id, summary_path, fit_path, summary, panel, candidate, unit, residual, models, loss))
        source_hashes[f"task_{task_id:03d}"] = {
            "summary_json": sha256(summary_path),
            "target_fit_long_csv": sha256(fit_path),
        }
    assert common is not None and common_moments is not None
    assert common_targets is not None and common_weights is not None
    return tasks, common_moments, common_targets, common_weights, {
        "contract": common,
        "task_files": source_hashes,
    }


def main() -> None:
    args = parse_args()
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise ContractError(f"Refusing to overwrite nonempty output: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    tasks, moments, targets, weights, provenance = load_panel(args.coordinate_results_dir)
    anchor = tasks[0]
    pairs = [(tasks[1 + 2 * index], tasks[2 + 2 * index]) for index in range(11)]
    jacobian, side_rows, backward, forward = central_jacobian(anchor, pairs)
    report = args.selected_report_dir.resolve()
    selected_summary_path = report / "summary.json"
    selected_fit_path = report / "best_target_fit.csv"
    selected_summary = read_json(selected_summary_path)
    selected_candidate = dict(selected_summary.get("best_candidate") or {})
    selected_label = str(selected_candidate.get("candidate") or "")
    if not selected_label.startswith("task_"):
        raise ContractError(f"Selected report has invalid candidate: {selected_label!r}")
    try:
        selected_task_id = int(selected_label.removeprefix("task_"))
    except ValueError as error:
        raise ContractError(
            f"Selected report has invalid candidate: {selected_label!r}"
        ) from error
    if not 1 <= selected_task_id <= len(tasks):
        raise ContractError(f"Selected task is outside the coordinate panel: {selected_task_id}")
    selected_parameters_path = (
        report.parent / f"task_{selected_task_id:03d}" / "parameter_table.csv"
    )
    selected_rows_raw = read_csv(selected_fit_path)
    selected_moments, selected_targets, selected_models, selected_weights, selected_residual, selected_loss = parse_target_rows(
        selected_rows_raw, selected_label
    )
    if selected_moments != moments or not np.array_equal(selected_targets, targets) or not np.array_equal(selected_weights, weights):
        raise ContractError("Selected report target contract differs from the coordinate panel")
    if not math.isclose(selected_loss, float(selected_candidate["transition_loss"]), rel_tol=1e-10, abs_tol=1e-9):
        raise ContractError("Selected report loss does not reconstruct")
    selected_task = tasks[selected_task_id - 1]
    if not np.array_equal(selected_models, selected_task.model) or not np.array_equal(
        selected_residual, selected_task.residual
    ):
        raise ContractError(
            f"Selected best_target_fit.csv is not numerically identical to {selected_label}"
        )

    total_loss = float(selected_residual @ selected_residual)
    fit_rows: list[dict[str, Any]] = []
    for index, row in enumerate(selected_rows_raw):
        contribution = selected_residual[index] ** 2
        fit_rows.append(
            {
                "rank_by_loss": 0,
                "moment": moments[index],
                "target": targets[index],
                "model": selected_models[index],
                "gap": selected_models[index] - targets[index],
                "standardized_gap": selected_residual[index],
                "weight": weights[index],
                "loss_contribution": contribution,
                "share_of_total_loss": contribution / total_loss,
            }
        )
    fit_rows.sort(key=lambda row: float(row["loss_contribution"]), reverse=True)
    for rank, row in enumerate(fit_rows, start=1):
        row["rank_by_loss"] = rank
    write_csv(outdir / "calibration_gap_decomposition.csv", fit_rows)

    domain = list(anchor.panel["domain"])
    parameter_rows_raw = read_csv(selected_parameters_path)
    parameter_lookup = {row["parameter"]: row for row in parameter_rows_raw}
    parameter_rows: list[dict[str, Any]] = []
    sensitivity_rows: list[dict[str, Any]] = []
    for dimension, domain_row in enumerate(domain):
        name = str(domain_row["name"])
        column = jacobian[:, dimension]
        top_index = int(np.argmax(np.abs(column)))
        selected_parameter = parameter_lookup[name]
        parameter_rows.append(
            {
                "dimension": dimension,
                "parameter": name,
                "selected_value": finite(selected_parameter["value"], f"{name} selected value"),
                "lower_bound": finite(domain_row["lower"], f"{name} lower bound"),
                "upper_bound": finite(domain_row["upper"], f"{name} upper bound"),
                "transform": str(domain_row["transform"]),
                "near_bound": str(selected_parameter["near_bound"]).lower() == "true",
                "anchor_loss": anchor.loss,
                "minus_loss": pairs[dimension][0].loss,
                "plus_loss": pairs[dimension][1].loss,
                "best_panel_side": "minus" if pairs[dimension][0].loss < pairs[dimension][1].loss else "plus",
                "best_panel_loss": min(pairs[dimension][0].loss, pairs[dimension][1].loss),
                "column_l2_norm": float(np.linalg.norm(column)),
                "most_sensitive_moment": moments[top_index],
                "most_sensitive_abs_weighted_derivative": abs(float(column[top_index])),
                "side_disagreement_ratio": side_rows[dimension]["side_disagreement_ratio"],
                "forward_backward_cosine": side_rows[dimension]["forward_backward_cosine"],
            }
        )
        for moment_index, moment in enumerate(moments):
            sensitivity_rows.append(
                {
                    "dimension": dimension,
                    "parameter": name,
                    "moment": moment,
                    "central_d_standardized_gap_d_unit": float(column[moment_index]),
                    "backward_d_standardized_gap_d_unit": float(backward[moment_index, dimension]),
                    "forward_d_standardized_gap_d_unit": float(forward[moment_index, dimension]),
                    "central_d_model_moment_d_unit": float(column[moment_index] / math.sqrt(weights[moment_index])),
                }
            )
    write_csv(outdir / "parameter_diagnostics.csv", parameter_rows)
    write_csv(outdir / "moment_parameter_sensitivities.csv", sensitivity_rows)
    write_csv(outdir / "coordinate_side_consistency.csv", side_rows)

    jacobian_rows = []
    for moment_index, moment in enumerate(moments):
        jacobian_rows.append(
            {
                "moment": moment,
                "anchor_standardized_gap": float(anchor.residual[moment_index]),
                **{
                    str(domain[dimension]["name"]): float(jacobian[moment_index, dimension])
                    for dimension in range(11)
                },
            }
        )
    write_csv(outdir / "weighted_local_jacobian.csv", jacobian_rows)

    u, singular, vt = np.linalg.svd(jacobian, full_matrices=False)
    del u
    singular_rows = [
        {
            "order": index + 1,
            "singular_value": float(value),
            "relative_to_largest": float(value / singular[0]),
        }
        for index, value in enumerate(singular)
    ]
    write_csv(outdir / "singular_values.csv", singular_rows)
    weak_rows: list[dict[str, Any]] = []
    for weak_order, singular_index in enumerate(range(10, max(-1, 7), -1), start=1):
        vector = vt[singular_index]
        for dimension, loading in enumerate(vector):
            weak_rows.append(
                {
                    "weak_order": weak_order,
                    "singular_value": float(singular[singular_index]),
                    "parameter": str(domain[dimension]["name"]),
                    "loading": float(loading),
                    "squared_loading_share": float(loading**2),
                }
            )
    write_csv(outdir / "weak_directions.csv", weak_rows)

    rank = rank_report(jacobian)
    dominant_misses = [row["moment"] for row in fit_rows[:3]]
    top_three_share = sum(float(row["share_of_total_loss"]) for row in fit_rows[:3])
    theta1 = parameter_lookup["theta1"]
    summary = {
        "status": "complete_read_only_diagnostic",
        "scope": "local finite-difference diagnosis; not a recalibration and not an identification proof away from the anchor",
        "selected_candidate": selected_label,
        "selected_loss": total_loss,
        "anchor_candidate": "task_001",
        "anchor_loss": anchor.loss,
        "selected_improvement_from_anchor": anchor.loss - total_loss,
        "dominant_calibration_misses": dominant_misses,
        "top_three_loss_share": top_three_share,
        "first_birth_rooms_shortfall": float(
            selected_models[moments.index("housing_increment_0to1")]
            - targets[moments.index("housing_increment_0to1")]
        ),
        "mean_rooms_excess": float(
            selected_models[moments.index("aggregate_mean_occupied_rooms_18_85")]
            - targets[moments.index("aggregate_mean_occupied_rooms_18_85")]
        ),
        "owner_rate_gap_pp": 100.0
        * float(selected_models[moments.index("own_rate")] - targets[moments.index("own_rate")]),
        "young_owner_rate_gap_pp": (
            100.0
            * float(
                selected_models[moments.index("own_rate_2534")]
                - targets[moments.index("own_rate_2534")]
            )
            if "own_rate_2534" in moments
            else None
        ),
        "theta1_near_lower_bound": str(theta1["near_bound"]).lower() == "true",
        "local_weighted_jacobian": rank,
        "rank_gate_at_relative_1e_3_passed": rank["relative_ranks"]["relative_0.001"] == 11,
        "largest_side_disagreement": max(side_rows, key=lambda row: float(row["side_disagreement_ratio"])),
        "interpretation": (
            "The largest fitted gaps are housing moments, while the active target set has no causal housing-cost-to-fertility moment. "
            "The weighted local Jacobian diagnoses numerical sensitivity of the calibration moments, but it does not identify the policy fertility elasticity."
        ),
        "input_contract": provenance["contract"],
        "input_hashes": {
            **provenance["task_files"],
            "selected_summary_json": sha256(selected_summary_path),
            "selected_best_target_fit_csv": sha256(selected_fit_path),
            "selected_parameter_table_csv": sha256(selected_parameters_path),
        },
    }
    write_json(outdir / "summary.json", summary)

    short_labels = {
        "housing_increment_0to1": "First-birth rooms",
        "prime30_55_parent_3plus_minus_1to2_mean_rooms": "3+ vs 1-2 rooms",
        "aggregate_mean_occupied_rooms_18_85": "Mean rooms",
        "own_rate": "Ownership",
        "own_rate_2534": "Young ownership (25--34)",
        "mean_age_first_birth": "Age at first birth",
        "share_first_births_age30plus": "First births 30+",
        "old_total_wealth_to_annual_income_p90_p50_7684": "Old wealth p90/p50",
        "own_family_gap": "Family ownership gap",
        "annual_bequest_flow_to_aggregate_wealth": "Bequest/wealth",
        "aggregate_wealth_to_annual_gross_labor_earnings": "Wealth/earnings",
        "childless_rate": "Childlessness",
        "tfr": "Completed fertility",
    }
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.2), constrained_layout=True)
    ordered = list(reversed(fit_rows))
    axes[0, 0].barh(
        [short_labels.get(str(row["moment"]), str(row["moment"])) for row in ordered],
        [float(row["loss_contribution"]) for row in ordered],
        color="#3b6c8e",
    )
    axes[0, 0].set_title("Where the calibration loss comes from")
    axes[0, 0].set_xlabel("Weighted loss contribution")
    axes[0, 0].grid(axis="x", alpha=0.2)

    indices = np.arange(len(moments))
    axes[0, 1].barh(
        indices,
        selected_residual,
        color=np.where(selected_residual >= 0.0, "#b0493d", "#2c7fb8"),
    )
    axes[0, 1].axvline(0.0, color="#555555", linewidth=0.8)
    axes[0, 1].set_yticks(indices, [short_labels.get(moment, moment) for moment in moments])
    axes[0, 1].set_title("Signed fit gaps")
    axes[0, 1].set_xlabel("Standardized gap (model minus target)")
    axes[0, 1].grid(axis="x", alpha=0.2)

    column_scale = np.maximum(np.linalg.norm(jacobian, axis=0), 1e-15)
    heat = jacobian / column_scale
    image = axes[1, 0].imshow(heat, aspect="auto", cmap="RdBu_r", vmin=-1.0, vmax=1.0)
    axes[1, 0].set_yticks(indices, [short_labels.get(moment, moment) for moment in moments])
    axes[1, 0].set_xticks(
        np.arange(11),
        [str(row["name"]).replace("_", "\n") for row in domain],
        rotation=45,
        ha="right",
        fontsize=8,
    )
    axes[1, 0].set_title("Local moment sensitivity by parameter")
    axes[1, 0].set_xlabel("Each column normalized to unit length")
    fig.colorbar(image, ax=axes[1, 0], fraction=0.03, pad=0.02)

    axes[1, 1].semilogy(np.arange(1, len(singular) + 1), singular / singular[0], marker="o", color="#6a51a3")
    axes[1, 1].axhline(1e-3, color="#b0493d", linestyle="--", linewidth=1.0, label="declared rank threshold")
    axes[1, 1].set_xticks(np.arange(1, len(singular) + 1))
    axes[1, 1].set_title("Local identification spectrum")
    axes[1, 1].set_xlabel("Singular-value order")
    axes[1, 1].set_ylabel("Relative singular value")
    axes[1, 1].grid(alpha=0.2)
    axes[1, 1].legend(frameon=False)
    for axis in axes.flat:
        axis.spines[["top", "right"]].set_visible(False)
    fig.savefig(outdir / "calibration_gap_identification.pdf", bbox_inches="tight")
    fig.savefig(outdir / "calibration_gap_identification.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    print(
        "E5F_CALIBRATION_PANEL_DIAGNOSTIC_COMPLETE "
        f"loss={total_loss:.12g} rank_1e-3={rank['relative_ranks']['relative_0.001']}/11 "
        f"output={outdir}",
        flush=True,
    )


if __name__ == "__main__":
    try:
        main()
    except ContractError as error:
        raise SystemExit(f"E5F_CALIBRATION_PANEL_DIAGNOSTIC_REJECTED {error}") from error
