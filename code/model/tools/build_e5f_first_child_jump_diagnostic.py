#!/usr/bin/env python3
"""Collect the fixed first-child housing-jump diagnostic.

This script performs no model solve and no selection.  It validates the six
predeclared external jump values, reconciles every twelve-row loss ledger, and
writes one compact diagnostic packet.  The exercise is a mechanism screen; it
does not turn the jump into an estimated parameter.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt


JUMPS = (0.0, 0.10, 0.20, 0.25, 0.35, 0.50)
KEY_MOMENTS = (
    "housing_increment_0to1",
    "prime30_55_parent_3plus_minus_1to2_mean_rooms",
    "aggregate_mean_occupied_rooms_18_85",
    "own_family_gap",
    "own_rate",
    "tfr",
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-center-sha256", required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-target-set", required=True)
    parser.add_argument("--expected-target-fingerprint", required=True)
    parser.add_argument("--expected-code-bundle-sha256", required=True)
    parser.add_argument("--expected-launcher-sha256", required=True)
    parser.add_argument("--expected-baseline-loss", type=float, required=True)
    return parser.parse_args()


def require_equal(actual: Any, expected: Any, label: str) -> None:
    if str(actual) != str(expected):
        raise RuntimeError(f"{label} mismatch: {actual!r} versus {expected!r}")


def read_target_fit(path: Path, expected_loss: float) -> list[dict[str, Any]]:
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    if len(rows) != 12 or len({row["moment"] for row in rows}) != 12:
        raise RuntimeError(f"{path}: target ledger must have twelve unique rows")
    contribution_sum = 0.0
    for row in rows:
        for name in ("target", "model", "gap", "weight", "loss_contribution"):
            value = float(row[name])
            if not math.isfinite(value):
                raise RuntimeError(f"{path}: nonfinite {name}")
        target = float(row["target"])
        model = float(row["model"])
        gap = float(row["gap"])
        weight = float(row["weight"])
        contribution = float(row["loss_contribution"])
        if not math.isclose(gap, model - target, rel_tol=0.0, abs_tol=2e-10):
            raise RuntimeError(f"{path}: gap identity failed for {row['moment']}")
        if not math.isclose(contribution, weight * gap * gap, rel_tol=0.0, abs_tol=2e-9):
            raise RuntimeError(f"{path}: loss identity failed for {row['moment']}")
        contribution_sum += contribution
    if not math.isclose(contribution_sum, expected_loss, rel_tol=0.0, abs_tol=2e-8):
        raise RuntimeError(f"{path}: target contributions do not sum to total loss")
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write an empty table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    root = args.results_root.resolve()
    out = args.output_dir.resolve()
    if out.exists():
        raise RuntimeError(f"refusing to overwrite output directory: {out}")

    frontier: list[dict[str, Any]] = []
    target_long: list[dict[str, Any]] = []
    input_hashes: dict[str, str] = {}
    for index, jump in enumerate(JUMPS, start=1):
        task = root / f"factor_{index:03d}"
        summary_path = task / "summary.json"
        fit_path = task / "target_fit_long.csv"
        contract_path = task / "jump_task_contract.json"
        if not summary_path.is_file() or not fit_path.is_file() or not contract_path.is_file():
            raise RuntimeError(f"factor {index}: incomplete artifacts")
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
        require_equal(contract.get("schema"), "e5f_first_child_jump_task_contract_v1", "task contract")
        require_equal(contract.get("status"), "complete", "task completion")
        require_equal(contract.get("mode"), "factor", "task mode")
        require_equal(contract.get("summary_sha256"), digest(summary_path), "task summary hash")
        require_equal(contract.get("center_sha256"), args.expected_center_sha256, "task center")
        require_equal(contract.get("source_sha256"), args.expected_source_sha256, "task source")
        require_equal(contract.get("target_set"), args.expected_target_set, "task target set")
        require_equal(
            contract.get("target_fingerprint"),
            args.expected_target_fingerprint,
            "task target fingerprint",
        )
        require_equal(
            contract.get("code_bundle_sha256"),
            args.expected_code_bundle_sha256,
            "task code bundle",
        )
        require_equal(
            contract.get("launcher_sha256"),
            args.expected_launcher_sha256,
            "task launcher",
        )
        if not math.isclose(
            float(contract.get("hbar_first_child_jump")), jump, rel_tol=0.0, abs_tol=1e-14
        ):
            raise RuntimeError(f"factor {index}: task-contract jump mismatch")
        require_equal(
            summary.get("status"),
            "complete_transition_calibration_panel_task",
            f"factor {index} status",
        )
        require_equal(summary.get("source_sha256"), args.expected_source_sha256, "source")
        require_equal(summary.get("target_set"), args.expected_target_set, "target set")
        require_equal(
            summary.get("target_fingerprint"),
            args.expected_target_fingerprint,
            "target fingerprint",
        )
        require_equal(
            summary["code_fingerprints"]["bundle_sha256"],
            args.expected_code_bundle_sha256,
            "code bundle",
        )
        panel = summary.get("panel_design", {})
        require_equal(panel.get("panel_size"), 1, "panel size")
        require_equal(panel.get("task_id"), 1, "panel task")
        require_equal(panel.get("design"), "anchor", "panel design")
        require_equal(panel.get("center_sha256"), args.expected_center_sha256, "center")
        profile = summary.get("model_profile", {})
        require_equal(profile.get("name"), "e5f-income-entry", "model profile")
        require_equal(
            profile.get("first_child_room_jump_status"),
            "externally fixed diagnostic; not included in the free-parameter count",
            "jump status",
        )
        if not math.isclose(
            float(profile.get("first_child_room_jump")), jump, rel_tol=0.0, abs_tol=1e-14
        ):
            raise RuntimeError(f"factor {index}: wrong first-child jump")

        best = summary["best_candidate"]
        loss = float(best["transition_loss"])
        if not math.isfinite(loss):
            raise RuntimeError(f"factor {index}: nonfinite loss")
        fit = read_target_fit(fit_path, loss)
        fit_by_name = {row["moment"]: row for row in fit}
        row: dict[str, Any] = {
            "factor_index": index,
            "hbar_first_child_jump": jump,
            "transition_loss": loss,
            "old_psi_child": float(summary["old_psi_child"]),
            "new_psi_child": float(best["new_psi_child"]),
        }
        for moment in KEY_MOMENTS:
            row[f"{moment}__target"] = float(fit_by_name[moment]["target"])
            row[f"{moment}__model"] = float(fit_by_name[moment]["model"])
            row[f"{moment}__loss"] = float(fit_by_name[moment]["loss_contribution"])
        frontier.append(row)
        for fit_row in fit:
            target_long.append(
                {
                    "factor_index": index,
                    "hbar_first_child_jump": jump,
                    **fit_row,
                }
            )
        input_hashes[f"factor_{index:03d}/summary.json"] = digest(summary_path)
        input_hashes[f"factor_{index:03d}/target_fit_long.csv"] = digest(fit_path)
        input_hashes[f"factor_{index:03d}/jump_task_contract.json"] = digest(contract_path)

    if not math.isclose(
        float(frontier[0]["transition_loss"]),
        float(args.expected_baseline_loss),
        rel_tol=0.0,
        abs_tol=2e-8,
    ):
        raise RuntimeError("zero-jump factor does not reproduce the certified baseline loss")

    out.mkdir(parents=True)
    write_csv(out / "jump_frontier.csv", frontier)
    write_csv(out / "target_fit_long.csv", target_long)

    jump_values = [float(row["hbar_first_child_jump"]) for row in frontier]
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), constrained_layout=True)
    panels = (
        ("transition_loss", "SMM loss", None),
        ("housing_increment_0to1", "First-birth rooms response", True),
        (
            "prime30_55_parent_3plus_minus_1to2_mean_rooms",
            "3+ minus 1–2 rooms",
            True,
        ),
        ("aggregate_mean_occupied_rooms_18_85", "Mean occupied rooms", True),
    )
    for axis, (field, label, show_target) in zip(axes.flat, panels, strict=True):
        if field == "transition_loss":
            values = [float(row[field]) for row in frontier]
            target = None
        else:
            values = [float(row[f"{field}__model"]) for row in frontier]
            target = float(frontier[0][f"{field}__target"])
        axis.plot(jump_values, values, marker="o", color="#1f4e79", linewidth=2)
        if show_target and target is not None:
            axis.axhline(target, color="#a61c00", linestyle="--", linewidth=1.5, label="Target")
            axis.legend(frameon=False)
        axis.set_xlabel("First-child floor jump")
        axis.set_ylabel(label)
        axis.grid(alpha=0.2)
    fig.savefig(out / "jump_frontier.png", dpi=180)
    fig.savefig(out / "jump_frontier.pdf")
    plt.close(fig)

    best = min(frontier, key=lambda row: float(row["transition_loss"]))
    summary = {
        "schema": "e5f_first_child_jump_diagnostic_v1",
        "status": "complete_fixed_external_jump_diagnostic",
        "interpretation": (
            "Mechanism screen only. The jump is fixed externally in each row and "
            "is not an estimated eleventh parameter."
        ),
        "jump_contract": list(JUMPS),
        "expected_baseline_loss": float(args.expected_baseline_loss),
        "best_diagnostic_row": best,
        "input_sha256": input_hashes,
        "outputs": {
            name: digest(out / name)
            for name in ("jump_frontier.csv", "target_fit_long.csv", "jump_frontier.png", "jump_frontier.pdf")
        },
    }
    (out / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (out / "README.md").write_text(
        "# First-child housing-jump diagnostic\n\n"
        "Each row holds the jump fixed, re-solves the old steady state, and reruns "
        "the complete 2007--2023 transition from the certified center. The jump is "
        "not yet a calibrated parameter. Zero must reproduce the certified baseline.\n\n"
        f"Best fixed diagnostic jump: {float(best['hbar_first_child_jump']):.2f}; "
        f"loss: {float(best['transition_loss']):.6f}.\n",
        encoding="utf-8",
    )
    print(
        "E5F_FIRST_CHILD_JUMP_DIAGNOSTIC_COMPLETE "
        f"best_jump={float(best['hbar_first_child_jump']):.2f} "
        f"best_loss={float(best['transition_loss']):.9f}"
    )


if __name__ == "__main__":
    main()
