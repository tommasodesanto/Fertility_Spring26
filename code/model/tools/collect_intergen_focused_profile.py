#!/usr/bin/env python3
"""Collect the focused relaxed h_bar_0 x chi profile and polish branches."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


def finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def strict(record: dict[str, Any]) -> bool:
    return bool(
        record.get("status") == "ok"
        and record.get("strict_converged") is True
        and finite(record.get("rank_loss"))
        and finite(record.get("market_residual"))
        and float(record["market_residual"]) <= 1e-4
    )


def records(path: Path) -> list[dict[str, Any]]:
    rows = []
    if not path.exists():
        return rows
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            payload = json.loads(line)
            if isinstance(payload, dict):
                payload["_line"] = line_number
                rows.append(payload)
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    root = args.results_root.resolve()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cells = []
    profile_root = root / "profile"
    for task_dir in sorted(profile_root.glob("task_*"), key=lambda p: int(p.name.split("_")[-1])):
        metadata_path = task_dir / "metadata.json"
        if not metadata_path.exists():
            continue
        metadata = json.loads(metadata_path.read_text())
        fixed = dict(metadata.get("fixed_theta", {}))
        all_rows = records(task_dir / "cases.jsonl")
        good = sorted((r for r in all_rows if strict(r)), key=lambda r: float(r["rank_loss"]))
        best = good[0] if good else None
        cells.append(
            {
                "task": task_dir.name,
                "h_bar_0": fixed.get("h_bar_0"),
                "chi": fixed.get("chi"),
                "cases": len(all_rows),
                "strict_cases": len(good),
                "infeasible_cases": sum(r.get("status") == "infeasible_theta" for r in all_rows),
                "best_loss": None if best is None else float(best["rank_loss"]),
                "best_residual": None if best is None else float(best["market_residual"]),
                "best_line": None if best is None else int(best["_line"]),
                "best_theta_json": None if best is None else json.dumps(best.get("theta", {}), sort_keys=True),
            }
        )
    cells.sort(key=lambda r: (float(r["h_bar_0"]), float(r["chi"])))
    write_csv(args.outdir / "profile_cells.csv", cells)

    h_values = sorted({float(r["h_bar_0"]) for r in cells})
    chi_values = sorted({float(r["chi"]) for r in cells})
    lookup = {(float(r["h_bar_0"]), float(r["chi"])): r for r in cells}
    matrix = ["# Profile Loss Matrix", "", "Each cell is the best strict repaired-objective loss after reoptimizing all other active coordinates.", ""]
    matrix.append("| $h_0$ / $\\chi$ | " + " | ".join(f"{x:.2f}" for x in chi_values) + " |")
    matrix.append("| --- | " + " | ".join("---:" for _ in chi_values) + " |")
    for h0 in h_values:
        values = []
        for chi in chi_values:
            row = lookup.get((h0, chi))
            value = None if row is None else row["best_loss"]
            values.append("--" if value is None else f"{float(value):.6f}")
        matrix.append(f"| {h0:.2f} | " + " | ".join(values) + " |")
    (args.outdir / "profile_loss_matrix.md").write_text("\n".join(matrix) + "\n")

    branches = []
    for branch_dir in sorted((root / "polish").glob("*")):
        if not branch_dir.is_dir():
            continue
        all_rows = []
        for path in branch_dir.rglob("cases.jsonl"):
            all_rows.extend(records(path))
        good = sorted((r for r in all_rows if strict(r)), key=lambda r: float(r["rank_loss"]))
        best = good[0] if good else None
        branches.append(
            {
                "branch": branch_dir.name,
                "cases": len(all_rows),
                "strict_cases": len(good),
                "infeasible_cases": sum(r.get("status") == "infeasible_theta" for r in all_rows),
                "best_loss": None if best is None else float(best["rank_loss"]),
                "best_residual": None if best is None else float(best["market_residual"]),
                "best_theta_json": None if best is None else json.dumps(best.get("theta", {}), sort_keys=True),
            }
        )
    write_csv(args.outdir / "polish_branches.csv", branches)

    finite_cells = [r for r in cells if finite(r.get("best_loss"))]
    finite_branches = [r for r in branches if finite(r.get("best_loss"))]
    best_cell = min(finite_cells, key=lambda r: float(r["best_loss"])) if finite_cells else None
    best_branch = min(finite_branches, key=lambda r: float(r["best_loss"])) if finite_branches else None
    summary = {
        "results_root": str(root),
        "profile_cell_count": len(cells),
        "profile_complete": len(cells) == 25,
        "best_profile_cell": best_cell,
        "polish_branch_count": len(branches),
        "best_polish_branch": best_branch,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    readme = ["# Focused relaxed calibration report", ""]
    readme.append(f"- Profile cells collected: **{len(cells)}/25**")
    if best_cell is not None:
        readme.append(
            f"- Best profile cell: loss **{float(best_cell['best_loss']):.8g}** at "
            f"$h_0={float(best_cell['h_bar_0']):.2f}$ and $\\chi={float(best_cell['chi']):.2f}$."
        )
    if best_branch is not None:
        readme.append(
            f"- Best free-polish branch: `{best_branch['branch']}`, loss "
            f"**{float(best_branch['best_loss']):.8g}**."
        )
    readme.extend(
        [
            "",
            "Artifacts: `profile_loss_matrix.md`, `profile_cells.csv`, "
            "`polish_branches.csv`, and `summary.json`.",
        ]
    )
    (args.outdir / "README.md").write_text("\n".join(readme) + "\n")


if __name__ == "__main__":
    main()
