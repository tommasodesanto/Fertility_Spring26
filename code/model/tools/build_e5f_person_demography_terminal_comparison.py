#!/usr/bin/env python3
"""Build a verified comparison of the corrected 1% and 2% terminal roots."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_BASELINE = (
    ROOT
    / "output/model/e5f_person_demography_terminal_root_rebated-tax1-baseline_20260826a_production/summary.json"
)
DEFAULT_REFORM = (
    ROOT
    / "output/model/e5f_person_demography_terminal_root_rebated-tax2-reform_20260826a_production/summary.json"
)
DEFAULT_OUTPUT = (
    ROOT / "output/model/e5f_person_demography_terminal_comparison_20260826a"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-summary", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--reform-summary", type=Path, default=DEFAULT_REFORM)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path: Path, expected_case: str) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "complete_corrected_terminal_root":
        raise RuntimeError(f"Incomplete terminal root: {path}")
    if payload.get("case") != expected_case:
        raise RuntimeError(f"Wrong terminal case in {path}")
    if not all((payload.get("gates") or {}).values()):
        raise RuntimeError(f"Terminal root gates failed: {path}")
    return payload


def main() -> None:
    args = parse_args()
    baseline_path = args.baseline_summary.resolve()
    reform_path = args.reform_summary.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    baseline = load(baseline_path, "rebated-tax1-baseline")
    reform = load(reform_path, "rebated-tax2-reform")
    comparable_hashes = ("person_demography", "funded_policy", "perfect_foresight", "solver")
    for name in comparable_hashes:
        if baseline["source_hashes"].get(name) != reform["source_hashes"].get(name):
            raise RuntimeError(f"Terminal roots use different {name} hashes")
    if baseline.get("terminal_demographic_closure") != reform.get(
        "terminal_demographic_closure"
    ):
        raise RuntimeError("Terminal roots use different demographic closures")
    if baseline.get("supply_contract") != reform.get("supply_contract"):
        raise RuntimeError("Terminal roots use different housing-supply contracts")

    metrics = [
        ("asset_price", "Asset price"),
        ("renter_price", "Rent-equivalent price"),
        ("equal_transfer_period_units", "Equal four-year rebate"),
        ("resident_persons_actual_scale", "Resident persons"),
        ("household_heads_model_units", "Household heads"),
        ("annual_births_per_head", "Annual births per head"),
        ("renewal_ratio", "Demographic renewal ratio"),
    ]
    rows: list[dict[str, Any]] = []
    for key, label in metrics:
        left = float(baseline[key])
        right = float(reform[key])
        rows.append(
            {
                "metric": key,
                "label": label,
                "rebated_tax_1_percent": left,
                "rebated_tax_2_percent": right,
                "reform_minus_baseline": right - left,
                "percent_change": 100.0 * (right / left - 1.0),
            }
        )
    with (output_dir / "terminal_comparison.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    labels = ["1% rebated", "2% rebated"]
    colors = ["#26547c", "#d1495b"]
    figure, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    axes[0, 0].bar(labels, [baseline["asset_price"], reform["asset_price"]], color=colors)
    axes[0, 0].set_title("Terminal asset price")
    axes[0, 1].bar(labels, [baseline["renter_price"], reform["renter_price"]], color=colors)
    axes[0, 1].set_title("Terminal rent-equivalent price")
    axes[1, 0].bar(
        labels,
        [baseline["resident_persons_actual_scale"] / 1e6, reform["resident_persons_actual_scale"] / 1e6],
        color=colors,
    )
    axes[1, 0].set(title="Terminal resident persons", ylabel="Millions")
    axes[1, 1].bar(
        labels,
        [baseline["annual_births_per_head"], reform["annual_births_per_head"]],
        color=colors,
    )
    axes[1, 1].set_title("Annual births per household head")
    for axis in axes.flat:
        axis.grid(axis="y", alpha=0.2)
    figure.savefig(output_dir / "terminal_comparison.png", dpi=200)
    figure.savefig(output_dir / "terminal_comparison.pdf")
    plt.close(figure)

    row_by_metric = {row["metric"]: row for row in rows}
    summary = {
        "status": "complete_verified_terminal_comparison",
        "promotion_status": "not_promoted",
        "interpretation": (
            "This is a comparison of conditional terminal equilibria under the "
            "same fixed terminal survival, migration, birth-sex-share, headship, "
            "housing-supply, and psi_child closures. It is not yet a filtered "
            "demographic forecast or a solved transition path."
        ),
        "baseline_summary": str(baseline_path),
        "baseline_summary_sha256": sha256(baseline_path),
        "reform_summary": str(reform_path),
        "reform_summary_sha256": sha256(reform_path),
        "resident_person_difference": row_by_metric[
            "resident_persons_actual_scale"
        ]["reform_minus_baseline"],
        "resident_person_percent_change": row_by_metric[
            "resident_persons_actual_scale"
        ]["percent_change"],
        "asset_price_percent_change": row_by_metric["asset_price"][
            "percent_change"
        ],
        "renter_price_percent_change": row_by_metric["renter_price"][
            "percent_change"
        ],
        "annual_births_per_head_percent_change": row_by_metric[
            "annual_births_per_head"
        ]["percent_change"],
        "terminal_demographic_closure": baseline["terminal_demographic_closure"],
        "supply_contract": baseline["supply_contract"],
        "shared_source_hashes": {
            name: baseline["source_hashes"][name] for name in comparable_hashes
        },
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
