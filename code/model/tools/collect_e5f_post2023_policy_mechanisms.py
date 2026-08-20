#!/usr/bin/env python3
"""Validate and collect the five post-2023 housing-policy mechanism cases."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import run_e5f_post2023_policy_mechanisms as policy


CASE_ORDER = tuple(policy.POLICIES)
COMPARE_FIELDS = (
    "asset_price",
    "housing_user_cost",
    "owner_rate",
    "dependent_child_owner_rate",
    "topcode_adjusted_births_per_adult",
    "population_index_2023",
    "housing_demand_per_adult",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"Expected an object in {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def validate_case(case_dir: Path, expected_case: str) -> tuple[dict[str, Any], list[dict[str, str]]]:
    summary_path = case_dir / "summary.json"
    manifest_path = case_dir / "manifest.json"
    path_path = case_dir / "policy_path.csv"
    for required in (summary_path, manifest_path, path_path):
        if not required.is_file() or required.stat().st_size == 0:
            raise RuntimeError(f"Missing case artifact: {required}")
    summary = read_json(summary_path)
    if summary.get("status") != "complete_post2023_policy_mechanism_case":
        raise RuntimeError(f"{expected_case} is incomplete")
    spec = summary.get("policy") or {}
    if spec.get("case") != expected_case:
        raise RuntimeError(f"Case mismatch at {case_dir}: {spec.get('case')}")
    expected = policy.POLICIES[expected_case]
    if not math.isclose(
        float(spec.get("supply_intercept_multiplier", math.nan)),
        expected.supply_multiplier,
        rel_tol=0.0,
        abs_tol=1e-15,
    ):
        raise RuntimeError(f"{expected_case} supply policy mismatch")
    if bool(spec.get("dependent_child_ltv95")) != expected.dependent_child_ltv95:
        raise RuntimeError(f"{expected_case} credit policy mismatch")
    if not math.isclose(
        float(spec.get("annual_property_tax_rate", math.nan)),
        expected.annual_property_tax_rate,
        rel_tol=0.0,
        abs_tol=1e-15,
    ):
        raise RuntimeError(f"{expected_case} property-tax policy mismatch")
    if spec.get("fiscal_treatment") != (
        "property-tax revenue is discarded; no rebate and no purchase grant"
    ):
        raise RuntimeError(f"{expected_case} fiscal treatment mismatch")
    manifest = read_json(manifest_path)
    if manifest.get("schema") != "e5f_post2023_policy_mechanism_manifest_v1":
        raise RuntimeError(f"{expected_case} manifest schema mismatch")
    for name, expected_hash in (manifest.get("artifacts") or {}).items():
        artifact = case_dir / str(name)
        if not artifact.is_file() or sha256_file(artifact) != expected_hash:
            raise RuntimeError(f"{expected_case} artifact hash mismatch: {name}")
    rows = read_csv(path_path)
    expected_dates = int(summary["post_2023_periods"]) + 1
    if len(rows) != expected_dates:
        raise RuntimeError(f"{expected_case} path length mismatch")
    years = [int(row["calendar_year"]) for row in rows]
    if years != list(range(2023, 2023 + 4 * expected_dates, 4)):
        raise RuntimeError(f"{expected_case} calendar mismatch: {years}")
    for row in rows:
        if row["policy_case"] != expected_case:
            raise RuntimeError(f"{expected_case} row label mismatch")
        if not math.isclose(
            float(row["annual_property_tax_rate"]),
            expected.annual_property_tax_rate,
            rel_tol=0.0,
            abs_tol=1e-15,
        ):
            raise RuntimeError(f"{expected_case} property tax is not active")
        if float(row["property_tax_lump_sum_transfer"]) != 0.0:
            raise RuntimeError(f"{expected_case} unexpectedly rebates tax revenue")
        if row["purchase_grant"].strip().lower() not in ("false", "0"):
            raise RuntimeError(f"{expected_case} unexpectedly activates a grant")
        for field in COMPARE_FIELDS + (
            "relative_market_residual",
            "mass_accounting_residual",
            "feasibility_frontier_projection_mass",
        ):
            if not math.isfinite(float(row[field])):
                raise RuntimeError(f"{expected_case} nonfinite {field}")
        if abs(float(row["relative_market_residual"])) > 2e-4:
            raise RuntimeError(f"{expected_case} market residual failed")
        if abs(float(row["mass_accounting_residual"])) > 2e-8:
            raise RuntimeError(f"{expected_case} mass residual failed")
        if float(row["feasibility_frontier_projection_mass"]) > 1e-6:
            raise RuntimeError(f"{expected_case} projection mass failed")
    return summary, rows


def build_effect_rows(paths: dict[str, list[dict[str, str]]]) -> list[dict[str, Any]]:
    baseline_by_year = {
        int(row["calendar_year"]): row for row in paths["baseline"]
    }
    rows: list[dict[str, Any]] = []
    for case in CASE_ORDER:
        for row in paths[case]:
            year = int(row["calendar_year"])
            base = baseline_by_year[year]
            item: dict[str, Any] = {
                "policy_case": case,
                "policy_label": policy.POLICIES[case].label,
                "calendar_year": year,
            }
            for field in COMPARE_FIELDS:
                value = float(row[field])
                base_value = float(base[field])
                item[field] = value
                item[f"baseline_{field}"] = base_value
                item[f"difference_{field}"] = value - base_value
                item[f"percent_difference_{field}"] = (
                    100.0 * (value / base_value - 1.0)
                    if abs(base_value) > 1e-15
                    else math.nan
                )
            rows.append(item)
    return rows


def make_figure(paths: dict[str, list[dict[str, str]]], outdir: Path) -> None:
    panels = (
        ("housing_user_cost", "Housing user cost"),
        ("dependent_child_owner_rate", "Ownership: dependent-child households"),
        ("topcode_adjusted_births_per_adult", "Adjusted births per adult"),
        ("population_index_2023", "Adult households (2023=1)"),
    )
    colors = {
        "baseline": "#222222",
        "supply-plus-20": "#2c7fb8",
        "dependent-child-ltv95": "#d95f0e",
        "combined": "#31a354",
        "property-tax-2pct-no-rebate": "#756bb1",
    }
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.0))
    for axis, (field, title) in zip(axes.flat, panels):
        for case in CASE_ORDER:
            rows = paths[case]
            base_2023 = float(paths["baseline"][0][field])
            values = [float(row[field]) for row in rows]
            if field == "housing_user_cost":
                values = [value / base_2023 for value in values]
                title_used = "Housing user cost (baseline 2023=1)"
            else:
                title_used = title
            axis.plot(
                [int(row["calendar_year"]) for row in rows],
                values,
                label=policy.POLICIES[case].label,
                color=colors[case],
                linewidth=2.0,
            )
        axis.set_title(title_used)
        axis.set_xlabel("Year")
        axis.grid(alpha=0.2)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0.0, 0.12, 1.0, 1.0))
    fig.savefig(outdir / "policy_mechanism_comparison.png", dpi=220)
    fig.savefig(outdir / "policy_mechanism_comparison.pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-post-2023-periods", type=int, required=True)
    args = parser.parse_args()
    case_root = args.case_root.resolve()
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    summaries: dict[str, dict[str, Any]] = {}
    paths: dict[str, list[dict[str, str]]] = {}
    for case in CASE_ORDER:
        summary, rows = validate_case(case_root / case, case)
        if int(summary["post_2023_periods"]) != int(args.expected_post_2023_periods):
            raise RuntimeError(f"{case} has the wrong horizon")
        summaries[case] = summary
        paths[case] = rows

    reference = summaries["baseline"]
    identity_fields = (
        ("reconstruction", "matched_2023_distribution_sha256"),
        ("reconstruction", "matched_2023_queue_sha256"),
    )
    for case in CASE_ORDER[1:]:
        if summaries[case]["input_hashes"] != reference["input_hashes"]:
            raise RuntimeError(f"{case} input contract differs from baseline")
        for parent, field in identity_fields:
            if summaries[case][parent][field] != reference[parent][field]:
                raise RuntimeError(f"{case} does not start from the same 2023 {field}")
    baseline_population = {
        int(row["calendar_year"]): float(row["population_index_2023"])
        for row in paths["baseline"]
    }
    for case in CASE_ORDER[1:]:
        for row in paths[case]:
            year = int(row["calendar_year"])
            if year >= 2043:
                continue
            gap = abs(
                float(row["population_index_2023"]) - baseline_population[year]
            )
            if gap > 2e-10:
                raise RuntimeError(
                    f"{case} changes adult population before the 20-year entry lag: "
                    f"year={year}, gap={gap:.3e}"
                )

    combined_rows: list[dict[str, Any]] = []
    for case in CASE_ORDER:
        combined_rows.extend(paths[case])
    effects = build_effect_rows(paths)
    write_csv(outdir / "policy_paths.csv", combined_rows)
    write_csv(outdir / "policy_effects_relative_to_baseline.csv", effects)
    make_figure(paths, outdir)

    terminal_year = 2023 + 4 * int(args.expected_post_2023_periods)
    summary = {
        "status": "complete_post2023_policy_mechanism_comparison",
        "scope": (
            "closed-national finite-horizon mechanism comparison from one common "
            "fitted 2023 state; property-tax revenue is discarded; no welfare or "
            "terminal-steady-state claim"
        ),
        "cases": list(CASE_ORDER),
        "post_2023_periods": int(args.expected_post_2023_periods),
        "terminal_year": terminal_year,
        "common_2023_distribution_sha256": reference["reconstruction"][
            "matched_2023_distribution_sha256"
        ],
        "common_2023_queue_sha256": reference["reconstruction"][
            "matched_2023_queue_sha256"
        ],
        "input_hashes": reference["input_hashes"],
        "case_summary_sha256": {
            case: sha256_file(case_root / case / "summary.json") for case in CASE_ORDER
        },
        "selected_effects": {
            case: {
                str(year): next(
                    row
                    for row in effects
                    if row["policy_case"] == case and int(row["calendar_year"]) == year
                )
                for year in sorted({2023, terminal_year} | ({2043} if terminal_year >= 2043 else set()))
            }
            for case in CASE_ORDER
        },
    }
    write_json(outdir / "summary.json", summary)
    readme = [
        "# Housing-policy mechanism comparison",
        "",
        "All five cases begin from the same fitted 2023 pre-fertility household distribution and inherited birth queue. The population closure is the closed national benchmark (`M=0`, `rho=1`).",
        "",
        "The policies are a permanent 20-percent shift in the housing-supply intercept, 95-percent financing for owned housing choices by households with dependent children, their combination, and an annual property-tax increase from 1 to 2 percent. The credit policy is not restricted to first-time buyers.",
        "",
        "This is a finite-horizon mechanism experiment. Property-tax revenue is discarded in every case: there is no rebate, purchase grant, or balanced-budget closure. It contains no welfare calculation, immigration response, or claim that the terminal year is a steady state. Because births enter the adult-household population after twenty years, population differences should not be attributed to the policy before 2043.",
        "",
        "See `policy_mechanism_comparison.pdf` and `policy_effects_relative_to_baseline.csv`.",
    ]
    (outdir / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    artifacts = (
        "policy_paths.csv",
        "policy_effects_relative_to_baseline.csv",
        "policy_mechanism_comparison.png",
        "policy_mechanism_comparison.pdf",
        "summary.json",
        "README.md",
    )
    write_json(
        outdir / "manifest.json",
        {
            "schema": "e5f_post2023_policy_mechanism_comparison_manifest_v1",
            "artifacts": {name: sha256_file(outdir / name) for name in artifacts},
        },
    )
    print(
        f"POLICY_COMPARISON_COMPLETE cases={len(CASE_ORDER)} "
        f"terminal_year={terminal_year}",
        flush=True,
    )


if __name__ == "__main__":
    main()
