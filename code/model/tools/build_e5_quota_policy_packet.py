#!/usr/bin/env python3
"""Validate raw E5 policy runs and build the quota-closure advisor packet."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTDIR = ROOT / "output/model/eqscale_seq_e5_policy_quota_closure_20260807"
OUTSIDE_SHARE = 0.169
MULTIPLIER = (1.0 - OUTSIDE_SHARE) / OUTSIDE_SHARE
PERIOD_YEARS = 4.0
ANNUAL_INTEREST_RATE = 0.04
ANNUAL_DEPRECIATION_RATE = 0.02
BASELINE_CASE = "rebated_tax1_baseline"
POLICIES = (
    ("rebated_tax2", "Rebated 2% tax"),
    ("rebated_tax2_grant0p4_Hge6", "Rebated 2% tax + purchase grant"),
)
ARMS = {
    "floor": {
        "label": "Floor, psi nonnegative winner",
        "calibration": ROOT
        / "output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806",
        "old_policy": ROOT
        / "output/model/intergen_e5f_child_room_floor_psinneg_policy_empirical_entry_20260806",
    },
    "tilt": {
        "label": "Tilt, psi nonnegative winner",
        "calibration": ROOT
        / "output/model/eqscale_seq_e5_maturation_repair_psinneg_extended_20260806",
        "old_policy": ROOT
        / "output/model/eqscale_seq_e5_maturation_repair_psinneg_policy_empirical_entry_20260806",
    },
    "chain6": {
        "label": "Maturation-repair chain-6 winner",
        "calibration": ROOT
        / "output/model/eqscale_seq_e5_maturation_repair_recalibration_20260805",
        "old_policy": ROOT
        / "output/model/eqscale_seq_e5_maturation_repair_policy_entry_20260806",
    },
}
POLICY_ARMS = ("floor", "tilt")
SCREENED_ARM = "chain6"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--raw-root", type=Path)
    parser.add_argument("--production-job-id", action="append", default=[])
    parser.add_argument("--smoke-job-id", action="append", default=[])
    parser.add_argument(
        "--reporting-only",
        action="store_true",
        help="Regenerate only README.md and tables/ from the existing policy summary.",
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def keyed(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["case"]: row for row in rows}


def number(row: dict[str, str], name: str) -> float:
    return float(row[name])


def assert_close(label: str, actual: float, expected: float, tolerance: float) -> None:
    if not math.isfinite(actual) or abs(actual - expected) > tolerance:
        raise RuntimeError(
            f"{label} failed: actual={actual:.16g}, expected={expected:.16g}, "
            f"tolerance={tolerance:.3g}"
        )


def renter_user_cost_rate(annual_property_tax_rate: float) -> float:
    """Per-period renter user cost under the active one-market mapping."""
    q = (1.0 + ANNUAL_INTEREST_RATE) ** PERIOD_YEARS - 1.0
    delta = 1.0 - (1.0 - ANNUAL_DEPRECIATION_RATE) ** PERIOD_YEARS
    return q + delta + float(annual_property_tax_rate) * PERIOD_YEARS


def copy_calibration_artifacts(outdir: Path) -> None:
    parameter_rows: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    for arm, contract in ARMS.items():
        calibration = Path(contract["calibration"])
        report = calibration / "report"
        arm_calibration = outdir / "calibration" / arm
        arm_calibration.mkdir(parents=True, exist_ok=True)
        for name in ("parameter_table_full.csv", "target_fit_full.csv"):
            shutil.copy2(report / name, arm_calibration / name)
        for row in read_csv(report / "parameter_table_full.csv"):
            parameter_rows.append({"arm": arm, **row})
        for row in read_csv(report / "target_fit_full.csv"):
            target_rows.append({"arm": arm, **row})
        shutil.copytree(
            calibration / "diagnostic_packet/standard",
            outdir / "diagnostics" / arm / "standard",
            dirs_exist_ok=True,
        )
    write_csv(outdir / "calibration/parameter_table_all_arms.csv", parameter_rows)
    write_csv(outdir / "calibration/target_fit_all_arms.csv", target_rows)


def build_plot(outdir: Path, rows: list[dict[str, Any]]) -> None:
    arms = list(POLICY_ARMS)
    colors = {"floor": "#17365D", "tilt": "#C45A3B", "chain6": "#2B6F68"}
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.8))
    positions = np.arange(len(POLICIES), dtype=float)
    width = 0.24
    for arm_index, arm in enumerate(arms):
        arm_rows = {row["case"]: row for row in rows if row["arm"] == arm}
        offset = (arm_index - 0.5) * width
        axes[0].bar(
            positions + offset,
            [arm_rows[case]["quota_population_change_percent"] for case, _ in POLICIES],
            width,
            color=colors[arm],
            label=arm,
        )
        axes[1].bar(
            positions + offset,
            [arm_rows[case]["tfr_change_percent"] for case, _ in POLICIES],
            width,
            color=colors[arm],
        )
        axes[2].bar(
            positions + offset,
            [arm_rows[case]["house_price_change_percent"] for case, _ in POLICIES],
            width,
            color=colors[arm],
        )
    short_labels = ["2% tax", "2% tax + grant"]
    for axis, title, ylabel in (
        (axes[0], "A. Household population", "Percent change"),
        (axes[1], "B. Completed fertility", "Percent change"),
        (axes[2], "C. House price and rent", "Percent change"),
    ):
        axis.axhline(0.0, color="#777777", linewidth=0.8)
        axis.set_xticks(positions, short_labels)
        axis.set_title(title)
        axis.set_ylabel(ylabel)
        axis.grid(axis="y", alpha=0.2)
    axes[0].legend(frameon=False)
    fig.suptitle("Experimental repaired-E5 funded policies: quota population closure")
    fig.tight_layout()
    figure_dir = outdir / "diagnostics/policy"
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / "quota_policy_comparison.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def build_readme(
    outdir: Path,
    wide_rows: list[dict[str, Any]],
    acceptance: dict[str, Any],
    production_job_ids: list[str],
    smoke_job_ids: list[str],
) -> None:
    lines = [
        "# Experimental repaired-E5 policy packet: quota population closure",
        "",
        "## Bottom line",
        "",
        (
            "This packet is an **experimental E5-repair arm**, not a promoted model. "
            "It reports long-run stationary comparisons, not forecasts. The 0.169 "
            "outside-origin entrant-share anchor is provisional and still needs a "
            "national ACS re-anchor."
        ),
        "",
        (
            "The change is deliberately small: the funded baseline now identifies a "
            "fixed retention rate and a fixed outside inflow. Policy population moves "
            "only because the policy changes mature city-born household formation. "
            "The unidentified logit entry response is bypassed and retained only as "
            "the final labeled sensitivity row."
        ),
        "",
        (
            "That is why the population effects fall from the old roughly 9–10 percent "
            "range to roughly 1–3 percent: the former was driven mostly by the external "
            "taste scale; the latter is the transparent renewal multiplier "
            f"(1−0.169)/0.169 = {MULTIPLIER:.3f}, modestly damped by housing feedback."
        ),
        "",
        (
            "The maturation-repair chain-6 winner is shown in `candidate_screening.csv` "
            "but has no corrected policy rows: the same 0.169 anchor implies "
            "Rbar = 1.0550, which violates the explicit Rbar <= 1 feasibility guard. "
            "It is excluded rather than assigned a fallback."
        ),
        "",
        "## Main results",
        "",
        "| Arm | Policy | ΔTFR (%) | Births/HH (%) | Household population (%) | Required immigration at fixed population (identical to percent of fertility gap closed) (%) | Price (%) | Rent (%) | Damping decomposition: arithmetic → realized (damping, pp) | Unidentified entry-response sensitivity (lambda=2) (%) |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in wide_rows:
        lines.append(
            "| {arm} | {policy} | {tfr:.3f} | {births:.3f} | {population:.3f} | "
            "{gap:.3f} | {price:.3f} | {rent:.3f} | {arithmetic:.3f} → "
            "{realized:.3f} ({damping:.3f}) | {logit:.3f} |".format(
                arm=row["arm"],
                policy=row["policy_label"],
                tfr=row["tfr_change_percent"],
                births=row["births_per_household_change_percent"],
                population=row["quota_population_change_percent"],
                gap=row["fertility_gap_closed_percent"],
                price=row["house_price_change_percent"],
                rent=row["rent_change_percent"],
                arithmetic=row["arithmetic_multiplier_times_dlnB_percent"],
                realized=row["realized_dln_population_percent"],
                damping=row["housing_feedback_damping_percentage_points"],
                logit=row["old_logit_population_change_percent"],
            )
        )
    lines.extend(
        [
            "",
            "The merged required-immigration figure is reported as the reduction in required immigration, so it is numerically identical to the percent of the fertility gap closed. Aggregate stationary births = population + births/HH (accounting identity); see `policy_summary.csv`. The difference between the arithmetic renewal response and the realized population response is the housing-feedback damping reported in the detailed tables.",
            "",
            "## Closure and interpretation",
            "",
            r"At baseline, \(\bar R=(1-s^{out})E_0/B_0\) and \(\bar M=s^{out}E_0\). In each counterfactual, \((\bar R,\bar M)\) are fixed and the driver solves \(S E_0=\bar R S B_0(\text{policy})+\bar M\) jointly with the house price and balanced-budget transfer. Households are the population unit. Renter-paid rent is \((q+\delta+\tau_H)p\), so it differs from the price change when the property-tax rate changes.",
            "",
            "The baseline retention rate is a derived accounting identity, not a directly estimated behavioral elasticity. The experiment assumes quota-rationed outside inflow and policy-invariant retention. A balanced-growth-path closure remains an explicit design question to reconsider after this meeting; it is not implemented here.",
            "",
            "## Acceptance tests and residuals",
            "",
            f"All {acceptance['test_count']} encoded checks passed. Quota baselines reproduce S=1 and the 0.169 share to {acceptance['baseline_tolerance']:.0e}; feasible-arm fixed-population CSVs are byte-identical across modes and to the prior packets, while the two chain-6 screening CSVs are byte-identical to each other; quota objects contain no outside value or taste-scale parameter; and all policy market/fiscal residuals are below {acceptance['residual_tolerance']:.1e}.",
            "",
            "- `policy_summary.csv`: one wide row per arm and policy.",
            "- `policy_summary_long.csv` and `tables/`: requested advisor-facing metric rows; the unidentified logit sensitivity is last.",
            "- `residuals.csv`: market and fiscal residuals for every baseline and policy under each closure, plus fixed-population decompositions.",
            "- `calibration/`: complete free/fixed parameter and target-fit tables for all three unchanged certified theta vectors.",
            "- `diagnostics/`: the unchanged standard solution plots for each certified arm plus the policy comparison plot.",
            "- `raw_runs/`: deterministic driver outputs and checkpoints from Torch.",
            "",
            "## Run record",
            "",
            f"- Local loose-tolerance smoke: {', '.join(smoke_job_ids) if smoke_job_ids else 'completed before Torch submission (local diagnostic)'}.",
            f"- Torch production array(s): {', '.join(production_job_ids) if production_job_ids else 'see metadata.json'}.",
            "- No calibration, target, policy, solver-numerics, ACS, half-life, CE-unit, promoted-model, or paper-LaTeX object was changed.",
            "",
        ]
    )
    (outdir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def rebuild_reporting_layout(outdir: Path) -> None:
    """Regenerate the advisor-facing layout without modifying model artifacts."""
    wide_rows: list[dict[str, Any]] = []
    for raw_row in read_csv(outdir / "policy_summary.csv"):
        row: dict[str, Any] = dict(raw_row)
        for key in (
            "tfr_change_percent",
            "births_per_household_change_percent",
            "quota_population_change_percent",
            "fertility_gap_closed_percent",
            "house_price_change_percent",
            "rent_change_percent",
            "arithmetic_multiplier_times_dlnB_percent",
            "realized_dln_population_percent",
            "housing_feedback_damping_percentage_points",
            "old_logit_population_change_percent",
        ):
            row[key] = float(raw_row[key])
        quota_rows = keyed(read_csv(outdir / "raw_runs" / row["arm"] / "quota" / "summary.csv"))
        baseline = quota_rows[BASELINE_CASE]
        policy = quota_rows[row["case"]]
        baseline_rent = renter_user_cost_rate(number(baseline, "annual_property_tax_rate")) * number(
            baseline, "price"
        )
        policy_rent = renter_user_cost_rate(number(policy, "annual_property_tax_rate")) * number(
            policy, "price"
        )
        row["rent_change_percent"] = 100.0 * (policy_rent / baseline_rent - 1.0)
        wide_rows.append(row)

    for row in wide_rows:
        metric_rows = (
            ("completed fertility change", row["tfr_change_percent"], "percent"),
            ("births per household change", row["births_per_household_change_percent"], "percent"),
            ("household population change under quota", row["quota_population_change_percent"], "percent"),
            (
                "required immigration at fixed population "
                "(identical to percent of fertility gap closed)",
                row["fertility_gap_closed_percent"],
                "percent reduction",
            ),
            ("house price change", row["house_price_change_percent"], "percent"),
            ("rent change", row["rent_change_percent"], "percent"),
            ("4.917 times fixed-population dlnB", row["arithmetic_multiplier_times_dlnB_percent"], "log-percent"),
            ("realized dln household population", row["realized_dln_population_percent"], "log-percent"),
            ("housing-feedback damping", row["housing_feedback_damping_percentage_points"], "percentage points"),
            (
                "unidentified entry-response sensitivity (lambda=2)",
                row["old_logit_population_change_percent"],
                "percent",
            ),
        )
        table_rows = [
            {
                "arm": row["arm"],
                "case": row["case"],
                "policy": row["policy_label"],
                "row_order": order,
                "metric": metric,
                "value": value,
                "unit": unit,
            }
            for order, (metric, value, unit) in enumerate(metric_rows, start=1)
        ]
        write_csv(outdir / "tables" / f"{row['arm']}_{row['case']}.csv", table_rows)

    acceptance = json.loads((outdir / "acceptance_tests.json").read_text(encoding="utf-8"))
    metadata = json.loads((outdir / "metadata.json").read_text(encoding="utf-8"))
    build_readme(
        outdir,
        wide_rows,
        acceptance,
        list(metadata.get("torch_production_job_ids", [])),
        list(metadata.get("torch_smoke_job_ids", [])),
    )


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    raw_root = (args.raw_root or (outdir / "raw_runs")).resolve()
    if args.reporting_only:
        rebuild_reporting_layout(outdir)
        return
    outdir.mkdir(parents=True, exist_ok=True)
    wide_rows: list[dict[str, Any]] = []
    long_rows: list[dict[str, Any]] = []
    residual_rows: list[dict[str, Any]] = []
    checks: list[dict[str, Any]] = []
    max_residual = 0.0
    residual_tolerance = 2.6e-5

    for arm in POLICY_ARMS:
        contract = ARMS[arm]
        quota_dir = raw_root / arm / "quota"
        logit_dir = raw_root / arm / "logit"
        quota_metadata = json.loads((quota_dir / "metadata.json").read_text(encoding="utf-8"))
        logit_metadata = json.loads((logit_dir / "metadata.json").read_text(encoding="utf-8"))
        if quota_metadata.get("status") != "complete" or quota_metadata.get("smoke"):
            raise RuntimeError(f"{arm}: quota run is not a completed production run")
        if logit_metadata.get("status") != "complete" or logit_metadata.get("smoke"):
            raise RuntimeError(f"{arm}: logit run is not a completed production run")
        if quota_metadata.get("closure_mode") != "quota":
            raise RuntimeError(f"{arm}: quota metadata has wrong closure mode")
        if logit_metadata.get("closure_mode") != "logit":
            raise RuntimeError(f"{arm}: logit metadata has wrong closure mode")

        quota_objects = json.loads(
            (quota_dir / "benchmark_quota_objects.json").read_text(encoding="utf-8")
        )
        for forbidden in ("outside_value", "entry_taste_scale"):
            if forbidden in quota_objects:
                raise RuntimeError(f"{arm}: quota closure leaked {forbidden}")
        checks.append({"arm": arm, "test": "quota_logit_objects_absent", "passed": True})

        quota_summary = keyed(read_csv(quota_dir / "summary.csv"))
        logit_summary = keyed(read_csv(logit_dir / "summary.csv"))
        fixed_quota = (quota_dir / "fixed_population_summary.csv").read_bytes()
        fixed_logit = (logit_dir / "fixed_population_summary.csv").read_bytes()
        fixed_prior = (Path(contract["old_policy"]) / "fixed_population_summary.csv").read_bytes()
        if fixed_quota != fixed_logit:
            raise RuntimeError(f"{arm}: quota/logit fixed-population rows are not byte-identical")
        if fixed_quota != fixed_prior:
            raise RuntimeError(f"{arm}: fixed-population rows differ from prior packet")
        checks.append({"arm": arm, "test": "fixed_population_bitwise", "passed": True})
        for name in ("summary.csv", "target_fit_long.csv"):
            if (logit_dir / name).read_bytes() != (
                Path(contract["old_policy"]) / name
            ).read_bytes():
                raise RuntimeError(f"{arm}: default logit {name} differs from prior packet")
        checks.append(
            {"arm": arm, "test": "default_logit_rows_bitwise", "passed": True}
        )

        baseline = quota_summary[BASELINE_CASE]
        assert_close(
            f"{arm} quota baseline scale",
            number(baseline, "population_scale"),
            1.0,
            1.0e-12,
        )
        assert_close(
            f"{arm} quota baseline outside share",
            number(baseline, "outside_origin_entrant_share"),
            OUTSIDE_SHARE,
            1.0e-12,
        )
        checks.append({"arm": arm, "test": "baseline_reproduction", "passed": True})

        fixed_demography = keyed(read_csv(quota_dir / "fixed_population_demography.csv"))
        baseline_demography = fixed_demography[BASELINE_CASE]
        entry_flow = float(quota_objects["baseline_entry_flow"])
        retention = float(quota_objects["retention_rate"])
        outside_inflow = float(quota_objects["outside_inflow"])
        baseline_mature = number(baseline_demography, "mature_cityborn_flow")

        for case, policy_label in POLICIES:
            policy = quota_summary[case]
            logit_policy = logit_summary[case]
            policy_demography = fixed_demography[case]
            mature = number(policy_demography, "mature_cityborn_flow")
            dln_mature = math.log(mature / baseline_mature)
            arithmetic = 100.0 * MULTIPLIER * dln_mature
            realized_dln_scale = 100.0 * math.log(number(policy, "population_scale"))
            damping = arithmetic - realized_dln_scale
            population_change = 100.0 * (number(policy, "population_scale") - 1.0)
            if abs(population_change) >= 6.0:
                raise RuntimeError(
                    f"{arm}/{case}: quota response resembles leaked logit response: "
                    f"{population_change:.4f}%"
                )
            if abs(arithmetic) > 1.0e-8:
                ratio = realized_dln_scale / arithmetic
                if not 0.5 <= ratio <= 1.2:
                    raise RuntimeError(
                        f"{arm}/{case}: quota renewal sanity ratio is {ratio:.4f}"
                    )
            required_inflow = entry_flow - retention * mature
            if required_inflow <= 0.0:
                raise RuntimeError(f"{arm}/{case}: fixed-population required inflow is nonpositive")
            immigration_change = 100.0 * (required_inflow / outside_inflow - 1.0)
            gap_closed = -immigration_change
            tfr_change = number(policy, "tfr") - number(baseline, "tfr")
            tfr_change_percent = 100.0 * (
                number(policy, "tfr") / number(baseline, "tfr") - 1.0
            )
            births_per_household_change = 100.0 * (
                number(policy, "normalized_births") / number(baseline, "normalized_births")
                - 1.0
            )
            total_births_change = 100.0 * (
                number(policy, "total_births") / number(baseline, "total_births") - 1.0
            )
            price_change = 100.0 * (
                number(policy, "price") / number(baseline, "price") - 1.0
            )
            rent_change = 100.0 * (
                renter_user_cost_rate(number(policy, "annual_property_tax_rate"))
                * number(policy, "price")
                / (
                    renter_user_cost_rate(number(baseline, "annual_property_tax_rate"))
                    * number(baseline, "price")
                )
                - 1.0
            )
            old_logit = number(logit_policy, "population_change_percent")
            wide = {
                "arm": arm,
                "arm_label": contract["label"],
                "case": case,
                "policy_label": policy_label,
                "tfr_change_children": tfr_change,
                "tfr_change_percent": tfr_change_percent,
                "births_per_household_change_percent": births_per_household_change,
                "stationary_total_births_change_percent": total_births_change,
                "quota_population_change_percent": population_change,
                "required_immigration_change_level": required_inflow - outside_inflow,
                "required_immigration_change_percent": immigration_change,
                "fertility_gap_closed_percent": gap_closed,
                "house_price_change_percent": price_change,
                "rent_change_percent": rent_change,
                "renewal_multiplier": MULTIPLIER,
                "fixed_population_mature_birth_flow_dln_percent": 100.0 * dln_mature,
                "arithmetic_multiplier_times_dlnB_percent": arithmetic,
                "realized_dln_population_percent": realized_dln_scale,
                "housing_feedback_damping_percentage_points": damping,
                "old_logit_population_change_percent": old_logit,
            }
            wide_rows.append(wide)
            metric_rows = (
                ("completed fertility change", tfr_change_percent, "percent"),
                ("births per household change", births_per_household_change, "percent"),
                ("household population change under quota", population_change, "percent"),
                (
                    "required immigration at fixed population "
                    "(identical to percent of fertility gap closed)",
                    gap_closed,
                    "percent reduction",
                ),
                ("house price change", price_change, "percent"),
                ("rent change", rent_change, "percent"),
                ("4.917 times fixed-population dlnB", arithmetic, "log-percent"),
                ("realized dln household population", realized_dln_scale, "log-percent"),
                ("housing-feedback damping", damping, "percentage points"),
                (
                    "unidentified entry-response sensitivity (lambda=2)",
                    old_logit,
                    "percent",
                ),
            )
            one_policy_rows = []
            for order, (metric, value, unit) in enumerate(metric_rows, start=1):
                metric_row = {
                    "arm": arm,
                    "case": case,
                    "policy": policy_label,
                    "row_order": order,
                    "metric": metric,
                    "value": value,
                    "unit": unit,
                }
                long_rows.append(metric_row)
                one_policy_rows.append(metric_row)
            write_csv(outdir / "tables" / f"{arm}_{case}.csv", one_policy_rows)
        checks.append({"arm": arm, "test": "renewal_sanity_and_no_logit_leak", "passed": True})

        for closure_name, rows in (("quota", quota_summary), ("logit", logit_summary)):
            for case, row in rows.items():
                market = abs(number(row, "market_residual"))
                fiscal = abs(number(row, "scaled_fiscal_residual"))
                max_residual = max(max_residual, market, fiscal)
                residual_rows.append(
                    {
                        "arm": arm,
                        "closure": closure_name,
                        "case": case,
                        "market_residual": number(row, "market_residual"),
                        "fiscal_residual": number(row, "scaled_fiscal_residual"),
                        "passes_2p6e_minus5": market <= residual_tolerance
                        and fiscal <= residual_tolerance,
                    }
                )
                if market > residual_tolerance or fiscal > residual_tolerance:
                    raise RuntimeError(
                        f"{arm}/{closure_name}/{case}: residual tolerance failed: "
                        f"market={market:.4g}, fiscal={fiscal:.4g}"
                    )
        fixed_rows = keyed(read_csv(quota_dir / "fixed_population_summary.csv"))
        for case, row in fixed_rows.items():
            market = abs(number(row, "market_residual"))
            fiscal = abs(number(row, "government_budget_residual"))
            max_residual = max(max_residual, market, fiscal)
            residual_rows.append(
                {
                    "arm": arm,
                    "closure": "fixed_population_decomposition",
                    "case": case,
                    "market_residual": number(row, "market_residual"),
                    "fiscal_residual": number(row, "government_budget_residual"),
                    "passes_2p6e_minus5": market <= residual_tolerance
                    and fiscal <= residual_tolerance,
                }
            )
            if market > residual_tolerance or fiscal > residual_tolerance:
                raise RuntimeError(f"{arm}/fixed/{case}: residual tolerance failed")
        checks.append({"arm": arm, "test": "market_and_fiscal_residuals", "passed": True})

    screened_contract = ARMS[SCREENED_ARM]
    screened_quota_dir = raw_root / SCREENED_ARM / "quota"
    screened_logit_dir = raw_root / SCREENED_ARM / "logit"
    screened_demography = keyed(
        read_csv(screened_quota_dir / "fixed_population_demography.csv")
    )
    screened_baseline = screened_demography[BASELINE_CASE]
    screened_entry = number(screened_baseline, "entry_flow")
    screened_mature = number(screened_baseline, "mature_cityborn_flow")
    screened_retention = (1.0 - OUTSIDE_SHARE) * screened_entry / screened_mature
    if screened_retention <= 1.0:
        raise RuntimeError(
            f"{SCREENED_ARM}: expected feasibility guard did not bind; "
            f"Rbar={screened_retention:.12g}"
        )
    screened_quota_fixed = (
        screened_quota_dir / "fixed_population_summary.csv"
    ).read_bytes()
    screened_logit_fixed = (
        screened_logit_dir / "fixed_population_summary.csv"
    ).read_bytes()
    if screened_quota_fixed != screened_logit_fixed:
        raise RuntimeError(
            f"{SCREENED_ARM}: quota/logit screening fixed rows are not byte-identical"
        )
    screening_rows = [
        {
            "arm": SCREENED_ARM,
            "arm_label": screened_contract["label"],
            "outside_origin_entrant_share": OUTSIDE_SHARE,
            "baseline_entry_flow_E0": screened_entry,
            "baseline_mature_cityborn_flow_B0": screened_mature,
            "implied_retention_rate_Rbar": screened_retention,
            "feasible_Rbar_le_1": False,
            "policy_rows_included": False,
            "reason": "explicit feasibility guard: Rbar exceeds one; no fallback",
        }
    ]
    write_csv(outdir / "candidate_screening.csv", screening_rows)
    checks.append(
        {
            "arm": SCREENED_ARM,
            "test": "infeasible_retention_fails_loudly",
            "passed": True,
            "retention_rate": screened_retention,
        }
    )
    checks.append(
        {
            "arm": SCREENED_ARM,
            "test": "screening_fixed_population_bitwise_across_modes",
            "passed": True,
        }
    )
    for case, row in keyed(
        read_csv(screened_quota_dir / "fixed_population_summary.csv")
    ).items():
        market = abs(number(row, "market_residual"))
        fiscal = abs(number(row, "government_budget_residual"))
        max_residual = max(max_residual, market, fiscal)
        residual_rows.append(
            {
                "arm": SCREENED_ARM,
                "closure": "fixed_population_screened_candidate",
                "case": case,
                "market_residual": number(row, "market_residual"),
                "fiscal_residual": number(row, "government_budget_residual"),
                "passes_2p6e_minus5": market <= residual_tolerance
                and fiscal <= residual_tolerance,
            }
        )
        if market > residual_tolerance or fiscal > residual_tolerance:
            raise RuntimeError(f"{SCREENED_ARM}/fixed/{case}: residual tolerance failed")

    write_csv(outdir / "policy_summary.csv", wide_rows)
    write_csv(outdir / "policy_summary_long.csv", long_rows)
    write_csv(outdir / "residuals.csv", residual_rows)
    copy_calibration_artifacts(outdir)
    build_plot(outdir, wide_rows)
    acceptance = {
        "status": "passed",
        "test_count": len(checks),
        "baseline_tolerance": 1.0e-12,
        "residual_tolerance": residual_tolerance,
        "maximum_absolute_reported_residual": max_residual,
        "checks": checks,
    }
    (outdir / "acceptance_tests.json").write_text(
        json.dumps(acceptance, indent=2) + "\n", encoding="utf-8"
    )
    metadata = {
        "status": "complete",
        "packet": "experimental repaired-E5 quota population closure",
        "interpretation": "long-run stationary household comparison; not a forecast",
        "promotion_status": "experimental; not promoted",
        "outside_origin_entrant_share": {
            "value": OUTSIDE_SHARE,
            "classification": "empirically normalized",
            "measure": "across-CBSA outside-origin entrant share",
            "qualification": "provisional; national ACS re-anchor pending",
        },
        "retention_rate": {"classification": "derived identity; candidate-specific"},
        "outside_inflow": {"classification": "derived identity; candidate-specific"},
        "entry_taste_scale": {"classification": "not used in quota production"},
        "outside_value": {"classification": "not used in quota production"},
        "logit": {
            "classification": "unidentified entry-response sensitivity only",
            "taste_scale": 2.0,
        },
        "policy_arms": {arm: ARMS[arm]["label"] for arm in POLICY_ARMS},
        "screened_out_arm": {
            "arm": SCREENED_ARM,
            "label": screened_contract["label"],
            "implied_retention_rate": screened_retention,
            "reason": "Rbar exceeds one at the 0.169 anchor; no corrected policy rows",
        },
        "policies": [case for case, _ in POLICIES],
        "raw_runs": str(raw_root),
        "torch_production_job_ids": args.production_job_id,
        "torch_smoke_job_ids": args.smoke_job_id,
        "acceptance": acceptance,
        "out_of_scope_unchanged": [
            "estimated parameters",
            "12-moment system",
            "solver numerics",
            "ACS re-anchoring pull",
            "half-life computation",
            "CE-unit sensitivity curves",
            "promoted-model files",
            "paper LaTeX",
        ],
    }
    (outdir / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    build_readme(
        outdir,
        wide_rows,
        acceptance,
        args.production_job_id,
        args.smoke_job_id,
    )
    print(f"wrote validated packet: {outdir}")


if __name__ == "__main__":
    main()
