#!/usr/bin/env python3
"""Certify and combine separately solved funded-policy PF horizons."""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for search_path in (MODEL_ROOT, TOOLS_ROOT):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

import run_e5f_perfect_foresight_rebated_property_tax as policy
import run_e5f_perfect_foresight_transition as pf


CASES = ("rebated-tax1-baseline", "rebated-tax2-reform")
EARLY_PATH_PERIODS = 5
EARLY_PATH_RELATIVE_TOLERANCE = 1.0e-3


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, action="append", required=True)
    parser.add_argument("--reform-dir", type=Path, action="append", required=True)
    parser.add_argument(
        "--baseline-stability-path",
        action="append",
        default=[],
        metavar="HORIZON=CSV",
    )
    parser.add_argument(
        "--reform-stability-path",
        action="append",
        default=[],
        metavar="HORIZON=CSV",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--early-path-periods", type=int, default=EARLY_PATH_PERIODS
    )
    parser.add_argument(
        "--diagnostic-unconverged",
        action="store_true",
        help=(
            "Write a clearly unpromoted comparison packet even when declared "
            "solution gates fail. Strict certification remains the default."
        ),
    )
    parser.add_argument(
        "--diagnostic-note",
        action="append",
        default=[],
        help="Repeatable limitation note stored in diagnostic output.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def load_csv(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"Refusing an empty path: {path}")
    return rows


def load_stability_path(specification: str) -> dict[str, Any]:
    horizon_text, separator, path_text = specification.partition("=")
    if not separator:
        raise ValueError(
            f"Stability path must have HORIZON=CSV form: {specification}"
        )
    horizon = int(horizon_text)
    path = Path(path_text).resolve()
    rows = load_csv(path)
    if len(rows) != horizon:
        raise RuntimeError(
            f"Horizon/path length mismatch at {path}: {horizon}, {len(rows)}"
        )
    if [int(row["period"]) for row in rows] != list(range(horizon)):
        raise RuntimeError(f"Noncontiguous path periods at {path}")
    return {
        "horizon": horizon,
        "rows": rows,
        "transition_path": path,
        "transition_path_sha256": pf.file_sha256(path),
        "source_kind": "atomic_checkpoint",
    }


def load_packet(path: Path, expected_case: str) -> dict[str, Any]:
    root = path.resolve()
    summary_path = root / "summary.json"
    summary = load_json(summary_path)
    if summary.get("status") != (
        "complete_unpromoted_perfect_foresight_funded_policy_diagnostic"
    ):
        raise RuntimeError(f"Incomplete funded-policy packet: {root}")
    if summary.get("cases") != [expected_case]:
        raise RuntimeError(
            f"Packet case mismatch at {root}: {summary.get('cases')} versus {expected_case}"
        )
    case_solutions = (summary.get("solutions") or {}).get(expected_case) or {}
    if len(case_solutions) != 1:
        raise RuntimeError(f"Each staged packet must contain one horizon: {root}")
    horizon = int(next(iter(case_solutions)))
    horizon_summary = dict(case_solutions[str(horizon)])
    transition_path = (
        root / expected_case / f"horizon_{horizon:03d}" / "transition_path.csv"
    )
    path_rows = load_csv(transition_path)
    if len(path_rows) != horizon:
        raise RuntimeError(
            f"Horizon/path length mismatch at {root}: {horizon}, {len(path_rows)}"
        )
    if [int(row["period"]) for row in path_rows] != list(range(horizon)):
        raise RuntimeError(f"Noncontiguous path periods at {root}")
    return {
        "root": root,
        "summary": summary,
        "summary_path": summary_path,
        "summary_sha256": pf.file_sha256(summary_path),
        "horizon": horizon,
        "horizon_summary": horizon_summary,
        "rows": path_rows,
        "transition_path": transition_path,
        "transition_path_sha256": pf.file_sha256(transition_path),
    }


def exact_common_contract(left: dict[str, Any], right: dict[str, Any]) -> None:
    fields = (
        "common_initial_state",
        "housing_supply",
        "psi_persistence_per_four_year_period",
        "terminal_psi_child",
        "shock_interpretation",
    )
    for field in fields:
        if left["summary"].get(field) != right["summary"].get(field):
            raise RuntimeError(f"Policy packets disagree on {field}")
    left_provenance = left["summary"].get("provenance") or {}
    right_provenance = right["summary"].get("provenance") or {}
    for field in (
        "selected_report_sha256",
        "selected_transition_sha256",
        "stationary_source_sha256",
        "perfect_foresight_driver_sha256",
        "driver_sha256",
    ):
        if left_provenance.get(field) != right_provenance.get(field):
            raise RuntimeError(f"Policy packets disagree on provenance field {field}")


def packet_solution_gates(packet: dict[str, Any]) -> dict[str, Any]:
    result = packet["horizon_summary"]
    case = packet["summary"]["cases"][0]
    terminal = (packet["summary"].get("terminal_steady_states") or {}).get(case) or {}
    return {
        "case": case,
        "horizon": packet["horizon"],
        "path_status": result.get("status"),
        "path_converged": bool(result.get("converged", False)),
        "maximum_market_residual": float(result.get("market_residual", math.inf)),
        "maximum_fiscal_residual": float(result.get("fiscal_residual", math.inf)),
        "terminal_fixed_point_passed": (
            (terminal.get("fixed_point_gates") or {}).get("status") == "passed"
        ),
        "terminal_path_convergence_passed": bool(
            (result.get("terminal_convergence") or {}).get("all_checks_pass", False)
        ),
    }


def path_arrays(rows: Sequence[dict[str, Any]]) -> dict[str, np.ndarray]:
    births_per_adult = np.array(
        [
            float(row["birth_children_topcode_adjusted"])
            / float(row["adult_population"])
            for row in rows
        ]
    )
    return {
        "asset_price": np.array([float(row["asset_price"]) for row in rows]),
        "renter_price": np.array([float(row["renter_price"]) for row in rows]),
        "equal_transfer": np.array(
            [float(row["equal_transfer_period_units"]) for row in rows]
        ),
        "births_per_adult": births_per_adult,
        "adult_population": np.array(
            [float(row["adult_population"]) for row in rows]
        ),
        "owner_rate": np.array([float(row["owner_rate"]) for row in rows]),
    }


def stability_rows(
    packets: dict[int, dict[str, Any]], *, label: str, comparison_periods: int
) -> list[dict[str, Any]]:
    horizons = sorted(packets)
    reference_horizon = horizons[-1]
    reference = path_arrays(packets[reference_horizon]["rows"])
    rows: list[dict[str, Any]] = []
    for horizon in horizons[:-1]:
        candidate = path_arrays(packets[horizon]["rows"])
        common = min(comparison_periods, horizon, reference_horizon)
        for variable in reference:
            left = candidate[variable][:common]
            right = reference[variable][:common]
            scale = np.maximum(np.abs(right), 1e-12)
            rows.append(
                {
                    "comparison": label,
                    "short_horizon": horizon,
                    "reference_horizon": reference_horizon,
                    "comparison_periods": common,
                    "variable": variable,
                    "maximum_absolute_gap": float(np.max(np.abs(left - right))),
                    "maximum_relative_gap": float(
                        np.max(np.abs(left - right) / scale)
                    ),
                }
            )
    return rows


def effect_stability_rows(
    effects: dict[int, list[dict[str, Any]]], *, comparison_periods: int
) -> list[dict[str, Any]]:
    horizons = sorted(effects)
    reference_horizon = horizons[-1]
    fields = (
        "asset_price_percent_change",
        "renter_price_percent_change",
        "equal_transfer_percent_change",
        "owner_rate_pp_change",
        "births_per_adult_percent_change",
        "adult_population_percent_change",
    )
    rows: list[dict[str, Any]] = []
    for horizon in horizons[:-1]:
        common = min(comparison_periods, horizon, reference_horizon)
        for field in fields:
            left = np.array(
                [float(row[field]) for row in effects[horizon][:common]]
            )
            right = np.array(
                [float(row[field]) for row in effects[reference_horizon][:common]]
            )
            # Effects are already expressed in percent or percentage points.
            # Near-zero effects therefore use one percentage point as the
            # scale; dividing by numerical noise would make stability
            # mechanically fail before the demographic lag begins.
            scale = np.maximum(np.abs(right), 1.0)
            rows.append(
                {
                    "comparison": "policy_effect",
                    "short_horizon": horizon,
                    "reference_horizon": reference_horizon,
                    "comparison_periods": common,
                    "variable": field,
                    "maximum_absolute_gap": float(np.max(np.abs(left - right))),
                    "maximum_relative_gap": float(
                        np.max(np.abs(left - right) / scale)
                    ),
                }
            )
    return rows


def make_figure(
    baseline_rows: Sequence[dict[str, Any]],
    reform_rows: Sequence[dict[str, Any]],
    effects: Sequence[dict[str, Any]],
    output_dir: Path,
    *,
    filename_stem: str,
) -> None:
    years = [int(row["calendar_year"]) for row in baseline_rows]
    figure, axes = plt.subplots(3, 2, figsize=(10.3, 9.0), constrained_layout=True)
    levels = (
        ("asset_price", "House price"),
        ("renter_price", "Rent"),
        ("equal_transfer_period_units", "Equal transfer"),
        ("adult_population", "Adult population"),
    )
    for axis, (field, title) in zip(axes.flat[:4], levels, strict=True):
        axis.plot(years, [float(row[field]) for row in baseline_rows], label="Rebated 1%")
        axis.plot(years, [float(row[field]) for row in reform_rows], label="Rebated 2%")
        axis.set_title(title)
    axes[0, 0].legend(frameon=False)
    axes[2, 0].plot(
        years,
        [float(row["births_per_adult_percent_change"]) for row in effects],
        color="#B4433C",
    )
    axes[2, 0].axhline(0.0, color="0.5", linewidth=1)
    axes[2, 0].set_title("Births per adult: reform effect (%)")
    axes[2, 1].plot(
        years,
        [float(row["adult_population_percent_change"]) for row in effects],
        color="#21476B",
    )
    axes[2, 1].axhline(0.0, color="0.5", linewidth=1)
    axes[2, 1].set_title("Adult population: reform effect (%)")
    for axis in axes.flat:
        axis.grid(alpha=0.2)
        axis.set_xlabel("Year")
    figure.savefig(output_dir / f"{filename_stem}.png", dpi=220)
    figure.savefig(output_dir / f"{filename_stem}.pdf")
    plt.close(figure)


def write_readme(
    output_dir: Path,
    *,
    status: str,
    gates: Sequence[dict[str, Any]],
    limitations: Sequence[str],
) -> None:
    lines = [
        "# Perfect-Foresight Rebated Property-Tax Comparison",
        "",
        f"Status: `{status}`.",
        "",
        "This packet compares the rebated 1% and rebated 2% property-tax paths. "
        "It does not alter or replace the certified unrebated 1% status-quo baseline.",
        "",
        "## Numerical gates",
        "",
        "| Case | Horizon | Market residual | Fiscal residual | Path converged | Terminal path converged |",
        "|---|---:|---:|---:|:---:|:---:|",
    ]
    for gate in gates:
        lines.append(
            "| {case} | {horizon} | {market:.8g} | {fiscal:.8g} | {path} | {terminal} |".format(
                case=gate["case"],
                horizon=gate["horizon"],
                market=gate["maximum_market_residual"],
                fiscal=gate["maximum_fiscal_residual"],
                path="yes" if gate["path_converged"] else "no",
                terminal="yes" if gate["terminal_path_convergence_passed"] else "no",
            )
        )
    if limitations:
        lines.extend(["", "## Limitations", ""])
        lines.extend(f"- {note}" for note in limitations)
    lines.extend(
        [
            "",
            "The CSV files contain the exact policy effects and horizon-stability checks; "
            "the figure is a visual diagnostic, not a certification of the full asymptotic transition.",
            "",
        ]
    )
    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    if args.early_path_periods <= 0:
        raise ValueError("--early-path-periods must be positive")
    if len(args.baseline_dir) != len(args.reform_dir):
        raise ValueError("Supply the same number of baseline and reform packets.")
    if len(args.baseline_stability_path) != len(args.reform_stability_path):
        raise ValueError(
            "Supply the same number of baseline and reform stability paths."
        )
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    packets: dict[str, dict[int, dict[str, Any]]] = {case: {} for case in CASES}
    for expected_case, directories in zip(
        CASES, (args.baseline_dir, args.reform_dir), strict=True
    ):
        for directory in directories:
            packet = load_packet(directory, expected_case)
            horizon = int(packet["horizon"])
            if horizon in packets[expected_case]:
                raise RuntimeError(f"Duplicate {expected_case} horizon {horizon}")
            packets[expected_case][horizon] = packet
    horizons = sorted(packets[CASES[0]])
    if horizons != sorted(packets[CASES[1]]):
        raise RuntimeError("Baseline and reform horizon sets differ.")
    if len(horizons) < 2 and not args.diagnostic_unconverged:
        raise RuntimeError("Certification requires at least two horizons.")

    for horizon in horizons:
        exact_common_contract(
            packets[CASES[0]][horizon], packets[CASES[1]][horizon]
        )
    reference = packets[CASES[0]][horizons[0]]
    for case in CASES:
        for horizon in horizons:
            exact_common_contract(reference, packets[case][horizon])

    stability_packets: dict[str, dict[int, dict[str, Any]]] = {
        case: dict(packets[case]) for case in CASES
    }
    stability_sources: dict[str, dict[str, dict[str, Any]]] = {
        case: {} for case in CASES
    }
    for case, specifications in zip(
        CASES,
        (args.baseline_stability_path, args.reform_stability_path),
        strict=True,
    ):
        for specification in specifications:
            item = load_stability_path(specification)
            horizon = int(item["horizon"])
            if horizon in stability_packets[case]:
                raise RuntimeError(f"Duplicate {case} stability horizon {horizon}")
            stability_packets[case][horizon] = item
            stability_sources[case][str(horizon)] = {
                "source_kind": item["source_kind"],
                "transition_path": str(item["transition_path"]),
                "transition_path_sha256": item["transition_path_sha256"],
            }
    comparison_horizons = sorted(stability_packets[CASES[0]])
    if comparison_horizons != sorted(stability_packets[CASES[1]]):
        raise RuntimeError("Baseline and reform stability horizons differ.")
    if len(comparison_horizons) < 2:
        raise RuntimeError("Horizon stability requires at least two paths.")

    gates = [
        packet_solution_gates(packets[case][horizon])
        for case in CASES
        for horizon in horizons
    ]
    acceptance_violations: list[str] = []
    for gate in gates:
        if not gate["path_converged"]:
            acceptance_violations.append(
                f"Unconverged housing/fiscal path: {gate}"
            )
        if gate["maximum_market_residual"] > 2e-4:
            acceptance_violations.append(f"Housing gate failed: {gate}")
        if gate["maximum_fiscal_residual"] > 2.5e-5:
            acceptance_violations.append(f"Fiscal gate failed: {gate}")
        if not gate["terminal_fixed_point_passed"]:
            acceptance_violations.append(f"Terminal fixed point failed: {gate}")

    longest = horizons[-1]
    for gate in gates:
        if gate["horizon"] == longest and not gate[
            "terminal_path_convergence_passed"
        ]:
            acceptance_violations.append(
                f"Longest path has not reached its terminal SS: {gate}"
            )
    if acceptance_violations and not args.diagnostic_unconverged:
        raise RuntimeError(acceptance_violations[0])

    effects = {
        horizon: policy.effect_rows(
            stability_packets[CASES[0]][horizon]["rows"],
            stability_packets[CASES[1]][horizon]["rows"],
        )
        for horizon in comparison_horizons
    }
    for horizon, rows in effects.items():
        pf.write_csv(output_dir / f"policy_effects_h{horizon:03d}.csv", rows)

    stability = []
    stability.extend(
        stability_rows(
            stability_packets[CASES[0]],
            label=CASES[0],
            comparison_periods=args.early_path_periods,
        )
    )
    stability.extend(
        stability_rows(
            stability_packets[CASES[1]],
            label=CASES[1],
            comparison_periods=args.early_path_periods,
        )
    )
    stability.extend(
        effect_stability_rows(
            effects, comparison_periods=args.early_path_periods
        )
    )
    pf.write_csv(output_dir / "horizon_stability.csv", stability)
    maximum_stability_gap = max(
        float(row["maximum_relative_gap"]) for row in stability
    )
    maximum_level_path_stability_gap = max(
        float(row["maximum_relative_gap"])
        for row in stability
        if row["comparison"] != "policy_effect"
    )
    maximum_policy_effect_stability_gap = max(
        float(row["maximum_relative_gap"])
        for row in stability
        if row["comparison"] == "policy_effect"
    )
    if maximum_stability_gap > EARLY_PATH_RELATIVE_TOLERANCE:
        acceptance_violations.append(
            "Early path/effect stability gate failed: "
            f"{maximum_stability_gap}"
        )
        if not args.diagnostic_unconverged:
            raise RuntimeError(acceptance_violations[-1])

    baseline_terminal = packets[CASES[0]][longest]["summary"][
        "terminal_steady_states"
    ][CASES[0]]
    reform_terminal = packets[CASES[1]][longest]["summary"][
        "terminal_steady_states"
    ][CASES[1]]
    terminal_effects = {
        "asset_price_percent_change": 100.0
        * (float(reform_terminal["asset_price"]) / float(baseline_terminal["asset_price"]) - 1.0),
        "renter_price_percent_change": 100.0
        * (float(reform_terminal["renter_price"]) / float(baseline_terminal["renter_price"]) - 1.0),
        "equal_transfer_percent_change": 100.0
        * (float(reform_terminal["equal_transfer"]) / float(baseline_terminal["equal_transfer"]) - 1.0),
        "births_per_adult_percent_change": 100.0
        * (float(reform_terminal["births_per_adult"]) / float(baseline_terminal["births_per_adult"]) - 1.0),
        "owner_rate_pp_change": 100.0
        * (float(reform_terminal["owner_rate"]) - float(baseline_terminal["owner_rate"])),
        "adult_population_percent_change": 100.0
        * (float(reform_terminal["adult_population"]) / float(baseline_terminal["adult_population"]) - 1.0),
    }
    diagnostic_mode = bool(args.diagnostic_unconverged)
    packet_status = (
        "complete_unpromoted_perfect_foresight_funded_policy_comparison_diagnostic"
        if diagnostic_mode
        else "certified_unpromoted_perfect_foresight_funded_policy_comparison"
    )
    figure_stem = (
        "diagnostic_policy_comparison"
        if diagnostic_mode
        else "certified_policy_comparison"
    )
    make_figure(
        packets[CASES[0]][longest]["rows"],
        packets[CASES[1]][longest]["rows"],
        effects[longest],
        output_dir,
        filename_stem=figure_stem,
    )
    early_figure_periods = min(args.early_path_periods, longest)
    make_figure(
        packets[CASES[0]][longest]["rows"][:early_figure_periods],
        packets[CASES[1]][longest]["rows"][:early_figure_periods],
        effects[longest][:early_figure_periods],
        output_dir,
        filename_stem=f"{figure_stem}_first_{early_figure_periods}_periods",
    )
    summary = {
        "status": packet_status,
        "promotion_status": "not_promoted",
        "diagnostic_unconverged": diagnostic_mode,
        "acceptance_status": (
            "passed" if not acceptance_violations else "failed"
        ),
        "acceptance_violations": acceptance_violations,
        "diagnostic_notes": list(args.diagnostic_note),
        "cases": list(CASES),
        "horizons": comparison_horizons,
        "solution_horizons": horizons,
        "longest_horizon": longest,
        "early_path_periods": args.early_path_periods,
        "early_path_relative_tolerance": EARLY_PATH_RELATIVE_TOLERANCE,
        "maximum_early_path_relative_gap": maximum_stability_gap,
        "maximum_level_path_relative_gap": maximum_level_path_stability_gap,
        "maximum_policy_effect_relative_gap": (
            maximum_policy_effect_stability_gap
        ),
        "solution_gates": gates,
        "impact_effects_2023": effects[longest][0],
        "terminal_steady_state_effects": terminal_effects,
        "source_packets": {
            case: {
                str(horizon): {
                    "root": str(packets[case][horizon]["root"]),
                    "summary": str(packets[case][horizon]["summary_path"]),
                    "summary_sha256": packets[case][horizon]["summary_sha256"],
                    "transition_path": str(
                        packets[case][horizon]["transition_path"]
                    ),
                    "transition_path_sha256": packets[case][horizon][
                        "transition_path_sha256"
                    ],
                }
                for horizon in horizons
            }
            for case in CASES
        },
        "stability_path_sources": stability_sources,
        "common_contract": {
            key: reference["summary"].get(key)
            for key in (
                "common_initial_state",
                "housing_supply",
                "psi_persistence_per_four_year_period",
                "terminal_psi_child",
                "shock_interpretation",
            )
        },
    }
    pf.write_json(output_dir / "certification.json", summary)
    write_readme(
        output_dir,
        status=packet_status,
        gates=gates,
        limitations=list(args.diagnostic_note),
    )
    print(json.dumps(pf.jsonable(summary), indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
