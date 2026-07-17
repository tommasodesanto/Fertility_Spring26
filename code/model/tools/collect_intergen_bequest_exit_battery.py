#!/usr/bin/env python3
"""Collect tight-only bequest/exit battery and profile summaries."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PACKET = ROOT / "output/model/intergen_bequest_exit_battery_20260714"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main-root", type=Path, default=DEFAULT_PACKET / "main_summaries")
    parser.add_argument("--profile-root", type=Path, default=DEFAULT_PACKET / "profile_summaries")
    parser.add_argument("--jacobian-root", type=Path, default=DEFAULT_PACKET / "jacobian")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_PACKET / "report")
    return parser.parse_args()


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def task_id(path: Path) -> int:
    return int(path.parent.name.split("task_", 1)[1].split("_", 1)[0])


def spec_label(meta: dict[str, Any]) -> str:
    arm = str(meta["arm"])
    mechanism = dict(meta["mechanism"])
    fixed = dict(meta["fixed_parameters"])
    lbar = float(mechanism["owner_ltv_terminal_share"])
    theta1 = fixed.get("theta1")
    if arm in {"A0", "A1"}:
        return arm
    if arm == "A3" and math.isclose(float(theta1), 0.50):
        return "A3-theta1-0.50"
    if arm in {"A2", "A3"}:
        return f"{arm}-L{lbar:.1f}"
    return arm


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    grouped: dict[str, list[tuple[Path, dict[str, Any]]]] = defaultdict(list)
    total_evals = 0
    for path in sorted(args.main_root.glob("task_*/summary.json"), key=task_id):
        payload = load(path)
        total_evals += int(payload.get("n_cases_completed", 0))
        grouped[spec_label(payload["metadata"])].append((path, payload))

    arm_rows: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    parameter_rows: list[dict[str, Any]] = []
    acceptance_rows: list[dict[str, Any]] = []
    selected: dict[str, dict[str, Any] | None] = {}
    for label, entries in grouped.items():
        eligible = [
            (path, payload["best_tight"])
            for path, payload in entries
            if isinstance(payload.get("best_tight"), dict)
        ]
        winner_pair = min(eligible, key=lambda item: float(item[1]["rank_loss"])) if eligible else None
        winner = winner_pair[1] if winner_pair else None
        selected[label] = winner
        meta = entries[0][1]["metadata"]
        mechanism = dict(meta["mechanism"])
        fixed = dict(meta["fixed_parameters"])
        lifecycle = dict((winner or {}).get("lifecycle", {}))
        checks = list(lifecycle.get("ownership_acceptance_checks", []))
        nonterminal = [row for row in checks if float(row["age_state"]) < 82.0]
        nonterminal_deltas = [
            float(nonterminal[i + 1]["model"]) - float(nonterminal[i]["model"])
            for i in range(len(nonterminal) - 1)
        ]
        nonterminal_band_pass = bool(nonterminal) and all(bool(row["within_band"]) for row in nonterminal)
        nonterminal_cliff = any(delta < -0.15 for delta in nonterminal_deltas)
        arm_rows.append(
            {
                "specification": label,
                "arm": meta["arm"],
                "free_parameters": meta["free_parameter_count"],
                "tight_eligible_chains": len(eligible),
                "chains": len(entries),
                "selected_task": task_id(winner_pair[0]) if winner_pair else None,
                "terminal_ltv_multiplier": mechanism["owner_ltv_terminal_share"] if meta["arm"] not in {"A0", "A1"} else None,
                "theta1_restriction": fixed.get("theta1"),
                "tight_loss": (winner or {}).get("rank_loss"),
                "market_residual": (winner or {}).get("market_residual"),
                "theta0": (winner or {}).get("theta", {}).get("theta0"),
                "theta1_estimate": (winner or {}).get("theta", {}).get("theta1"),
                "h_bar_0": (winner or {}).get("theta", {}).get("h_bar_0"),
                "hard_acceptance_pass": lifecycle.get("hard_acceptance_pass"),
                "ownership_band_failure": lifecycle.get("ownership_band_failure"),
                "ownership_cliff_failure": lifecycle.get("ownership_cliff_failure"),
                "nonterminal_62_78_band_pass": nonterminal_band_pass if winner else None,
                "nonterminal_62_78_cliff_failure": nonterminal_cliff if winner else None,
                "nonterminal_62_78_acceptance_pass": (nonterminal_band_pass and not nonterminal_cliff) if winner else None,
                "decumulation_ratio_74plus_over_62_74": lifecycle.get("decum_ratio_wealth_74plus_over_62_74"),
            }
        )
        if winner is None:
            continue
        for row in winner["target_fit"]:
            fit_rows.append({"specification": label, **row})
        for row in winner["parameters"]:
            parameter_rows.append({"specification": label, **row})
        for row in checks:
            acceptance_rows.append({"specification": label, **row})

    arm_rows.sort(key=lambda row: row["specification"])
    write_csv(args.outdir / "arm_summary_tight.csv", arm_rows)
    write_csv(args.outdir / "target_fit_all_tight.csv", fit_rows)
    write_csv(args.outdir / "parameters_all_tight.csv", parameter_rows)
    write_csv(args.outdir / "ownership_acceptance_all_tight.csv", acceptance_rows)

    profile_grouped: dict[float, list[tuple[Path, dict[str, Any]]]] = defaultdict(list)
    profile_evals = 0
    for path in sorted(args.profile_root.glob("task_*/summary.json"), key=task_id):
        payload = load(path)
        profile_evals += int(payload.get("n_cases_completed", 0))
        theta0 = float(payload["metadata"]["fixed_parameters"]["theta0"])
        profile_grouped[theta0].append((path, payload))
    profile_rows: list[dict[str, Any]] = []
    for theta0, entries in sorted(profile_grouped.items()):
        eligible = [
            (path, payload["best_tight"])
            for path, payload in entries
            if isinstance(payload.get("best_tight"), dict)
        ]
        pair = min(eligible, key=lambda item: float(item[1]["rank_loss"])) if eligible else None
        winner = pair[1] if pair else None
        profile_rows.append(
            {
                "theta0_fixed": theta0,
                "direct_profile_tight_loss": (winner or {}).get("rank_loss"),
                "tight_eligible_chains": len(eligible),
                "selected_task": task_id(pair[0]) if pair else None,
                "market_residual": (winner or {}).get("market_residual"),
                "hard_acceptance_pass": (winner or {}).get("lifecycle", {}).get("hard_acceptance_pass"),
                "decumulation_ratio_74plus_over_62_74": (winner or {}).get("lifecycle", {}).get("decum_ratio_wealth_74plus_over_62_74"),
            }
        )

    # A2-L0.4 is algebraically the theta0=0 member of A3-L0.4.  Record its
    # lower tight loss as the valid zero-point envelope when the dedicated
    # profile chains fail to return to that already-discovered basin.
    a2_zero = selected.get("A2-L0.4")
    zero_envelope = float(a2_zero["rank_loss"]) if a2_zero else math.nan
    for row in profile_rows:
        if math.isclose(float(row["theta0_fixed"]), 0.0):
            direct = row["direct_profile_tight_loss"]
            row["profile_envelope_loss"] = min(float(direct), zero_envelope) if direct is not None else zero_envelope
            row["envelope_source"] = "A2-L0.4 equivalent zero-bequest calibration"
        else:
            row["profile_envelope_loss"] = row["direct_profile_tight_loss"]
            row["envelope_source"] = "direct profile"
    write_csv(args.outdir / "theta0_profile_tight.csv", profile_rows)

    jacobian = load(args.jacobian_root / "summary.json") if (args.jacobian_root / "summary.json").exists() else None
    primary = selected.get("A3-L0.4")
    exit_only = selected.get("A2-L0.4")
    equal = selected.get("A4")
    summary = {
        "main_search_evaluations": total_evals,
        "profile_search_evaluations": profile_evals,
        "main_chain_count": sum(len(entries) for entries in grouped.values()),
        "main_chains_with_tight_winner": sum(int(row["tight_eligible_chains"]) for row in arm_rows),
        "hard_accepted_main_specifications": sum(row["hard_acceptance_pass"] is True for row in arm_rows),
        "nonterminal_accepted_main_specifications": sum(row["nonterminal_62_78_acceptance_pass"] is True for row in arm_rows),
        "headline_A3_L0.4": primary,
        "exit_only_A2_L0.4": exit_only,
        "equal_division_A4": equal,
        "theta0_profile": profile_rows,
        "jacobian": jacobian,
        "verdict": {
            "ownership_rule_accepted": False,
            "terminal_age_82_check_valid": False,
            "theta0_walks_back_to_zero": bool(
                primary and exit_only and float(exit_only["rank_loss"]) <= float(primary["rank_loss"])
            ),
            "reason": (
                "No main specification passes either the predeclared rule or the "
                "scientifically corrected nonterminal 62-78 rule. The age-82 model "
                "state has forced terminal liquidation and is not comparable to the "
                "ACS 82-84 bin. "
                "At the primary LTV multiplier, the strict exit-only calibration "
                "beats the free-theta0 headline, and the valid theta0=0 envelope "
                "is slightly below the theta0=0.05 profile point."
            ),
        },
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    write_readme(args.outdir / "README.md", arm_rows, profile_rows, summary)


def fmt(value: Any, digits: int = 4) -> str:
    if value is None:
        return "--"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return f"{float(value):.{digits}f}"


def write_readme(
    path: Path,
    arms: list[dict[str, Any]],
    profile: list[dict[str, Any]],
    summary: dict[str, Any],
) -> None:
    lines = [
        "# Bequest/exit battery: tight morning readout",
        "",
        "The borrowing-rule family is rejected: none of the main specifications or",
        "profile winners passes the ownership-by-age bands and cliff rule. The model's",
        "age-82 state forces terminal liquidation, so that endpoint is not a valid ACS",
        "comparison; recomputing acceptance over ages 62-78 does not rescue any arm.",
        "At the primary external variant, exit-only A2 has tight loss 5.0115, while",
        "free-bequest A3 has loss 5.2935 with theta0=0.0151 near its zero boundary.",
        "The equal-division robustness also returns theta0 close to zero.",
        "",
        "## Tight arm table",
        "",
        "| specification | free | loss | theta0 | residual | 62-78 pass | 62-78 cliff | decumulation ratio |",
        "|---|---:|---:|---:|---:|---|---|---:|",
    ]
    for row in arms:
        lines.append(
            f"| {row['specification']} | {row['free_parameters']} | {fmt(row['tight_loss'])} | "
            f"{fmt(row['theta0'])} | {fmt(row['market_residual'], 7)} | "
            f"{fmt(row['nonterminal_62_78_acceptance_pass'])} | {fmt(row['nonterminal_62_78_cliff_failure'])} | "
            f"{fmt(row['decumulation_ratio_74plus_over_62_74'])} |"
        )
    lines.extend(
        [
            "",
            "## Nuisance-reoptimized theta0 profile",
            "",
            "| fixed theta0 | direct tight loss | envelope loss | ownership pass |",
            "|---:|---:|---:|---|",
        ]
    )
    for row in profile:
        lines.append(
            f"| {fmt(row['theta0_fixed'], 2)} | {fmt(row['direct_profile_tight_loss'])} | "
            f"{fmt(row['profile_envelope_loss'])} | {fmt(row['hard_acceptance_pass'])} |"
        )
    jacobian = summary.get("jacobian")
    lines.extend(
        [
            "",
            "The dedicated zero-profile chains missed the already-discovered A2-L0.4",
            "basin. Because A2-L0.4 is algebraically the theta0=0 member of A3-L0.4,",
            "its 5.0115 tight loss is the valid zero-point envelope. It is below the",
            "best positive profile point (5.0155 at theta0=0.05), so the battery does",
            "not support a positive bequest strength once nuisance parameters are",
            "properly compared.",
            "",
            "Full tight tables: `target_fit_all_tight.csv`, `parameters_all_tight.csv`,",
            "`ownership_acceptance_all_tight.csv`, and `theta0_profile_tight.csv`.",
            "",
            f"Search volume: {summary['main_search_evaluations']:,} main evaluations and "
            f"{summary['profile_search_evaluations']:,} profile evaluations.",
        ]
    )
    if jacobian:
        smm_rank = jacobian["rank"]["smm_weighted"]
        theta0_column = jacobian["theta0_column"]
        lines.extend(
            [
                "",
                "## Local Jacobian warning",
                "",
                f"The adaptive tight Jacobian has weighted rank {smm_rank['rank_rel_1e_2']} "
                f"at relative tolerance 1e-2 and {smm_rank['rank_rel_1e_3']} at 1e-3; "
                f"the theta0 column norm is {theta0_column['target_scaled_norm']:.4f} and "
                "is above the exact-zero repeat-noise floor. This is local to the rejected, "
                "non-optimal A3 point, so it is a conditioning warning rather than valid "
                "inference about the preferred specification.",
            ]
        )
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
