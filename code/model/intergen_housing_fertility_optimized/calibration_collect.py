#!/usr/bin/env python3
"""Collect strict winners from the isolated six-hour M5 continuation."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

from .m5_profile import (
    M5_LOSS,
    M5_OPTIMIZED_STRICT_ROOT_LOSS,
    M5_PROFILE_NAME,
    m5_target_system,
)
from .new_moment_profile import NEW_MOMENT_PROFILE_NAME, new_moment_target_system


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile", choices=("m5", "new-moments"), default="m5")
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def repeat_is_exact(records: list[dict[str, Any]]) -> bool:
    return bool(
        len(records) == 2
        and all(record.get("strict_converged") for record in records)
        and records[0].get("rank_loss") == records[1].get("rank_loss")
        and records[0].get("price") == records[1].get("price")
        and records[0].get("moments") == records[1].get("moments")
    )


def main() -> None:
    args = parse_args()
    target_system = m5_target_system() if args.profile == "m5" else new_moment_target_system()
    profile_name = M5_PROFILE_NAME if args.profile == "m5" else NEW_MOMENT_PROFILE_NAME
    args.output.mkdir(parents=True, exist_ok=True)
    chains: list[dict[str, Any]] = []
    eligible: list[tuple[Path, dict[str, Any], dict[str, Any]]] = []
    for summary_path in sorted(args.run_root.glob("chain_*/summary.json")):
        summary = json.loads(summary_path.read_text())
        metadata = dict(summary.get("metadata") or {})
        tight_records = list(summary.get("tight_repeats") or [])
        exact = repeat_is_exact(tight_records)
        tight = tight_records[0] if exact else None
        row = {
            "chain": summary_path.parent.name,
            "method": metadata.get("method"),
            "seed": metadata.get("seed"),
            "start_mix": metadata.get("start_mix"),
            "initial_step": metadata.get("initial_step"),
            "status": summary.get("status"),
            "cases_completed": summary.get("n_cases_completed"),
            "converged_cases": summary.get("n_converged"),
            "infeasible_cases": summary.get("n_infeasible"),
            "failed_cases": summary.get("n_failed"),
            "elapsed_sec": summary.get("elapsed_sec"),
            "best_search_loss": (summary.get("best_search") or {}).get("rank_loss"),
            "strict_loss": tight.get("rank_loss") if tight else None,
            "strict_residual": tight.get("market_residual") if tight else None,
            "repeat_exact": exact,
            "eligible": bool(exact and summary.get("eligible")),
        }
        chains.append(row)
        if row["eligible"]:
            eligible.append((summary_path.parent, summary, tight))

    if not eligible:
        write_csv(args.output / "all_chains.csv", chains)
        raise RuntimeError("no chain has two exact strict winner repeats")
    selected_dir, selected_summary, selected = min(
        eligible, key=lambda item: float(item[2]["rank_loss"])
    )
    selected_loss = float(selected["rank_loss"])
    canonical_improvement = float(M5_LOSS) - selected_loss if args.profile == "m5" else None
    evaluator_improvement = (
        float(M5_OPTIMIZED_STRICT_ROOT_LOSS) - selected_loss if args.profile == "m5" else None
    )
    fit = list(selected["target_fit"])
    parameters = list(selected["parameters"])
    if len(fit) != target_system.count:
        raise RuntimeError("selected strict winner does not contain the complete target-fit table")
    if not math.isclose(
        sum(float(row["loss_contribution"]) for row in fit),
        selected_loss,
        rel_tol=1e-12,
        abs_tol=1e-12,
    ):
        raise RuntimeError("selected target-fit contributions do not sum to the strict loss")
    if sum(row.get("role") == "estimated" for row in parameters) != 14:
        raise RuntimeError("selected parameter table does not contain all 14 free parameters")

    write_csv(args.output / "all_chains.csv", chains)
    write_csv(args.output / "selected_target_fit.csv", fit)
    write_csv(args.output / "selected_parameters.csv", parameters)
    selected_payload = {
        "status": (
            "beats_canonical_m5"
            if args.profile == "m5" and canonical_improvement > 0.0
            else "improves_optimized_evaluator_only"
            if args.profile == "m5" and evaluator_improvement > 0.0
            else "selected_lowest_strict_loss"
            if args.profile == "new-moments"
            else "no_strict_improvement"
        ),
        "profile": profile_name,
        "target_system": target_system.name,
        "target_fingerprint": target_system.fingerprint,
        "free_parameter_count": 14,
        "target_count": target_system.count,
        "external_restrictions": {"theta_n": 0.0},
        "canonical_m5_loss": M5_LOSS,
        "optimized_strict_root_m5_loss": M5_OPTIMIZED_STRICT_ROOT_LOSS,
        "selected_loss": selected_loss,
        "loss_improvement_vs_canonical_m5": canonical_improvement,
        "loss_improvement_vs_optimized_evaluator": evaluator_improvement,
        "relative_improvement_vs_canonical_m5": (
            canonical_improvement / float(M5_LOSS) if canonical_improvement is not None else None
        ),
        "selected_chain": selected_dir.name,
        "eligible_chain_count": len(eligible),
        "completed_chain_count": len(chains),
        "selected": selected,
        "selected_chain_summary": selected_summary,
    }
    (args.output / "results.json").write_text(
        json.dumps(selected_payload, indent=2, sort_keys=True) + "\n"
    )
    report = [
        "# Calibration search collection",
        "",
        f"Profile: `{profile_name}`; target system: `{target_system.name}`.",
        "",
        f"- Completed chains: `{len(chains)}`",
        f"- Eligible chains with two exact strict repeats: `{len(eligible)}`",
        f"- Selected chain: `{selected_dir.name}`",
        f"- Selected strict loss: `{selected_loss:.12g}`",
        f"- Verdict: **{selected_payload['status'].replace('_', ' ').upper()}**",
        "",
        "The complete target-fit, parameter/bounds, and chain tables are adjacent to this report.",
        "No result is eligible on its loose search loss alone.",
        "",
    ]
    if args.profile == "m5":
        report[7:7] = [
            f"- Canonical M5 loss: `{M5_LOSS:.12g}`",
            f"- Optimized strict-root M5 baseline: `{M5_OPTIMIZED_STRICT_ROOT_LOSS:.12g}`",
            f"- Improvement versus canonical M5: `{canonical_improvement:.12g}`",
            f"- Improvement versus optimized-evaluator M5: `{evaluator_improvement:.12g}`",
        ]
    (args.output / "REPORT.md").write_text("\n".join(report))
    print(json.dumps({key: value for key, value in selected_payload.items() if key not in {"selected", "selected_chain_summary"}}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
