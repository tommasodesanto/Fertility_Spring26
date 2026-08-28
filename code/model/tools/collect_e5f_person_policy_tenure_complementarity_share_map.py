#!/usr/bin/env python3
"""Audit endpoint tenure shares and test full-path complementarity feasibility."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable

import matplotlib.pyplot as plt
import numpy as np


MARKET_TOLERANCE = 2.0e-4
FISCAL_TOLERANCE = 2.5e-5
AFFINITY_TOLERANCE = 5.0e-14
AFFINE_COLUMNS = (
    "housing_demand",
    "resident_persons",
    "household_heads",
    "birth_children",
    "birth_children_topcode_adjusted",
    "next_resident_persons",
    "next_household_heads",
    "property_tax_revenue",
    "equal_transfer_outlays",
    "government_budget_residual",
    "owner_rate",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--share-zero-dir", type=Path, required=True)
    parser.add_argument("--share-half-dir", type=Path, required=True)
    parser.add_argument("--share-one-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if [int(row["period"]) for row in rows] != list(range(len(rows))):
        raise RuntimeError(f"Noncontiguous transition path: {path}")
    return rows


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write an empty table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def signed_market_residual(row: dict[str, str]) -> float:
    supply = float(row["housing_supply"])
    if not math.isfinite(supply) or supply <= 0.0:
        raise RuntimeError("Housing supply must be finite and positive")
    residual = (supply - float(row["housing_demand"])) / supply
    stored = float(row["relative_market_residual"])
    if abs(abs(residual) - stored) > 2.0e-15:
        raise RuntimeError("Stored market residual is not the absolute signed residual")
    return residual


def allowed_interval(y_zero: float, y_one: float, tolerance: float) -> tuple[float, float] | None:
    slope = float(y_one) - float(y_zero)
    if abs(slope) < 1.0e-30:
        return (0.0, 1.0) if abs(float(y_zero)) <= float(tolerance) else None
    roots = (
        (-float(tolerance) - float(y_zero)) / slope,
        (float(tolerance) - float(y_zero)) / slope,
    )
    lower = max(0.0, min(roots))
    upper = min(1.0, max(roots))
    return (lower, upper) if lower <= upper else None


def validate_packet(
    directory: Path,
    *,
    expected_share: float,
    expected_initial_hash: str,
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, str]]]:
    summary_path = directory / "summary.json"
    path_file = directory / "transition_path.csv"
    oracle_path = directory / "tenure_complementarity_oracle_summary.json"
    summary = load_json(summary_path)
    oracle = load_json(oracle_path)
    rows = load_rows(path_file)
    failures: list[str] = []
    if oracle.get("initial_path_sha256") != expected_initial_hash:
        failures.append("initial_path_hash")
    if float((oracle.get("transformation") or {}).get("owner_share", math.nan)) != float(expected_share):
        failures.append("owner_share")
    if oracle.get("production_summary_sha256") != file_sha256(summary_path):
        failures.append("summary_hash")
    if oracle.get("production_transition_path_sha256") != file_sha256(path_file):
        failures.append("path_hash")
    if summary.get("history_reproduction_status") != "passed":
        failures.append("history")
    if not all((summary.get("accounting_gates") or {}).values()):
        failures.append("accounting")
    if not bool((summary.get("terminal_convergence") or {}).get("all_checks_pass")):
        failures.append("terminal_path")
    if not all(((summary.get("terminal_root") or {}).get("gates") or {}).values()):
        failures.append("terminal_root")
    if max(abs(float(row["feasibility_frontier_projection_mass"])) for row in rows) != 0.0:
        failures.append("feasibility")
    if failures:
        raise RuntimeError(f"Packet {directory} failed: {failures}")
    return summary, oracle, rows


def main() -> None:
    args = parse_args()
    manifest = load_json(args.manifest.resolve())
    fixed_hash = str(manifest["fixed_path_sha256"])
    packets: dict[float, tuple[dict[str, Any], dict[str, Any], list[dict[str, str]]]] = {}
    for share, directory in (
        (0.0, args.share_zero_dir.resolve()),
        (0.5, args.share_half_dir.resolve()),
        (1.0, args.share_one_dir.resolve()),
    ):
        packets[share] = validate_packet(
            directory,
            expected_share=share,
            expected_initial_hash=fixed_hash,
        )

    zero_rows = packets[0.0][2]
    half_rows = packets[0.5][2]
    one_rows = packets[1.0][2]
    identities = [
        [(row["period"], row["calendar_year"]) for row in rows]
        for rows in (zero_rows, half_rows, one_rows)
    ]
    if identities[0] != identities[1] or identities[0] != identities[2]:
        raise RuntimeError("Endpoint and midpoint paths are not aligned")

    affinity: dict[str, float] = {}
    for column in AFFINE_COLUMNS:
        affinity[column] = max(
            abs(
                float(half_rows[index][column])
                - 0.5
                * (float(zero_rows[index][column]) + float(one_rows[index][column]))
            )
            for index in range(len(zero_rows))
        )
    signed_market_affinity = max(
        abs(
            signed_market_residual(half_rows[index])
            - 0.5
            * (
                signed_market_residual(zero_rows[index])
                + signed_market_residual(one_rows[index])
            )
        )
        for index in range(len(zero_rows))
    )
    affinity["signed_market_residual"] = signed_market_affinity
    if affinity["signed_market_residual"] > AFFINITY_TOLERANCE:
        raise RuntimeError("The signed market residual is not affine in the atom share")
    if affinity["government_budget_residual"] > AFFINITY_TOLERANCE:
        raise RuntimeError("The fiscal residual is not affine in the atom share")

    constraint_rows: list[dict[str, Any]] = []
    intersection_lower = 0.0
    intersection_upper = 1.0
    definitions: tuple[tuple[str, float, Callable[[dict[str, str]], float]], ...] = (
        ("market", MARKET_TOLERANCE, signed_market_residual),
        ("fiscal", FISCAL_TOLERANCE, lambda row: float(row["government_budget_residual"])),
    )
    for name, tolerance, measure in definitions:
        for index, zero_row in enumerate(zero_rows):
            zero_value = measure(zero_row)
            one_value = measure(one_rows[index])
            interval = allowed_interval(zero_value, one_value, tolerance)
            if interval is None:
                lower = None
                upper = None
                status = "empty_at_this_date"
                intersection_lower = 1.0
                intersection_upper = 0.0
            else:
                lower, upper = interval
                status = "nonempty_at_this_date"
                intersection_lower = max(intersection_lower, lower)
                intersection_upper = min(intersection_upper, upper)
            constraint_rows.append(
                {
                    "object": name,
                    "period": int(zero_row["period"]),
                    "calendar_year": int(zero_row["calendar_year"]),
                    "tolerance": tolerance,
                    "share_zero_residual": zero_value,
                    "share_half_residual": measure(half_rows[index]),
                    "share_one_residual": one_value,
                    "allowed_share_lower": lower,
                    "allowed_share_upper": upper,
                    "status": status,
                }
            )

    output_dir = args.output_dir.resolve()
    existing = sorted(path.name for path in output_dir.iterdir()) if output_dir.exists() else []
    allowed_existing = {args.manifest.resolve().name}
    unexpected = [name for name in existing if name not in allowed_existing]
    if unexpected:
        raise FileExistsError(f"Refusing to overwrite existing collector output: {unexpected}")
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "date_share_constraints.csv", constraint_rows)

    shares = np.linspace(0.0, 1.0, 401)
    market_maxima: list[float] = []
    fiscal_maxima: list[float] = []
    for share in shares:
        market_maxima.append(
            max(
                abs(
                    (1.0 - share) * signed_market_residual(zero_rows[index])
                    + share * signed_market_residual(one_rows[index])
                )
                for index in range(len(zero_rows))
            )
        )
        fiscal_maxima.append(
            max(
                abs(
                    (1.0 - share) * float(zero_rows[index]["government_budget_residual"])
                    + share * float(one_rows[index]["government_budget_residual"])
                )
                for index in range(len(zero_rows))
            )
        )

    critical = [row for row in constraint_rows if row["status"] == "empty_at_this_date"]
    figure, axes = plt.subplots(1, 2, figsize=(12.0, 4.5))
    axes[0].plot(shares, np.asarray(market_maxima) / MARKET_TOLERANCE, label="Housing market")
    axes[0].plot(shares, np.asarray(fiscal_maxima) / FISCAL_TOLERANCE, label="Fiscal")
    axes[0].axhline(1.0, color="black", linestyle="--", linewidth=1.0, label="Tolerance")
    axes[0].set(xlabel="Owner share at the indifferent atom", ylabel="Maximum residual / tolerance")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.2)

    for row in critical:
        values = (1.0 - shares) * float(row["share_zero_residual"]) + shares * float(row["share_one_residual"])
        axes[1].plot(
            shares,
            values / float(row["tolerance"]),
            label=f"{row['object']} {row['calendar_year']}",
        )
    axes[1].axhline(1.0, color="black", linestyle="--", linewidth=1.0)
    axes[1].axhline(-1.0, color="black", linestyle="--", linewidth=1.0)
    axes[1].set(xlabel="Owner share at the indifferent atom", ylabel="Signed residual / tolerance")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(alpha=0.2)
    figure.suptitle("Can mixing at the 2351 tenure atom clear the fixed transition path?")
    figure.tight_layout()
    figure.savefig(output_dir / "share_feasibility.png", dpi=180)
    figure.savefig(output_dir / "share_feasibility.pdf")
    plt.close(figure)

    summary = {
        "status": "complete_no_feasible_share" if intersection_lower > intersection_upper else "complete_feasible_share_interval",
        "interpretation": (
            "Endpoint and midpoint packets test whether a zero-gap atom-specific tenure share "
            "can satisfy every fixed-path market and fiscal gate."
        ),
        "manifest": str(args.manifest.resolve()),
        "manifest_sha256": file_sha256(args.manifest.resolve()),
        "fixed_path_sha256": fixed_hash,
        "packet_hashes": {
            str(share): {
                "summary": file_sha256(directory.resolve() / "summary.json"),
                "transition_path": file_sha256(directory.resolve() / "transition_path.csv"),
            }
            for share, directory in (
                (0.0, args.share_zero_dir),
                (0.5, args.share_half_dir),
                (1.0, args.share_one_dir),
            )
        },
        "affinity_maximum_absolute_errors": affinity,
        "affinity_tolerance_for_residuals": AFFINITY_TOLERANCE,
        "market_tolerance": MARKET_TOLERANCE,
        "fiscal_tolerance": FISCAL_TOLERANCE,
        "feasible_share_interval": (
            [intersection_lower, intersection_upper]
            if intersection_lower <= intersection_upper
            else None
        ),
        "empty_single_date_constraints": critical,
        "minimum_market_violation_ratio": float(np.min(np.asarray(market_maxima) / MARKET_TOLERANCE)),
        "minimum_fiscal_violation_ratio": float(np.min(np.asarray(fiscal_maxima) / FISCAL_TOLERANCE)),
    }
    write_json(output_dir / "summary.json", summary)


if __name__ == "__main__":
    main()
