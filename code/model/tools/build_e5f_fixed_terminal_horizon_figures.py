#!/usr/bin/env python3
"""Plot the supplemental fixed-terminal 10-period original-queue transition.

This reader consumes the immutable JSON artifacts written by the experiment
controller.  It performs no solving, interpolation, convergence assessment,
or fertility remeasurement.  In particular, the fertility series is the
existing ``period_tfr_topcode_adjusted`` diagnostic.
"""

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
import numpy as np


FERTILITY_FIELD = "period_tfr_topcode_adjusted"
QUANTITY_FIELDS = (
    (FERTILITY_FIELD, "Fertility rate (top-code adjusted period diagnostic)"),
    ("adult_population", "Household population"),
    ("asset_price", "House price"),
    ("renter_price", "Rent"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, required=True)
    return parser.parse_args()


def load(path: Path, required: bool = True) -> Any:
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def nested(record: Any, key: str) -> Any:
    """Find a named quantity in a row, diagnostics, or endpoint payload."""
    if not isinstance(record, dict):
        return None
    if key in record:
        return record[key]
    for child_key in ("quantities", "fertility", "diagnostics", "payload", "endpoint_reference", "final"):
        child = record.get(child_key)
        if isinstance(child, dict) and key in child:
            return child[key]
    return None


def finite_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def row_status(row: dict[str, Any]) -> str:
    raw = row.get("status", row.get("root_status", row.get("mapping_status", "complete")))
    if isinstance(raw, dict):
        raw = raw.get("label", raw.get("status", "incomplete"))
    return str(raw)


def reference_row(reference: Any, terminal: bool = False) -> dict[str, Any]:
    if not isinstance(reference, dict):
        return {}
    # stationary_reference commonly stores initial/terminal side by side.
    side = reference.get("terminal" if terminal else "initial")
    if isinstance(side, dict):
        return side
    return reference


def source_hash(case: Path) -> str:
    digest = hashlib.sha256()
    for name in ("rows.json", "fertility.json", "root_receipt.json", "terminal_distance.json", "stationary_reference.json"):
        path = case / name
        if path.exists():
            digest.update(name.encode())
            digest.update(path.read_bytes())
    return digest.hexdigest()


def safe_limits(values: list[float | None]) -> tuple[float, float] | None:
    data = np.asarray([v for v in values if v is not None], dtype=float)
    if data.size == 0:
        return None
    lo, hi = float(np.min(data)), float(np.max(data))
    span = hi - lo
    scale = max(abs(lo), abs(hi), 1.0)
    # Avoid turning machine-level roundoff into a visually large movement.
    if span <= 1e-5 * scale:
        pad = max(1e-3 * scale, 1e-6)
    else:
        pad = 0.08 * span
    return lo - pad, hi + pad


def main() -> None:
    args = parse_args()
    case = args.case_dir.resolve()
    rows_raw = load(case / "rows.json")
    fertility_raw = load(case / "fertility.json")
    receipt = load(case / "root_receipt.json", required=False) or {}
    distance_raw = load(case / "terminal_distance.json", required=False)
    references = load(case / "stationary_reference.json", required=False) or {}
    if not isinstance(rows_raw, list) or not rows_raw:
        raise ValueError("rows.json must be a nonempty list")
    if not isinstance(fertility_raw, list):
        raise ValueError("fertility.json must be a list of calendar-year diagnostics")
    if any(not isinstance(row, dict) for row in rows_raw):
        raise ValueError("Every quantity row must be a record")
    rows = [dict(row) for row in rows_raw]
    fertility_by_year = {
        float(item["calendar_year"]): item
        for item in fertility_raw
        if isinstance(item, dict) and "calendar_year" in item
    }
    years = [finite_float(row.get("calendar_year")) for row in rows]
    if any(year is None for year in years):
        raise ValueError("each row requires a finite calendar_year")
    years_f = [float(year) for year in years]
    if len(fertility_by_year) != len(fertility_raw) or set(fertility_by_year) != set(years_f):
        raise ValueError("Fertility and quantity dates must match uniquely")
    origin = years_f[0]
    x = [year - origin for year in years_f]
    initial = reference_row(references, terminal=False)
    terminal = reference_row(references, terminal=True)
    verified_endpoint = bool(references.get("stationary_endpoint_verified", False))
    root_converged = bool(receipt.get("finite_horizon_market_fiscal_converged", receipt.get("converged", False)))
    status_text = str(receipt.get("status", receipt.get("status_label", "INCOMPLETE / status unspecified")))
    complete_rows = [row_status(row).lower() in {"complete", "valid", "accepted", "converged"} for row in rows]

    series: dict[str, list[float | None]] = {}
    for field, _ in QUANTITY_FIELDS:
        values = []
        for row in rows:
            if field == FERTILITY_FIELD:
                diag = fertility_by_year.get(float(row["calendar_year"]), {})
                values.append(finite_float(nested(diag, FERTILITY_FIELD)))
            else:
                values.append(finite_float(nested(row, field)))
        series[field] = values
        if any(value is None for value in values):
            raise ValueError("Missing or nonfinite required series: " + field)

    # Optional endpoint-distance series can be row-level or a list of dated records.
    distance_by_year: dict[float, float | None] = {}
    if isinstance(distance_raw, list):
        for item in distance_raw:
            if isinstance(item, dict) and "calendar_year" in item:
                distance_by_year[float(item["calendar_year"])] = finite_float(item.get("distance", item.get("endpoint_distance")))
    elif isinstance(distance_raw, dict):
        values = distance_raw.get("rows", distance_raw.get("diagnostics"))
        if isinstance(values, list):
            for item in values:
                if isinstance(item, dict) and "calendar_year" in item:
                    distance_by_year[float(item["calendar_year"])] = finite_float(item.get("distance", item.get("endpoint_distance")))
        elif finite_float(distance_raw.get("distance", distance_raw.get("endpoint_distance"))) is not None:
            distance_by_year[years_f[-1]] = finite_float(distance_raw.get("distance", distance_raw.get("endpoint_distance")))
    distance_values = [distance_by_year.get(year, finite_float(nested(row, "endpoint_distance"))) for year, row in zip(years_f, rows, strict=True)]
    if any(v is not None for v in distance_values):
        series["endpoint_distance"] = distance_values

    fields = list(QUANTITY_FIELDS)
    if "endpoint_distance" in series:
        fields.append(("endpoint_distance", "Distance to stationary endpoint"))
    terminal_gaps = {}
    if isinstance(distance_raw, dict):
        for key, label in [("population_relative_gap", "Population"),
                           ("distribution_relative_l1", "Full distribution (L1)"),
                           ("queue_relative_max", "Entrant queue (max)")]:
            value = finite_float(distance_raw.get(key))
            if value is not None:
                terminal_gaps[label] = 100.0 * value
    panel_count = len(fields) + bool(terminal_gaps)
    fig, axes = plt.subplots(2, 3 if panel_count > 4 else 2, figsize=(13.0, 7.2), squeeze=False)
    axes_flat = axes.ravel()
    qa_rows: list[dict[str, Any]] = []
    for i, (field, label) in enumerate(fields):
        ax = axes_flat[i]
        ys = series[field]
        plotted = np.asarray([np.nan if value is None else value for value in ys], dtype=float)
        ax.plot(x, plotted, marker="o", linewidth=1.8, color="tab:blue", label="10-period transition")
        ref = finite_float(nested(initial, field))
        if ref is not None:
            ax.axhline(ref, color="0.35", linestyle="--", linewidth=1.1, label="Initial SS reference")
        tref = finite_float(nested(terminal, field))
        if verified_endpoint and tref is not None:
            ax.axhline(tref, color="tab:red", linestyle=":", linewidth=1.2, label="Verified stationary endpoint")
        missing = [xv for xv, yv in zip(x, ys, strict=True) if yv is None]
        if missing:
            ax.scatter(missing, [0.0] * len(missing), marker="x", color="crimson", zorder=4, label="Missing/incomplete")
        ax.set_title(label)
        ax.set_xlabel("Years since transition start")
        ax.grid(alpha=0.25)
        limits = safe_limits(ys + ([ref] if ref is not None else []) + ([tref] if verified_endpoint and tref is not None else []))
        if limits:
            ax.set_ylim(*limits)
        ax.legend(fontsize=7, frameon=False)
    if terminal_gaps:
        ax = axes_flat[len(fields)]
        ax.barh(list(terminal_gaps), list(terminal_gaps.values()), color="tab:orange")
        ax.set_title("Deviation from new SS after 40 years")
        ax.set_xlabel("Percent of terminal reference")
        ax.grid(axis="x", alpha=0.25)
    for ax in axes_flat[panel_count:]:
        ax.remove()
    fig.suptitle("Supplemental original-queue transition: fixed terminal horizon", fontsize=13)
    incomplete = not root_converged or not all(complete_rows)
    warning = "UNCONVERGED ROOT — diagnostic path only" if incomplete else "Finite root converged; assess endpoint gaps. Fertility alone does not establish convergence."
    fig.text(0.5, 0.012, f"{warning}  |  root status: {status_text}", ha="center", fontsize=8, color="crimson" if incomplete else "0.25")
    fig.tight_layout(rect=(0, 0.035, 1, 0.94))
    png, pdf = case / "fixed_terminal_horizon.png", case / "fixed_terminal_horizon.pdf"
    fig.savefig(png, dpi=180)
    fig.savefig(pdf)
    plt.close(fig)

    for idx, row in enumerate(rows):
        out = {"calendar_year": years_f[idx], "years_since_start": x[idx], "status": row_status(row)}
        out.update({field: series[field][idx] for field, _ in fields})
        qa_rows.append(out)
    with (case / "fixed_terminal_horizon_raw.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(qa_rows[0]))
        writer.writeheader()
        writer.writerows(qa_rows)
    qa = {
        "sourceSHA256": source_hash(case),
        "fertility_field": FERTILITY_FIELD,
        "fertility_definition": "period_fertility_diagnostics['period_tfr_topcode_adjusted']",
        "number_rows": len(rows),
        "complete_rows": int(sum(complete_rows)),
        "incomplete_rows": int(len(rows) - sum(complete_rows)),
        "verified_stationary_endpoint_reference_used": bool(verified_endpoint),
        "endpoint_distance_available": bool(terminal_gaps) or "endpoint_distance" in series,
        "terminal_gaps_percent": terminal_gaps,
        "finite_root_converged": root_converged,
        "graph_paths": {"png": str(png), "pdf": str(pdf)},
        "raw_csv": str(case / "fixed_terminal_horizon_raw.csv"),
        "series": series,
    }
    (case / "fixed_terminal_horizon_qa.json").write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
