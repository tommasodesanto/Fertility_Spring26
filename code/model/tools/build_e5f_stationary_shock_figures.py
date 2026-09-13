#!/usr/bin/env python3
"""Build a six-panel diagnostic IRF from native stationary-shock rows.

The utility deliberately treats its input as a plotting contract: it does not
solve, interpolate, or assess convergence of the supplied equilibrium rows.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


FIELDS = (
    ("birth_children", "Births per four-year period"),
    ("adult_population", "Household heads"),
    ("housing_demand/adult_population", "Rooms per household"),
    ("asset_price", "House price"),
    ("renter_price", "Rent"),
    ("pension_period", "Pension per retiree"),
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case-dir", type=Path, required=True)
    p.add_argument("--fertility-rate", action="store_true",
                   help="Use native age-specific fertility diagnostics, never a birth-count proxy.")
    p.add_argument("--include-pre-shock", action="store_true",
                   help="Show the initial steady state and the impact jump explicitly.")
    return p.parse_args()


def read_json(path: Path, required: bool = True) -> Any:
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def source_hash(case: Path) -> str:
    h = hashlib.sha256()
    for name in ("rows.json", "irf_contract.json", "initial_reference.json", "terminal_reference.json",
                 "fertility.json", "stationary_reference.json"):
        path = case / name
        if path.exists():
            h.update(name.encode())
            h.update(path.read_bytes())
    return h.hexdigest()


def value(row: dict[str, Any], field: str) -> float | None:
    if field == "housing_demand/adult_population":
        if field in row:
            return None if row[field] is None else float(row[field])
        den = row.get("adult_population")
        if den is None or float(den) <= 0:
            raise ValueError("adult_population must be positive for every plotted housing ratio")
        num = row.get("housing_demand")
        return None if num is None else float(num) / float(den)
    if field == "pension_period":
        raw = row.get("pension_period", row.get("pension"))
    else:
        raw = row.get(field)
    return None if raw is None else float(raw)


def x_values(rows: list[dict[str, Any]], contract: dict[str, Any]) -> list[float]:
    if all("period" in row for row in rows):
        return [4.0 * float(row["period"]) for row in rows]
    years = [float(row["calendar_year"]) for row in rows]
    origin = contract.get("shock_calendar_year", years[0])
    return [year - float(origin) for year in years]


def main() -> None:
    args = parse_args()
    case = args.case_dir.resolve()
    rows = read_json(case / "rows.json")
    contract = read_json(case / "irf_contract.json")
    if not isinstance(rows, list) or not rows:
        raise ValueError("rows.json must contain a nonempty list")
    if not isinstance(contract, dict):
        raise ValueError("irf_contract.json must contain an object")
    rows = [dict(row) for row in rows]
    xs = x_values(rows, contract)
    initial = read_json(case / "initial_reference.json", required=False) or {}
    terminal = read_json(case / "terminal_reference.json", required=False)
    terminal_ok = isinstance(terminal, dict) and terminal.get("verified") is True
    fields = list(FIELDS)
    if args.fertility_rate:
        diagnostics = read_json(case / "fertility.json")
        references = read_json(case / "stationary_reference.json")
        by_year = {float(item["calendar_year"]): item for item in diagnostics}
        if (len(by_year) != len(diagnostics)
                or len(rows) != len(diagnostics)
                or set(by_year) != {float(row["calendar_year"]) for row in rows}):
            raise ValueError("Native fertility observations must match every dated row uniquely")
        for row in rows:
            d = by_year[float(row["calendar_year"])]
            flows = np.asarray(d["birth_flow_topcode_adjusted"], dtype=float)
            masses = np.asarray(d["age_mass"], dtype=float)
            if (flows.shape != masses.shape or not np.isfinite(flows).all()
                    or not np.isfinite(masses).all() or np.any(masses < 0)):
                raise ValueError("Invalid native fertility age cells")
            rate = float(np.divide(flows, masses, out=np.zeros_like(flows),
                                   where=masses > 1e-15).sum())
            if (abs(float(flows.sum()) - float(row["birth_children_topcode_adjusted"])) > 2e-10
                    or abs(float(masses.sum()) - float(row["adult_population"])) > 2e-10
                    or abs(rate - float(d["period_tfr_topcode_adjusted"])) > 2e-10):
                raise ValueError("Native fertility measurement does not reproduce its dated row")
            row["period_tfr_topcode_adjusted"] = rate
        initial = dict(references["initial"]["quantities"])
        terminal = dict(references["terminal"]["quantities"])
        for target, source in ((initial, references["initial"]),
                               (terminal, references["terminal"])):
            target["period_tfr_topcode_adjusted"] = float(source["fertility"]["period_tfr_topcode_adjusted"])
        terminal_ok = references.get("stationary_endpoint_verified") is True
        if not terminal_ok:
            raise ValueError("Fertility figure requires a verified stationary reference")
        fields[0] = ("period_tfr_topcode_adjusted", "Period fertility rate")

    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.4), constrained_layout=False)
    axes = axes.ravel()
    qa_series: dict[str, dict[str, list[float | None]]] = {}
    for ax, (field, label) in zip(axes, fields, strict=True):
        ys = [value(row, field) for row in rows]
        qa_series[field] = {"x": xs, "y": ys}
        plotted = np.array([np.nan if y is None else y for y in ys], dtype=float)
        ref = value(initial, field) if initial else None
        if args.include_pre_shock:
            if ref is None:
                raise ValueError("Missing pre-shock reference: " + field)
            ax.plot([-4.0, 0.0] + xs, [ref, ref] + plotted.tolist(),
                    marker="o", linewidth=1.8, label="Path")
        else:
            ax.plot(xs, plotted, marker="o", linewidth=1.8, label="IRF")
        if ref is not None:
            ax.axhline(ref, color="0.35", linestyle="--", linewidth=1.1, label="Initial reference")
        if terminal_ok:
            tref = value(terminal, field)
            if tref is not None:
                ax.plot([xs[-1]], [tref], marker="D", color="tab:red", markersize=4.5, label="Terminal equilibrium")
        missing = [x for x, y in zip(xs, ys, strict=True) if y is None]
        if missing:
            ax.scatter(missing, [0.0] * len(missing), marker="x", color="crimson", zorder=4, label="Missing data")
            ax.text(0.02, 0.04, "Missing data marked ×", transform=ax.transAxes, color="crimson", fontsize=8)
        ax.set_title(label)
        ax.set_xlabel("Years since shock")
        ax.grid(alpha=0.25)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(fontsize=7, frameon=False)
    fig.suptitle(str(contract.get("label", "Stationary preference shock IRF")), fontsize=13)
    status = str(contract.get("status_label", "status unspecified"))
    description=str(contract.get("shock_description", "Permanent preference shock."))
    fig.text(0.5, 0.012, "Diagnostic under original household-entry rule. " + description + " " + status,
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.035, 1, 0.94))
    stem = "irf_fertility" if args.fertility_rate else ("irf_preshock" if args.include_pre_shock else "irf")
    png = case / (stem + ".png")
    pdf = case / (stem + ".pdf")
    fig.savefig(png, dpi=180)
    fig.savefig(pdf)
    plt.close(fig)

    qa = {
        "sourceSHA256": source_hash(case),
        "number_rows": len(rows),
        "series": qa_series,
        "graph_paths": {"png": str(png), "pdf": str(pdf)},
        "contract_label": contract.get("label"),
        "status_label": contract.get("status_label"),
        "terminal_reference_used": terminal_ok,
        "native_period_fertility_used": args.fertility_rate,
        "pre_shock_point_included": args.include_pre_shock,
    }
    qa_name = "qa.json" if stem == "irf" else stem + "_qa.json"
    (case / qa_name).write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
