#!/usr/bin/env python3
"""Build compact supplemental figures and slide tables from authenticated receipts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

LABEL = "2007 stationary reference — block0506, September 28 verified export"
BLUE = "#2C6E9B"
ORANGE = "#D57932"


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def valid_receipt(path: Path) -> dict:
    receipt = read_json(path)
    require(receipt.get("status") == "passed", f"Receipt did not pass: {path}")
    require(receipt.get("reference_label") == LABEL, f"Reference label differs: {path}")
    return receipt


def fit_value(path: Path | None, moment: str) -> float | None:
    if path is None or not path.exists():
        return None
    for row in read_csv(path):
        if row.get("moment") == moment:
            value = row.get("model", "")
            if value:
                return float(value)
    return None


def as_float(value, label: str) -> float:
    require(value is not None and value != "", f"Missing value for {label}")
    if isinstance(value, (list, tuple)):
        require(len(value) == 1, f"Expected one reported market for {label}")
        value = value[0]
    result = float(value)
    require(math.isfinite(result), f"Nonfinite value for {label}")
    return result


def set_metadata(paths: list[Path]) -> dict:
    unique = sorted({p.resolve() for p in paths} | {Path(__file__).resolve()})
    return {
        "reference_label": LABEL,
        "figures_are_supplemental": True,
        "sources": [{"path": str(path), "sha256": sha(path)} for path in unique],
    }


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def tex_escape(value: str) -> str:
    return (str(value).replace("\\", r"\textbackslash{}").replace("&", r"\&")
            .replace("%", r"\%").replace("$", r"\$").replace("#", r"\#")
            .replace("_", r"\_").replace("{", r"\{").replace("}", r"\}"))


def write_booktabs(path: Path, headers: list[str], rows: list[list[str]], caption: str,
                   label: str, note: str) -> None:
    lines = [r"\begin{table}[!htbp]", r"\centering", r"\caption{" + tex_escape(caption) + "}",
             r"\label{" + label + "}", r"\begin{tabular}{" + "l" + "r" * (len(headers)-1) + "}",
             r"\toprule", " & ".join(tex_escape(x) for x in headers) + " " + chr(92) * 2, r"\midrule"]
    for row in rows:
        lines.append(" & ".join(tex_escape(x) for x in row) + " " + chr(92) * 2)
    lines += [r"\bottomrule", r"\end{tabular}", r"\par\smallskip{\scriptsize " + tex_escape(note) + "}",
              r"\end{table}", ""]
    with path.open("x", encoding="utf-8") as stream:
        stream.write("\n".join(lines))


def validate_price_data(rows: list[dict[str, str]], completed: dict) -> dict:
    require(completed.get("status") == "passed", "Elasticity completed receipt must pass")
    require(completed.get("reference_label") == LABEL, "Elasticity reference label differs")
    fields = {"regime", "scope", "price_factor", "outcome", "value"}
    require(bool(rows) and fields.issubset(rows[0]), "Elasticity comparison schema differs")
    by_key = {}
    factors = (.98, .99, 1., 1.01, 1.02)
    for row in rows:
        regime, scope = row["regime"], row["scope"]
        require(regime in ("reference", "credit") and scope in ("impact", "cohort"),
                "Unexpected regime or scope in elasticity data")
        key = (regime, scope, row["outcome"], as_float(row["price_factor"], "price factor"))
        require(key not in by_key, f"Duplicate price response row: {key}")
        by_key[key] = as_float(row["value"], str(key))
    impact_key = "births_per_household"
    for regime in ("reference", "credit"):
        for factor in factors:
            require((regime, "impact", impact_key, factor) in by_key,
                    f"Missing immediate-birth response: {regime} {factor}")
            require((regime, "cohort", "completed_fertility", factor) in by_key,
                    f"Missing cohort-fertility response: {regime} {factor}")
    return by_key


def figure_price_response(path_png: Path, path_pdf: Path, by_key: dict) -> None:
    import matplotlib.pyplot as plt

    factors = (.98, .99, 1., 1.01, 1.02)
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), sharex=True)
    for ax, scope, outcome, title, ylabel in (
        (axes[0], "impact", "births_per_household", "Immediate births", "Percent change from own-regime baseline"),
        (axes[1], "cohort", "completed_fertility", "Completed fertility", "Percent change from own-regime baseline"),
    ):
        for regime, color, label in (("reference", BLUE, "Reference credit"), ("credit", ORANGE, "Solvency-only credit")):
            baseline = by_key[(regime, scope, outcome, 1.)]
            x = [(f - 1.) * 100 for f in factors]
            y = [(by_key[(regime, scope, outcome, f)] / baseline - 1.) * 100 for f in factors]
            ax.plot(x, y, color=color, marker="o", linewidth=1.8, markersize=4, label=label)
        ax.axhline(0, color="#777777", linewidth=.75)
        ax.axvline(0, color="#777777", linewidth=.75)
        ax.set_title(title)
        ax.set_xlabel("Housing asset price change (%)")
        ax.set_ylabel(ylabel)
        ax.set_xlim(-2, 2)
        ax.grid(axis="y", color="#dddddd", linewidth=.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].legend(frameon=False, loc="best")
    fig.suptitle("Housing costs and fertility", y=.99, fontsize=12)
    fig.text(.5, .035, "Supplemental figure · Prescribed price; not GE. Price and mapped rent move together.",
             ha="center", fontsize=7)
    fig.text(.5, .012, LABEL, ha="center", fontsize=6.7)
    fig.tight_layout(rect=(0, .10, 1, .94))
    fig.savefig(path_png, dpi=220, bbox_inches="tight")
    fig.savefig(path_pdf, bbox_inches="tight")
    plt.close(fig)


def figure_supply(path_png: Path, path_pdf: Path, rows: list[dict[str, str]]) -> None:
    import matplotlib.pyplot as plt

    by_case = {row["case"]: row for row in rows}
    need = ("frozen_baseline", "credit_fixed_baseline_stock", "credit_original_elastic_stock")
    require(all(case in by_case for case in need), "Supply comparison is missing an endpoint")
    base = by_case["frozen_baseline"]
    cases = (("credit_fixed_baseline_stock", "Credit, fixed baseline stock"),
             ("credit_original_elastic_stock", "Credit, original elastic supply"))
    metrics = ("Household population", "Physical housing supply", "House price / implied rent")
    endpoints = [
        [100*(as_float(by_case[c]["population"], "population") / as_float(base["population"], "baseline population")-1)
         for c, _ in cases],
        [100*(as_float(by_case[c]["absolute_supply"], "supply") / as_float(base["absolute_supply"], "baseline supply")-1)
         for c, _ in cases],
        [100*(as_float(by_case[c]["price"], "price") / as_float(base["price"], "baseline price")-1)
         for c, _ in cases],
    ]
    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    centers = range(len(metrics)); width = .34
    for i, (case, label) in enumerate(cases):
        vals = [column[i] for column in endpoints]
        positions = [center + (i-.5)*width for center in centers]
        bars = ax.bar(positions, vals, width, color=(BLUE, ORANGE)[i], label=label)
        for bar, value in zip(bars, vals):
            ax.annotate(f"{value:.3f}", (bar.get_x() + bar.get_width()/2, value),
                        xytext=(0, 3 if value >= 0 else -9), textcoords="offset points",
                        ha="center", va="bottom" if value >= 0 else "top", fontsize=7)
    ax.axhline(0, color="#777777", linewidth=.8)
    ax.set_xticks(list(centers))
    ax.set_xticklabels(metrics)
    ax.set_ylabel("Change from frozen baseline (%)")
    ax.margins(y=.16)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=.5)
    ax.legend(frameon=False, loc="best")
    fig.suptitle("Borrowing and housing supply", y=.98, fontsize=12)
    fig.text(.5, .035, "Supplemental figure · Stationary endpoints, not a transition; population in household units.",
             ha="center", fontsize=7)
    fig.text(.5, .012, LABEL, ha="center", fontsize=6.7)
    fig.tight_layout(rect=(0, .10, 1, .92))
    fig.savefig(path_png, dpi=220, bbox_inches="tight")
    fig.savefig(path_pdf, bbox_inches="tight")
    plt.close(fig)


def build_borrowing(receipts: list[dict], fits: list[Path | None], out: Path) -> None:
    names = ["Reference", "Matched grid", "Credit (fixed prices)", "Credit (GE)"]
    q0 = as_float(receipts[0].get("price"), "reference price")
    cols = ["frozen_reference", "matched_grid_baseline_fixed_prices",
            "solvency_only_fixed_prices", "solvency_only_ge"]
    # Index values are computed from each receipt; only the GE column is an endogenous price.
    table_values = {
        "Price index (reference = 100)": [100*as_float(r.get("price"), "price")/q0 for r in receipts],
        "Stationary household population index (reference = 100)": [100., None, None,
            100*as_float(receipts[3].get("closure", {}).get("population_scale"), "GE population")],
        "Completed fertility": [as_float(r.get("completed_fertility"), "completed fertility") for r in receipts],
        "Mean age at first birth (years)": [],
        "Homeownership (%)": [100*as_float(r.get("cohort_summary", {}).get("ownership_rate"), "ownership") for r in receipts],
        "Housing services per household (rooms)": [as_float(r.get("cohort_summary", {}).get("rooms_per_household"), "rooms") for r in receipts],
    }
    for i, (receipt, fit) in enumerate(zip(receipts, fits)):
        value = receipt.get("mean_first_birth_age")
        if value is None:
            value = fit_value(fit, "nchs_mean_age")
        require(value is not None, f"Missing nchs_mean_age for {names[i]}; provide its fit CSV")
        table_values["Mean age at first birth (years)"].append(float(value))
    csv_rows = []
    tex_rows = []
    for label, values in table_values.items():
        csv_rows.append({"outcome": label, **{col: "" if value is None else value
                                               for col, value in zip(cols, values)}})
        tex_rows.append([label] + ["--" if value is None else f"{value:.3f}" for value in values])
    write_csv(out / "borrowing_comparison.csv", ["outcome"] + cols, csv_rows)
    write_booktabs(out / "borrowing_comparison.tex", ["Outcome"] + names, tex_rows,
                   "Supplemental borrowing-regime comparison", "tab:borrowing-comparison",
                   "Matched-grid and partial-equilibrium columns report conditional-cohort outcomes. "
                   "Stationary population is defined only for the reference and general-equilibrium columns. " + LABEL)


def build_elasticity_table(rows: list[dict[str, str]], out: Path) -> None:
    required = {"regime", "scope", "outcome", "step", "central_log_elasticity"}
    require(bool(rows) and required.issubset(rows[0]), "Elasticities schema differs")
    wanted = {"first_births": "First births", "second_births": "Second births",
              "births_per_household": "Total births", "completed_fertility": "Completed fertility"}
    selected = [r for r in rows if r["outcome"] in wanted]
    lookup = {}
    for row in selected:
        scope = row["scope"]
        if row["outcome"] == "completed_fertility":
            require(scope == "cohort", "Completed fertility must be cohort scope")
        key = (row["regime"], scope, row["outcome"], as_float(row["step"], "step"))
        lookup[key] = row.get("central_log_elasticity", "")
    cols = ["Reference ±1%", "Reference ±2%", "Credit ±1%", "Credit ±2%"]
    csv_fields = ["outcome", "scope"] + [f"r{i}" for i in range(4)]
    csv_rows, tex_rows = [], []
    row_specs = (("births_per_household", "Total births (impact)", "impact"),
                 ("first_births", "First births (impact)", "impact"),
                 ("second_births", "Second births (impact)", "impact"),
                 ("completed_fertility", "Completed fertility (cohort)", "cohort"))
    for outcome, label, scope in row_specs:
        values = []
        for regime in ("reference", "credit"):
            for step in (.01, .02):
                row = lookup.get((regime, scope, outcome, step))
                require(row not in (None, ""), f"Missing elasticity {regime}/{scope}/{outcome}/{step}")
                values.append(float(row))
        csv_rows.append({"outcome": label, "scope": scope,
                         **{f"r{i}": v for i, v in enumerate(values)}})
        tex_rows.append([label] + ["--" if v is None else f"{v:.3f}" for v in values])
    write_csv(out / "local_elasticities.csv", csv_fields, csv_rows)
    write_booktabs(out / "local_elasticities.tex", ["Response"] + cols, tex_rows,
                   "Supplemental local price elasticities; impact births and cohort completed fertility",
                   "tab:local-elasticities",
                   "Impact responses use the frozen inherited distribution; completed fertility is cohort-based. " + LABEL)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-receipt", type=Path, required=True)
    parser.add_argument("--grid-receipt", type=Path, required=True)
    parser.add_argument("--credit-receipt", type=Path, required=True)
    parser.add_argument("--ge-receipt", type=Path, required=True)
    parser.add_argument("--reference-fit", type=Path, required=True)
    parser.add_argument("--grid-fit", type=Path, required=True)
    parser.add_argument("--credit-fit", type=Path, required=True)
    parser.add_argument("--ge-fit", type=Path, required=True)
    parser.add_argument("--price-comparison", type=Path)
    parser.add_argument("--elasticities", type=Path)
    parser.add_argument("--elasticity-completed", type=Path)
    parser.add_argument("--supply-comparison", type=Path)
    parser.add_argument("--supply-verification", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output = args.output.resolve()
    require(not args.output.exists(), "Refusing to overwrite an existing output directory")
    receipt_paths = [args.reference_receipt, args.grid_receipt, args.credit_receipt, args.ge_receipt]
    receipts = [valid_receipt(p) for p in receipt_paths]
    require(receipts[0].get("case") == "control" and receipts[1].get("case") == "grid_control" and
            receipts[2].get("case") == "credit" and receipts[3].get("case") == "selected_repeat",
            "Borrowing comparison receipt roles differ")
    for path in (args.reference_fit, args.grid_fit, args.credit_fit, args.ge_fit):
        if path is not None:
            require(path.exists(), f"Fit file not found: {path}")
    price_inputs = (args.price_comparison, args.elasticities, args.elasticity_completed)
    have_price_inputs = all(path is not None for path in price_inputs)
    require(have_price_inputs or all(path is None for path in price_inputs),
            "Pass all three elasticity inputs together, or omit all three")
    price_by_key = None
    elasticity_rows = None
    if have_price_inputs:
        completed = read_json(args.elasticity_completed)
        require(completed.get("comparison_sha256") == sha(args.price_comparison),
                "Price comparison SHA-256 differs from the passed completion receipt")
        require(completed.get("elasticities_sha256") == sha(args.elasticities),
                "Elasticities SHA-256 differs from the passed completion receipt")
        price_by_key = validate_price_data(read_csv(args.price_comparison), completed)
        elasticity_rows = read_csv(args.elasticities)
    supply_inputs = (args.supply_comparison, args.supply_verification)
    have_supply_inputs = all(path is not None for path in supply_inputs)
    require(have_supply_inputs or all(path is None for path in supply_inputs),
            "Pass both supply comparison and verification, or omit both")
    supply_rows = None
    if have_supply_inputs:
        supply_verification = read_json(args.supply_verification)
        require(supply_verification.get("status") == "passed" and supply_verification.get("model_solves") == 0,
                "Supply verification must pass with zero model solves")
        require(supply_verification.get("reference_label") == LABEL, "Supply reference label differs")
        supply_rows = read_csv(args.supply_comparison)
        require(bool(supply_rows) and {"case", "price", "absolute_supply", "population"}.issubset(supply_rows[0]),
                "Supply comparison schema differs")
    args.output.mkdir(parents=True, exist_ok=False)
    if have_supply_inputs:
        figure_supply(args.output/"supplemental_housing_supply.png",
                      args.output/"supplemental_housing_supply.pdf", supply_rows)
    build_borrowing(receipts, [args.reference_fit, args.grid_fit, args.credit_fit, args.ge_fit], args.output)
    if have_price_inputs:
        figure_price_response(args.output/"supplemental_price_response.png",
                              args.output/"supplemental_price_response.pdf", price_by_key)
        build_elasticity_table(elasticity_rows, args.output)
    dependencies = receipt_paths + [p for p in (args.reference_fit, args.grid_fit, args.credit_fit, args.ge_fit) if p]
    dependencies += [p for p in price_inputs if p is not None]
    dependencies += [p for p in supply_inputs if p is not None]
    runner = Path(__file__).with_name("run_existing_render.sh")
    if runner.exists():
        dependencies.append(runner)
    manifest = set_metadata(dependencies)
    manifest.update({"tables": ["borrowing_comparison.csv", "borrowing_comparison.tex"]})
    if have_supply_inputs:
        manifest["housing_supply"] = {"png": "supplemental_housing_supply.png", "pdf": "supplemental_housing_supply.pdf"}
    if have_price_inputs:
        manifest["price_response"] = {"png": "supplemental_price_response.png", "pdf": "supplemental_price_response.pdf"}
        manifest["tables"] += ["local_elasticities.csv", "local_elasticities.tex"]
    with (args.output/"manifest.json").open("x", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, ensure_ascii=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
