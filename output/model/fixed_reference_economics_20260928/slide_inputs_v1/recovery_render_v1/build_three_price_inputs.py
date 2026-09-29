#!/usr/bin/env python3
"""Render the approved three-price prescribed-price comparison and elasticity table."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

LABEL = "2007 stationary reference — block0506, September 28 verified export"
FACTORS = (0.99, 1.0, 1.01)
REGIMES = {"reference": "Baseline borrowing limits", "credit": "Lifetime repayment only"}
BLUE = "#2C6E9B"
ORANGE = "#D57932"
TABLE_ROWS = (
    ("births_per_household", "Total births (impact)", "impact"),
    ("first_births", "First births (impact)", "impact"),
    ("second_births", "Second births (impact)", "impact"),
    ("completed_fertility", "Completed fertility (cohort)", "cohort"),
)


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def number(value, description: str) -> float:
    require(value is not None and value != "", f"Missing value: {description}")
    result = float(value)
    require(math.isfinite(result), f"Nonfinite value: {description}")
    return result


def validate_completion(completed: dict, comparison_path: Path, elasticities_path: Path) -> None:
    require(completed.get("status") == "passed", "Three-price completion status must be passed")
    require(completed.get("complete_three_price") is True, "Three-price completion flag is missing")
    require(completed.get("complete_five_price") is False, "Five-price completion must be explicitly false")
    require(completed.get("reference_label") == LABEL, "Completion reference label differs")
    require(completed.get("factors") == [0.99, 1.0, 1.01], "Completion factors must be exactly .99, 1, 1.01")
    require(completed.get("step_sizes") == [0.01], "Completion step size must be exactly .01")
    require(completed.get("comparison_sha256") == sha(comparison_path),
            "comparison.csv SHA-256 differs from completion receipt")
    require(completed.get("elasticities_sha256") == sha(elasticities_path),
            "elasticities.csv SHA-256 differs from completion receipt")


def validate_comparison(rows: list[dict[str, str]]) -> dict:
    required = {"regime", "scope", "price_factor", "outcome", "value", "unit",
                "prescribed_price", "mapped_rent"}
    require(bool(rows) and required.issubset(rows[0]), "comparison.csv schema differs")
    by_key = {}
    factors = set(FACTORS)
    for row in rows:
        regime, scope = row["regime"], row["scope"]
        require(regime in REGIMES, f"Unexpected regime: {regime}")
        require(scope in ("impact", "cohort"), f"Unexpected response scope: {scope}")
        factor = number(row["price_factor"], "price_factor")
        require(factor in factors, f"Unexpected price factor {factor}; this renderer supports only .99, 1, 1.01")
        key = (regime, scope, row["outcome"], factor)
        require(key not in by_key, f"Duplicate comparison row: {key}")
        by_key[key] = number(row["value"], str(key))
        number(row["prescribed_price"], "prescribed_price")
        number(row["mapped_rent"], "mapped_rent")
    regimes = {row["regime"] for row in rows}
    scopes = {row["scope"] for row in rows}
    require(regimes == set(REGIMES), "Both reference and credit regimes are required")
    require(scopes == {"impact", "cohort"}, "Both impact and cohort scopes are required")
    for regime in REGIMES:
        for scope in ("impact", "cohort"):
            present = {k[3] for k in by_key if k[0] == regime and k[1] == scope}
            require(present == factors, f"Incomplete three-price coverage for {regime}/{scope}")
        for factor in FACTORS:
            require((regime, "impact", "births_per_household", factor) in by_key,
                    f"Missing immediate births at {regime}/{factor}")
            require((regime, "cohort", "completed_fertility", factor) in by_key,
                    f"Missing cohort completed fertility at {regime}/{factor}")
    return by_key


def validate_elasticities(rows: list[dict[str, str]]) -> dict:
    required = {"regime", "scope", "outcome", "step", "central_log_elasticity",
                "lower_one_sided_log_elasticity", "upper_one_sided_log_elasticity"}
    require(bool(rows) and required.issubset(rows[0]), "elasticities.csv schema differs")
    lookup = {}
    for row in rows:
        regime, scope = row["regime"], row["scope"]
        require(regime in REGIMES, f"Unexpected elasticity regime: {regime}")
        require(scope in ("impact", "cohort"), f"Unexpected elasticity scope: {scope}")
        step = number(row["step"], "elasticity step")
        require(step == 0.01, f"Unexpected elasticity step: {step}")
        key = (regime, scope, row["outcome"], step)
        require(key not in lookup, f"Duplicate elasticity row: {key}")
        lookup[key] = {
            "central": number(row["central_log_elasticity"], str(key) + "/central"),
            "lower": number(row["lower_one_sided_log_elasticity"], str(key) + "/lower"),
            "upper": number(row["upper_one_sided_log_elasticity"], str(key) + "/upper"),
        }
    for regime in REGIMES:
        for outcome, _, scope in TABLE_ROWS:
            require((regime, scope, outcome, 0.01) in lookup,
                    f"Missing required elasticity: {regime}/{scope}/{outcome}")
    return lookup


def verify_elasticities_from_comparison(elasticities: dict, comparison: dict) -> None:
    def log_slope(left: float, right: float, factor_left: float, factor_right: float) -> float:
        require(left > 0 and right > 0, "Displayed elasticities require positive comparison outcomes")
        return (math.log(right) - math.log(left)) / (math.log(factor_right) - math.log(factor_left))

    for regime in REGIMES:
        for outcome, _, scope in TABLE_ROWS:
            key = (regime, scope, outcome, 0.01)
            for factor in FACTORS:
                require((regime, scope, outcome, factor) in comparison,
                        f"Missing comparison value for {regime}/{scope}/{outcome}/{factor}")
            observed = elasticities[key]
            expected = {
                "central": log_slope(comparison[(regime, scope, outcome, 0.99)],
                                     comparison[(regime, scope, outcome, 1.01)], 0.99, 1.01),
                "lower": log_slope(comparison[(regime, scope, outcome, 0.99)],
                                   comparison[(regime, scope, outcome, 1.0)], 0.99, 1.0),
                "upper": log_slope(comparison[(regime, scope, outcome, 1.0)],
                                   comparison[(regime, scope, outcome, 1.01)], 1.0, 1.01),
            }
            for statistic, value in expected.items():
                actual = observed[statistic]
                require(abs(actual-value) <= 1e-10 * max(1.0, abs(value)),
                        f"{statistic} elasticity disagrees with comparison data for {key}")


def tex_escape(value: str) -> str:
    return (str(value).replace("\\", r"\textbackslash{}").replace("&", r"\&")
            .replace("%", r"\%").replace("$", r"\$").replace("#", r"\#")
            .replace("_", r"\_").replace("{", r"\{").replace("}", r"\}"))


def write_table(out: Path, elasticities: dict, test_only: bool) -> None:
    regimes = list(REGIMES)
    rows_csv, rows_tex = [], []
    for outcome, response, scope in TABLE_ROWS:
        values = [elasticities[(regime, scope, outcome, 0.01)]["central"] for regime in regimes]
        rows_csv.append({"response": response, "scope": scope,
                         "baseline_borrowing_limits": values[0],
                         "lifetime_repayment_only": values[1], "reference_label": LABEL,
                         "synthetic_test_only": test_only})
        rows_tex.append([response, *(f"{value:.3f}" for value in values)])
    with (out / "local_elasticities.csv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows_csv[0]))
        writer.writeheader()
        writer.writerows(rows_csv)
    caption = "Supplemental central 1% price elasticities" + (" (synthetic test only)" if test_only else "")
    note = (r"All preference parameters, including $\psi_{\mathrm{child}}$, are fixed. Impact uses the same inherited households. "
            r"The cohort distribution has unit household mass; completed fertility is free to change. "
            "Not general equilibrium or a transition. " + LABEL)
    spec = "lrr"
    lines = [r"\begin{table}[!htbp]", r"\centering", r"\caption{" + tex_escape(caption) + "}",
             r"\label{tab:recovery-price-elasticities}", r"\begin{tabular}{" + spec + "}",
             r"\toprule", "Response & Baseline limits & Lifetime repayment only " + chr(92)*2,
             r"\midrule"]
    lines.extend(" & ".join(tex_escape(cell) for cell in row) + " " + chr(92)*2 for row in rows_tex)
    lines += [r"\bottomrule", r"\end{tabular}", r"\par\smallskip{\scriptsize " + note + "}",
              r"\end{table}", ""]
    with (out / "local_elasticities.tex").open("x", encoding="utf-8") as stream:
        stream.write("\n".join(lines))


def render_figure(out: Path, values: dict, test_only: bool) -> None:
    import matplotlib.pyplot as plt

    factors = FACTORS
    x = [(factor - 1.0) * 100 for factor in factors]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), sharex=True)
    panels = ((axes[0], "impact", "births_per_household", "Immediate births"),
              (axes[1], "cohort", "completed_fertility", "Completed fertility"))
    for ax, scope, outcome, title in panels:
        for regime, color in zip(REGIMES, (BLUE, ORANGE)):
            y0 = values[(regime, scope, outcome, 1.0)]
            require(y0 > 0, f"Baseline level must be positive for {regime}/{scope}/{outcome}")
            y = [(values[(regime, scope, outcome, factor)] / y0 - 1.0) * 100 for factor in factors]
            ax.plot(x, y, color=color, marker="o", markersize=4, linewidth=1.8, label=REGIMES[regime])
        ax.axhline(0, color="#777777", linewidth=.75)
        ax.axvline(0, color="#777777", linewidth=.75)
        ax.set_title(title)
        ax.set_xlabel("Housing asset price change (%)")
        ax.set_ylabel("Change from own-regime baseline (%)")
        ax.set_xlim(-1, 1)
        ax.set_xticks([-1, 0, 1])
        ax.grid(axis="y", color="#dddddd", linewidth=.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].legend(frameon=False, loc="best")
    title = "Housing costs and fertility" + (" · SYNTHETIC TEST ONLY" if test_only else "")
    fig.suptitle(title, y=.99, fontsize=12)
    footer_1 = "Supplemental figure · House price and mapped rent move together; all preference parameters are fixed."
    footer_2 = ("Impact uses the same inherited households. Cohort distribution has unit household mass; "
                "completed fertility is free to change. Not general equilibrium or a transition.")
    if test_only:
        footer_1 = "SYNTHETIC TEST ONLY · " + footer_1
    fig.text(.5, .042, footer_1, ha="center", fontsize=6.8)
    fig.text(.5, .022, footer_2, ha="center", fontsize=6.8)
    fig.text(.5, .012, LABEL, ha="center", fontsize=6.8)
    fig.tight_layout(rect=(0, .10, 1, .94))
    fig.savefig(out / "price_response.png", dpi=220, bbox_inches="tight")
    fig.savefig(out / "price_response.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--elasticities", type=Path, required=True)
    parser.add_argument("--completed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--test-only", action="store_true", help="Brand output as synthetic validation, never a result")
    args = parser.parse_args()
    args.output = args.output.resolve()
    require(not args.output.exists(), "Refusing to overwrite an existing output directory")
    completed = read_json(args.completed)
    validate_completion(completed, args.comparison, args.elasticities)
    comparison = validate_comparison(read_csv(args.comparison))
    elasticities = validate_elasticities(read_csv(args.elasticities))
    verify_elasticities_from_comparison(elasticities, comparison)
    args.output.mkdir(parents=True, exist_ok=False)
    render_figure(args.output, comparison, args.test_only)
    write_table(args.output, elasticities, args.test_only)
    inputs = [args.comparison, args.elasticities, args.completed, Path(__file__).resolve()]
    manifest = {
        "reference_label": LABEL,
        "status": "synthetic_test_only" if args.test_only else "passed",
        "complete_three_price": True,
        "complete_five_price": False,
        "factors": list(FACTORS),
        "step_sizes": [0.01],
        "figures_are_supplemental": True,
        "synthetic_test_only": args.test_only,
        "outputs": ["price_response.png", "price_response.pdf", "local_elasticities.csv", "local_elasticities.tex"],
        "sources": [{"path": str(path.resolve()), "sha256": sha(path)} for path in inputs],
    }
    with (args.output / "manifest.json").open("x", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, ensure_ascii=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
