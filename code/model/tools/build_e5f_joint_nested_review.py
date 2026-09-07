#!/usr/bin/env python3
"""Build an evidence-driven simultaneous-choice calibration readout.

The selected directory must contain ``summary.json``, ``target_fit_long.csv``,
``parameter_table.csv``, and ``case_receipt.json``.  Policy evidence is optional;
when supplied, every completed case must share the selected summary hash and pass
its recorded date-level checks.  Narrative interpretation is deliberately supplied
by the lead in a small text file or JSON object, rather than hard-wired here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import PageBreak, Paragraph, SimpleDocTemplate, Spacer

import build_e5f_overnight_numerical_report as layout


EXPECTED_POLICIES = {
    "baseline", "supply-plus-20", "dependent-child-ltv95",
    "property-tax-2pct-no-rebate",
}
POLICY_LABELS = {
    "baseline": "Baseline", "supply-plus-20": "Housing supply +20%",
    "dependent-child-ltv95": "Dependent-child LTV 95%",
    "property-tax-2pct-no-rebate": "Property tax doubled",
}
MOMENT_LABELS = {**layout.MOMENT_LABELS,
    "own_family_gap": "Parent/nonparent ownership gap, ages 30-55",
    "own_rate": "Ownership, ages 30-55",
    "old_total_wealth_to_annual_income_p90_p50_7684": "Wealth/income p90/p50, ages 76-84",
}
PARAMETER_LABELS = {
    "kappa_fert": "First-birth taste dispersion", "kappa_fert_continuation": "Later-birth taste dispersion",
    "beta_annual": "Annual discount factor", "tenure_choice_kappa": "Outer taste scale",
    "joint_nest_lambda": "Nest dissimilarity", "chi": "Owner housing-service premium",
    "H0": "Housing-supply normalization", "theta0": "Bequest utility weight",
    "theta1": "Bequest wealth shift", "hbar_child_rooms": "Per-child room floor",
    "first_birth_fixed_cost": "First-birth utility cost",
    "hbar_first_child_jump": "First-child room-floor jump",
    "psi_child_change_2023": "Child-value change, 2007-23",
    "psi_child_2007": "Child value, 2007", "psi_child_2023": "Child value, 2023",
    "housing_supply_elasticity": "Housing-supply elasticity",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"Expected JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def finite(value: Any, label: str) -> float:
    try:
        answer = float(value)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"Non-numeric {label}: {value!r}") from error
    if not math.isfinite(answer):
        raise RuntimeError(f"Non-finite {label}: {value!r}")
    return answer


def fmt(value: Any) -> str:
    return f"{finite(value, 'reported value'):.6g}"


def bool_value(value: Any) -> bool:
    return str(value).strip().lower() == "true"


def load_narrative(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    raw = path.read_text(encoding="utf-8").strip()
    if not raw:
        return {}
    if path.suffix.lower() == ".json":
        value = json.loads(raw)
        if not isinstance(value, dict) or any(not isinstance(v, str) for v in value.values()):
            raise RuntimeError("Narrative JSON must be an object of string fields")
        return value
    return {"interpretation": raw}


def validate_reference(fit_path, parameter_path, selected_fit):
    if fit_path is None and parameter_path is None:
        return None
    if fit_path is None or parameter_path is None:
        raise RuntimeError("A benchmark comparison requires its complete fit and parameter tables")
    fit, parameters = read_csv(fit_path), read_csv(parameter_path)
    targets = {row["moment"]: row for row in selected_fit}
    if len(fit) != 12 or {row["moment"] for row in fit} != set(targets):
        raise RuntimeError("Benchmark moment set differs")
    loss = 0.
    for row in fit:
        for key in ("target", "weight"):
            if finite(row[key], key) != finite(targets[row["moment"]][key], key):
                raise RuntimeError("Benchmark target or weight differs")
        gap = finite(row["model"], "benchmark model") - float(row["target"])
        if not math.isclose(gap, finite(row["gap"], "benchmark gap"), rel_tol=0., abs_tol=1e-10):
            raise RuntimeError("Benchmark gap does not reproduce")
        contribution = float(row["weight"]) * gap**2
        if not math.isclose(contribution, finite(row["loss_contribution"], "benchmark contribution"), rel_tol=0., abs_tol=1e-8):
            raise RuntimeError("Benchmark objective does not reproduce")
        loss += contribution
    if (len(parameters) != 15 or len({r["parameter"] for r in parameters}) != 15
            or sum(bool_value(r["is_free_parameter"]) for r in parameters) != 11):
        raise RuntimeError("Retained benchmark requires all eleven free and four restricted entries")
    for row in parameters:
        value = finite(row["value"], "benchmark parameter")
        if bool_value(row["is_free_parameter"]):
            if not finite(row["lower_bound"], "bound") <= value <= finite(row["upper_bound"], "bound"):
                raise RuntimeError("Benchmark parameter lies outside its reported bounds")
    return dict(fit=fit, parameters=parameters, loss=loss,
                fit_sha256=sha(fit_path), parameters_sha256=sha(parameter_path))


def validate_arrays_and_budget(arrays, budget, label):
    if finite(budget.get("budget_excess_mass"), label + " budget") > 2e-10:
        raise RuntimeError(f"Budget gate failed: {label}")
    if finite(arrays.get("occupied_negative_steps"), label + " value drops") != 0:
        raise RuntimeError(f"Occupied value gate failed: {label}")
    probabilities = arrays.get("probabilities")
    if not isinstance(probabilities, dict) or not probabilities:
        raise RuntimeError(f"Missing probability diagnostics: {label}")
    for name, bounds in probabilities.items():
        low, high = finite(bounds.get("minimum"), name), finite(bounds.get("maximum"), name)
        if finite(bounds.get("nonfinite"), name) != 0 or not 0 <= low <= high <= 1:
            raise RuntimeError(f"Probability gate failed: {label}, {name}")


def validate_selected(selected: Path) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, Any], dict[str, Any], int]:
    for name in ("summary.json", "target_fit_long.csv", "parameter_table.csv", "case_receipt.json"):
        if not (selected / name).is_file():
            raise RuntimeError(f"Selected candidate is missing {name}")
    fit, parameters = read_csv(selected / "target_fit_long.csv"), read_csv(selected / "parameter_table.csv")
    summary, receipt = read_json(selected / "summary.json"), read_json(selected / "case_receipt.json")
    if len(fit) != 12 or {r.get("moment") for r in fit} != set(layout.MOMENT_LABELS):
        raise RuntimeError(f"Selected candidate requires exactly 12 target rows, found {len(fit)}")
    free = [row for row in parameters if bool_value(row.get("is_free_parameter"))]
    if len(free) != 11 or len(parameters) != 14:
        raise RuntimeError("Selected candidate requires 11 free and 3 fixed/derived parameter rows")
    if len({r.get("parameter") for r in parameters}) != 14:
        raise RuntimeError("Duplicate parameter rows")
    domain = {row["name"]: row for row in summary["panel_design"]["domain"]}
    if {row["parameter"] for row in free} != set(domain):
        raise RuntimeError("Free parameter table differs from selected search domain")
    for row in free:
        spec = domain[row["parameter"]]
        value, lower, upper = [finite(row[k], k) for k in ("value", "lower_bound", "upper_bound")]
        if lower != spec["lower"] or upper != spec["upper"] or not lower <= value <= upper:
            raise RuntimeError("Parameter value or bounds differ from selected search domain")
        if row.get("near_bound", "").lower() not in ("true", "false"):
            raise RuntimeError("Missing parameter bound status")
    for row in parameters:
        finite(row.get("value"), row.get("parameter"))
    if receipt.get("status") != "complete" or receipt.get("standard_graph_count") != 17:
        raise RuntimeError("Selected candidate receipt is incomplete")
    loss = finite(summary.get("best_candidate", {}).get("transition_loss", summary.get("best_candidate", {}).get("loss")), "selected loss")
    contribution_sum = 0.0
    for row in fit:
        for column in ("target", "model", "gap", "weight", "loss_contribution"):
            finite(row.get(column), f"fit {row.get('moment')} {column}")
        gap = finite(row["model"], "model") - finite(row["target"], "target")
        if not math.isclose(gap, finite(row["gap"], "gap"), abs_tol=1e-10):
            raise RuntimeError(f"Target gap mismatch: {row.get('moment')}")
        if finite(row["weight"], "weight") <= 0:
            raise RuntimeError("Target weight must remain positive")
        contribution = gap * gap * finite(row["weight"], "weight")
        if not math.isclose(contribution, finite(row["loss_contribution"], "loss contribution"), abs_tol=1e-8):
            raise RuntimeError(f"Target contribution mismatch: {row.get('moment')}")
        contribution_sum += contribution
    if not math.isclose(contribution_sum, loss, abs_tol=1e-7):
        raise RuntimeError("Selected loss does not equal the full target-fit table")
    artifacts = receipt.get("artifact_sha256")
    if not isinstance(artifacts, dict) or not artifacts:
        raise RuntimeError("Selected receipt has no artifact hashes")
    graphs = sorted((selected / "standard_diagnostics").glob("*.png"))
    required = {"summary.json", "target_fit_long.csv", "parameter_table.csv", "dated_state.pkl.gz",
                *[str(g.relative_to(selected)) for g in graphs]}
    if len(graphs) != 17 or not required.issubset(artifacts):
        raise RuntimeError("Selected receipt omits required source or graph artifacts")
    checked = 0
    for relative, expected in artifacts.items():
        artifact = selected / relative
        if not artifact.is_file() or sha(artifact) != expected:
            raise RuntimeError(f"Selected artifact hash mismatch: {artifact}")
        checked += 1
    gates = receipt.get("gates", {})
    limits = dict(market=2e-4, mass=2e-10, population=2e-10,
                  stationary_measurement=2e-8, childless_measurement=2e-10)
    if set(gates) != set(limits):
        raise RuntimeError("Selected receipt omits required numerical gates")
    for name, limit in limits.items():
        gate = gates[name]
        value = finite(gate.get("value"), name)
        if gate.get("passed") is not True or finite(gate.get("limit"), name) != limit or not 0 <= value <= limit:
            raise RuntimeError(f"Selected numerical gate failed: {name}")
    if finite(receipt.get("loss"), "receipt loss") != loss:
        raise RuntimeError("Selected receipt has a different objective")
    validate_arrays_and_budget(receipt["policy_array_diagnostic"], receipt["budget_diagnostic"], "selected")
    return fit, parameters, summary, receipt, checked


def validate_policy(policy_root: Path | None, selected_summary_sha: str, summary: dict) -> tuple[dict[str, Any] | None, list[dict[str, Any]], str]:
    if policy_root is None:
        return None, [], "Not supplied: no policy result is represented as passed."
    receipt_path = policy_root / "equilibrium_receipt.json"
    if not receipt_path.is_file():
        raise RuntimeError("Policy directory is missing equilibrium_receipt.json")
    overall = read_json(receipt_path)
    if overall.get("status") != "complete" or overall.get("failures") != {}:
        raise RuntimeError("Policy experiment is incomplete; report failure separately")
    if (overall.get("scientific_bundle") != summary["code_fingerprints"]["bundle_sha256"]
            or overall.get("target_fingerprint") != summary["target_fingerprint"]):
        raise RuntimeError("Policy source or target fingerprint differs from selected candidate")
    handoff_path = policy_root / "inherited_state_verification.json"
    if not handoff_path.is_file() or sha(handoff_path) != overall.get("inherited_state_verification_sha256"):
        raise RuntimeError("Policy inherited-state verification hash mismatch")
    handoff = read_json(handoff_path)
    if (handoff.get("status") != "exact_feasibility_replay" or
            handoff.get("source_summary_sha256") != selected_summary_sha or handoff.get("birth_queue_source_year") != 2019):
        raise RuntimeError("Policy inherited-state verification differs from selected history")
    if overall.get("smoke") not in (True, False):
        raise RuntimeError("Missing policy horizon classification")
    expected_years = [2023, 2027] if overall["smoke"] else list(range(2023, 2064, 4))
    cases = overall.get("cases")
    if not isinstance(cases, dict) or set(cases) != EXPECTED_POLICIES:
        raise RuntimeError("Policy receipt must contain the four standard policy cases")
    if overall.get("selected_summary_sha256") != selected_summary_sha:
        raise RuntimeError("Policy evidence belongs to a different selected candidate")
    checks: list[dict[str, Any]] = []
    all_years: set[int] = set()
    expected_dates: int | None = None
    for name in sorted(EXPECTED_POLICIES):
        case_dir, case = policy_root / name, cases[name]
        local = read_json(case_dir / "receipt.json")
        if local != case or case.get("status") != "complete":
            raise RuntimeError(f"Policy case receipt mismatch or incomplete: {name}")
        if case.get("source_summary_sha256") != selected_summary_sha:
            raise RuntimeError(f"Policy case has mismatched selected source: {name}")
        dates = int(case.get("dates", 0))
        if expected_dates is None:
            expected_dates = dates
        if dates not in (2, 11) or dates != expected_dates:
            raise RuntimeError("Policy cases must consistently contain either 2 or 11 dates")
        gates = case.get("gates", {})
        if finite(gates.get("maximum_market_residual"), "policy market residual") > 2e-4 or finite(gates.get("maximum_mass_residual"), "policy mass residual") > 2e-10 or int(gates.get("maximum_nonfinite_count", 0)) != 0:
            raise RuntimeError(f"Policy gates failed: {name}")
        rows = read_csv(case_dir / "policy_path.csv")
        years = [int(finite(r.get("calendar_year"), "policy calendar year")) for r in rows]
        if len(rows) != dates or years != expected_years:
            raise RuntimeError(f"Policy path has inconsistent dates: {name}")
        all_years.update(years)
        for row in rows:
            for field in ("birth_children_topcode_adjusted", "owner_rate", "housing_demand_per_adult",
                          "asset_price", "topcode_adjusted_births_per_adult", "population_index_2023"):
                finite(row.get(field), f"{name} dated {field}")
        for year in years:
            date_dir = case_dir / f"date_{year}"
            budget, arrays = read_json(date_dir / "budget_summary.json"), read_json(date_dir / "policy_array_summary.json")
            validate_arrays_and_budget(arrays, budget, f"{name}, {year}")
            graphs = list((date_dir / "standard_diagnostics").glob("*.png"))
            if len(graphs) != 17:
                raise RuntimeError(f"Policy date lacks the 17 standard graphs: {name}, {year}")
            checks.append({"policy": name, "year": year, "budget_excess_mass": budget["budget_excess_mass"], "standard_graphs": 17})
    effects_path = policy_root / "policy_effects.csv"
    if not effects_path.is_file():
        raise RuntimeError("Policy effects table is missing")
    if effects_path.is_file():
        effects = read_csv(effects_path)
        paths = {name: {int(finite(row["calendar_year"], "year")): row for row in read_csv(policy_root / name / "policy_path.csv")} for name in EXPECTED_POLICIES}
        expected = {(name, year) for name in EXPECTED_POLICIES - {"baseline"} for year in (min(all_years), max(all_years))}
        if len(effects) != len(expected) or {(r.get("policy"), int(finite(r.get("year"), "effect year"))) for r in effects} != expected:
            raise RuntimeError("Policy effect rows are incomplete or duplicated")
        for row in effects:
            base, policy = paths["baseline"][int(row["year"])], paths[row["policy"]][int(row["year"])]
            recomputed = {
                "births_percent": 100 * (finite(policy["birth_children_topcode_adjusted"], "births") / finite(base["birth_children_topcode_adjusted"], "baseline births") - 1),
                "ownership_pp": 100 * (finite(policy["owner_rate"], "ownership") - finite(base["owner_rate"], "baseline ownership")),
                "rooms_percent": 100 * (finite(policy["housing_demand_per_adult"], "rooms") / finite(base["housing_demand_per_adult"], "baseline rooms") - 1),
            }
            if any(not math.isclose(finite(row[key], key), value, abs_tol=1e-10) for key, value in recomputed.items()):
                raise RuntimeError("Recorded policy effect does not reproduce paths")
    status = "Complete 44-date policy paths verified." if expected_dates == 11 else "Complete two-date smoke verified; this is not a 44-date full path."
    return overall, checks, status


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--selected-dir", type=Path, required=True, help="Complete selected candidate directory.")
    parser.add_argument("--reference-fit", type=Path, help="Optional complete retained-benchmark target table; requires its parameters.")
    parser.add_argument("--reference-parameters", type=Path, help="Complete retained-benchmark free and restricted parameters.")
    parser.add_argument("--output", type=Path, required=True, help="New PDF path; refuses to overwrite.")
    parser.add_argument("--policy-results", type=Path, help="Optional root containing four policy cases and receipt.")
    parser.add_argument("--search-verification", type=Path, help="Optional completed search/final verification JSON; hashed in receipt.")
    parser.add_argument("--narrative", type=Path, help="Optional lead-supplied .txt or JSON strings: summary, interpretation, outstanding.")
    parser.add_argument("--fixture-label", default="", help="Required label for smoke/preview use, e.g. 'SMOKE PREVIEW - NOT A SEARCH'.")
    args = parser.parse_args()
    selected, output = args.selected_dir.resolve(), args.output.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    fit, parameters, summary, _, selected_artifacts = validate_selected(selected)
    reference = validate_reference(args.reference_fit, args.reference_parameters, fit)
    selected_summary_sha = sha(selected / "summary.json")
    search = read_json(args.search_verification.resolve()) if args.search_verification else None
    if search:
        # This report summarizes the independently checked original controller
        # receipt; a generic status JSON is not final calibration verification.
        if (search.get("status") != "pass" or search.get("exact_repeats") != 2
                or search.get("exact_standard_pngs") != 17):
            raise RuntimeError("Two exact final repetitions are required")
        origin = search.get("selected", {})
        if finite(origin.get("loss"), "repeated selected loss") != float(summary["best_candidate"]["transition_loss"]):
            raise RuntimeError("Final repetitions concern another selected loss")
        if origin.get("unit_vector") != summary["panel_design"]["unit_vector"]:
            raise RuntimeError("Final repetitions concern another parameter vector")
        if origin.get("plan_sha256") != read_json(selected / "case_receipt.json")["plan_sha256"]:
            raise RuntimeError("Final repetitions concern another source plan")
    narrative = load_narrative(args.narrative.resolve() if args.narrative else None)
    policy, policy_checks, policy_status = validate_policy(args.policy_results.resolve() if args.policy_results else None, selected_summary_sha, summary)
    loss = finite(summary["best_candidate"].get("transition_loss", summary["best_candidate"].get("loss")), "loss")
    fixture = args.fixture_label.strip()
    if ("smoke" in {part.lower() for part in selected.parts} or (policy and policy["smoke"])) and not fixture:
        raise RuntimeError("Smoke inputs require --fixture-label; they cannot be presented as an unlabeled calibration readout")

    fonts = Path("/System/Library/Fonts/Supplemental")
    pdfmetrics.registerFont(TTFont("Review", str(fonts / "Arial.ttf")))
    pdfmetrics.registerFont(TTFont("ReviewBold", str(fonts / "Arial Bold.ttf")))
    pdfmetrics.registerFontFamily("Review", normal="Review", bold="ReviewBold")
    styles = getSampleStyleSheet()
    for name, size, leading, face, color in (("RTitle", 22, 26, "ReviewBold", layout.BLUE), ("RHead", 15, 18, "ReviewBold", layout.BLUE), ("RBody", 9.5, 13, "Review", colors.black), ("RSmall", 7.5, 9.5, "Review", colors.black), ("RTable", 8.5, 11, "Review", colors.black), ("RTableHead", 8.5, 11, "ReviewBold", colors.white)):
        styles.add(ParagraphStyle(name, fontName=face, fontSize=size, leading=leading, textColor=color, spaceAfter=7))
    story: list[Any] = []
    def add(text: str, style: str = "RBody") -> None:
        story.append(Paragraph(text, styles[style]))
    def heading(text: str) -> None:
        story.append(PageBreak()); add(text, "RHead")
    def table(rows: list[list[Any]], widths: list[float]) -> None:
        body = [[Paragraph(escape(str(cell)).replace("\n", "<br/>"), styles["RTableHead" if i == 0 else "RTable"]) for cell in row] for i, row in enumerate(rows)]
        story.append(layout.table(body, widths, font_size=8.5))

    add("Simultaneous-choice calibration readout", "RTitle")
    add((fixture + " | " if fixture else "") + "Experimental simultaneous tenure and birth-attempt choice", "RSmall")
    if narrative.get("summary"): add(escape(narrative["summary"]))
    add(f"<b>Selected objective: {loss:.3f}.</b> All twelve target rows and eleven free parameters appear in the following tables.")
    if reference:
        add(f"<b>Retained benchmark objective: {reference['loss']:.3f}.</b> The target and weight rows match exactly. Its full fit and parameter tables appear in the appendix; the shock specification differs.")
    add("<b>Reproduction:</b> " + ("Two final histories reproduce the selected fit and all seventeen standard graphs exactly." if search else "Final repetitions of a searched candidate are not yet certified."))
    add(f"<b>Policy status:</b> {escape(policy_status)}")
    if narrative.get("interpretation"): add("<b>Interpretation:</b> " + escape(narrative["interpretation"]))
    if narrative.get("outstanding"): add("<b>Outstanding issues:</b> " + escape(narrative["outstanding"]))
    add("<b>Economic specification.</b> Tenure and birth attempts are chosen jointly after the taste shocks are observed. Tenure is committed before conception success; housing size, consumption and saving can adapt within tenure afterward. The inner taste scale is λκ, with 0 &lt; λ ≤ 1. This restriction is experimental.", "RSmall")
    add("<b>Policy closure.</b> Each date clears housing markets with current prices treated as permanent. After 2023, outside entry is zero, retention is one, and the inherited four-slot birth queue converts births to new households at 1/2.1. Supply elasticity is fixed at 0.63. Tax revenue is discarded. These are finite paths under temporary equilibrium.", "RSmall")

    heading("Complete target fit")
    add("Shares and share gaps are fractions: 0.01 equals one percentage point. The objective is the sum of weight × gap²; weights need not equal empirical precision.", "RSmall")
    rows = [["Moment", "Target", "Model", "Gap", "Weight", "Loss contribution"]]
    rows += [[MOMENT_LABELS.get(row["moment"], row["moment"]), *[fmt(row[x]) for x in ("target", "model", "gap", "weight", "loss_contribution")]] for row in fit]
    table(rows, [180, 60, 60, 62, 70, 76])
    add(f"Full-table objective check: sum of contributions = {loss:.9f}.", "RSmall")
    add("For the ownership gap, the active model parent group has at least one child at home; the nonparent group has no previous birth. The code retains a legacy new-parent variable name.", "RSmall")

    heading("Free parameters, bounds, and fixed/derived entries")
    free = [row for row in parameters if bool_value(row.get("is_free_parameter"))]
    fixed = [row for row in parameters if not bool_value(row.get("is_free_parameter"))]
    table([["Free parameter", "Value", "Lower", "Upper", "Near bound?"]] + [[PARAMETER_LABELS.get(row["parameter"], row["parameter"]), fmt(row["value"]), fmt(row["lower_bound"]), fmt(row["upper_bound"]), row.get("near_bound", "")] for row in free], [220, 80, 70, 70, 82])
    add("Near-bound flags use the existing convention: within 2% of the full physical parameter range, rather than distance in the transformed search coordinate.", "RSmall")
    story.append(Spacer(1, 10))
    table([["Fixed/derived entry", "Value", "Recorded status"]] + [[PARAMETER_LABELS.get(row["parameter"], row["parameter"]), fmt(row["value"]), {"psi_child_2007": "Normalized to old completed fertility 2.1", "psi_child_2023": "Old level plus fitted 2007–23 change", "housing_supply_elasticity": "Externally fixed"}[row["parameter"]]] for row in fixed], [180, 90, 252])

    heading("Policy outcomes and verification coverage")
    if policy is None:
        add("No policy directory was supplied. Policy outcomes are unavailable and are not represented as passed.")
    else:
        effects = read_csv(args.policy_results.resolve() / "policy_effects.csv") if (args.policy_results.resolve() / "policy_effects.csv").is_file() else []
        shown_years = {min(int(row["year"]) for row in effects), max(int(row["year"]) for row in effects)}
        table([["Policy", "Year", "Births (%)", "Ownership (pp)", "Rooms (%)"]] + [[POLICY_LABELS.get(row["policy"], row["policy"]), row["year"], fmt(row["births_percent"]), fmt(row["ownership_pp"]), fmt(row["rooms_percent"])] for row in effects if int(row["year"]) in shown_years], [190, 55, 85, 100, 85])
        add("Each effect compares the policy with the baseline at the same date and inherited starting population. The table shows the first and last computed dates; complete dated results are retained in each policy's policy_path.csv.", "RSmall")
        add(policy_status, "RSmall")
    if policy_checks:
        add(f"All {len(policy_checks)} dated packets contain the unchanged seventeen graphs and pass the occupied-value, budget and probability checks. Maximum budget-violating mass: {max(x['budget_excess_mass'] for x in policy_checks):.3g}, below 2e-10.", "RSmall")

    if reference:
        heading("Retained benchmark: complete target fit")
        table([["Moment", "Target", "Model", "Gap", "Weight", "Loss contribution"]] + [
            [MOMENT_LABELS[row["moment"]], *[fmt(row[k]) for k in ("target", "model", "gap", "weight", "loss_contribution")]]
            for row in reference["fit"]], [180, 60, 60, 62, 70, 76])
        add(f"Complete-table objective: {reference['loss']:.9f}. This is the retained production comparison; experimental search does not promote a replacement.", "RSmall")
        heading("Retained benchmark: all parameters and restrictions")
        table([["Parameter", "Value", "Lower", "Upper", "Free?", "Near bound?"]] + [
            [{"tenure_choice_kappa": "Tenure taste dispersion"}.get(row["parameter"], PARAMETER_LABELS.get(row["parameter"], row["parameter"])), fmt(row["value"]),
             fmt(row["lower_bound"]) if bool_value(row["is_free_parameter"]) else "-",
             fmt(row["upper_bound"]) if bool_value(row["is_free_parameter"]) else "-",
             row["is_free_parameter"], row["near_bound"]] for row in reference["parameters"]], [190, 70, 65, 65, 57, 75])
        add("The old child-value intercept is normalized to completed fertility 2.1; its 2023 level follows from the estimated change. Housing-supply elasticity and the tenure taste scale are externally fixed. The original model estimates separate first- and later-birth taste scales.", "RSmall")

    graph_dir = selected / "standard_diagnostics"
    graphs = sorted(graph_dir.glob("*.png"))
    if graphs:
        heading("Standard diagnostic graphs")
        add("The established standard graph set is reproduced without redesign. These are selected-candidate diagnostics, not additional policy conclusions.", "RSmall")
        for index, graph in enumerate(graphs):
            add(escape(graph.stem.replace("_", " ")), "RSmall")
            story.append(layout.image_fit(graph, 515, 300))
            if index % 2 == 1 and index != len(graphs) - 1: story.append(PageBreak())

    def footer(canvas: Any, doc: Any) -> None:
        canvas.saveState(); canvas.setFont("Review", 7); canvas.setFillColor(layout.MID_GREY)
        canvas.drawString(36, 20, "Simultaneous-choice experiment | Discussion copy | Production unchanged")
        canvas.drawRightString(A4[0] - 36, 20, f"Page {doc.page}"); canvas.restoreState()
    output.parent.mkdir(parents=True, exist_ok=True)
    SimpleDocTemplate(str(output), pagesize=A4, rightMargin=36, leftMargin=36, topMargin=34, bottomMargin=34, title="Simultaneous-choice calibration readout", author="Research discussion draft").build(story, onFirstPage=footer, onLaterPages=footer)
    verification = {"status": "numerical_source_checks_passed_visual_review_pending", "pdf": str(output), "pdf_sha256": sha(output), "selected_dir": str(selected), "selected_summary_sha256": selected_summary_sha, "selected_artifacts_checked": selected_artifacts, "reference": reference, "fit_rows": len(fit), "fit_cells": fit, "parameter_cells": parameters, "selected_receipt_sha256": sha(selected / "case_receipt.json"), "policy_receipt_sha256": sha(args.policy_results.resolve() / "equilibrium_receipt.json") if policy else None, "free_parameters": len(free), "fixed_or_derived_parameters": len(fixed), "loss": loss, "fixture_label": fixture or None, "policy_status": policy_status, "policy_date_checks": policy_checks, "search_verification_sha256": sha(args.search_verification.resolve()) if args.search_verification else None, "narrative_sha256": sha(args.narrative.resolve()) if args.narrative else None, "builder_sha256": sha(Path(__file__).resolve()), "production_promoted": False}
    verification_path = output.parent / f"{output.stem}_verification.json"
    verification_path.write_text(json.dumps(verification, indent=2) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
