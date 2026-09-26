#!/usr/bin/env python3
"""Torch-only saved-packet comparison and export for the four utility arms.

No model solve, production launch, retry, or source mutation is performed here.
The scientific signature and exact comparison descend from the reviewed frozen
nightpair_20260925_v1/render_pair.py; imports of that file are deliberately
avoided because they also import its old runner and PDF dependencies.  Each
checkpoint hash is authenticated separately: gzip bytes need not be identical.
Native classes/observers/plots are loaded only by the pinned runner's setup.

Run the PDF skill's artifact-operation marker before the first actual export.
Rendered PNGs are retained for lead visual review; automatic rendering is not
represented as a completed visual inspection.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import pickle
import shutil
import tempfile
from collections.abc import Mapping
from xml.sax.saxutils import escape


STANDARD_NAMES = (
    "fertility_by_age", "ownership_by_age", "housing_market", "housing_prices",
    "market_clearing_by_market", "market_clearing_residuals", "owner_rungs", "tenure_services",
    "fertility_policy_by_age_income_state", "housing_by_age_income_state",
    "ownership_by_age_income_state", "liquid_wealth_by_age_income_state", "income_state_outcomes",
    "policy_childless_renter_age30", "wealth_dist_childless_renter_age30",
    "policy_childless_renter_age42", "wealth_dist_childless_renter_age42",
)
FIXED_NAMES = (
    "theta1", "psi_child", "payroll_tax", "pension_period", "housing_supply_elasticity",
    "tenure_choice_kappa", "alpha_cons", "sigma", "selling_cost", "financed_share",
    "annual_depreciation", "period_depreciation", "annual_property_tax", "period_property_tax",
    "income_process", "entrant_conversion_factor", "adult_entry_birth_to_household_conversion",
    "child_benefit_exponent", "utility_reference_rent", "pension_to_gross_worker_earnings",
)
TARGET_FIELDS = ("moment", "target", "model", "gap", "weight", "loss_contribution")
PARAMETER_FIELDS = ("parameter", "estimate", "lower", "upper", "near_bound", "status")
ENV_PIN = "EXPECTED_UTILITY_COMPARISON_CONTRACT_SHA256"
TIMING_FIELDS = {
    "chosen_solve_seconds", "objective_stationary_solve_seconds", "complete_objective_wall_seconds",
}


def read_json(path):
    return json.loads(Path(path).read_text())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_text(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"


def write_new_json(path, value):
    """Publish an entire file atomically; never replace an existing output."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(_json_text(value))
        stream.flush()
        os.fsync(stream.fileno())
    try:
        os.link(temporary, path)
    finally:
        temporary.unlink()


def case_directory(path):
    path = Path(path).resolve()
    return path if (path / "receipt.json").is_file() else path / "case"


def _number(value, label):
    if isinstance(value, bool):
        raise RuntimeError(label + " must be numeric")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(label + " must be numeric") from exc
    if not math.isfinite(result):
        raise RuntimeError(label + " must be finite")
    return result


def same_numerics(a, b, tolerance=0.):
    """Reviewed nightpair semantics: rtol=0, explicit atol, equal_nan=True.

    The four-arm caller only accepts zero tolerance. Missing dictionary keys or
    changed array shapes fail. NaN equality is retained for native value arrays;
    the target/parameter/receipt finite-value checks are separate.
    """
    import numpy as np
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(
            same_numerics(a[k], b[k], tolerance) for k in a)
    if isinstance(a, (list, tuple, np.ndarray)):
        try:
            left, right = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
            return left.shape == right.shape and bool(np.allclose(
                left, right, rtol=0, atol=tolerance, equal_nan=True))
        except (TypeError, ValueError):
            return list(a) == list(b)
    if isinstance(a, (int, float, np.number)) and isinstance(b, (int, float, np.number)):
        return bool(np.isclose(a, b, rtol=0, atol=tolerance, equal_nan=True))
    return a == b


def read_table(path, fields, key):
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or set(reader.fieldnames) != set(fields):
            raise RuntimeError(f"{path}: incomplete or unexpected table columns")
        rows = list(reader)
    if any(set(row) != set(fields) or any(v is None for v in row.values()) for row in rows):
        raise RuntimeError(f"{path}: malformed table row")
    names = [row[key] for row in rows]
    if any(not name for name in names) or len(set(names)) != len(names):
        raise RuntimeError(f"{path}: missing or duplicate {key}")
    return rows


def validate_tables(fits, parameters, objective, receipt, free_count):
    """Validate complete economic tables against the pinned target/restrictions."""
    targets = {row["restriction_id"]: row for row in objective["target_rows"]}
    if len(targets) != 13 or len(fits) != 13 or {r["moment"] for r in fits} != set(targets):
        raise RuntimeError("require all 13 unique pinned target rows")
    weighted = 0
    for row in fits:
        source = targets[row["moment"]]
        values = {k: _number(row[k], k) for k in ("target", "model", "gap")}
        if values["target"] != float(source["target"]):
            raise RuntimeError("target value differs from pinned objective")
        if values["gap"] != values["model"] - values["target"]:
            raise RuntimeError("target gap differs from model minus target")
        if source["actual_weight"] is None:
            if row["moment"] != "initial_normalization" or row["weight"] or row["loss_contribution"]:
                raise RuntimeError("only completed-fertility normalization may be unweighted")
        else:
            weight = _number(row["weight"], "weight")
            contribution = _number(row["loss_contribution"], "loss contribution")
            if weight <= 0 or contribution < 0 or weight != float(source["actual_weight"]):
                raise RuntimeError("target weight/contribution differs or is invalid")
            # The frozen objective computes squared gaps; allow no changed score
            # in repeated tables. This arithmetic check only admits roundoff.
            expected = weight * values["gap"] ** 2
            if not math.isclose(contribution, expected, rel_tol=2e-14, abs_tol=1e-14):
                raise RuntimeError("loss contribution differs from weighted squared gap")
            weighted += 1
    if weighted != 12:
        raise RuntimeError("require exactly 12 positively weighted target rows")
    if sum(float(r["loss_contribution"]) for r in fits if r["loss_contribution"]) != receipt["loss"]:
        raise RuntimeError("receipt loss differs from complete target table")
    restrictions = {row["parameter"]: row for row in objective["parameter_restrictions"]}
    names = [row["parameter"] for row in parameters]
    if (len(restrictions) != free_count or len(names) != free_count + len(FIXED_NAMES)
            or set(names) != set(restrictions) | set(FIXED_NAMES) or len(set(names)) != len(names)):
        raise RuntimeError("require all free and fixed parameter/restriction rows")
    if set(receipt["point"]) != set(restrictions):
        raise RuntimeError("point differs from the complete free-coordinate set")
    for row in parameters:
        name, value = row["parameter"], _number(row["estimate"], "parameter estimate")
        if not row["status"]:
            raise RuntimeError("parameter restriction status is absent")
        if name in restrictions:
            bound = restrictions[name]
            low, high = float(bound["lower"]), float(bound["upper"])
            if (float(row["lower"]) != low or float(row["upper"]) != high
                    or not low <= value <= high or value != receipt["point"][name]):
                raise RuntimeError("free parameter value or restriction differs")
            near = min(value - low, high - value) <= .01 * (high - low)
            if row["near_bound"] != str(near):
                raise RuntimeError("near-bound flag differs from frozen 1%-of-range screen")
        elif row["lower"] or row["upper"] or row["near_bound"]:
            raise RuntimeError("fixed/derived row unexpectedly carries search bounds")


def _load_contract(contract):
    import run_e5f_utility_comparison as runner
    if isinstance(contract, Mapping):
        result = dict(contract)
        # Preparation places budget.json beside contract.json. Recheck the
        # original file so a mutated in-memory mapping cannot bypass its pin.
        path = Path(result["files"]["budget"]["path"]).parent / "contract.json"
        if runner.verified_contract(path) != result:
            raise RuntimeError("in-memory contract differs from its pinned original")
    else:
        result = runner.verified_contract(Path(contract))
    expected = os.environ.get(ENV_PIN, "")
    if len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected):
        raise RuntimeError("reviewed comparison-contract SHA256 environment required")
    own = result["files"].get(Path(__file__).name)
    if not own or Path(own["path"]).resolve() != Path(__file__).resolve() or sha256(__file__) != own["sha256"]:
        raise RuntimeError("collector source is not pinned to the executing file")
    required = result["required_readout"]
    if (required["numerical_tolerance"] != 0 or required["repetitions"] != 2
            or required["target_rows"] != 13 or required["weighted_rows"] != 12
            or required["standard_graph_count"] != 17 or not required["compare_original_to_each_repeat"]):
        raise RuntimeError("comparison/readout contract differs from the reviewed requirements")
    return result, runner


def _receipt_science(receipt):
    # Runtime durations and per-file gzip hashes are provenance, not numerical
    # equality objects. Every other scientific receipt field is compared.
    value = {k: v for k, v in receipt.items() if k not in TIMING_FIELDS | {"case_checkpoint_sha256"}}
    value["normalization"] = {k: v for k, v in receipt["normalization"].items()
                              if k != "stationary_solve_seconds"}
    return value


def scientific_checkpoint(case, contract, arm, runtime, tax):
    """Authenticate before unpickling, then extract the frozen scientific signature."""
    import numpy as np
    case = case_directory(case)
    receipt = read_json(case / "receipt.json")
    expected = {
        "status": "verified_experimental_point", "utility_comparison_arm": arm,
        "comparison_contract_sha256": os.environ[ENV_PIN],
        "target_system_sha256": contract["arms"][arm]["objective"]["sha256"],
        "source_manifest_sha256": contract["parent_source_inventory"]["sha256"],
        "selected_checkpoint_sha256": tax.CHECKPOINT_SHA,
        "free_count": contract["arms"][arm]["free_count"], "weighted_count": 12, "display_count": 13,
        "benefit_exponent": contract["arms"][arm]["benefit_exponent"],
        "utility_reference_rent": contract["reference_rent"],
    }
    for key, value in expected.items():
        if receipt.get(key) != value:
            raise RuntimeError(f"{case}: receipt {key} differs from pinned contract")
    checkpoint = case / "initial_state.pkl.gz"
    if receipt["case_checkpoint_sha256"] != sha256(checkpoint):
        raise RuntimeError(f"{case}: native checkpoint hash differs")
    objective = read_json(contract["arms"][arm]["objective"]["path"])
    fits = read_table(case / "target_fit.csv", TARGET_FIELDS, "moment")
    params = read_table(case / "parameters.csv", PARAMETER_FIELDS, "parameter")
    validate_tables(fits, params, objective, receipt, expected["free_count"])
    with gzip.open(checkpoint, "rb") as stream:
        packet = pickle.load(stream)
    P, solution, evaluation = packet["parameters"], packet["solution"], packet["evaluation"]
    if (P.adult_entry_clock != "split_birth_vintage" or P.utility_comparison_arm != arm
            or P.utility_child_benefit_exponent != expected["benefit_exponent"]
            or P.utility_reference_rent != expected["utility_reference_rent"]):
        raise RuntimeError("saved parameter object differs from entry/utility contract")
    actual = tax.actual_parameters(P)
    actual.update(theta1=P.theta1, psi_child=P.psi_child, payroll_tax=P.tau_pay,
        pension_period=P.pension, housing_supply_elasticity=P.xi_supply[0],
        tenure_choice_kappa=P.tenure_choice_kappa, alpha_cons=P.alpha_cons, sigma=P.sigma,
        selling_cost=P.psi, financed_share=P.phi[0], annual_depreciation=receipt["annual_depreciation"],
        period_depreciation=P.delta, annual_property_tax=receipt["annual_property_tax"],
        period_property_tax=P.tau_H, income_process=15, entrant_conversion_factor=P.entrant_conversion_factor,
        adult_entry_birth_to_household_conversion=1/2.1,
        child_benefit_exponent=P.utility_child_benefit_exponent,
        utility_reference_rent=P.utility_reference_rent,
        pension_to_gross_worker_earnings=contract["pension_ratio"])
    if any(float(row["estimate"]) != float(actual[row["parameter"]]) for row in params):
        raise RuntimeError("parameter table differs from saved native parameter object")
    if float(P.psi_child) != receipt["normalization"]["psi_child"]:
        raise RuntimeError("saved psi differs from normalization receipt")
    signature = dict(price=receipt["price"], native_price=np.asarray(solution.p_eq),
        evaluated_price=np.asarray(evaluation.policy.price),
        native_V=np.asarray(solution.V), native_g=np.asarray(solution.g),
        V=np.asarray(evaluation.policy.V), g=np.asarray(evaluation.g_current),
        stationary_g_pre=np.asarray(packet["stationary_g_pre"]),
        moments=runtime["chain"].extract_moments(solution, P), psi=float(P.psi_child),
        normalization={k: v for k, v in receipt["normalization"].items() if k != "stationary_solve_seconds"},
        loss=receipt["loss"], receipt=_receipt_science(receipt),
        target_table=fits, parameter_table=params)
    pins = {name: sha256(case / name) for name in
            ("receipt.json", "initial_state.pkl.gz", "target_fit.csv", "parameters.csv")}
    return dict(signature=signature, packet=packet, receipt=receipt, fits=fits, parameters=params,
                case=case, pins=pins)


def _verify_loaded(contract, arm, original_case, other_cases, required_count, runtime, tax):
    if isinstance(required_count, bool) or required_count not in (1, 2):
        raise RuntimeError("required_count must be 1 for smoke or 2 for selected repetitions")
    if len(other_cases) != required_count:
        raise RuntimeError(f"require exactly {required_count} comparison cases, received {len(other_cases)}")
    paths = [case_directory(original_case), *[case_directory(path) for path in other_cases]]
    if len(set(paths)) != len(paths):
        raise RuntimeError("original and repetition paths must be distinct")
    selected = scientific_checkpoint(paths[0], contract, arm, runtime, tax)
    repeats = []
    for index, path in enumerate(paths[1:]):
        other = scientific_checkpoint(path, contract, arm, runtime, tax)
        # Each comparison is against the ORIGINAL, never merely repeat 1 vs 2.
        for field, value in selected["signature"].items():
            if not same_numerics(value, other["signature"][field], 0.):
                raise RuntimeError(f"original vs repeat {index}: scientific field {field} differs at zero tolerance")
        repeats.append(dict(case=str(path), pins=other["pins"], exact_match_to_original=True))
    result = dict(schema="utility_comparison_repetition_v1", status="exact_comparison_passed",
        arm=arm, comparison_contract_sha256=os.environ[ENV_PIN], original_case=str(paths[0]),
        original_pins=selected["pins"], compared_cases=repeats, required_count=required_count,
        repetitions_verified=len(repeats), numerical_absolute_tolerance=0., numerical_relative_tolerance=0.,
        exact_repeat_claim=required_count == 2, compared_fields=list(selected["signature"]),
        native_solve_count=0)
    return result, selected


def verify_cases(contract, arm, original_case, other_cases, required_count):
    """Return an exact receipt or raise; one smoke counterpart or BOTH repeats.

    A mapping must come from runner.verified_contract; a path is verified here.
    Setup only imports the frozen native code and observers, without a GE solve.
    """
    contract, runner = _load_contract(contract)
    with tempfile.TemporaryDirectory(prefix="utility_compare_") as temporary:
        _, _, _, tax, _, _, runtime, _ = runner.setup(contract, arm, Path(temporary) / "runtime")
        result, _ = _verify_loaded(contract, arm, original_case, other_cases, required_count, runtime, tax)
    return result


def verify_repetitions(contract, arm, selected_case_path, repeat_case_paths):
    return verify_cases(contract, arm, selected_case_path, repeat_case_paths, required_count=2)


def parameter_pages(rows, page_size=10):
    if page_size <= 0:
        raise ValueError("page_size must be positive")
    return [rows[start:start + page_size] for start in range(0, len(rows), page_size)]


def fmt(value):
    if value is None or value == "":
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number and abs(number) < .0005:
        return f"{number:.3e}"
    return f"{number:,.3f}"


def _graphs(case, packet, runtime, output):
    source = case / "standard_diagnostics"
    expected = {name + ".png" for name in STANDARD_NAMES}
    existing = {p.name for p in source.glob("*.png")}
    if existing and existing != expected:
        raise RuntimeError("saved standard diagnostic gallery is partial or changed")
    if existing:
        destination = output / "standard_diagnostics"
        destination.mkdir()
        for name in sorted(expected):
            shutil.copyfile(source / name, destination / name)
    else:
        runtime["audit"].standard_diagnostics(packet, output, validate_production_young=False)
    if {p.name for p in (output / "standard_diagnostics").glob("*.png")} != expected:
        raise RuntimeError("native exporter did not produce the unchanged 17 standard diagnostic names")
    return {name: sha256(output / "standard_diagnostics" / name) for name in sorted(expected)}


def render_report(contract, arm, selected, comparison, output, *, layout_fixture=False):
    """Report-only PDF and every-page PNGs. Requires later visual inspection."""
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle, Image
    import pymupdf as fitz
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="SmallCell", fontName="Helvetica", fontSize=8, leading=10))
    styles["BodyText"].fontSize, styles["BodyText"].leading = 10, 14
    cell = lambda value: Paragraph(escape(str(value)), styles["SmallCell"])
    title = ("LAYOUT FIXTURE - NO NEW CALIBRATION" if layout_fixture else
             "Overnight utility comparison: " + arm.replace("_", " "))
    story = [Paragraph(title, styles["Title"])]
    count = contract["arms"][arm]["free_count"]
    lines = [
        f"Selected weighted loss: {fmt(selected['receipt']['loss'])}. The complete fit has 13 rows, "
        f"12 positively weighted moments and {count} free structural coordinates. Psi is normalized to completed fertility 2.1.",
        "Both independently executed repetitions match the original selected saved solution at zero tolerance. "
        "This verifies repetition; it does not establish identification, optimizer convergence or grid robustness.",
        "The comparison uses a common experimental material-utility normalization at the fixed reference rent. "
        "The mild concavity exponent 0.86 is a fixed sensitivity, not an externally estimated fact.",
        "The first-birth room response uses an unmatched PSID proxy. The model sums positive estates from all households; "
        "the SCF bequest target is child-directed. These measurement limitations remain in all arms.",
        "The adopted pension-to-gross-working-earnings ratio sets baseline payroll tax from baseline demographics. "
        "Equal retiree pensions balance PAYGO; transitions are outside this comparison.",
        "Entry wealth, estate recipients and mortality mapping remain inherited. Adult entry splits birth vintages "
        "equally at 16 and 20 years, with the division by 2.1 applied once.",
        "Exact source values, targets, weights, restrictions, near-bound flags and file fingerprints remain in the accompanying CSV and JSON files.",
    ]
    if layout_fixture:
        lines = [
            "LAYOUT FIXTURE - NO NEW CALIBRATION. This report tests pagination and image placement only. "
            "It is not a result from any new utility arm and makes no claim about scientific repetition.",
            "The 13 target rows, previously reported parameter rows and all 17 standard figures are copied from "
            "the supplied saved historical case. The old housing-floor row is omitted from this layout example; "
            "five visibly marked, unestimated placeholder rows exercise the 29-row share-arm layout.",
            "Empty estimates in placeholder rows are intentional. Bounds shown there test formatting; they "
            "do not report estimates or authorize a run. No native model is loaded or solved by this fixture.",
            "The future comparison must disclose that the reference-rent normalization is experimental, that "
            "0.86 is a fixed curvature sensitivity, that the PSID observer is an unmatched proxy, and that "
            "positive all-estates differs from the child-directed SCF target.",
            "All pages are rendered to PNG for visual inspection. Successful file creation or text checks "
            "alone do not certify that visual review is complete.",
        ]
    for line in lines:
        story += [Spacer(1, 10), Paragraph(escape(line), styles["BodyText"])]
    def add_table(title, headers, rows, widths):
        story.extend([PageBreak(), Paragraph(title, styles["Heading1"]), Spacer(1, 10)])
        table = Table([[cell(x) for x in headers]] + [[cell(x) for x in row] for row in rows],
                      colWidths=widths, repeatRows=1, hAlign="LEFT")
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dfebf4")),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f3f5f7")]),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("TOPPADDING", (0, 0), (-1, -1), 6), ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ]))
        story.append(table)
    add_table("Complete target fit", ["Moment", "Target", "Model", "Gap", "Weight", "Loss"],
        [[r["moment"]] + [fmt(r[k]) for k in TARGET_FIELDS[1:]] for r in selected["fits"]],
        [164, 108, 108, 108, 108, 132])
    pages = parameter_pages(selected["parameters"])
    for number, rows in enumerate(pages, 1):
        add_table(f"Parameters and restrictions ({number}/{len(pages)})",
            ["Parameter", "Estimate", "Lower", "Upper", "Near bound", "Status"],
            [[r["parameter"], fmt(r["estimate"]), fmt(r["lower"]), fmt(r["upper"]),
              r["near_bound"], r["status"]] for r in rows], [192, 103, 68, 68, 62, 235])
    for name in STANDARD_NAMES:
        story.append(PageBreak())
        picture = Image(str(output / "standard_diagnostics" / (name + ".png")))
        scale = min(728 / picture.imageWidth, 530 / picture.imageHeight)
        picture.drawWidth, picture.drawHeight = picture.imageWidth * scale, picture.imageHeight * scale
        story.append(picture)
    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 8)
        canvas.drawString(32, 18, "LAYOUT FIXTURE - NO NEW CALIBRATION" if layout_fixture else
                          "Experimental four-arm recalibration | saved-packet export")
        canvas.drawRightString(760, 18, str(doc.page))
        canvas.restoreState()
    path = output / ("layout_fixture.pdf" if layout_fixture else f"utility_{arm}_review.pdf")
    # Six-point internal frame padding leaves the intended 728-point content width.
    doc = SimpleDocTemplate(str(path), pagesize=(792, 612), leftMargin=26, rightMargin=26,
                           topMargin=30, bottomMargin=30)
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    qa = output / "pdf_qa"
    qa.mkdir()
    with fitz.open(path) as rendered:
        text = "\n".join(page.get_text() for page in rendered)
        compact_text = "".join(text.split())
        for row in selected["fits"]:
            if row["moment"] not in compact_text:
                raise RuntimeError("PDF omitted target row " + row["moment"])
        for row in selected["parameters"]:
            if row["parameter"] not in compact_text:
                raise RuntimeError("PDF omitted parameter row " + row["parameter"])
        for index, page in enumerate(rendered):
            page.get_pixmap(matrix=fitz.Matrix(1.25, 1.25), alpha=False).save(qa / f"page_{index + 1:02d}.png")
        page_count = len(rendered)
    return dict(pdf_name=path.name, pdf_sha256=sha256(path), pdf_pages=page_count,
                rendered_pages=page_count, visual_review_status="pending_lead_image_inspection")


def render_layout_fixture(target_csv, parameter_csv, diagnostics_dir, output_dir):
    """Lead-invoked Torch fixture from an old 25-row floor-case report packet.

    This is intentionally separate from collect_arm and scientific verification.
    It accepts historical tables solely for layout, uses a new output directory,
    and does not import the runner or native model. Invoke the PDF operation
    marker before calling; inspect the resulting PNGs before approving layout.
    """
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError("layout fixture output already exists")
    fits = read_table(target_csv, TARGET_FIELDS, "moment")
    original = read_table(parameter_csv, PARAMETER_FIELDS, "parameter")
    if len(fits) != 13 or len(original) != 25 or "h_P" not in {r["parameter"] for r in original}:
        raise RuntimeError("layout fixture requires the supplied historical 13-target/25-parameter floor packet")
    source = Path(diagnostics_dir)
    if {p.stem for p in source.glob("*.png")} != set(STANDARD_NAMES):
        raise RuntimeError("layout fixture requires exactly the existing 17 standard PNG names")
    parameters = [dict(row) for row in original if row["parameter"] != "h_P"]
    for name in ("delta_alpha_jump", "delta_alpha", "child_benefit_exponent",
                 "utility_reference_rent", "pension_to_gross_worker_earnings"):
        parameters.append(dict(parameter=name, estimate="", lower="0.0" if name.startswith("delta_") else "",
            upper="0.25" if name.startswith("delta_") else "", near_bound="",
            status="LAYOUT PLACEHOLDER; no estimate"))
    if len(parameters) != 29:
        raise RuntimeError("layout fixture parameter count differs")
    output.mkdir(parents=True, exist_ok=False)
    destination = output / "standard_diagnostics"
    destination.mkdir()
    for name in STANDARD_NAMES:
        shutil.copyfile(source / (name + ".png"), destination / (name + ".png"))
    selected = dict(fits=fits, parameters=parameters, receipt=dict(loss=0.))
    contract = dict(arms=dict(shares_linear=dict(free_count=9)))
    report = render_report(contract, "shares_linear", selected, {}, output, layout_fixture=True)
    receipt = dict(status="layout_fixture_rendered_pending_visual_review", native_solve_count=0,
        scientific_repeat_claim=False, layout_fixture=True, target_rows=13, parameter_rows=29,
        source_target_csv=dict(path=str(Path(target_csv).resolve()), sha256=sha256(target_csv)),
        source_parameter_csv=dict(path=str(Path(parameter_csv).resolve()), sha256=sha256(parameter_csv)),
        **report)
    write_new_json(output / "layout_fixture_receipt.json", receipt)
    return receipt


def collect_arm(contract_path, arm, run_root, output_dir):
    """Export a new directory atomically; publish explicit failure if repeats fail.

    Controller contract: run_root/arm/selected.json contains original_case_output
    and repeat_case_outputs (parent directories or their case subdirectories).
    No incomplete result receives a verified PDF. No output is overwritten.
    """
    output = Path(output_dir).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    reservation = output.parent / ("." + output.name + ".reservation")
    with reservation.open("x"):
        pass
    try:
        if output.exists():
            raise FileExistsError("collection output already exists")
        with tempfile.TemporaryDirectory(prefix="." + output.name + ".", dir=output.parent) as temp:
            staging = Path(temp) / "export"
            staging.mkdir()
            arm_root = Path(run_root) / arm
            summary = arm_root / "controller_summary.json"
            if summary.is_file():
                shutil.copyfile(summary, staging / summary.name)
            try:
                contract, runner = _load_contract(contract_path)
                selection = read_json(arm_root / "selected.json")
                if selection["arm"] != arm or selection["contract_sha256"] != os.environ[ENV_PIN]:
                    raise RuntimeError("selected receipt differs from arm/contract")
                with tempfile.TemporaryDirectory(prefix="utility_export_runtime_") as scratch:
                    _, _, _, tax, _, _, runtime, _ = runner.setup(contract, arm, Path(scratch) / "runtime")
                    comparison, selected = _verify_loaded(contract, arm, selection["original_case_output"],
                        selection["repeat_case_outputs"], 2, runtime, tax)
                    if selection["point"] != selected["receipt"]["point"] or selection["loss"] != selected["receipt"]["loss"]:
                        raise RuntimeError("selected point/loss differs from original receipt")
                    graphs = _graphs(selected["case"], selected["packet"], runtime, staging)
                    for name in ("receipt.json", "target_fit.csv", "parameters.csv"):
                        shutil.copyfile(selected["case"] / name, staging / name)
                    shutil.copyfile(arm_root / "selected.json", staging / "selected.json")
                    write_new_json(staging / "repetition_verification.json", comparison)
                    report = render_report(contract, arm, selected, comparison, staging)
                result = dict(status="export_complete_pending_visual_review", arm=arm,
                    comparison_contract_sha256=os.environ[ENV_PIN], exact_repeat_claim=True,
                    repetitions_verified=2, target_rows=13, weighted_rows=12,
                    free_parameters=contract["arms"][arm]["free_count"], parameter_rows=len(selected["parameters"]),
                    standard_graph_count=17, graph_sha256=graphs, native_solve_count=0,
                    selected_case=str(selected["case"]), loss=selected["receipt"]["loss"], **report)
            except Exception as exc:
                # No retry and no partial comparison becomes a certification.
                # A failed rendering may leave intermediate files, explicitly
                # non-deliverable under this receipt.
                result = dict(status="comparison_incomplete", arm=arm, exact_repeat_claim=False,
                    error_type=type(exc).__name__, error=str(exc), automatic_retry=False,
                    output_certified=False, native_solve_count=0)
                (staging / "failure.txt").write_text(
                    "Collection incomplete; no final comparison is certified.\n" + str(exc) + "\n")
            write_new_json(staging / "collection_receipt.json", result)
            os.rename(staging, output)
        return result
    finally:
        reservation.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("compare", "smoke", "collect"), required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--arm", required=True)
    parser.add_argument("--original", type=Path)
    parser.add_argument("--repeat", type=Path, action="append", default=[])
    parser.add_argument("--required-count", type=int, choices=(1, 2))
    parser.add_argument("--run-root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "collect":
        if args.run_root is None:
            parser.error("collect requires --run-root")
        result = collect_arm(args.contract, args.arm, args.run_root, args.output)
    else:
        count = 1 if args.stage == "smoke" else args.required_count
        if args.original is None or count is None:
            parser.error("comparison requires --original and --required-count (smoke uses 1)")
        if args.stage == "smoke" and args.required_count not in (None, 1):
            parser.error("smoke requires exactly one comparison case")
        result = verify_cases(args.contract, args.arm, args.original, args.repeat, count)
        write_new_json(args.output, result)
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 2 if result["status"] == "comparison_incomplete" else 0


if __name__ == "__main__":
    raise SystemExit(main())
