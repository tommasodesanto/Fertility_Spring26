"""Compare two verified saved stationary-GE reports; never solve the model."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
HERE = Path(__file__).resolve().parent
NEW_COLLECTION = ROOT / "output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1/collection"
OLD_COLLECTION = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection"
NEW = NEW_COLLECTION / "binary/selected_root"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: object) -> str:
    if value is None or value == "":
        return "—"
    if isinstance(value, float):
        return f"{value:.9g}"
    return str(value)


def markdown(rows: list[dict[str, object]], columns: list[str]) -> str:
    header = "| " + " | ".join(columns) + " |"
    rule = "|" + "|".join("---" for _ in columns) + "|"
    body = ["| " + " | ".join(fmt(row[c]).replace("|", "\\|") for c in columns) + " |" for row in rows]
    return "\n".join([header, rule, *body])


def numeric(row: dict[str, str], key: str) -> float | None:
    return float(row[key]) if row[key] else None


def main() -> None:
    newer = read_csv(NEW / "target_fit_new_contract.csv")
    older = read_csv(OLD_COLLECTION / "alternative_target_fit.csv")
    new_params = read_csv(NEW / "parameters_estate_a.csv")
    old_params = read_csv(OLD_COLLECTION / "alternative_parameters.csv")
    new_receipt = json.loads((NEW_COLLECTION / "final_collection.json").read_text())
    old_receipt = json.loads((OLD_COLLECTION / "collection.json").read_text())
    new_winner = new_receipt["best_by_arm"]["binary"]
    old_winner = old_receipt["winners"]["alternative"]
    assert new_winner["chain"] == 1 and old_winner["chain"] == 13
    assert len(newer) == len(older) == 14
    assert len(new_params) == len(old_params) == 31
    assert [r["moment"] for r in newer] == [r["moment"] for r in older]
    assert [r["parameter"] for r in new_params] == [r["parameter"] for r in old_params]

    moments = []
    for n, o in zip(newer, older):
        moment = n["moment"]
        assert n["role"] == o["role"] and n["weight"] == o["weight"]
        if moment != "wealth_earnings":
            assert n["target"] == o["target"]
        weight = numeric(n, "weight")
        nmodel, omodel = float(n["model"]), float(o["model"])
        ntarget, otarget = float(n["target"]), float(o["target"])
        new_native = numeric(n, "loss_contribution")
        old_native = numeric(o, "loss_contribution")
        if weight is not None:
            assert abs(new_native - weight * (nmodel - ntarget) ** 2) < 1e-8
            assert abs(old_native - weight * (omodel - otarget) ** 2) < 1e-8
        moments.append({
            "moment": moment, "role": n["role"], "new_target": ntarget,
            "old_target": otarget, "weight": weight,
            "new_model": nmodel, "new_gap_new_target": nmodel - ntarget,
            "new_loss_new_target": new_native,
            "new_gap_old_target": nmodel - otarget,
            "new_loss_old_target": None if weight is None else weight * (nmodel - otarget) ** 2,
            "old_model": omodel, "old_gap_old_target": omodel - otarget,
            "old_loss_old_target": old_native,
            "old_gap_new_target": omodel - ntarget,
            "old_loss_new_target": None if weight is None else weight * (omodel - ntarget) ** 2,
        })
    native_new = sum(r["new_loss_new_target"] or 0 for r in moments)
    native_old = sum(r["old_loss_old_target"] or 0 for r in moments)
    assert abs(native_new - new_winner["native_loss"]) < 1e-9
    assert abs(native_old - old_winner["native_loss"]) < 1e-9
    write_csv(HERE / "moments_full.csv", moments)

    bounds = json.loads((NEW_COLLECTION / "binary/provenance/start_contract.json").read_text())["bounds"]
    parameters = []
    for n, o in zip(new_params, old_params):
        name = n["parameter"]
        if name in bounds:
            assert [float(n["lower"]), float(n["upper"])] == bounds[name]
        parameters.append({
            "parameter": name,
            "role": "estimated/free" if name in bounds else ("derived" if name in ("H0", "child_benefit_CRRA_coefficient", "payroll_tax", "pension_period", "period_depreciation", "period_property_tax") else "fixed/retained"),
            "new_estimate": n["estimate"], "new_lower": n["lower"],
            "new_upper": n["upper"], "new_native_near_bound": n["near_bound"],
            "new_native_status": n["status"],
            "old_estimate": o["estimate"], "old_lower": o["lower"],
            "old_upper": o["upper"], "old_native_near_bound": o["near_bound"],
            "old_native_status": o["status"],
        })
    assert sum(r["role"] == "estimated/free" for r in parameters) == 10
    write_csv(HERE / "parameters_full.csv", parameters)
    cross_new_old = sum(r["new_loss_old_target"] or 0 for r in moments)
    cross_old_new = sum(r["old_loss_new_target"] or 0 for r in moments)
    manifest = {}
    for path in (NEW / "target_fit_new_contract.csv", NEW / "parameters_estate_a.csv",
                 OLD_COLLECTION / "alternative_target_fit.csv", OLD_COLLECTION / "alternative_parameters.csv",
                 NEW_COLLECTION / "final_collection.json", OLD_COLLECTION / "collection.json"):
        manifest[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    (HERE / "source_hashes.json").write_text(json.dumps(manifest, indent=2) + "\n")

    p = {r["parameter"]: r for r in parameters}
    lines = [
        "# One-birth Estate-A recovery versus working revised-timing anchor",
        "",
        "Reproduce from saved reports only: `python3 build_comparison.py`. This performs no model solve. Both selected points passed fresh native and exact-repeat checks; neither has an optimizer-convergence, grid-convergence, or identification-rank certificate. The Estate-A point is experimental and is not adopted as the paper baseline.",
        "",
        "## Native and cross-contract scores",
        "",
        "| Saved model moments | New wealth target 4.458387 | Old wealth target 6.926584 |",
        "|---|---:|---:|",
        f"| One-birth Estate-A chain 1 | {native_new:.9f} (native) | {cross_new_old:.9f} (arithmetic rescore) |",
        f"| Working revised-timing chain 13 | {cross_old_new:.9f} (arithmetic rescore) | {native_old:.9f} (native) |",
        "",
        "An arithmetic rescore replaces the target in `weight × (saved model moment − target)²`. It does not re-estimate parameters or solve either economy. The two rows differ in the estate utility and death-flow definition, bequest measurement, wealth target, beta search bound, and all ten estimated parameter values. Thus neither the native losses across columns nor cross-scored losses identify a causal Estate-A effect or establish which economic specification is better.",
        "",
        "The working chain-13 anchor uses the author-adopted post-interest transaction timing and the old PSID wealth/earnings target. Both points retain the same consumption/housing flow-utility form: constant consumption share `α=0.733`, physical parent housing floor `h_P>0`, and equivalence scale `e(m)=((2+0.7m)/2)^0.7`. The separate normalized-CES, parent-dependent-share, `h_P=0` experiment is not either point here. The experimental chain-1 point retains the one-birth menu, soft financing with financed share 0.8, and the other fixed economic inputs. The old estate values its housing leg at gross `Ph′`; Estate A instead uses terminal estate wealth `W = b′ + (1−ψ)Ph′` with selling cost `ψ = 0.06` and no extra interest on `b′`. The net housing leg enters utility, native death-flow accounting, and the empirical bequest observer. The new aggregate wealth/earnings target is 4.45838713455674 versus 6.92658379107299, with unchanged numerical weight 7.595098472533724. The lower search bound for annual beta changes from 0.94 to 0.93. The recovered point re-estimates all ten free coordinates. Entry wealth/income distributions, earnings, transfers/floors, housing/financing rules, other targets and weights are reported unchanged in the recovery status and source receipt; they were not individually replayed for this comparison. SCF wealth scope and estate recipient/creditor mapping remain provisional.",
        "",
        "## Complete moment comparison",
        "",
        "`new` denotes Estate A; `old` denotes the working anchor. Empty weight/loss is the normalization row; zero-weight rows are validation moments. Full precision and both cross-contract gap/loss columns are in [moments_full.csv](moments_full.csv).",
        "",
        markdown(moments, ["moment", "role", "new_target", "old_target", "weight", "new_model", "new_gap_new_target", "new_loss_new_target", "old_model", "old_gap_old_target", "old_loss_old_target"]),
        "",
        "## Complete parameter comparison",
        "",
        "The ten `estimated/free` bounds come from the Estate-A start contract and old native parameter table. `H0` is derived for household scale `N0=1`; its displayed interval is advisory, not an optimization bound. `Near` is each native report's flag, not a new threshold calculation. Full native status wording is in [parameters_full.csv](parameters_full.csv).",
        "",
        markdown(parameters, ["parameter", "role", "new_estimate", "new_lower", "new_upper", "new_native_near_bound", "old_estimate", "old_lower", "old_upper", "old_native_near_bound"]),
        "",
        "## Provenance and limits",
        "",
        f"New selected chain 1: array {new_receipt['job_id']}, inventory `{new_receipt['inventory_sha256']}`, target fingerprint `{new_receipt['target_fingerprint']}`, weight fingerprint `{new_receipt['weight_fingerprint']}`. Old selected chain 13: 48-chain collection, target fingerprint `{old_receipt['target_fingerprint']}`, weight fingerprint `{old_receipt['weight_fingerprint']}`. [Input hashes](source_hashes.json) pin the exact six saved files used here.",
        "",
        "Sources: `CALIBRATION_STATUS.md` (working anchor and Estate-A recovery sections); `output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1/collection/RESULTS.md` and `binary/provenance/start_contract.json`; `output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/collection.json`. These are local saved reports, not a new numerical verification.",
    ]
    (HERE / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
