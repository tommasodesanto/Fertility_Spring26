#!/usr/bin/env python3
"""Compare three saved selected points without solving or altering targets."""

import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
BASE = ROOT / "output/model/fixed_reference_economics_20260928"
TIMING = BASE / "soft_timing_calibration_20261002_v1/collection"
WEALTH = BASE / "alternative_wealth_local_20261003_v1/collection"
SOURCES = {
    "original_timing_old_target": (TIMING / "original_target_fit.csv", TIMING / "original_parameters.csv"),
    "alternative_timing_old_target": (TIMING / "alternative_target_fit.csv", TIMING / "alternative_parameters.csv"),
    "alternative_timing_new_wealth_target": (WEALTH / "winner_target_fit.csv", WEALTH / "winner_parameters.csv"),
}


def records(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write(name, rows, columns):
    with (HERE / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, columns)
        writer.writeheader()
        writer.writerows(rows)


def number(value):
    return float(value) if value != "" else None


def main():
    timing_receipt = json.loads((TIMING / "collection.json").read_text())
    wealth_receipt = json.loads((WEALTH / "verification.json").read_text())
    assert timing_receipt["status"] == "complete"
    assert timing_receipt["counts"] == {"no_admissible_candidate": 2, "verified": 46}
    assert timing_receipt["winners"]["original"]["chain"] == 15
    assert timing_receipt["winners"]["alternative"]["chain"] == 13
    assert wealth_receipt["status"] == "all_ten_selected_numerically_verified"
    assert wealth_receipt["count"] == 10 and wealth_receipt["winner_chain"] == 2
    assert wealth_receipt["optimizer_convergence_certified"] is False
    assert timing_receipt["target_fingerprint"] != wealth_receipt["target_fingerprint"]

    fits = {arm: records(paths[0]) for arm, paths in SOURCES.items()}
    parameters = {arm: records(paths[1]) for arm, paths in SOURCES.items()}
    names = [row["moment"] for row in next(iter(fits.values()))]
    assert len(names) == 14 and len(set(names)) == 14
    assert all([row["moment"] for row in fit] == names for fit in fits.values())
    assert all(len(p) == 31 for p in parameters.values())
    old_target = number(fits["original_timing_old_target"][names.index("wealth_earnings")]["target"])
    new_target = number(fits["alternative_timing_new_wealth_target"][names.index("wealth_earnings")]["target"])
    weight = number(fits["original_timing_old_target"][names.index("wealth_earnings")]["weight"])
    assert abs(old_target - 6.92658379107299) < 1e-12
    assert abs(new_target - 4.45838713455674) < 1e-12
    assert abs(weight - 7.595098472533724) < 1e-12

    rows = []
    summary = {}
    for arm, fit in fits.items():
        native_loss = sum(number(r["loss_contribution"]) or 0 for r in fit)
        old_score = 0.0
        new_score = 0.0
        for row in fit:
            moment = row["moment"]
            target, model = number(row["target"]), number(row["model"])
            row_weight = number(row["weight"])
            contribution = number(row["loss_contribution"])
            if row_weight is None:
                assert contribution is None
            else:
                assert abs((contribution or 0) - row_weight * (model - target) ** 2) < 1e-8
            if moment != "wealth_earnings":
                ref = fits["original_timing_old_target"][names.index(moment)]
                assert target == number(ref["target"])
                assert row_weight == number(ref["weight"])
            old_contrib = (row_weight or 0) * (model - (old_target if moment == "wealth_earnings" else target)) ** 2
            new_contrib = (row_weight or 0) * (model - (new_target if moment == "wealth_earnings" else target)) ** 2
            old_score += old_contrib
            new_score += new_contrib
            rows.append({"arm": arm, **row, "old_contract_target": old_target if moment == "wealth_earnings" else target,
                         "old_contract_contribution": old_contrib,
                         "new_contract_target": new_target if moment == "wealth_earnings" else target,
                         "new_contract_contribution": new_contrib})
        expected = old_score if arm != "alternative_timing_new_wealth_target" else new_score
        assert abs(native_loss - expected) < 1e-8
        receipt_loss = (wealth_receipt["winner_native_loss"] if arm == "alternative_timing_new_wealth_target"
                        else timing_receipt["winners"]["original" if arm == "original_timing_old_target"
                                                        else "alternative"]["native_loss"])
        assert abs(native_loss - receipt_loss) < 1e-8
        summary[arm] = {"native_loss": native_loss, "common_old_target_score": old_score,
                        "common_new_target_score": new_score,
                        "wealth_model": number(fit[names.index("wealth_earnings")]["model"]),
                        "scored_moments": sum(r["role"] == "scored" for r in fit),
                        "free_parameters": sum("free in" in r["status"] for r in parameters[arm])}
        assert summary[arm]["scored_moments"] == summary[arm]["free_parameters"] == 10

    write("complete_moment_comparison.csv", rows, list(rows[0]))
    param_rows = []
    for arm, params in parameters.items():
        for row in params:
            param_rows.append({"arm": arm, **row})
    write("complete_parameter_comparison.csv", param_rows, list(param_rows[0]))
    digest = {}
    for arm, paths in SOURCES.items():
        digest[arm] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    for path in [TIMING / "collection.json", WEALTH / "verification.json",
                 ROOT / "output/model/wealth_numerator_match_20261002/wealth_numerator_match_results.csv"]:
        digest[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    (HERE / "source_fingerprints.json").write_text(json.dumps(digest, indent=2) + "\n")
    (HERE / "score_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
