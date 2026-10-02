"""Validate and tabulate the collected native strict-purchase result; no model solve."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
TARGET_FIELDS = ("moment", "role", "target", "baseline_model", "strict_model",
                 "baseline_gap", "strict_gap", "weight", "baseline_loss_contribution",
                 "strict_loss_contribution")
PARAMETER_FIELDS = ("parameter", "baseline_estimate", "strict_estimate", "lower", "upper",
                    "baseline_near_bound", "strict_near_bound", "baseline_status", "strict_status")


def close(a: object, b: object, *, tol: float = 1e-9) -> bool:
    x, y = float(a), float(b)
    return math.isfinite(x) and math.isfinite(y) and math.isclose(x, y, rel_tol=tol, abs_tol=tol)


def compare_optional_number(a: str, b: str, label: str) -> None:
    assert (a == "") == (b == ""), label
    if a != "":
        assert close(a, b, tol=1e-11), label


def fmt(value: object) -> str:
    return "" if value is None else str(value)


def table(rows: list[dict[str, str]], columns: tuple[str, ...]) -> str:
    return ("| " + " | ".join(columns) + " |\n"
            + "| " + " | ".join("---" for _ in columns) + " |\n"
            + "".join("| " + " | ".join(row[c].replace("|", "\\|") for c in columns) + " |\n"
                      for row in rows))


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--completed", type=Path, default=HERE / "collected/run/completed.json")
    parser.add_argument("--out", type=Path, default=HERE / "readout")
    args = parser.parse_args()
    assert args.completed.is_file(), f"Collected completed.json missing: {args.completed}"
    assert not args.out.exists(), f"Refusing existing readout directory: {args.out}"
    incumbent = json.loads((HERE / "incumbent.json").read_text())
    manifest = json.loads((HERE / "manifest.json").read_text())
    result = json.loads(args.completed.read_text())
    assert result["status"] == "fixed_coordinate_strict_purchase_experiment_passed"
    assert result["experimental_not_adopted"] is True
    assert result["all_ten_coordinates_fixed"] is True
    assert result["source_manifest_sha256"] == hashlib.sha256((HERE / "manifest.json").read_bytes()).hexdigest()
    assert incumbent["chain"] == manifest["reference_chain"] == 2
    assert incumbent["case"] == manifest["reference_case"] == "0064_nm"
    assert result["target_fingerprint"] == incumbent["target_fingerprint"] == manifest["target_fingerprint"]
    assert result["weight_fingerprint"] == incumbent["weight_fingerprint"] == manifest["weight_fingerprint"]
    assert len(incumbent["parameters"]) == 10
    assert len(incumbent["postcheck_target_fit"]) == len(result["target_fit"]) == 14
    assert len(incumbent["postcheck_parameters"]) == len(result["parameters"]) == 31
    assert "repeat" in result  # Native ROOT/REPEAT comparison completed before the packet was written.

    target_rows: list[dict[str, str]] = []
    for baseline, strict in zip(incumbent["postcheck_target_fit"], result["target_fit"]):
        moment = baseline["moment"]
        assert strict["moment"] == moment
        assert strict["role"] == baseline["role"], moment
        compare_optional_number(baseline["target"], strict["target"], moment + " target")
        compare_optional_number(baseline["weight"], strict["weight"], moment + " weight")
        for label, row in (("baseline", baseline), ("strict", strict)):
            if row["gap"] != "":
                assert close(float(row["model"]) - float(row["target"]), row["gap"], tol=1e-8), (moment, label, "gap")
            if row["loss_contribution"] != "":
                assert close(float(row["weight"]) * float(row["gap"]) ** 2,
                             row["loss_contribution"], tol=1e-7), (moment, label, "loss contribution")
        target_rows.append(dict(zip(TARGET_FIELDS, (
            moment, fmt(baseline["role"]), fmt(baseline["target"]), fmt(baseline["model"]),
            fmt(strict["model"]), fmt(baseline["gap"]), fmt(strict["gap"]), fmt(baseline["weight"]),
            fmt(baseline["loss_contribution"]), fmt(strict["loss_contribution"])))))
    baseline_loss = sum(float(row["loss_contribution"]) for row in incumbent["postcheck_target_fit"]
                        if row["loss_contribution"] != "")
    strict_loss = sum(float(row["loss_contribution"]) for row in result["target_fit"]
                      if row["loss_contribution"] != "")
    assert close(baseline_loss, incumbent["base_loss"], tol=1e-7)
    assert close(strict_loss, result["loss"], tol=1e-7)

    parameter_rows: list[dict[str, str]] = []
    strict_by_name = {row["parameter"]: row for row in result["parameters"]}
    assert len(strict_by_name) == 31
    assert list(strict_by_name) == [row["parameter"] for row in incumbent["postcheck_parameters"]]
    assert set(incumbent["parameters"]).issubset(strict_by_name)
    for baseline in incumbent["postcheck_parameters"]:
        name = baseline["parameter"]
        strict = strict_by_name[name]
        strict_status = ("fixed at reference estimate for this experiment"
                         if name in incumbent["parameters"] else fmt(strict["status"]))
        for bound in ("lower", "upper"):
            compare_optional_number(fmt(baseline[bound]), fmt(strict[bound]), name + " " + bound)
        if name in incumbent["parameters"]:
            assert close(strict["estimate"], incumbent["parameters"][name], tol=1e-11), name
            assert close(baseline["estimate"], incumbent["parameters"][name], tol=1e-11), name
        parameter_rows.append(dict(zip(PARAMETER_FIELDS, (
            name, fmt(baseline["estimate"]), fmt(strict["estimate"]), fmt(baseline["lower"]),
            fmt(baseline["upper"]), fmt(baseline["near_bound"]), fmt(strict["near_bound"]),
            fmt(baseline["status"]), strict_status))))

    summary = {
        "reference": "normalized_v2_chain2_0064_nm",
        "experimental_change": manifest["experimental_change"],
        "baseline_loss": baseline_loss,
        "strict_loss": strict_loss,
        "strict_minus_baseline_loss": strict_loss - baseline_loss,
        "fixed_parameter_count": 10,
        "target_rows": len(target_rows),
        "parameter_rows": len(parameter_rows),
        "target_fingerprint": result["target_fingerprint"],
        "weight_fingerprint": result["weight_fingerprint"],
        "completed_json": str(args.completed),
        "root_repeat_comparison": result["repeat"],
    }
    readme = ("# Strict purchase origination readout\n\n"
              "Experimental fixed-coordinate comparison with normalized v2 chain 2, case `0064_nm`. "
              "Only purchase origination eligibility excludes current income. The ten calibrated coordinates "
              "and target/weight fingerprints are unchanged. Initial population remains \\(N_0=1\\); "
              "the fertility-price root normalizes completed fertility to 2.1. The equilibrium price may "
              "change and \\(H_0\\) is derived separately from housing demand, so neither is held fixed.\n\n"
              f"Scored loss: baseline {baseline_loss:.12g}; strict {strict_loss:.12g}; "
              f"strict minus baseline {strict_loss-baseline_loss:+.12g}.\n\n"
              "## Full target fit\n\n" + table(target_rows, TARGET_FIELDS) + "\n"
              "## All reported parameters\n\n" + table(parameter_rows, PARAMETER_FIELDS))
    args.out.mkdir(parents=True)
    write_csv(args.out / "target_comparison.csv", TARGET_FIELDS, target_rows)
    write_csv(args.out / "parameter_comparison.csv", PARAMETER_FIELDS, parameter_rows)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (args.out / "README.md").write_text(readme)
    print(f"Wrote {len(target_rows)} target rows and {len(parameter_rows)} parameter rows to {args.out}")


if __name__ == "__main__":
    main()
