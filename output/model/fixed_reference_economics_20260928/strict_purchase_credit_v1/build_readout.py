"""Check the collected three-case result and render full fixed-price comparisons."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRICE = 0.7152515073815459
H0 = 6.778473404808042
TARGET_COLUMNS = ("moment", "role", "target", "strict80_model", "strict90_model",
                  "strict80_gap", "strict90_gap", "weight", "strict80_loss_contribution",
                  "strict90_loss_contribution", "model_relative_change_percent")
PARAM_COLUMNS = ("parameter", "strict80_estimate", "strict90_estimate", "lower", "upper",
                 "strict80_near_bound", "strict90_near_bound", "strict80_status", "strict90_status")
FIXED = "fixed at strict-80 reference estimate for this experiment"
POLICY = "externally fixed policy input for this experiment"


def near(a: object, b: object, tol: float = 1e-10) -> bool:
    x, y = float(a), float(b)
    return math.isfinite(x) and math.isfinite(y) and math.isclose(x, y, abs_tol=tol, rel_tol=tol)


def optional_equal(a: str, b: str, label: str, tol: float = 1e-10) -> None:
    assert (a == "") == (b == ""), label
    if a != "":
        assert near(a, b, tol), label


def value(x: object) -> str:
    return "" if x is None else str(x)


def table(rows: list[dict[str, str]], columns: tuple[str, ...]) -> str:
    return ("| " + " | ".join(columns) + " |\n"
            + "| " + " | ".join("---" for _ in columns) + " |\n"
            + "".join("| " + " | ".join(row[k].replace("|", "\\|") for k in columns) + " |\n"
                      for row in rows))


def csv_write(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def check_target(rows: list[dict[str, str]], strict_reference: list[dict[str, str]]) -> float:
    assert len(rows) == len(strict_reference) == 14
    for row, ref in zip(rows, strict_reference):
        assert row["moment"] == ref["moment"] and row["role"] == ref["role"]
        optional_equal(row["target"], ref["target"], row["moment"] + " target")
        optional_equal(row["weight"], ref["weight"], row["moment"] + " weight")
        if row["gap"] != "":
            assert near(float(row["model"]) - float(row["target"]), row["gap"], 1e-8)
        if row["loss_contribution"] != "":
            assert near(float(row["weight"]) * float(row["gap"]) ** 2,
                        row["loss_contribution"], 1e-7), row["moment"]
    return sum(float(row["loss_contribution"] or 0) for row in rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--completed", type=Path, default=HERE / "collected/run/completed.json")
    parser.add_argument("--out", type=Path, default=HERE / "readout")
    args = parser.parse_args()
    assert args.completed.is_file(), f"Collected result missing: {args.completed}"
    assert not args.out.exists(), f"Refusing existing readout: {args.out}"
    result = json.loads(args.completed.read_text())
    strict = json.loads((HERE / "input/strict80/collected/run/completed.json").read_text())
    incumbent = json.loads((HERE / "input/strict80/incumbent.json").read_text())
    manifest = json.loads((HERE / "input/strict80/manifest.json").read_text())
    assert result["status"] == "fixed_price_diagnostic_passed"
    assert result["total_lifecycle_solves"] == 3
    assert result["repeat"]["status"] == "exact_full_ge_repeat_passed"
    assert len(result["repeat"]["standard_plot_hashes"]) == 17
    assert result["target_fingerprint"] == strict["target_fingerprint"] == manifest["target_fingerprint"]
    assert result["weight_fingerprint"] == strict["weight_fingerprint"] == manifest["weight_fingerprint"]
    assert near(result["fixed_price"], PRICE) and near(result["fixed_H0"], H0)
    assert near(result["fixed_population"], 1.0)
    assert len(incumbent["parameters"]) == 10 and strict["all_ten_coordinates_fixed"] is True
    cases = result["cases"]
    assert len(cases) == 3
    assert [case["label"] for case in cases] == ["strict80_replay", "strict90", "strict90_repeat"]
    assert [float(case["phi"]) for case in cases] == [0.8, 0.9, 0.9]
    for case in cases:
        assert case["status"] == "passed_fixed_price_diagnostic" and case["lifecycle_solves"] == 1
        assert near(case["price"], PRICE) and near(case["H0"], H0)
        assert len(case["parameters"]) == 31
        total = check_target(case["target_fit"], strict["target_fit"])
        assert near(total, case["loss"], 1e-7)
        assert math.isfinite(float(case["renewal_residual"]))
        assert math.isfinite(float(case["housing_residual"]))
    baseline, policy, repeat = cases
    assert len(baseline["target_fit"]) == len(strict["target_fit"])
    for left, right in zip(baseline["target_fit"], strict["target_fit"]):
        for key in ("target", "model", "gap", "weight", "loss_contribution"):
            optional_equal(left[key], right[key], left["moment"] + " strict80 replay " + key)
    assert near(baseline["loss"], strict["loss"], 1e-8)
    assert policy["target_fit"] == repeat["target_fit"]
    assert policy["parameters"] == repeat["parameters"]

    target_rows: list[dict[str, str]] = []
    for old, new in zip(baseline["target_fit"], policy["target_fit"]):
        old_model, new_model = float(old["model"]), float(new["model"])
        relative = "" if old_model == 0 else str(100 * (new_model - old_model) / old_model)
        target_rows.append(dict(zip(TARGET_COLUMNS, (
            old["moment"], old["role"], old["target"], old["model"], new["model"],
            old["gap"], new["gap"], old["weight"], old["loss_contribution"],
            new["loss_contribution"], relative))))

    parameter_rows: list[dict[str, str]] = []
    old_params, new_params = baseline["parameters"], policy["parameters"]
    assert [row["parameter"] for row in old_params] == [row["parameter"] for row in new_params]
    assert len(old_params) == len(new_params) == 31
    for old, new in zip(old_params, new_params):
        name = old["parameter"]
        optional_equal(old["lower"], new["lower"], name + " lower")
        optional_equal(old["upper"], new["upper"], name + " upper")
        if name in incumbent["parameters"]:
            assert near(old["estimate"], incumbent["parameters"][name])
            assert near(new["estimate"], incumbent["parameters"][name])
            old_status = new_status = FIXED
        elif name == "financed_share":
            assert near(old["estimate"], 0.8) and near(new["estimate"], 0.9)
            old_status = new_status = POLICY
        else:
            assert near(old["estimate"], new["estimate"]), name
            old_status, new_status = old["status"], new["status"]
        parameter_rows.append(dict(zip(PARAM_COLUMNS, (
            name, old["estimate"], new["estimate"], old["lower"], old["upper"],
            value(old["near_bound"]), value(new["near_bound"]), old_status, new_status))))

    by_moment = {row["moment"]: row for row in target_rows}
    def metric(name: str) -> dict[str, float]:
        row = by_moment[name]
        a, b = float(row["strict80_model"]), float(row["strict90_model"])
        return dict(strict80=a, strict90=b, change=b-a,
                    relative_change_percent=(100 * (b-a)/a if a != 0 else None))
    summary = dict(status="validated_fixed_price_diagnostic", cases=3,
        lifecycle_solves=3, fixed_price=PRICE, fixed_H0=H0, fixed_N0=1.0,
        financed_share=[0.8, 0.9],
        scored_loss=dict(strict80=baseline["loss"], strict90=policy["loss"],
                         change=policy["loss"]-baseline["loss"],
                         relative_change_percent=100*(policy["loss"]-baseline["loss"])/baseline["loss"]),
        completed_fertility=metric("initial_normalization"),
        early_fertility=metric("early_fertility"),
        mean_age_first_birth=metric("nchs_mean_age"),
        ownership_30_55=metric("ownership_30_55"),
        recent_parent_ownership=metric("recent_parent_ownership"),
        mean_rooms=metric("mean_rooms"),
        renewal_residual=dict(strict80=baseline["renewal_residual"], strict90=policy["renewal_residual"]),
        housing_market_residual=dict(strict80=baseline["housing_residual"], strict90=policy["housing_residual"]),
        target_fingerprint=result["target_fingerprint"], weight_fingerprint=result["weight_fingerprint"],
        exact_90_repeat=result["repeat"]["status"])
    readme = ("# Fixed-price credit diagnostic readout\n\n"
        "Strict wealth-only purchase origination, with the financed share changed from 80% to 90% "
        "for both buyers and owner-stayers. The price, housing-supply coefficient, initial population, "
        "and ten calibrated coordinates are fixed. This is a partial-equilibrium diagnostic; birth "
        "replacement and housing residuals are measured, not cleared. The 80% replay matches the "
        "completed strict-80 target fit; the 90% repeat matches the native tables, closure and 17 plots exactly.\n\n"
        f"Scored loss: {baseline['loss']:.12g} at 80%; {policy['loss']:.12g} at 90%; "
        f"change {policy['loss']-baseline['loss']:+.12g}.\n\n"
        f"Birth replacement residual: {baseline['renewal_residual']:.12g} to "
        f"{policy['renewal_residual']:.12g}. Housing-market residual: "
        f"{baseline['housing_residual']:.12g} to {policy['housing_residual']:.12g}.\n\n"
        "## Full target comparison\n\n" + table(target_rows, TARGET_COLUMNS) + "\n"
        "## All reported parameters\n\n" + table(parameter_rows, PARAM_COLUMNS))
    args.out.mkdir(parents=True)
    csv_write(args.out / "target_comparison.csv", TARGET_COLUMNS, target_rows)
    csv_write(args.out / "parameter_comparison.csv", PARAM_COLUMNS, parameter_rows)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (args.out / "README.md").write_text(readme)
    print(f"Wrote {len(target_rows)} target rows and {len(parameter_rows)} parameter rows")


if __name__ == "__main__":
    main()
