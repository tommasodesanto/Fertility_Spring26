#!/usr/bin/env python3
"""Torch-only synthetic tests; every generated artifact is marked TEST ONLY."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
BUILDER = HERE / "build_three_price_inputs.py"
LABEL = "2007 stationary reference — block0506, September 28 verified export"
FACTORS = (0.99, 1.0, 1.01)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fixtures(folder: Path) -> tuple[Path, Path, Path, list[dict], list[dict]]:
    folder.mkdir(parents=True, exist_ok=True)
    comp = folder / "synthetic_TEST_ONLY_comparison.csv"
    elast = folder / "synthetic_TEST_ONLY_elasticities.csv"
    done = folder / "synthetic_TEST_ONLY_completed.json"
    comparison = []
    for regime in ("reference", "credit"):
        for scope in ("impact", "cohort"):
            outcomes = ("births_per_household", "first_births", "second_births") if scope == "impact" else (
                "births_per_household", "first_births", "second_births", "completed_fertility")
            for outcome in outcomes:
                base = 0.12 if outcome == "births_per_household" else (
                    0.05 if outcome == "first_births" else 0.04 if outcome == "second_births" else 2.1)
                elasticity = 0.25 if regime == "reference" else -0.15
                for factor in FACTORS:
                    value = base * factor**elasticity
                    comparison.append({"regime": regime, "scope": scope, "price_factor": factor,
                        "outcome": outcome, "value": value, "unit": "synthetic test value",
                        "prescribed_price": factor, "mapped_rent": factor})
    elasticities = []
    for regime in ("reference", "credit"):
        for outcome, scope in (("births_per_household", "impact"), ("first_births", "impact"),
                               ("second_births", "impact"), ("completed_fertility", "cohort")):
            elasticities.append({"regime": regime, "scope": scope, "outcome": outcome, "step": 0.01,
                "central_log_elasticity": 0.25 if regime == "reference" else -0.15,
                "lower_one_sided_log_elasticity": 0.25 if regime == "reference" else -0.15,
                "upper_one_sided_log_elasticity": 0.25 if regime == "reference" else -0.15})
    write_csv(comp, comparison)
    write_csv(elast, elasticities)
    completed = {"status": "passed", "reference_label": LABEL, "complete_three_price": True,
        "complete_five_price": False, "factors": list(FACTORS), "step_sizes": [0.01],
        "comparison_sha256": sha(comp), "elasticities_sha256": sha(elast)}
    done.write_text(json.dumps(completed, indent=2) + "\n", encoding="utf-8")
    return comp, elast, done, comparison, elasticities


def invoke(comp: Path, elast: Path, done: Path, out: Path, success: bool) -> subprocess.CompletedProcess:
    command = [sys.executable, str(BUILDER), "--comparison", str(comp), "--elasticities", str(elast),
               "--completed", str(done), "--output", str(out), "--test-only"]
    result = subprocess.run(command, text=True, capture_output=True, env={**os.environ, "MPLBACKEND": "Agg"})
    assert (result.returncode == 0) == success, result.stderr[-1500:]
    return result


def update_completion(comp: Path, elast: Path, done: Path, **updates) -> None:
    value = json.loads(done.read_text(encoding="utf-8"))
    value.update(comparison_sha256=sha(comp), elasticities_sha256=sha(elast), **updates)
    done.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def expect_rejection(root: Path, name: str, comp_rows: list[dict], elast_rows: list[dict],
                     completion_updates: dict | None = None) -> None:
    folder = root / name
    comp, elast, done, _, _ = fixtures(folder)
    write_csv(comp, comp_rows)
    write_csv(elast, elast_rows)
    update_completion(comp, elast, done, **(completion_updates or {}))
    out = folder / "rejected_output"
    invoke(comp, elast, done, out, success=False)
    assert not out.exists(), f"Rejected case unexpectedly wrote output: {name}"


def main() -> None:
    root = Path(sys.argv[1]).resolve()
    assert not root.exists(), f"Refusing to overwrite test folder: {root}"
    comp, elast, done, comp_rows, elast_rows = fixtures(root / "positive_TEST_ONLY")
    positive_out = root / "positive_TEST_ONLY" / "rendered_TEST_ONLY"
    invoke(comp, elast, done, positive_out, success=True)
    manifest = json.loads((positive_out / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["synthetic_test_only"] is True
    assert all((positive_out / name).is_file() and (positive_out / name).stat().st_size > 0
               for name in ("price_response.png", "price_response.pdf", "local_elasticities.csv", "local_elasticities.tex"))
    assert LABEL in (positive_out / "local_elasticities.csv").read_text(encoding="utf-8")

    missing = [r for r in comp_rows if not (r["regime"] == "reference" and r["scope"] == "impact" and
             r["outcome"] == "births_per_household" and float(r["price_factor"]) == 0.99)]
    expect_rejection(root, "missing_price_TEST_ONLY", missing, elast_rows)

    wrong_regime = [{**r, "regime": "unrecognized"} if r["regime"] == "credit" else r for r in comp_rows]
    expect_rejection(root, "wrong_regime_TEST_ONLY", wrong_regime, elast_rows)

    mismatched_elasticities = [dict(r) for r in elast_rows]
    mismatched_elasticities[0]["central_log_elasticity"] += 0.1
    expect_rejection(root, "wrong_slope_TEST_ONLY", comp_rows, mismatched_elasticities)

    folder = root / "wrong_hash_TEST_ONLY"
    c2, e2, d2, _, _ = fixtures(folder)
    payload = json.loads(d2.read_text(encoding="utf-8")); payload["comparison_sha256"] = "0" * 64
    d2.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    invoke(c2, e2, d2, folder / "rejected_output", success=False)
    assert not (folder / "rejected_output").exists()

    folder = root / "wrong_elasticity_hash_TEST_ONLY"
    c3, e3, d3, _, _ = fixtures(folder)
    payload = json.loads(d3.read_text(encoding="utf-8")); payload["elasticities_sha256"] = "0" * 64
    d3.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    invoke(c3, e3, d3, folder / "rejected_output", success=False)
    assert not (folder / "rejected_output").exists()

    expect_rejection(root, "five_price_claim_TEST_ONLY", comp_rows, elast_rows,
                     {"complete_five_price": True, "factors": [0.98, 0.99, 1.0, 1.01, 1.02]})

    report = {"status": "passed", "synthetic_test_only": True,
              "cases": ["three-price chart/table emitted with TEST ONLY marking", "missing factor rejected",
                        "wrong regime label rejected", "comparison-elasticity mismatch rejected",
                        "wrong comparison hash rejected",
                        "wrong elasticity hash rejected", "five-price claim rejected"]}
    (root / "synthetic_TEST_ONLY_test_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
