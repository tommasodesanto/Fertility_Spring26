#!/usr/bin/env python3
"""Audit recomputation of the July 10 overnight objective from STORED moment vectors.

Verifies, without any model solve:
  1. The active target/weight vector (candidate_replacement_post_audit_v1 +
     aggregate_mean_occupied_rooms_18_85 @ weight 6) reconstructed from
     calibration.py matches the targets/weights stored in the report fit arrays.
  2. Each stored contribution equals weight * (model - target)^2 and each
     stored gap equals model - target (exact float equality and 1e-12 tol).
  3. The sum of contributions equals the stored loss for both candidates
     (current_bound_best 14.780020699972585, housing_relaxed_best
     13.90346531862862).
  4. diagnostic_loss() applied to the stored fit moments (plus the stored
     market residual) reproduces the stored loss exactly, confirming the
     objective is sum_i w_i (m_i - t_i)^2 in raw units with no normalization
     and no residual penalty at these residuals.
  5. Per-moment loss shares for both candidates.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
sys.path.insert(0, str(ROOT / "code" / "model"))

from intergen_housing_fertility.calibration import diagnostic_loss, get_target_set  # noqa: E402

REPORT = ROOT / "output" / "model" / "combined_recalibration" / "overnight_20260710_report"
EXTRA_TARGET = {"aggregate_mean_occupied_rooms_18_85": 5.779970481941968}
EXTRA_WEIGHT = {"aggregate_mean_occupied_rooms_18_85": 6.0}
STORED_LOSS = {
    "current_bound_best": 14.780020699972585,
    "housing_relaxed_best": 13.90346531862862,
}

targets, weights = get_target_set("candidate_replacement_post_audit_v1")
targets.update(EXTRA_TARGET)
weights.update(EXTRA_WEIGHT)
assert set(targets) == set(weights)
print(f"active target system: {len(targets)} moments")

failures = []
for name in ("current_bound_best", "housing_relaxed_best"):
    rec = json.loads((REPORT / f"{name}.json").read_text())
    fit = rec["fit"]
    print(f"\n=== {name} ===")
    print(f"stored loss           : {rec['loss']!r}")
    print(f"stored contribution_sum: {rec['contribution_sum']!r}")
    print(f"stored market_residual: {rec['market_residual']!r}")
    fit_names = {row["moment"] for row in fit}
    if fit_names != set(targets):
        failures.append(f"{name}: fit moment set != active target set: "
                        f"missing={set(targets)-fit_names} extra={fit_names-set(targets)}")
    total = 0.0
    rows = []
    for row in sorted(fit, key=lambda r: r["moment"]):
        m = row["moment"]
        t_code, w_code = targets.get(m), weights.get(m)
        if t_code is not None and t_code != row["target"]:
            failures.append(f"{name}/{m}: stored target {row['target']!r} != code {t_code!r}")
        if w_code is not None and w_code != row["weight"]:
            failures.append(f"{name}/{m}: stored weight {row['weight']!r} != code {w_code!r}")
        gap = row["model"] - row["target"]
        contrib = row["weight"] * gap * gap
        if gap != row["gap"] and abs(gap - row["gap"]) > 1e-12:
            failures.append(f"{name}/{m}: gap mismatch stored {row['gap']!r} recomputed {gap!r}")
        if contrib != row["loss_contribution"] and abs(contrib - row["loss_contribution"]) > 1e-12:
            failures.append(f"{name}/{m}: contribution mismatch stored "
                            f"{row['loss_contribution']!r} recomputed {contrib!r}")
        total += contrib
        rows.append((contrib, m, gap, row["weight"]))
    exact = total == rec["loss"] == rec["contribution_sum"] == STORED_LOSS[name]
    print(f"recomputed sum        : {total!r}")
    print(f"matches stored loss exactly: {exact} (|diff|={abs(total - STORED_LOSS[name]):.3e})")
    if not exact and abs(total - STORED_LOSS[name]) > 1e-12:
        failures.append(f"{name}: recomputed sum {total!r} != stored loss {STORED_LOSS[name]!r}")

    # diagnostic_loss cross-check from a synthetic moments dict.
    moments = {row["moment"]: row["model"] for row in fit}
    moments["market_residual"] = rec["market_residual"]
    dl = diagnostic_loss(moments, targets=targets, weights=weights)
    print(f"diagnostic_loss(stored moments): {dl!r}  == stored: {dl == STORED_LOSS[name]}")
    if dl != STORED_LOSS[name] and abs(dl - STORED_LOSS[name]) > 1e-12:
        failures.append(f"{name}: diagnostic_loss {dl!r} != stored {STORED_LOSS[name]!r}")
    # residual-penalty branch check: penalty applies only above 5e-3
    moments_bad = dict(moments, market_residual=6e-3)
    assert abs(diagnostic_loss(moments_bad, targets=targets, weights=weights) - (dl + 100.0)) < 1e-9

    print(f"{'moment':62s} {'weight':>7s} {'gap':>12s} {'contrib':>10s} {'share':>7s}")
    for contrib, m, gap, w in sorted(rows, reverse=True):
        print(f"{m:62s} {w:7.1f} {gap:12.6f} {contrib:10.6f} {contrib/total:7.2%}")

print("\n=== RESULT ===")
if failures:
    print(f"{len(failures)} FAILURES:")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("ALL CHECKS PASSED: targets/weights match calibration.py; every stored "
      "contribution = w*(model-target)^2; sums equal stored losses exactly; "
      "diagnostic_loss reproduces both losses; +100 penalty branch confirmed "
      "inactive at these residuals.")
