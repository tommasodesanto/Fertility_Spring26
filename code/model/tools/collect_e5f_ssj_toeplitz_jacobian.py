"""Summarize the block-Toeplitz initial-Jacobian experiment from mirrored artifacts.

Reads ``derivative/derivative_receipt.json``, ``derivative/jacobian.json`` and,
when present, ``comparison.json`` from the local packet, plus the reference
ten-period root receipt.  Writes ``summary.md`` and
``reference_path_prediction_check.json`` next to them.  No model imports.

Usage:
    python collect_e5f_ssj_toeplitz_jacobian.py --packet <dir> --reference <root_receipt.json>
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def prediction_check(reference, jacobian, horizon, slope):
    """Relative error of predicted residual changes along the reference Broyden path."""
    default = np.diag(np.concatenate([np.full(horizon, -slope), np.full(2 * horizon, -200.0)]))
    rows = []
    history = reference["history"]
    for a, b in zip(history[:-1], history[1:]):
        if b.get("phase") != "iterate":
            continue
        dx = np.log(np.asarray(b["prices"], dtype=float)) - np.log(np.asarray(a["prices"], dtype=float))
        dy = np.asarray(b["residual"], dtype=float) - np.asarray(a["residual"], dtype=float)
        def err(J, sl=slice(None)):
            return float(np.linalg.norm(dy[sl] - (J @ dx)[sl]) / max(np.linalg.norm(dy[sl]), 1e-300))
        rows.append(dict(evaluation=b["evaluation"], diagonal_default=err(default), measured_toeplitz=err(jacobian),
                         measured_housing_block=err(jacobian, slice(0, horizon)),
                         measured_rebate_block=err(jacobian, slice(2 * horizon, 3 * horizon))))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--packet", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--derivative-from", type=Path, default=None,
                    help="packet whose derivative stage was reused (its receipt and Jacobian are read)")
    args = ap.parse_args()
    source = args.derivative_from or args.packet
    receipt = json.loads((source / "derivative" / "derivative_receipt.json").read_text())
    jacobian = np.asarray(json.loads((source / "derivative" / "jacobian.json").read_text())["jacobian"])
    reference = json.loads(args.reference.read_text())
    horizon = int(receipt["horizon"])
    check = prediction_check(reference, jacobian, horizon, 1.63)
    (args.packet / "reference_path_prediction_check.json").write_text(json.dumps(check, indent=2) + "\n")
    lines = ["# Block-Toeplitz initial Jacobian experiment (ten dates)", "",
             (f"Derivative stage reused from `{source}`." if args.derivative_from else "Derivative stage measured in this packet."), "",
             f"Derivative stage: {receipt['mappings']} native mappings, {receipt['elapsed_seconds']:.0f} s "
             f"(per mapping: {', '.join(f'{s:.0f}' for s in receipt['mapping_seconds'])} s).",
             f"Baseline scaled residual {receipt['baseline_max_abs']:.4e} (gate {receipt['baseline_gate']}); "
             f"stationary drift max {max(receipt['stationary_drift'].values()):.2e} (limit {receipt['drift_limit']}).",
             f"Step {receipt['step']} in log coordinates at date {receipt['perturbed_date']}; measured lags {receipt['measured_lags']}.",
             "No fake-news derivatives were constructed; lags outside the window are zero.", "",
             "## Lag profiles (d residual_t / d log u_s, lag = t - s)", "",
             "| block | " + " | ".join(str(k) for k in receipt["measured_lags"]) + " |",
             "|---|" + "---|" * len(receipt["measured_lags"])]
    for name, values in receipt["lag_profiles"].items():
        lines.append(f"| {name} | " + " | ".join(f"{v:.3g}" for v in values) + " |")
    lines += ["", "## Prediction of residual changes along the reference Broyden path (relative error)", "",
              "| evaluation | diagonal default | measured Toeplitz | housing block | rebate block |", "|---|---|---|---|---|"]
    for r in check:
        lines.append(f"| {r['evaluation']} | {r['diagonal_default']:.3f} | {r['measured_toeplitz']:.3f} | "
                     f"{r['measured_housing_block']:.3f} | {r['measured_rebate_block']:.3f} |")
    comparison = args.packet / "comparison.json"
    if comparison.exists():
        c = json.loads(comparison.read_text())
        lines += ["", "## Root comparison (identical start, controls, endpoint, eight mappings)", "",
                  f"Reference: status `{c['reference_status']}`, converged {c['reference_converged']}, best score {c['reference_best_score']:.4g}, {c['reference_elapsed_seconds']:.0f} s.",
                  f"Measured initial Jacobian: status `{c['new_status']}`, converged {c['new_converged']}, best score {c['new_best_score']:.4g}, {c['new_elapsed_seconds']:.0f} s.", "",
                  "| eval | phase | ref score | ref housing | ref PAYGO | ref rebate | ref s | new score | new housing | new PAYGO | new rebate | new s | note |",
                  "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        new = {r["evaluation"]: r for r in c["new_history"]}
        for r in c["reference_history"]:
            n = new.get(r["evaluation"])
            right = (f"{n['score']:.3e} | {n['max_housing']:.2e} | {n['max_paygo']:.2e} | {n['max_rebate']:.2e} | {n['evaluation_seconds']:.0f} | {n.get('safeguard') or n.get('reset_reason') or ''}"
                     if n else " | | | | | ")
            lines.append(f"| {r['evaluation']} | {r['phase']} | {r['score']:.3e} | {r['max_housing']:.2e} | {r['max_paygo']:.2e} | {r['max_rebate']:.2e} | {r['evaluation_seconds']:.0f} | {right} |")
    else:
        lines += ["", "Root comparison not yet collected."]
    (args.packet / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
