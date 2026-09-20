"""Deterministic resolution audit for the persistent/transitory income grid.

This tool performs only finite Markov algebra.  It does not simulate annual
panels, estimate a process, solve a household model, or change a target.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from statistics import NormalDist
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
CANDIDATE = ROOT / "output/model/native_financing_diagnostic_20260919/earnings_candidate/candidate.json"
MC_RECEIPT = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/income_aggregation_v1/results/full/receipt.json"
OUTPUT = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/income_grid_resolution_v1"
PERSISTENT_NODES = (5, 7, 9, 15, 25)
TRANSITORY_NODES = (3, 5)
LAGS = (0, 1, 2, 4)
QUANTILES = (0.01, 0.05, 0.50, 0.95, 0.99)
GATE_TOL = 1e-10


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def json_fingerprint(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _iid_lognormal_rule(nodes: int, log_sd: float) -> tuple[np.ndarray, np.ndarray]:
    hermite_nodes, weights = np.polynomial.hermite.hermgauss(int(nodes))
    weights = weights / math.sqrt(math.pi)
    log_var = float(log_sd) ** 2
    levels = np.exp(-0.5 * log_var + math.sqrt(2.0 * log_var) * hermite_nodes)
    levels /= float(weights @ levels)
    return levels, weights


def annual_level_covariance(v_p: float, v_e: float, rho: float, lag: int) -> float:
    distance = abs(int(lag))
    return float(np.exp(v_p * rho**distance + (v_e if distance == 0 else 0.0)) - 1.0)


def block_average_level_covariance(
    v_p: float, v_e: float, rho: float, block_size: int, block_lag: int
) -> float:
    size = int(block_size)
    if size < 1:
        raise ValueError("block_size must be positive")
    return float(
        sum(
            annual_level_covariance(v_p, v_e, rho, size * int(block_lag) + j - i)
            for i in range(size)
            for j in range(size)
        )
        / (size * size)
    )


def _load_inputs(candidate_path: Path, mc_receipt_path: Path) -> dict[str, Any]:
    candidate = json.loads(candidate_path.read_text())
    annual = candidate["annual_coefficients_recovered_from_nested_fitted_covariances"]
    mapping = candidate["four_year_diagnostic_mapping"]
    period_years = float(candidate["period_years"])
    if period_years != 4.0:
        raise ValueError("grid audit is pinned to the four-year model period")
    mc = json.loads(mc_receipt_path.read_text())
    if mc.get("status") != "completed" or not mc.get("validity", {}).get("all_exact_level_moment_checks_pass", False):
        raise ValueError("full income aggregation MC receipt is not a completed valid reference")
    return {
        "candidate": candidate,
        "rho_annual": float(annual["rho_annual"]),
        "persistent_variance": float(annual["persistent_variance"]),
        "transitory_variance": float(annual["transitory_variance"]),
        "rho_period": float(mapping["rho_period"]),
        "transitory_log_variance_period": float(mapping["transitory_log_variance_period"]),
        "period_years": 4,
        "mc": mc,
    }


def _constructor_payload(inputs: dict[str, Any]) -> dict[str, Any]:
    model_path = ROOT / "code/model"
    tools_path = ROOT / "code/model/tools"
    for path in (str(model_path), str(tools_path)):
        if path not in sys.path:
            sys.path.insert(0, path)
    from build_persistent_transitory_income_candidate import build_persistent_transitory_income_candidate
    from intergen_eqscale_seq_optimized.local_panel import income_process_fingerprint

    eta_sd = math.sqrt((1.0 - inputs["rho_annual"] ** 2) * inputs["persistent_variance"])
    overrides, metadata = build_persistent_transitory_income_candidate(
        rho_annual=inputs["rho_annual"],
        persistent_innovation_sd_annual=eta_sd,
        transitory_log_sd_period=math.sqrt(inputs["transitory_log_variance_period"]),
        period_years=4.0,
        persistent_states=5,
    )
    return {"overrides": overrides, "metadata": metadata, "fingerprint": income_process_fingerprint(overrides)}


def _grid_payload(inputs: dict[str, Any], persistent_nodes: int, transitory_nodes: int) -> dict[str, Any]:
    model_path = ROOT / "code/model"
    if str(model_path) not in sys.path:
        sys.path.insert(0, str(model_path))
    from intergen_eqscale_seq_optimized.local_panel import income_process_overrides, income_process_fingerprint

    eta_sd = math.sqrt((1.0 - inputs["rho_annual"] ** 2) * inputs["persistent_variance"])
    persistent = income_process_overrides(
        int(persistent_nodes), "rouwenhorst", eta_sd, inputs["rho_annual"]
    )
    eps_grid, eps_weights = _iid_lognormal_rule(transitory_nodes, math.sqrt(inputs["transitory_log_variance_period"]))
    z_grid = np.multiply.outer(np.asarray(persistent["z_grid"]), eps_grid).reshape(-1)
    z_weights = np.multiply.outer(np.asarray(persistent["z_weights"]), eps_weights).reshape(-1)
    pi_eps = np.broadcast_to(eps_weights, (transitory_nodes, transitory_nodes)).copy()
    pi_z = np.kron(np.asarray(persistent["Pi_z"]), pi_eps)
    overrides = {
        "z_grid": z_grid,
        "z_weights": z_weights,
        "Pi_z": pi_z,
        "income_shock_persistence": float(persistent["income_shock_persistence"]),
    }
    return {"overrides": overrides, "fingerprint": income_process_fingerprint(overrides)}


def _moments(payload: dict[str, Any]) -> dict[str, Any]:
    z = np.asarray(payload["z_grid"], dtype=float)
    w = np.asarray(payload["z_weights"], dtype=float)
    pi = np.asarray(payload["Pi_z"], dtype=float)
    w = w / w.sum()
    mean = float(w @ z)
    centered = z - mean
    logz = np.log(z)
    centered_log = logz - float(w @ logz)
    rows: dict[str, dict[str, float]] = {}
    for lag in LAGS:
        powered = np.linalg.matrix_power(pi, int(lag))
        rows[str(lag)] = {
            "level_covariance": float((w * centered) @ powered @ centered),
            "log_covariance": float((w * centered_log) @ powered @ centered_log),
        }
    return {"mean": mean, "moments": rows}


def _weighted_quantile(z: np.ndarray, weights: np.ndarray, q: float) -> float:
    order = np.argsort(z)
    sorted_z = np.asarray(z)[order]
    sorted_w = np.asarray(weights)[order]
    index = int(np.searchsorted(np.cumsum(sorted_w), q, side="left"))
    return float(sorted_z[min(index, sorted_z.size - 1)])


def _quantiles(payload: dict[str, Any]) -> dict[str, float]:
    return {
        str(q): _weighted_quantile(payload["z_grid"], payload["z_weights"], q)
        for q in QUANTILES
    }


def _continuous_endpoint_quantiles(inputs: dict[str, Any]) -> dict[str, float]:
    variance = inputs["persistent_variance"] + inputs["transitory_log_variance_period"]
    normal = NormalDist()
    return {str(q): float(math.exp(-0.5 * variance + math.sqrt(variance) * normal.inv_cdf(q))) for q in QUANTILES}


def _gates(payload: dict[str, Any], moments: dict[str, Any], inputs: dict[str, Any]) -> dict[str, Any]:
    z = np.asarray(payload["z_grid"], dtype=float)
    w = np.asarray(payload["z_weights"], dtype=float)
    pi = np.asarray(payload["Pi_z"], dtype=float)
    expected = {
        "mean": 1.0,
        "log_covariance_0": inputs["persistent_variance"] + inputs["transitory_log_variance_period"],
        "log_covariance_1": inputs["persistent_variance"] * inputs["rho_period"],
        "log_covariance_2": inputs["persistent_variance"] * inputs["rho_period"] ** 2,
    }
    actual = {
        "mean": moments["mean"],
        "log_covariance_0": moments["moments"]["0"]["log_covariance"],
        "log_covariance_1": moments["moments"]["1"]["log_covariance"],
        "log_covariance_2": moments["moments"]["2"]["log_covariance"],
    }
    gaps = {key: abs(actual[key] - expected[key]) for key in expected}
    return {
        "row_sum_max_gap": float(np.max(np.abs(pi.sum(axis=1) - 1.0))),
        "stationary_weight_max_gap": float(np.max(np.abs(w @ pi - w))),
        "mean_max_gap": float(gaps["mean"]),
        "log_moment_max_gap": float(max(gaps[k] for k in ("log_covariance_0", "log_covariance_1", "log_covariance_2"))),
        "passed": bool(
            np.all(z > 0.0)
            and np.isclose(w.sum(), 1.0, atol=GATE_TOL, rtol=0.0)
            and np.isclose(np.max(np.abs(pi.sum(axis=1) - 1.0)), 0.0, atol=GATE_TOL, rtol=0.0)
            and np.isclose(np.max(np.abs(w @ pi - w)), 0.0, atol=GATE_TOL, rtol=0.0)
            and all(gap <= GATE_TOL for gap in gaps.values())
        ),
    }


def _write_plots(output: Path, rows: list[dict[str, Any]], mc: dict[str, Any], continuous: dict[str, Any]) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [f"{r['persistent_states']}×{r['transitory_states']}" for r in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
    for lag in (0, 1, 2, 4):
        axes[0].plot(x, [r["moments"]["moments"][str(lag)]["level_covariance"] for r in rows], marker="o", label=f"lag {lag}")
    axes[0].axhline(mc["exact_block_average_level_covariances"]["0"], color="k", ls=":", lw=0.8, label="block exact lag 0")
    axes[0].set_xticks(x, labels, rotation=45, ha="right")
    axes[0].set_ylabel("level covariance")
    axes[0].set_title("Finite-grid level moments")
    axes[0].legend(frameon=False, fontsize=8)

    qlabels = [str(q) for q in QUANTILES]
    for q in qlabels:
        axes[1].plot(x, [r["quantiles"][q] for r in rows], marker="o", label=f"q={q}")
    axes[1].axhline(continuous["quantiles"]["0.99"], color="k", ls=":", lw=0.8, label="endpoint q=.99")
    axes[1].set_xticks(x, labels, rotation=45, ha="right")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("level quantile")
    axes[1].set_title("Endpoint grid tails")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    path = output / "grid_resolution_moments_quantiles.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    for q in ("0.01", "0.5", "0.99"):
        ax.plot(x, [r["quantiles"][q] for r in rows], marker="o", label=f"grid q={q}")
        ax.axhline(mc["monte_carlo_quantiles"][q]["estimate"], ls="--", lw=0.8, label=f"MC block q={q}")
    ax.set_xticks(x, labels, rotation=45, ha="right")
    ax.set_yscale("log")
    ax.set_ylabel("level quantile")
    ax.set_title("Endpoint grid versus block-average MC")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout()
    path2 = output / "grid_resolution_tail_comparison.png"
    fig.savefig(path2, dpi=150)
    plt.close(fig)
    return [path.name, path2.name]


def run(output: Path = OUTPUT, candidate_path: Path = CANDIDATE, mc_receipt_path: Path = MC_RECEIPT) -> dict[str, Any]:
    if output.exists():
        if not output.is_dir() or any(output.iterdir()):
            raise FileExistsError(f"refusing to overwrite nonempty output directory: {output}")
    else:
        output.mkdir(parents=True)
    candidate_path = candidate_path.resolve()
    mc_receipt_path = mc_receipt_path.resolve()
    constructor_path = ROOT / "code/model/tools/build_persistent_transitory_income_candidate.py"
    local_panel_path = ROOT / "code/model/intergen_eqscale_seq_optimized/local_panel.py"
    source_path = ROOT / "code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727/md_autocovariance_fit.csv"
    paths = {
        "driver": Path(__file__).resolve(),
        "candidate_json": candidate_path,
        "constructor": constructor_path,
        "local_panel": local_panel_path,
        "source_autocovariance_fit": source_path,
        "mc_reference_receipt": mc_receipt_path,
    }
    hashes_before = {name: sha256(path) for name, path in paths.items()}
    inputs = _load_inputs(candidate_path, mc_receipt_path)
    baseline = _constructor_payload(inputs)
    expected_baseline_fp = inputs["mc"]["payload_contract"]["rebuilt_15_state_payload_fingerprint"]
    baseline_reproduced = baseline["fingerprint"] == expected_baseline_fp
    if not baseline_reproduced:
        raise ValueError("5x3 production payload fingerprint does not reproduce the MC reference")
    base_moments = _moments(baseline["overrides"])
    base_gates = _gates(baseline["overrides"], base_moments, inputs)
    if not base_gates["passed"]:
        raise ValueError("5x3 production payload fails stationary/log-moment gates")
    mc = inputs["mc"]
    mc_quantiles = mc["monte_carlo"]["block_quantiles"]
    mc_bundle = {
        "monte_carlo_quantiles": mc_quantiles,
        "exact_block_average_level_covariances": mc["exact_block_average_level_covariances"],
    }
    continuous = {
        "quantiles": _continuous_endpoint_quantiles(inputs),
        "moments": {
            str(lag): {
                "level_covariance": float(
                    np.exp(
                        inputs["persistent_variance"] + inputs["transitory_log_variance_period"]
                        if lag == 0 else inputs["persistent_variance"] * inputs["rho_period"] ** lag
                    ) - 1.0
                ),
                "log_covariance": float(
                    inputs["persistent_variance"] + inputs["transitory_log_variance_period"]
                    if lag == 0 else inputs["persistent_variance"] * inputs["rho_period"] ** lag
                ),
            }
            for lag in LAGS
        },
    }
    rows: list[dict[str, Any]] = []
    for persistent_nodes in PERSISTENT_NODES:
        for transitory_nodes in TRANSITORY_NODES:
            payload = baseline["overrides"] if (persistent_nodes, transitory_nodes) == (5, 3) else _grid_payload(inputs, persistent_nodes, transitory_nodes)["overrides"]
            moments = _moments(payload)
            gates = _gates(payload, moments, inputs)
            quantiles = _quantiles(payload)
            rows.append({
                "persistent_states": persistent_nodes,
                "transitory_states": transitory_nodes,
                "joint_states": int(np.asarray(payload["z_grid"]).size),
                "moments": moments,
                "gates": gates,
                "quantiles": quantiles,
                "quantile_gap_vs_mc_block": {
                    q: quantiles[q] - float(mc_quantiles[q]["estimate"]) for q in quantiles
                },
                "quantile_gap_vs_continuous_endpoint": {
                    q: quantiles[q] - continuous["quantiles"][q] for q in quantiles
                },
            })
    all_gates = all(row["gates"]["passed"] for row in rows)
    hashes_after = {name: sha256(path) for name, path in paths.items()}
    plots = _write_plots(output, rows, mc_bundle, continuous)
    plots_present = all((output / name).is_file() for name in plots)
    result = {
        "status": "completed" if all_gates and baseline_reproduced and plots_present and hashes_before == hashes_after else "failed",
        "diagnostic": "deterministic_income_grid_resolution",
        "input_contract": {
            "persistent_nodes": list(PERSISTENT_NODES),
            "transitory_nodes": list(TRANSITORY_NODES),
            "annual_parameters": {
                "rho_annual": inputs["rho_annual"],
                "persistent_log_variance": inputs["persistent_variance"],
                "transitory_log_variance": inputs["transitory_variance"],
            },
            "period_years": inputs["period_years"],
            "rho_period": inputs["rho_period"],
            "transitory_log_variance_period": inputs["transitory_log_variance_period"],
            "no_household_solves": True,
            "no_new_monte_carlo": True,
        },
        "source_files": {name: str(path) for name, path in paths.items()},
        "source_hashes_before": hashes_before,
        "source_hashes_after": hashes_after,
        "source_hashes_unchanged": hashes_before == hashes_after,
        "baseline_5x3_reproduction": {
            "production_payload_fingerprint": baseline["fingerprint"],
            "mc_reference_payload_fingerprint": expected_baseline_fp,
            "exact_reproduction": baseline_reproduced,
            "stationary_and_log_moment_gates": base_gates,
        },
        "continuous_endpoint_reference": continuous,
        "mc_block_average_reference": mc_bundle,
        "rows": rows,
        "plots": {"files": plots, "all_present": plots_present},
    }
    write_json(output / "receipt.json", result)
    rows_for_csv = []
    for row in rows:
        flat = {
            "persistent_states": row["persistent_states"],
            "transitory_states": row["transitory_states"],
            "joint_states": row["joint_states"],
            "mean": row["moments"]["mean"],
            "gates_passed": row["gates"]["passed"],
        }
        for lag in LAGS:
            flat[f"level_cov_lag{lag}"] = row["moments"]["moments"][str(lag)]["level_covariance"]
            flat[f"log_cov_lag{lag}"] = row["moments"]["moments"][str(lag)]["log_covariance"]
        for q in QUANTILES:
            flat[f"q{q:g}"] = row["quantiles"][str(q)]
        rows_for_csv.append(flat)
    import csv
    with (output / "grid_moments.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows_for_csv[0]))
        writer.writeheader()
        writer.writerows(rows_for_csv)
    (output / "README.md").write_text(
        "# Income grid-resolution diagnostic\n\n"
        "This deterministic algebra audit compares 5/7/9/15/25 persistent states crossed with 3/5 iid-transitory quadrature nodes at the retained annual parameters. It does not run the household model, simulate new panels, estimate a process, or adopt a process. The 5x3 row reproduces the existing constructor and the prior full MC payload fingerprint.\n\n"
        "Grid `level_cov_lag*` and `log_cov_lag*` are endpoint Markov moments. The MC reference is the exact annual four-year block-average experiment already collected; its tail quantiles are therefore not directly comparable to endpoint-grid quantiles. The 25x5 variance happens to sit near the block-average variance through opposing discretization errors; refinement converges toward the continuous endpoint target, not the block-average target. The second plot displays that distinction explicitly.\n\n"
        f"Receipt status: **{result['status']}**. See `receipt.json` and `grid_moments.csv`; plots are `{', '.join(plots)}`.\n"
    )
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--candidate-json", type=Path, default=CANDIDATE)
    parser.add_argument("--mc-receipt", type=Path, default=MC_RECEIPT)
    args = parser.parse_args(argv)
    result = run(args.output, args.candidate_json, args.mc_receipt)
    return 0 if result["status"] == "completed" else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
