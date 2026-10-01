"""Read saved q0 arrays only; no engine imports, solves, or file writes."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import resource
import time

import numpy as np

PRE_HASH = "caf760cb04c6a941cb6890510a832d5c0b61b2744f4fa4c16b89e1c8f3317903"
BUNDLE_HASH = "8189b1682de32608cf5f4236c7dec7f103d11135a0b84eb1f961411017a54fab"
CASES = ("00_reference_p1.00", "04_lifetime_repayment_only_p1.00")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def recover(prob, inclusive, k, weights=None):
    require(prob.shape[:-1] == inclusive.shape and prob.shape[-1] == 4, "Probability/value shape")
    require(np.isfinite(prob).all() and ((prob >= 0) & (prob <= 1)).all(), "Invalid probability")
    require(np.all(prob[..., 2:] == 0), "Unexpected first-birth actions")
    sums = prob[..., :2].sum(axis=-1)
    bad_sums = (sums != 0) & (np.abs(sums - 1) > 1e-12)
    if weights is None:
        require(not bad_sums.any(), "Action probabilities do not sum to one")
    else:
        require(weights.shape == sums.shape, "Probability-sum weight shape")
        require(np.all(weights[bad_sums] == 0), "Action probabilities do not sum to one on occupied matched states")
    valid = ((prob[..., 0] > 0) & (prob[..., 0] < 1)
             & (prob[..., 1] > 0) & (prob[..., 1] < 1)
             & np.isfinite(inclusive) & (inclusive > -1e9)
             & (np.abs(sums - 1) <= 1e-12))
    actions = np.full(prob.shape[:-1] + (2,), np.nan)
    actions[valid] = inclusive[valid, None] + k * np.log(prob[valid, :2])
    gap = np.full(inclusive.shape, np.nan)
    gap[valid] = k * (np.log(prob[..., 1][valid]) - np.log(prob[..., 0][valid]))
    require(np.isfinite(actions[valid]).all() and np.isfinite(gap[valid]).all(), "Nonfinite recovered value")
    return valid, actions, gap, sums


def reduce(weights, probabilities, inclusive_values, k, ages, eligible):
    require(k > 0 and np.isfinite(k), "Invalid common logit scale")
    require(np.isfinite(weights).all() and np.all(weights >= 0), "Invalid pre-fertility mass")
    require(weights.shape == inclusive_values[0].shape == inclusive_values[1].shape, "Matched-state shapes")
    require(weights.shape[3] == len(ages) == len(eligible), "Age shapes")
    recovered = [recover(p, f, k, weights) for p, f in zip(probabilities, inclusive_values)]
    sum_diagnostics = []
    for _, _, _, sums in recovered:
        nonzero = sums != 0
        error = np.abs(sums-1)
        bad = nonzero & (error > 1e-12)
        sum_diagnostics.append(dict(global_nonzero_sum_count=int(nonzero.sum()),
                                    global_bad_sum_count=int(bad.sum()),
                                    global_maximum_nonzero_sum_error=float(error[nonzero].max(initial=0)),
                                    bad_sum_baseline_pre_mass=float(weights[bad].sum()),
                                    occupied_bad_sum_count=int((bad & (weights > 0)).sum()),
                                    occupied_maximum_nonzero_sum_error=float(error[nonzero & (weights > 0)].max(initial=0))))
    rows = []
    for j, age in enumerate(ages):
        w = weights[:, :, :, j, :]
        total = float(w.sum())
        v0, a0, g0, s0 = recovered[0]
        v1, a1, g1, s1 = recovered[1]
        mask = v0[:, :, :, j, :] & v1[:, :, :, j, :] & bool(eligible[j])
        included = float(w[mask].sum())
        row = dict(age=float(age), childless_pre_mass=total, recovered_mass=included,
                   excluded_mass=total-included, excluded_fraction=None if total == 0 else 1-included/total,
                   eligible_fecund_age=bool(eligible[j]),
                   baseline_unrecoverable_mass=float(w[~v0[:, :, :, j, :]].sum()),
                   expanded_unrecoverable_mass=float(w[~v1[:, :, :, j, :]].sum()))
        if included > 0:
            wn = w[mask] / included
            dp = probabilities[1][:, :, :, j, :, 1][mask] - probabilities[0][:, :, :, j, :, 1][mask]
            dg = g1[:, :, :, j, :][mask] - g0[:, :, :, j, :][mask]
            da = a1[:, :, :, j, :, :][mask] - a0[:, :, :, j, :, :][mask]
            row.update(delta_attempt=float(wn @ dp), delta_gap=float(wn @ dg),
                       delta_wait_value=float(wn @ da[:, 0]), delta_try_value=float(wn @ da[:, 1]),
                       positive_gap_fraction=float(wn[dg > 0].sum()), negative_gap_fraction=float(wn[dg < 0].sum()),
                       zero_gap_fraction=float(wn[dg == 0].sum()),
                       positive_attempt_fraction=float(wn[dp > 0].sum()), negative_attempt_fraction=float(wn[dp < 0].sum()))
        rows.append(row)
    total = float(weights.sum())
    included = sum(r["recovered_mass"] for r in rows)
    return dict(rows=rows, childless_pre_mass=total, recovered_mass=included,
                excluded_mass=total-included, weighted_excluded_fraction=None if total == 0 else 1-included/total,
                probability_sum_tolerance=1e-12,
                probability_sum_scope="Every global violation must have exactly zero baseline PRE childless mass; recovery mask unchanged",
                probability_sum_diagnostics=sum_diagnostics,
                maximum_nonzero_probability_sum_error=max(float(np.max(np.abs(x[3][x[3] > 0]-1), initial=0)) for x in recovered))


def main(root, bundle):
    started = time.monotonic()
    require(sha(bundle) == BUNDLE_HASH, "Frozen bundle metadata changed")
    P = json.loads(bundle.read_text())["parameters"]
    require(P["sequential_births"] and not P["readiness_gate_enabled"] and not P["fertility_nest_choice"], "Unsupported choice/readiness architecture")
    require(P["J"] == 17 and P["age_start"] == 18 and P["da"] == 4, "Age clock differs")
    ages = P["age_start"] + np.arange(P["J"]) * P["da"]
    # Exact frozen get_fecundity_by_age specialization; no model import.
    omega1 = float(P["fecundity_omega1"])
    if omega1 == 0:
        fec = np.ones(P["J"])
    else:
        fec = np.clip(1-omega1*np.exp(float(P["fecundity_omega2"])*(ages-P["age_start"])), 0, 1)
        decay = float(P["fecundity_terminal_decay"])
        require(decay >= 0 and np.isfinite(decay), "Invalid fecundity decay")
        if decay > 0:
            fec *= np.exp(-decay*np.maximum(ages-float(P["fecundity_tail_start_age"]), 0))
        fec[ages >= float(P["fecundity_terminal_age"])] = 0
    eligible = ((np.arange(P["J"])+1 >= P["A_f_start"]) & (np.arange(P["J"])+1 <= P["A_f_end"]) & (fec > 0))
    with np.load(root / "q0_reference_inherited_states.npz", allow_pickle=False) as a:
        pre = a["g_pre"]
    require(hashlib.sha256(pre.tobytes()).hexdigest() == PRE_HASH, "Baseline PRE identity differs")
    require(pre.shape == (120, 6, 1, 17, 9, 4, 4), "Baseline PRE shape")
    require(np.isfinite(pre).all() and np.all(pre >= 0), "Invalid full PRE distribution")
    require(np.all(pre[..., 0, 1:] == 0), "Unexpected childless readiness-state mass")
    weights = pre[..., 0, 0].copy()
    full_mass = float(pre.sum())
    del pre
    probs, values, scales, receipts = [], [], [], []
    for case in CASES:
        folder = root / case
        receipt = json.loads((folder / "receipt.json").read_text())
        require(receipt["price_factor"] == 1.0, "Not q0")
        require(sha(folder / "parameters.csv") == receipt["parameters_sha256"], "Parameter receipt changed")
        with (folder / "parameters.csv").open() as f:
            rows = list(csv.DictReader(f))
        scales.append(float(next(r["estimate"] for r in rows if r["parameter"] == "kappa_fert")))
        with np.load(folder / "solution_arrays.npz", allow_pickle=False) as arrays:
            probs.append(arrays["fert_probs"])
            values.append(arrays["fert_value"])
            grid = arrays["b_grid"]
            price = arrays["p_eq"]
        require(grid.shape == (120,) and grid[0] == -12 and grid[-1] == 3000, "Grid identity")
        if receipts:
            require(np.array_equal(grid, prior_grid) and np.array_equal(price, prior_price), "Matched grid/price differs")
            require(receipt["parameters_sha256"] == receipts[0]["parameters_sha256"] and receipt["source_binding_sha256"] == receipts[0]["source_binding_sha256"], "Matched reference differs")
        prior_grid, prior_price = grid, price
        receipts.append(receipt)
    require(scales[0] == scales[1] == 0.17614503485474298, "Common kappa differs")
    result = reduce(weights, probs, values, scales[0], ages, eligible)
    result.update(status="saved_data_only_completed", model_calls=0,
                  script_sha256=globals().get("EXECUTED_SOURCE_SHA256") or sha(__file__),
                  source_cases=list(CASES), source_receipts=receipts, baseline_pre_sha256=PRE_HASH,
                  frozen_bundle_sha256=BUNDLE_HASH, common_kappa=scales[0], fecundity_by_age=fec.tolist(),
                  full_pre_mass=full_mass, elapsed_seconds=time.monotonic()-started,
                  maximum_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  interpretation="Matched baseline PRE childless states; conditional summaries on joint interior recoverable support. Wait/try inclusive-action values, not raw success-branch values. No policy aggregation or causal attribution.")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, required=True)
    args = parser.parse_args()
    main(args.root, args.bundle)
