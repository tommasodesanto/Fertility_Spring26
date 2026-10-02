"""Selected quarter fit: one-date fixed-price financing response.

This driver makes at most two native household Bellman calls. It never solves
markets or advances the distribution. --preflight performs the same source and
array setup but stops before either Bellman call.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PACKET = HERE.parents[1]
BOOT = PACKET / "local_runtime/bootstrap.py"
COMPLETED = PACKET / "local_runtime/runs/local10_v1/chain54/postcheck/completed.json"
PRICE_RECEIPT = PACKET / "independent_diagnosis/value_decomposition_receipt.json"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def boot_overlay():
    # This is the frozen selected-postcheck read-only overlay. Its final branch
    # launches a search, so execute only the setup prefix as in quarter_overlap.
    prefix = BOOT.read_text().split("if '--preflight-context' in sys.argv:", 1)[0]
    if "DIGESTS=" not in prefix or "MAPPING=" not in prefix:
        raise RuntimeError("Frozen overlay structure changed")
    exec(compile(prefix, str(BOOT), "exec"),
         {"__file__": str(BOOT), "__name__": "quarter_fixedprice_boot"})


def age_flow(g, fert, fec, P):
    age_rows = []
    for j in range(P.J):
        age = float(P.age_start + P.da * j)
        if not (P.A_f_start <= j + 1 <= P.A_f_end):
            continue
        # g_pre is indexed by wealth, tenure, location, age, income,
        # children ever born, children at home. fert is attempt probability.
        mass = np.asarray(g[:, :, :, j, :, 0, 0])
        prob = np.asarray(fert[:, :, :, j, :, 1])
        if mass.shape != prob.shape:
            raise RuntimeError("Birth risk set and attempt policy shapes differ")
        if np.any(prob < -1e-12) or np.any(prob > 1 + 1e-12):
            raise RuntimeError("Invalid first-birth attempt probability")
        by_tenure = np.sum(mass * prob * fec[j], axis=(0, 2, 3))
        risk_by_tenure = np.sum(mass, axis=(0, 2, 3))
        age_rows.append(dict(age=age, flow_by_tenure=by_tenure.tolist(),
                             risk_by_tenure=risk_by_tenure.tolist(),
                             flow=float(by_tenure.sum()), risk=float(risk_by_tenure.sum())))
    flow = float(sum(row["flow"] for row in age_rows))
    risk = float(sum(row["risk"] for row in age_rows))
    return dict(flow=flow, fertile_childless_risk_mass=risk,
                fertile_childless_hazard=flow / risk, by_age=age_rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--smoke", action="store_true",
                        help="exercise both solve branches with saved fertility array, no Bellman calls")
    parser.add_argument("--price-factorial", action="store_true",
                        help="phi 80/100 at accepted temporary-path date-zero price and rent")
    parser.add_argument("--price-shapley", action="store_true",
                        help="phi 100 crossed asset-price and rent cells")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(int(os.environ.get(k, "1")) != 1 for k in
           ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")):
        raise RuntimeError("Exactly one numerical thread required")
    boot_overlay()
    sys.path.insert(0, str(PACKET / "mechanism"))
    import selected_runtime

    selected, contract, root, repeat = selected_runtime.authenticate_selected("quarter", COMPLETED)
    rt = selected_runtime.construct("quarter", COMPLETED,
                                    Path(tempfile.mkdtemp(prefix="quarter_fixedprice_auth_")) / "runtime")
    from refactor_lab.engine.parameters import get_fecundity_by_age
    stage = repeat / "stage/solution_arrays.npz"
    with np.load(stage, allow_pickle=False) as z:
        V, g, fert = (z[k] for k in ("V", "distribution.g_pre", "fert_probs"))
        b, price = (z[k] for k in ("b_grid", "p_eq"))
        shared = {k: z["shared." + k] for k in
                  ("cb_flat", "hb_flat", "psi_flat", "gb_flat", "alpha_flat", "escale_flat", "phi_choice")}
    P80 = rt.P
    SD80 = rt.model.precompute_shared(P80, b)
    if not np.array_equal(b, rt.grid) or not np.allclose(price, [rt.reference_price], rtol=0, atol=1e-12):
        raise RuntimeError("Saved grid or price differs from authenticated fit")
    for key, old in shared.items():
        if not np.allclose(np.asarray(getattr(SD80, key)), old, rtol=0, atol=1e-12):
            raise RuntimeError("Authenticated shared array differs: " + key)
    if not np.array_equal(np.asarray(P80.phi), np.full_like(np.asarray(P80.phi), 0.8)):
        raise RuntimeError("Selected baseline financing is not 80%")
    fec = get_fecundity_by_age(P80)
    saved = age_flow(g, fert, fec, P80)
    expected = 0.04974551614189717
    if abs(saved["flow"] - expected) > 1e-10:
        raise RuntimeError("Saved first-birth flow differs from native accepted control")
    receipt = dict(status="preflight_passed", selected_loss=selected["selected"]["loss"],
                   selected_receipt_sha256=sha(COMPLETED), stage_sha256=sha(stage),
                   bootstrap_sha256=sha(BOOT), saved_baseline=saved, price=float(price[0]),
                   source_engine_pins_sha256=sha(PACKET / "engine_pins.json"),
                   source_integration_pins_sha256=sha(PACKET / "source_pins.json"),
                   bellman_calls=0, elapsed_seconds=time.time() - started)
    write(out / "progress.json", receipt)
    if args.preflight:
        return

    def solve(P, SD, label, current_price, rent):
        if args.smoke:
            new_fert = fert.copy()
        else:
            objects = rt.model.solve_bellman_full_markov_income(
                np.array([rent]), np.array([current_price]), P, b, SD, continuation_V=V)
            new_fert = np.asarray(objects[7])
        result = age_flow(g, new_fert, fec, P)
        write(out / "progress.json", dict(status=label + "_complete",
                                          bellman_calls=0 if args.smoke else (1 if label == "control" else 2),
                                          elapsed_seconds=time.time() - started,
                                          flow=result["flow"]))
        return result, new_fert, rent

    P100 = copy.deepcopy(P80)
    P100.phi = np.full_like(np.asarray(P80.phi), 1.0)
    SD100 = rt.model.precompute_shared(P100, b)
    if args.price_factorial or args.price_shapley:
        source = json.loads(PRICE_RECEIPT.read_text())
        cases = [c for c in source["cases"] if c["arm"] == "quarter" and c["horizon"] == 48]
        if len(cases) != 1:
            raise RuntimeError("Unique accepted quarter h48 price comparison unavailable")
        case = cases[0]
        if (not source["checks"]["same_inherited_g_pre"]
                or not source["checks"]["first_birth_flow_reproduced_from_saved_state"]
                or case["policy_sha256"] != "3e7eb543c73d8d840b54ccad3d6a1b9eec820fc556e2e4102249bbc0a750f5b5"):
            raise RuntimeError("Saved current-price receipt identity differs")
        base_result_path = HERE / "local_run/production/result.json"
        baseline = json.loads(base_result_path.read_text())
        p0 = float(price[0])
        r0 = float(P80.user_cost_rate * p0)
        if (baseline.get("status") != "computed" or baseline.get("bellman_calls") != 2
                or baseline["source_receipt"]["selected_receipt_sha256"] != sha(COMPLETED)
                or abs(baseline["price"] - p0) > 1e-12
                or abs(baseline["rent"] - r0) > 1e-12
                or abs(baseline["control"]["flow"] - saved["flow"]) > 1e-10):
            raise RuntimeError("Saved baseline factorial corner failed source/price check")
        current_price = p0 * (1 + float(case["asset_price_percent_change"]) / 100)
        rent = r0 * (1 + float(case["rent_percent_change"]) / 100)
        if not (abs(current_price - 0.6764414785119185) < 1e-11
                and abs(rent - 0.12403942214723254) < 1e-11):
            raise RuntimeError("Derived accepted date-zero price or rent differs")
        if args.price_shapley:
            factorial_path = HERE / "local_run/factorial/result.json"
            factorial = json.loads(factorial_path.read_text())
            if (factorial.get("status") != "computed" or factorial.get("bellman_calls") != 2
                    or factorial["baseline_result_sha256"] != sha(base_result_path)
                    or factorial["price_receipt_sha256"] != sha(PRICE_RECEIPT)
                    or abs(factorial["p1"] - current_price) > 1e-12
                    or abs(factorial["r1"] - rent) > 1e-12
                    or abs(factorial["baseline_p0_r0_phi100"]["flow"] - baseline["policy"]["flow"]) > 1e-12):
                raise RuntimeError("Saved price-factorial corners failed authentication")
            asset_only, _, _ = solve(P100, SD100, "asset_only", current_price, r0)
            rent_only, _, _ = solve(P100, SD100, "rent_only", p0, rent)
            f00 = baseline["policy"]["flow"]
            f10 = asset_only["flow"]
            f01 = rent_only["flow"]
            f11 = factorial["policy_p1_r1_phi100"]["flow"]
            asset_shapley = 0.5 * ((f10 - f00) + (f11 - f01))
            rent_shapley = 0.5 * ((f01 - f00) + (f11 - f10))
            if abs((asset_shapley + rent_shapley) - (f11 - f00)) > 1e-14:
                raise RuntimeError("Price Shapley contributions do not add up")
            result = dict(status="smoke_passed_no_bellman_calls" if args.smoke else "computed",
                          definition="Two crossed current-price cells at phi100, saved g_pre and V80",
                          baseline_p0_r0_phi100=baseline["policy"],
                          asset_only_p1_r0_phi100=asset_only,
                          rent_only_p0_r1_phi100=rent_only,
                          joint_p1_r1_phi100=factorial["policy_p1_r1_phi100"],
                          raw_asset_effect_at_r0=f10-f00,
                          raw_rent_effect_at_p0=f01-f00,
                          raw_asset_effect_at_r1=f11-f01,
                          raw_rent_effect_at_p1=f11-f10,
                          asset_shapley=asset_shapley, rent_shapley=rent_shapley,
                          joint_current_price_effect=f11-f00,
                          price_receipt_sha256=sha(PRICE_RECEIPT),
                          baseline_result_sha256=sha(base_result_path),
                          factorial_result_sha256=sha(factorial_path),
                          accepted_policy_packet_sha256=case["policy_sha256"],
                          p0=p0, r0=r0, p1=current_price, r1=rent,
                          source_receipt=receipt, bellman_calls=0 if args.smoke else 2,
                          elapsed_seconds=time.time()-started)
            write(out / "result.json", result)
            write(out / "progress.json", dict(status="complete", bellman_calls=0 if args.smoke else 2,
                                              elapsed_seconds=time.time()-started))
            return
        control, _, _ = solve(P80, SD80, "price_control", current_price, rent)
        policy, _, _ = solve(P100, SD100, "price_policy", current_price, rent)
        result = dict(status="smoke_passed_no_bellman_calls" if args.smoke else "computed",
                      definition="Two-by-two current-price and financed-share factorial with fixed g_pre and V80",
                      baseline_p0_r0_phi80=baseline["control"],
                      baseline_p0_r0_phi100=baseline["policy"],
                      policy_p1_r1_phi80=control, policy_p1_r1_phi100=policy,
                      p0=p0, r0=r0, p1=current_price, r1=rent,
                      price_receipt_sha256=sha(PRICE_RECEIPT),
                      accepted_policy_packet_sha256=case["policy_sha256"],
                      accepted_control_packet_sha256=case["control_sha256"],
                      baseline_result_sha256=sha(base_result_path),
                      current_price_effect_at_phi80=control["flow"]-baseline["control"]["flow"],
                      current_price_effect_at_phi100=policy["flow"]-baseline["policy"]["flow"],
                      financing_effect_at_p0_r0=baseline["policy"]["flow"]-baseline["control"]["flow"],
                      financing_effect_at_p1_r1=policy["flow"]-control["flow"],
                      accepted_temporary_ge_policy_flow=0.0495313221352156,
                      transition_remainder_at_policy_price=0.0495313221352156-policy["flow"],
                      source_receipt=receipt,
                      bellman_calls=0 if args.smoke else 2,
                      elapsed_seconds=time.time()-started)
        write(out / "result.json", result)
        write(out / "progress.json", dict(status="complete", bellman_calls=0 if args.smoke else 2,
                                          elapsed_seconds=time.time()-started))
        return

    rent = float(P80.user_cost_rate * price[0])
    control, fert_control, _ = solve(P80, SD80, "control", float(price[0]), rent)
    max_prob_gap = float(np.max(np.abs(fert_control - fert)))
    if max_prob_gap > 1e-6 or abs(control["flow"] - saved["flow"]) > 1e-8:
        write(out / "result.json", dict(status="control_reproduction_failed",
                                       max_attempt_probability_gap=max_prob_gap,
                                       saved=saved, control=control, bellman_calls=1))
        raise RuntimeError("Native one-date control fails saved baseline reproduction")
    policy, _, _ = solve(P100, SD100, "policy", float(price[0]), rent)
    ge_flow = 0.0495313221352156
    result = dict(status="computed", definition="One date phi 0.8 to 1 at saved prices and saved next-age values",
                  economic_change="experimental: current-date financed share only, all four owner-state entries including stayers",
                  control=control, policy=policy, saved=saved, rent=rent, price=float(price[0]),
                  fixed_price_absolute_change=policy["flow"]-control["flow"],
                  fixed_price_relative_change_percent=100*(policy["flow"]/control["flow"]-1),
                  accepted_temporary_ge_date0_flow=ge_flow,
                  fixed_price_minus_ge_policy_flow=policy["flow"]-ge_flow,
                  max_control_attempt_probability_gap=max_prob_gap,
                  bellman_calls=0 if args.smoke else 2, elapsed_seconds=time.time()-started,
                  source_receipt=receipt)
    if args.smoke:
        result["status"] = "smoke_passed_no_bellman_calls"
    write(out / "result.json", result)
    write(out / "progress.json", dict(status="complete", bellman_calls=0 if args.smoke else 2,
                                      elapsed_seconds=time.time()-started))


if __name__ == "__main__":
    main()
