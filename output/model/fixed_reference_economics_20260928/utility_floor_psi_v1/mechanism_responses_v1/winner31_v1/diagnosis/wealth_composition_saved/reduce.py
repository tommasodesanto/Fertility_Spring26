"""Winner31 q0 saved-data accounting; no model solves or policy extrapolation."""
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
from types import SimpleNamespace

import numpy as np

PRE_HASH = "c6589a65ad74e624579c7abe2e4069ab4e07712f299496deefda6ec7a6d5e48a"
FIELDS = ("V", "c_pol", "hR_pol", "bp_pol", "tenure_choice", "tenure_probs",
          "loc_probs", "fert_probs", "fert_value", "fert2_probs",
          "g_beginning_distribution", "entry_by_loc", "p_eq", "b_grid")


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def forbidden(*args, **kwargs):
    raise RuntimeError("Household, lifecycle and GE solving forbidden")


def save_csv(path, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def origin_flows(native_apply, pre, policy, P):
    """Native birth operator on each origin family cell; keep wealth and g axes."""
    out = {}
    for n in range(3):
        for m in range(n + 1):
            masked = np.zeros_like(pre)
            masked[..., n, m] = pre[..., n, m]
            post, births, _ = native_apply(masked, policy.fert_probs, P, policy.fert2_probs)
            flow = post[..., n + 1:, :].sum(axis=(-2, -1))
            exposure = masked[..., n, m]
            require(flow.shape == exposure.shape, "Flow/exposure shape")
            require(float(flow.min()) >= -2e-12, "Negative native flow")
            require(abs(float(flow.sum()) - float(births)) < 2e-10, "Native origin flow")
            out[n, m] = (exposure, flow)
    return out


def main(packet, root, out):
    started = time.monotonic()
    require(not out.exists(), "Refuse existing output")
    out.mkdir(parents=True)
    sys.path.insert(0, str(packet))
    import fixed_price_responses as driver
    auth = driver.authenticate_candidate(out / "runtime_auth")
    cal = auth["context"]["prepared"].rt["primitive"].pf.calendar
    transition = auth["context"]["prepared"].rt["primitive"].pf.transition
    require(cal.apply_fertility is transition.apply_sequential_fertility, "Native birth operator")
    blocked = []
    for module in (auth["solver"], auth["context"]["prepared"].rt["model"], cal):
        for name in dir(module):
            if name.startswith("solve") and callable(getattr(module, name)):
                setattr(module, name, forbidden)
                blocked.append(module.__name__ + "." + name)
    with np.load(root / "q0_reference_inherited_states.npz", allow_pickle=False) as a:
        reference_pre = a["g_pre"]
    require(hashlib.sha256(reference_pre.tobytes()).hexdigest() == PRE_HASH, "Reference PRE hash")
    require(reference_pre.shape == (120, 6, 1, 17, 9, 4, 4), "Reference shape")
    base = root / "00_reference_p1.00"
    credit = root / "04_lifetime_repayment_only_p1.00"
    base_receipt = json.loads((base / "receipt.json").read_text())
    receipt = json.loads((credit / "receipt.json").read_text())
    require(base_receipt["price_factor"] == receipt["price_factor"] == 1.0, "q0 factor")
    require(base_receipt["source_binding_sha256"] == receipt["source_binding_sha256"], "Source binding")
    require(base_receipt["parameters_sha256"] == receipt["parameters_sha256"], "31-field parameter table")
    require(sha(base / "parameters.csv") == base_receipt["parameters_sha256"], "Base parameter CSV")
    require(sha(credit / "parameters.csv") == receipt["parameters_sha256"], "Credit parameter CSV")
    P = copy.deepcopy(auth["natural"])
    require(P.native_solvency_credit and not P.native_due_stayer_credit and P.unsecured_credit_limit is None,
            "Expected credit treatment")
    grid = auth["grid"]
    price = np.asarray([driver.BINDING["candidate_price"]])
    with np.load(credit / "solution_arrays.npz", allow_pickle=False) as a:
        data = {name: a[name] for name in FIELDS}
        for name in ("bp_pol_stay", "c_pol_stay"):
            if name in a.files:
                data[name] = a[name]
    require(np.array_equal(data.pop("b_grid"), grid), "Grid equality")
    require(np.array_equal(data["p_eq"], price), "Price equality")
    solution = SimpleNamespace(**data)
    P._fert2_probs = solution.fert2_probs.copy()
    shared = auth["solver"].precompute_shared(P, grid)
    policy = cal.policy_from_solution(solution, price, P, grid, shared)
    credit_pre, reconstruction = cal.reconstruct_stationary_pre_fertility(solution, policy, P, grid, shared)
    require(reconstruction["stationary_post_fertility_nesting_l1"] <= 5e-9, "Credit POST nesting")
    require(reconstruction["stationary_feasibility_projection_mass"] == 0, "Credit projection")
    require(abs(float(reference_pre.sum()) - 1) < 2e-10 and abs(float(credit_pre.sum()) - 1) < 2e-10,
            "PRE mass")
    f0 = origin_flows(cal.apply_fertility, reference_pre, policy, P)
    f1 = origin_flows(cal.apply_fertility, credit_pre, policy, P)
    rows = []
    overlap_state_max_rate_gap = 0.0
    for n in range(3):
        for m in range(n + 1):
            exposure0, flow0 = f0[n, m]
            exposure1, flow1 = f1[n, m]
            both = (exposure0 > 0) & (exposure1 > 0)
            if np.any(both):
                gap = np.abs(flow0[both] / exposure0[both] - flow1[both] / exposure1[both])
                overlap_state_max_rate_gap = max(overlap_state_max_rate_gap, float(gap.max()))
            for ten in range(exposure0.shape[1]):
                for loc in range(exposure0.shape[2]):
                    for age in range(exposure0.shape[3]):
                        for income in range(exposure0.shape[4]):
                            x0 = exposure0[:, ten, loc, age, income]
                            x1 = exposure1[:, ten, loc, age, income]
                            y0 = flow0[:, ten, loc, age, income]
                            y1 = flow1[:, ten, loc, age, income]
                            M0 = float(x0.sum()); M1 = float(x1.sum())
                            B0 = float(y0.sum()); B1 = float(y1.sum())
                            if M0 == 0 and M1 == 0:
                                continue
                            row = dict(age=18 + 4 * age, income_state=income, inherited_tenure=ten,
                                       location=loc, children_ever_born_before=n, children_at_home_before=m,
                                       reference_mass=M0, credit_mass=M1,
                                       reference_creditpolicy_births=B0, credit_creditpolicy_births=B1,
                                       reference_conditional_birth_rate=None, credit_conditional_birth_rate=None,
                                       shared_group=M0 > 0 and M1 > 0,
                                       wealth_term=0.0, group_mass_term=0.0, unmatched_term=0.0,
                                       total_change=B1 - B0)
                            if M0 > 0 and M1 > 0:
                                u0 = B0 / M0; u1 = B1 / M1
                                row["reference_conditional_birth_rate"] = u0
                                row["credit_conditional_birth_rate"] = u1
                                row["wealth_term"] = 0.5 * (M0 + M1) * (u1 - u0)
                                row["group_mass_term"] = 0.5 * (u1 + u0) * (M1 - M0)
                            else:
                                row["unmatched_term"] = B1 - B0
                            require(abs(row["wealth_term"] + row["group_mass_term"] + row["unmatched_term"] - row["total_change"]) < 3e-14,
                                    "Group add-up")
                            rows.append(row)
    require(overlap_state_max_rate_gap < 2e-10, "Fixed credit policy state birth rates changed")
    save_csv(out / "groups.csv", rows)
    age_rows = []
    for age in sorted({r["age"] for r in rows}):
        for name, sel in (("all_births", lambda r: True), ("first_births", lambda r: r["children_ever_born_before"] == 0)):
            selected = [r for r in rows if r["age"] == age and sel(r)]
            age_rows.append(dict(age=age, outcome=name,
                                 reference_creditpolicy_births=sum(r["reference_creditpolicy_births"] for r in selected),
                                 credit_creditpolicy_births=sum(r["credit_creditpolicy_births"] for r in selected),
                                 wealth_term=sum(r["wealth_term"] for r in selected),
                                 group_mass_term=sum(r["group_mass_term"] for r in selected),
                                 unmatched_term=sum(r["unmatched_term"] for r in selected),
                                 total_change=sum(r["total_change"] for r in selected),
                                 shared_groups=sum(r["shared_group"] for r in selected),
                                 unmatched_groups=sum(not r["shared_group"] for r in selected)))
    save_csv(out / "age.csv", age_rows)
    baseline_saved = json.loads((base / "closure.json").read_text())["baseline_state_impact"]
    credit_saved = json.loads((credit / "closure.json").read_text())
    targets = {"all_births": (float(credit_saved["baseline_state_impact"]["births"]),
                              float(credit_saved["cohort_summary"]["births"])),
               "first_births": (float(credit_saved["baseline_state_impact"]["first_births"]),
                                float(credit_saved["cohort_summary"]["first_births"]))}
    totals = {}
    for name, expected in targets.items():
        selected = [r for r in rows if name == "all_births" or r["children_ever_born_before"] == 0]
        observed = (sum(r["reference_creditpolicy_births"] for r in selected),
                    sum(r["credit_creditpolicy_births"] for r in selected))
        require(max(abs(a - b) for a, b in zip(observed, expected)) < 2e-10, "Native saved-flow reconciliation " + name)
        parts = {key: sum(r[key] for r in selected) for key in ("wealth_term", "group_mass_term", "unmatched_term", "total_change")}
        require(abs(sum(parts[k] for k in ("wealth_term", "group_mass_term", "unmatched_term")) - parts["total_change"]) < 2e-12,
                "Global add-up " + name)
        require(abs(parts["total_change"] - (expected[1] - expected[0])) < 2e-10,
                "Saved composition delta " + name)
        totals[name] = dict(reference_creditpolicy_births=observed[0], credit_creditpolicy_births=observed[1], **parts)
    require(abs(float(baseline_saved["births"]) - 0.1152538456218217) < 2e-12, "Reference winner q0")
    receipt_out = dict(status="passed_saved_data_only", model_solves=0, bellman_calls=0, ge_calls=0,
                       baseline_pre_sha256=PRE_HASH,
                       credit_pre_sha256=hashlib.sha256(credit_pre.tobytes()).hexdigest(),
                       source_binding_sha256=receipt["source_binding_sha256"],
                       parameters_sha256=receipt["parameters_sha256"],
                       script_sha256=sha(__file__),
                       support_status=receipt["support_status"],
                       post_nesting_l1=reconstruction["stationary_post_fertility_nesting_l1"],
                       feasibility_projection_mass=reconstruction["stationary_feasibility_projection_mass"],
                       pre_mass=[float(reference_pre.sum()), float(credit_pre.sum())],
                       birth_risk_mass=[float(reference_pre[..., :3, :].sum()), float(credit_pre[..., :3, :].sum())],
                       g="age x income x inherited tenure/housing x location x children ever born x children at home",
                       shared_groups=sum(r["shared_group"] for r in rows),
                       unmatched_groups=sum(not r["shared_group"] for r in rows),
                       overlap_state_max_rate_gap=overlap_state_max_rate_gap,
                       totals=totals, blocked_solver_entries=blocked,
                       elapsed_seconds=time.monotonic() - started,
                       max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                       interpretation="Symmetric two-order accounting under fixed credit policy; not causal wealth intervention")
    (out / "receipt.json").write_text(json.dumps(receipt_out, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": receipt_out["status"], "seconds": receipt_out["elapsed_seconds"],
                      "groups": len(rows), "totals": totals}, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.packet, args.root, args.out)
