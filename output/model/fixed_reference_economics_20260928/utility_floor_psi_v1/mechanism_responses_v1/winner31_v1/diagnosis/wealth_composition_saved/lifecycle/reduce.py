"""Saved-policy never-parent lifecycle readout; zero model solves."""
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


def require(ok, text):
    if not ok:
        raise ValueError(text)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def forbidden(*args, **kwargs):
    raise RuntimeError("Model solve forbidden in saved-policy readout")


def load_case(folder, P, auth, driver, cal, price):
    receipt = json.loads((folder / "receipt.json").read_text())
    require(receipt["price_factor"] == 1.0, "q0 factor")
    require(driver.sha(folder / "parameters.csv") == receipt["parameters_sha256"], "Parameter table")
    with np.load(folder / "solution_arrays.npz", allow_pickle=False) as arrays:
        data = {name: arrays[name] for name in FIELDS}
        for name in ("bp_pol_stay", "c_pol_stay"):
            if name in arrays.files:
                data[name] = arrays[name]
    require(np.array_equal(data.pop("b_grid"), auth["grid"]), "Grid")
    require(np.array_equal(data["p_eq"], price), "Price")
    sol = SimpleNamespace(**data)
    P._fert2_probs = sol.fert2_probs.copy()
    shared = auth["solver"].precompute_shared(P, auth["grid"])
    policy = cal.policy_from_solution(sol, price, P, auth["grid"], shared)
    pre, reconstruction = cal.reconstruct_stationary_pre_fertility(sol, policy, P, auth["grid"], shared)
    require(reconstruction["stationary_post_fertility_nesting_l1"] <= 5e-9, "POST nesting")
    require(reconstruction["stationary_feasibility_projection_mass"] == 0, "Projection")
    return policy, shared, pre, receipt, reconstruction


def cells(cal, grid, price, P, shared, policy, source_pre, scenario, inherited_tenure):
    mask = np.zeros_like(source_pre)
    tenure_slice = 0 if inherited_tenure == "renter" else slice(1, None)
    mask[:, tenure_slice, :, :, :, 0, 0] = source_pre[:, tenure_slice, :, :, :, 0, 0]
    counter = cal.SolveCounter()
    ev = cal.evaluate_period(price, mask, P, grid, shared, counter, supplied_policy=policy)
    require(counter.total == 0, "Saved policy unexpectedly solved")
    require(ev.feasibility_projection_mass == 0, "Masked distribution projected")
    require(abs(float(ev.g_current.sum()) - float(mask.sum())) < 2e-10, "Current mass")
    require(abs(float(ev.g_post_fertility.sum()) - float(mask.sum())) < 2e-10, "Post-fertility mass")
    out = []
    for j in range(7):
        pre_j = mask[:, :, :, j, :, :, :]
        current_j = ev.g_current[:, :, :, j, :, :, :]
        mass = float(pre_j.sum())
        if mass == 0:
            continue
        require(abs(float(current_j.sum()) - mass) < 2e-10, "Age mass")
        b = grid.reshape(-1, 1, 1, 1, 1, 1)
        next_b = policy.bp_pol[:, :, :, j, :, :, :]
        cons = policy.c_pol[:, :, :, j, :, :, :]
        renter = current_j[:, 0, :, :, :, :]
        owner = current_j[:, 1:, :, :, :, :]
        renter_next = next_b[:, 0, :, :, :, :]
        owner_next = next_b[:, 1:, :, :, :, :]
        birth_flow = float(ev.g_post_fertility[:, :, :, j, :, 1:, :].sum())
        pre_financial = float((pre_j * b).sum()) / mass
        current_financial = float((current_j * b).sum()) / mass
        next_financial = float((current_j * next_b).sum()) / mass
        current_renter_mass = float(renter.sum())
        current_owner_mass = float(owner.sum())
        row = dict(scenario=scenario, inherited_tenure=inherited_tenure, age=18 + 4*j,
                   origin_never_parent_mass=mass, first_birth_flow=birth_flow,
                   first_birth_rate=birth_flow/mass,
                   mean_inherited_financial_assets=pre_financial,
                   inherited_negative_financial_asset_share=float(pre_j[grid < 0].sum())/mass,
                   mean_current_post_transaction_financial_assets=current_financial,
                   mean_chosen_next_financial_assets=next_financial,
                   mean_next_minus_post_transaction_financial_assets=next_financial-current_financial,
                   mean_nonhousing_consumption=float((current_j * cons).sum())/mass,
                   destination_renter_mass=current_renter_mass,
                   destination_owner_mass=current_owner_mass,
                   destination_owner_rate=current_owner_mass/mass,
                   destination_renter_rate=current_renter_mass/mass,
                   renter_next_negative_asset_mass=float(renter[renter_next < 0].sum()),
                   owner_next_negative_asset_mass=float(owner[owner_next < 0].sum()),
                   renter_next_negative_asset_rate=None if current_renter_mass == 0 else float(renter[renter_next < 0].sum())/current_renter_mass,
                   owner_next_negative_asset_rate=None if current_owner_mass == 0 else float(owner[owner_next < 0].sum())/current_owner_mass,
                   renter_next_at_zero_mass=float(renter[np.abs(renter_next) <= 1e-10].sum()),
                   owner_next_at_grid_min_mass=float(owner[np.abs(owner_next-grid[0]) <= 1e-10].sum()))
        out.append(row)
    require(abs(sum(row["first_birth_flow"] for row in out)-float(ev.births)) < 2e-10,
            "Never-parent birth flow")
    return out


def main(packet, root, out):
    started = time.monotonic()
    require(not out.exists(), "Refuse existing output")
    out.mkdir(parents=True)
    sys.path.insert(0, str(packet))
    import fixed_price_responses as driver
    auth = driver.authenticate_candidate(out / "runtime_auth")
    cal = auth["context"]["prepared"].rt["primitive"].pf.calendar
    blocked = []
    for module in (auth["solver"], auth["context"]["prepared"].rt["model"], cal):
        for name in dir(module):
            if name.startswith("solve") and callable(getattr(module, name)):
                setattr(module, name, forbidden)
                blocked.append(module.__name__ + "." + name)
    with np.load(root / "q0_reference_inherited_states.npz", allow_pickle=False) as arrays:
        baseline_pre = arrays["g_pre"]
    require(hashlib.sha256(baseline_pre.tobytes()).hexdigest() == PRE_HASH, "Baseline PRE")
    price = np.asarray([driver.BINDING["candidate_price"]])
    p0 = copy.deepcopy(auth["P"])
    pc = copy.deepcopy(auth["natural"])
    base = root / "00_reference_p1.00"
    credit = root / "04_lifetime_repayment_only_p1.00"
    pol0, sd0, replay_pre, rec0, recon0 = load_case(base, p0, auth, driver, cal, price)
    polc, sdc, credit_pre, recc, reconc = load_case(credit, pc, auth, driver, cal, price)
    require(np.array_equal(replay_pre, baseline_pre), "Baseline PRE exact replay")
    require(rec0["source_binding_sha256"] == recc["source_binding_sha256"], "Source binding")
    require(rec0["parameters_sha256"] == recc["parameters_sha256"], "31 parameters")
    require(float(baseline_pre.sum()) > 0.999999999 and float(credit_pre.sum()) > 0.999999999,
            "PRE mass")
    rows = []
    for scenario, P, sd, policy, pre in (
        ("reference_policy_reference_PRE", p0, sd0, pol0, baseline_pre),
        ("credit_policy_reference_PRE", pc, sdc, polc, baseline_pre),
        ("credit_policy_credit_PRE", pc, sdc, polc, credit_pre)):
        for tenure in ("renter", "owner"):
            rows.extend(cells(cal, auth["grid"], price, P, sd, policy, pre, scenario, tenure))
    expected = {
        "reference_policy_reference_PRE": float(json.loads((base / "closure.json").read_text())["cohort_summary"]["first_births"]),
        "credit_policy_reference_PRE": float(json.loads((credit / "closure.json").read_text())["baseline_state_impact"]["first_births"]),
        "credit_policy_credit_PRE": float(json.loads((credit / "closure.json").read_text())["cohort_summary"]["first_births"]),
    }
    by_scenario = {name:sum(r["first_birth_flow"] for r in rows if r["scenario"] == name)
                   for name in expected}
    for key in expected:
        require(abs(by_scenario[key] - expected[key]) < 2e-10, "Saved first-birth reconciliation " + key)
    with (out / "lifecycle.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    receipt = dict(status="passed_saved_policy_lifecycle", model_solves=0,
                   baseline_pre_sha256=PRE_HASH,
                   credit_pre_sha256=hashlib.sha256(credit_pre.tobytes()).hexdigest(),
                   source_binding_sha256=recc["source_binding_sha256"],
                   parameters_sha256=recc["parameters_sha256"],
                   source_sha256=sha(__file__),
                   baseline_post_nesting_l1=recon0["stationary_post_fertility_nesting_l1"],
                   credit_post_nesting_l1=reconc["stationary_post_fertility_nesting_l1"],
                   projection_mass=[recon0["stationary_feasibility_projection_mass"],reconc["stationary_feasibility_projection_mass"]],
                   pre_mass=[float(baseline_pre.sum()),float(credit_pre.sum())],
                   first_births=by_scenario, expected_first_births=expected,
                   age_range="18-42, four-year cells", rows=len(rows),
                   negative_renter_next_assets_interpretation="unsecured net financial borrowing for destination renters; owner negative assets combine secured and other debt",
                   binding_natural_credit_floor_mass="unavailable: saved solution does not retain state-specific lifetime-solvency floor",
                   solver_traps=blocked, elapsed_seconds=time.monotonic()-started,
                   max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status":receipt["status"],"elapsed":receipt["elapsed_seconds"],"rows":len(rows),"first_births":by_scenario}))


if __name__ == "__main__":
    p=argparse.ArgumentParser(); p.add_argument("--packet",type=Path,required=True); p.add_argument("--root",type=Path,required=True); p.add_argument("--out",type=Path,required=True)
    a=p.parse_args(); main(a.packet,a.root,a.out)
