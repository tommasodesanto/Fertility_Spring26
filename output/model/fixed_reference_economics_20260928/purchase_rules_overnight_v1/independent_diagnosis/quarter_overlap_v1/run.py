"""Selected quarter fit: one-period renter-cap shadow and matched birth margin.

No equilibrium or Bellman recursion is solved. Native renter kernels reoptimize
current saving using the saved selected continuation value at each fertile age.
"""
from __future__ import annotations

import hashlib
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PACKET = HERE.parents[1]
ROOT = PACKET.parents[3]
BOOT = PACKET / "local_runtime/bootstrap.py"
COMPLETED = PACKET / "local_runtime/runs/local10_v1/chain54/postcheck/completed.json"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def summarize(x, w):
    total = float(w.sum())
    return dict(mass=total, mean=float(np.sum(x * w) / total) if total else None,
                positive_share=float(np.sum(w * (x > 1e-8)) / total) if total else None)


def weighted_quantile(x, w, q):
    active = w > 0
    if not np.any(active):
        return None
    xv, wv = x[active], w[active]
    order = np.argsort(xv)
    return float(xv[order][np.searchsorted(np.cumsum(wv[order]), q * wv.sum())])


def group_distribution(x, w):
    total = float(w.sum())
    return {"mean": float(np.sum(x*w)/total) if total else None,
            "p10": weighted_quantile(x, w, .1),
            "median": weighted_quantile(x, w, .5),
            "p90": weighted_quantile(x, w, .9)}


def main():
    # Replay the exact read-only overlay used by the local selected postcheck.
    # Its final dispatch launches a calibration chain, so execute only setup.
    prefix = BOOT.read_text().split("if '--preflight-context' in sys.argv:", 1)[0]
    if "DIGESTS=" not in prefix or "MAPPING=" not in prefix:
        raise RuntimeError("Frozen overlay structure changed")
    exec(compile(prefix, str(BOOT), "exec"), {"__file__": str(BOOT), "__name__": "quarter_overlap_boot"})
    sys.path.insert(0, str(PACKET / "mechanism"))
    sys.path.insert(0, str(PACKET / "buyer_diagnostics"))
    import selected_runtime
    from financial_access import matched_first_birth_access

    _, contract, root, repeat = selected_runtime.authenticate_selected("quarter", COMPLETED)
    rt = selected_runtime.construct("quarter", COMPLETED, Path(tempfile.mkdtemp(prefix="quarter_overlap_auth_")) / "runtime")
    from refactor_lab.engine.parameters import get_fecundity_by_age, parent_age_maturation_active, independent_child_maturation_active
    from refactor_lab.engine.household import apply_child_aging, apply_child_aging_exempt, renter_borrowing_floor, fixed_unsecured_credit_active
    from refactor_lab.engine.kernels import full_renter_block_kernel
    from refactor_lab.engine.shared import income_at_state, income_transition_values
    from refactor_lab.engine.utils import flat_nc

    stage = repeat / "stage/solution_arrays.npz"
    with np.load(stage, allow_pickle=False) as z:
        V, g, fert = (z[k] for k in ("V", "distribution.g_pre", "fert_probs"))
        saved_h = z["hR_pol"]
        b = z["b_grid"]
        price = z["p_eq"]
        saved_shared = {k: z["shared." + k] for k in ("cb_flat", "hb_flat", "psi_flat", "gb_flat", "alpha_flat", "escale_flat", "phi_choice")}
    P = rt.P
    SD = rt.model.precompute_shared(P, b)
    assert np.array_equal(b, rt.grid)
    assert np.allclose(price, [rt.reference_price], rtol=0, atol=1e-12)
    for k, old in saved_shared.items():
        if not np.allclose(np.asarray(getattr(SD, k)), old, rtol=0, atol=1e-12):
            raise RuntimeError("Authenticated shared array differs: " + k)
    if P.n_parity != 4 or P.n_child_states != 4 or P.I != 1 or P.hR_max != 6.0:
        raise RuntimeError("Unexpected selected quarter state or rental cap")
    fec = get_fecundity_by_age(P)
    access = matched_first_birth_access(P, SD, b, price, g, fert, rule="quarter")
    bridged = ~access["owner_feasible_at_80"].any(-1) & access["owner_feasible_at_100"].any(-1)
    at80 = access["owner_feasible_at_80"].any(-1)
    never_feasible = ~access["owner_feasible_at_100"].any(-1)
    pi = fec[None, None, :, None]
    p = fert[:, 0, :, :, :, 1]
    mass = g[:, 0, :, :, :, 0, 0]
    birth = mass * pi * p
    responsive = mass * pi * p * (1-p) / float(P.kappa_fert)
    success_responsive = responsive * pi
    if not np.allclose(birth, access["weights"], rtol=0, atol=1e-15):
        raise RuntimeError(f"First-birth origin weights differ: shape {birth.shape} / {access['weights'].shape}; max {np.max(np.abs(birth-access['weights']))}")

    _, _, Pi_z = income_transition_values(P)
    cap = np.zeros_like(birth)
    loose = np.zeros_like(birth)
    wait_cap = np.zeros_like(birth)
    wait_loose = np.zeros_like(birth)
    relaxed_parent_rooms = np.zeros_like(birth)
    relaxed_wait_rooms = np.zeros_like(birth)
    capped_parent_rooms = np.zeros_like(birth)
    capped_wait_rooms = np.zeros_like(birth)
    sweep_caps = (6.0, 5.5, 5.0, 4.5)
    sweep_parent_value = {c: np.zeros_like(birth) for c in sweep_caps}
    sweep_wait_value = {c: np.zeros_like(birth) for c in sweep_caps}
    sweep_parent_rooms = {c: np.zeros_like(birth) for c in sweep_caps}
    sweep_wait_rooms = {c: np.zeros_like(birth) for c in sweep_caps}
    baseline_h_max_error = 0.0
    r = float(P.user_cost_rate * price[0])
    cb = np.ascontiguousarray(SD.cb_flat.reshape(-1))
    hb = np.ascontiguousarray(SD.hb_flat.reshape(-1))
    psi = np.ascontiguousarray(SD.psi_flat.reshape(-1))
    gb = np.ascontiguousarray(SD.gb_flat.reshape(-1))
    alpha_v = np.ascontiguousarray(SD.alpha_flat.reshape(-1))
    esc = np.ascontiguousarray(SD.escale_flat.reshape(-1))
    n_c = SD.nc
    age_cells = []
    for j in range(P.J):
        if not (P.A_f_start <= j+1 <= P.A_f_end) or j == P.J-1:
            continue
        age_cells.append(float(P.age_start + P.da*j))
        for zz, z in enumerate(P.z_grid):
            nxt = np.zeros((len(b), 1+P.n_house, P.I, P.n_parity, P.n_child_states))
            for znext, prob in enumerate(Pi_z[zz]):
                nxt += float(prob) * V[:, :, :, j+1, znext]
            if bool(getattr(P, "use_age_survival", False)):
                # These selected fertile ages are below any terminal age.
                # Reconstruct the exact Bellman death branch when needed.
                if float(P.survival_probs[j]) < 1.0:
                    raise RuntimeError("Fertile-age survival needs explicit bequest reconstruction")
            standard = apply_child_aging(nxt, P, len(b), 1+P.n_house, P.I,
                                         P.n_parity, P.n_child_states, age_index=j)
            exempt = (apply_child_aging_exempt(nxt, P, len(b), 1+P.n_house, P.I,
                                              P.n_parity, P.n_child_states, age_index=j)
                      if parent_age_maturation_active(P) and independent_child_maturation_active(P)
                      else standard)
            y = float(income_at_state(P, 0, j, float(z)))
            rv = np.ascontiguousarray(P.R_gross*b + y)
            rvt = np.ascontiguousarray(P.R_gross*np.maximum(b, 0) + y)
            floor = np.maximum(renter_borrowing_floor(P, b, j), b[0])
            # The native current-period renter kernel is called with the same
            # parameters and continuation as the selected Bellman equation.
            vals = {}
            for label, continuation in (("wait", standard), ("parent", exempt)):
                vcr = np.ascontiguousarray(flat_nc(continuation[:, 0, 0], len(b), n_c))
                for max_rooms in (*sweep_caps, 100.0):
                    value, _, _, rooms = full_renter_block_kernel(
                        rv, rvt, vcr, np.zeros((len(b), n_c)), 0, b,
                        cb, hb, psi, gb, alpha_v, esc, r, max_rooms,
                        P.c_min, P.c_bar_0, P.h_bar_0, P.alpha_cons,
                        1.0-P.sigma, P.beta, float(P.debt_taper_weights[j+1]),
                        float(P.debt_caps[j+1]), (3-5**.5)/2, (5**.5-1)/2, 1e-3,
                        0, np.zeros(n_c), 0, 0, 0., 0., 6.,
                        bool(getattr(P, "native_exact_allocation_output", False)),
                        None, float(floor[0]) if fixed_unsecured_credit_active(P) else -np.inf,
                    )
                    vals[(label, max_rooms)] = (value, rooms)
            # Family cells use Fortran flattening: nn + npar*cs.
            childless_col, parent_col = 0, 1+P.n_parity
            cap[:, 0, j, zz] = vals[("parent", 6.)][0][:, parent_col]
            loose[:, 0, j, zz] = vals[("parent", 100.)][0][:, parent_col]
            relaxed_parent_rooms[:, 0, j, zz] = vals[("parent", 100.)][1][:, parent_col]
            capped_parent_rooms[:, 0, j, zz] = vals[("parent", 6.)][1][:, parent_col]
            wait_cap[:, 0, j, zz] = vals[("wait", 6.)][0][:, childless_col]
            wait_loose[:, 0, j, zz] = vals[("wait", 100.)][0][:, childless_col]
            relaxed_wait_rooms[:, 0, j, zz] = vals[("wait", 100.)][1][:, childless_col]
            capped_wait_rooms[:, 0, j, zz] = vals[("wait", 6.)][1][:, childless_col]
            for sweep_cap in sweep_caps:
                sweep_parent_value[sweep_cap][:, 0, j, zz] = vals[("parent", sweep_cap)][0][:, parent_col]
                sweep_wait_value[sweep_cap][:, 0, j, zz] = vals[("wait", sweep_cap)][0][:, childless_col]
                sweep_parent_rooms[sweep_cap][:, 0, j, zz] = vals[("parent", sweep_cap)][1][:, parent_col]
                sweep_wait_rooms[sweep_cap][:, 0, j, zz] = vals[("wait", sweep_cap)][1][:, childless_col]
            h_saved = saved_h[:, 0, 0, j, zz]
            h_current = vals[("wait", 6.)][1][:, childless_col]
            h_parent_saved = saved_h[:, 0, 0, j, zz, 1, 1]
            h_parent_current = vals[("parent", 6.)][1][:, parent_col]
            if parent_age_maturation_active(P) and independent_child_maturation_active(P):
                # Saved parent policy uses standard continuation outside the
                # successful-birth exempt branch, so only wait matches here.
                pass
            baseline_h_max_error = max(baseline_h_max_error,
                                       float(np.max(np.abs(h_saved[:, 0, 0]-h_current))))
    parent_shadow = loose-cap
    wait_shadow = wait_loose-wait_cap
    differential = parent_shadow-wait_shadow
    if not np.array_equal(sweep_parent_value[6.0], cap) or not np.array_equal(sweep_wait_value[6.0], wait_cap):
        raise RuntimeError("Cap-six value replay changed")
    if not np.array_equal(sweep_parent_rooms[6.0], capped_parent_rooms) or not np.array_equal(sweep_wait_rooms[6.0], capped_wait_rooms):
        raise RuntimeError("Cap-six room replay changed")
    occupied = responsive > 0
    if np.min(parent_shadow[occupied]) < -1e-6 or np.min(wait_shadow[occupied]) < -1e-6:
        raise RuntimeError("Relaxing cap lowered a renter branch value")
    if baseline_h_max_error > 0.01:
        raise RuntimeError(f"Native baseline renter policy reconstruction differs: {baseline_h_max_error}")
    relaxed_parent_max = float(np.max(relaxed_parent_rooms[occupied]))
    relaxed_wait_max = float(np.max(relaxed_wait_rooms[occupied]))
    if max(relaxed_parent_max, relaxed_wait_max) >= 99.99:
        raise RuntimeError("Relaxed rental cap 100 still binds at a responsive state")
    sections = {}
    ages = np.broadcast_to(np.asarray(P.age_start + P.da*np.arange(P.J))[None,None,:,None], birth.shape)
    wealth = np.broadcast_to(b[:,None,None,None], birth.shape)
    income = np.broadcast_to(np.asarray([[income_at_state(P, 0, j, float(z))
                                        for z in P.z_grid] for j in range(P.J)])[None,None,:,:], birth.shape)
    owner_sizes = np.asarray(P.H_own, dtype=float)
    newly_feasible_products = access["owner_feasible_at_100"] & ~access["owner_feasible_at_80"]
    def detailed(mask):
        w_birth = birth*mask
        w_response = responsive*mask
        fertile = (ages >= min(age_cells)) & (ages <= max(age_cells))
        w_mass = mass*mask*fertile
        total_birth = float(w_birth.sum())
        total_response = float(w_response.sum())
        total_mass = float(w_mass.sum())
        age_rows = []
        for age in age_cells:
            age_mask = ages == age
            age_rows.append({"age": age,
                             "mass_share": float(np.sum(w_mass*age_mask)/total_mass) if total_mass else None,
                             "birth_flow_share": float(np.sum(w_birth*age_mask)/total_birth) if total_birth else None,
                             "responsiveness_share": float(np.sum(w_response*age_mask)/total_response) if total_response else None})
        room_rows = {}
        for name, rooms in (("parent_cap6", capped_parent_rooms), ("parent_relaxed100", relaxed_parent_rooms),
                            ("wait_cap6", capped_wait_rooms), ("wait_relaxed100", relaxed_wait_rooms)):
            room_rows[name] = {"response_weighted_mean": float(np.sum(w_response*rooms)/total_response) if total_response else None,
                               "birth_weighted_mean": float(np.sum(w_birth*rooms)/total_birth) if total_birth else None,
                               "at_six_grid_cell_count": int(np.count_nonzero((w_response>0)&(rooms>=6-1e-5))),
                               "at_six_responsiveness_share": float(np.sum(w_response*(rooms>=6-1e-5))/total_response) if total_response else None}
        product_rows = []
        for k, size in enumerate(owner_sizes):
            feasible = newly_feasible_products[...,k] & mask
            product_rows.append({"owner_product_size": float(size),
                                 "positive_birth_cell_count": int(np.count_nonzero(feasible & (birth>0))),
                                 "share_of_group_birth_flow_with_new_access_to_product": float(np.sum(w_birth*feasible)/total_birth) if total_birth else None})
        return {"state_cell_count_positive_response": int(np.count_nonzero(w_response>0)),
                "pre_fertility_renter_mass": total_mass,
                "share_of_all_fertile_childless_renter_mass": total_mass/float(np.sum(mass*fertile)),
                "age_distribution_within_group": age_rows,
                "income": group_distribution(income, w_response),
                "wealth": group_distribution(wealth, w_response),
                "share_negative_wealth_response_weighted": float(np.sum(w_response*(wealth<0))/total_response) if total_response else None,
                "share_zero_or_negative_wealth_response_weighted": float(np.sum(w_response*(wealth<=0))/total_response) if total_response else None,
                "renter_rooms": room_rows,
                "max_parent_cap_shadow_positive_response": float(np.max(parent_shadow[w_response>0])) if total_response else None,
                "max_parent_minus_wait_shadow_positive_response": float(np.max(differential[w_response>0])) if total_response else None,
                "newly_feasible_owner_products": product_rows}
    for name, mask in (("all", np.ones_like(bridged, bool)), ("bridged", bridged),
                       ("feasible_80", at80), ("never_feasible_100", never_feasible),
                       ("ineligible_80", ~at80)):
        sections[name] = dict(
            birth_flow=float(np.sum(birth*mask)), responsiveness_mass=float(np.sum(responsive*mask)),
            birth_share=float(np.sum(birth*mask)/birth.sum()),
            responsiveness_share=float(np.sum(responsive*mask)/responsive.sum()),
            successful_birth_value_responsiveness_share=float(np.sum(success_responsive*mask)/success_responsive.sum()),
            parent_cap_shadow=summarize(parent_shadow, responsive*mask),
            wait_cap_shadow=summarize(wait_shadow, responsive*mask),
            parent_minus_wait_cap_shadow=summarize(differential, responsive*mask),
            parent_minus_wait_cap_shadow_success_weighted=summarize(differential, success_responsive*mask),
            parent_minus_wait_above_tenth_kappa_share=float(np.sum(responsive*mask*(differential>0.1*P.kappa_fert))/np.sum(responsive*mask)) if np.sum(responsive*mask)>0 else None,
            decomposition=detailed(mask),
        )
    if not np.allclose(sections["bridged"]["birth_flow"]+sections["feasible_80"]["birth_flow"], sections["all"]["birth_flow"], rtol=0, atol=1e-12):
        raise RuntimeError("Financial-access birth groups do not exhaust renter-origin births")
    if not np.allclose(sections["bridged"]["responsiveness_mass"]+sections["feasible_80"]["responsiveness_mass"], sections["all"]["responsiveness_mass"], rtol=0, atol=1e-12):
        raise RuntimeError("Financial-access responsiveness groups do not exhaust renter margin")
    if abs(sum(sections[name]["decomposition"]["share_of_all_fertile_childless_renter_mass"] for name in ("bridged", "feasible_80", "never_feasible_100"))-1) > 1e-10:
        raise RuntimeError("Financial-access groups do not exhaust fertile renter mass")
    for name in ("all", "bridged", "feasible_80"):
        for key in ("mass_share", "birth_flow_share", "responsiveness_share"):
            if abs(sum(row[key] for row in sections[name]["decomposition"]["age_distribution_within_group"])-1) > 1e-10:
                raise RuntimeError(f"Age decomposition fails for {name}/{key}")
    cap_sweep = {}
    for name, mask in (("bridged", bridged), ("feasible_80", at80)):
        w = success_responsive*mask
        total = float(w.sum())
        cap_sweep[name] = []
        for sweep_cap in sweep_caps:
            ps = loose-sweep_parent_value[sweep_cap]
            ws = wait_loose-sweep_wait_value[sweep_cap]
            diff = ps-ws
            if np.min(ps[w>0]) < -1e-6 or np.min(ws[w>0]) < -1e-6:
                raise RuntimeError(f"Lower cap improved renter branch value at {name}/{sweep_cap}")
            cap_sweep[name].append({"cap_rooms": sweep_cap,
                                    "parent_minus_wait_cap_shadow_success_weighted": float(np.sum(w*diff)/total),
                                    "parent_rooms_success_weighted": float(np.sum(w*sweep_parent_rooms[sweep_cap])/total),
                                    "wait_rooms_success_weighted": float(np.sum(w*sweep_wait_rooms[sweep_cap])/total),
                                    "parent_rooms_slack100_success_weighted": float(np.sum(w*relaxed_parent_rooms)/total),
                                    "wait_rooms_slack100_success_weighted": float(np.sum(w*relaxed_wait_rooms)/total)})
        if cap_sweep[name][0]["parent_minus_wait_cap_shadow_success_weighted"] != sections[name]["parent_minus_wait_cap_shadow_success_weighted"]["mean"]:
            raise RuntimeError("Cap-six success-weighted shadow changed")
    result = dict(status="computed", definition="One-period rental cap 6 to 100, reoptimizing saving at baseline prices and saved next-age value; owner menu and continuation fixed", rule="quarter", cap_rooms=6.0,
                  relaxed_rooms=100.0, age_cells=age_cells, sections=sections,
                  fixed_continuation_cap_sweep=cap_sweep,
                  exact_80_to_100_financial_access=access["share_ineligible_at_80_eligible_at_100"],
                  relaxed_max_rooms_positive_responsiveness=dict(parent=relaxed_parent_max, wait=relaxed_wait_max),
                  fixed_unsecured_credit_active=bool(fixed_unsecured_credit_active(P)),
                  baseline_renter_policy_max_room_error=baseline_h_max_error,
                  source=dict(selected_completed=str(COMPLETED), completed_sha256=sha(COMPLETED),
                              stage=str(stage), stage_sha256=sha(stage), bootstrap=str(BOOT), bootstrap_sha256=sha(BOOT),
                              root_parameters_sha256=sha(root/"parameters.csv"),
                              contract_sha256=sha(COMPLETED.parent/"input_contract.json")),
                  caution="Differential renter-only value shadow is not the change in the full fertility attempt/wait values, nor a causal policy response.")
    (HERE/"result.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"status":result["status"],"sections":sections,"baseline_h_max_error":baseline_h_max_error}, indent=2))


if __name__ == "__main__":
    main()
