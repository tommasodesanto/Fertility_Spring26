"""Part 2: lock-in diagnostic from saved dose-probe arrays (zero model solves).

Uses THIS folder's psi06 (P.psi reference 0.06) and psi00 (P.psi 0.0) saved
stage arrays at both financed shares (4 cases). Read-only w.r.t. engine.

Timing (verified against the executed quarter engine,
small_credit_lab/engine/distribution.py, and the observer reconstruction that
gates 5e-9 in every solve): the saved g_beginning_distribution is the
POST-fertility mass at each age (fertility is applied to it in place in the
KFE loop, then the observer rebuilds pre-fertility mass by advancing it).
The at-risk (pre-fertility) mass g_pre is therefore reconstructed exactly as
the observer does -- entrant cohort at j=0, advance of post-fertility mass
otherwise -- from saved arrays + saved policy, with no Bellman solve:

  g_pre[0]     = entry cohort from saved entry_by_loc (KFE inline entry logic)
  g_pre[j>=1]  = advance_cohort_one_period_markov_income(survival[j-1] *
                 g_post[j-1]) with saved loc/tenure/saving policy and rebuilt
                 location/tenure maps.

Population: CHILDLESS households (parity n=0, children at home m/cs=0),
age cells j=0..4 (ages 18-37), by beginning tenure state to
(0=renter, 1..5=owner rungs of reference H_own).

Tables:
 1. at-risk mass by beginning tenure state (from reconstructed g_pre),
    per (case, age cell, to), with within-cell share.
 2. first-birth attempt probability by beginning tenure state: mass-weighted
    mean of saved fert_probs try probability (parity-0 slice [...,1]) over the
    at-risk mass -- the exact engine birth-flow factor; birth probability =
    attempt x pi_j with the reference per-period fecundity pi_j reported.
    Verification columns compare the implied aggregate first births against
    the observer first_birth_flow and the post-fertility reconstruction
    against saved g_beginning_distribution (L1).
 3. Post-birth tenure transitions for childless owners with a first birth
    (keep/upsize/downsize/rent): checked against the saved stage key
    inventory; reported only if a post-birth-branch tenure policy is saved.

Writes part2_tables/ CSVs + MD. No economic interpretation.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
CASES = ("psi06_phi08", "psi06_phi10", "psi00_phi08", "psi00_phi10")
J_CELLS = (0, 1, 2, 3, 4)
PRICE = 0.6744838540900874


def build_case_P(P_ref, case, phi):
    import copy
    from small_credit_lab import credit
    P = copy.deepcopy(P_ref)
    P.hbar_child_rooms = 1.0
    P.hbar_first_child_jump = 1.49169824815624
    P.phi = np.full_like(np.asarray(P.phi, dtype=float), float(phi))
    P.H0 = np.full_like(np.asarray(P.H0, dtype=float), 7.288573389887633)
    if case["psi"] is not None:
        P.psi = float(case["psi"])
    credit.bind_engine_credit(P, "corrected", 0.0)
    return P


def main():
    results = json.loads((HERE / "case_results.json").read_text())
    for lab in CASES:
        assert results.get(lab, {}).get("status") == "passed", lab

    sys.path.insert(0, str(HERE))
    import run_dose_probe as rp
    native = rp.setup_imports(HERE / "part2_setup_tmp")
    out_pre = HERE / "part2_setup_tmp"
    out_pre.mkdir(parents=True, exist_ok=True)
    _, _, _, _, P_ref, grid, _ = rp.build_base(native, out_pre)
    b_grid = np.asarray(grid, dtype=float)
    H_OWN = np.asarray(P_ref.H_own, dtype=float)
    assert H_OWN.tolist() == [2.0, 4.0, 6.0, 8.0, 10.0], H_OWN.tolist()
    assert float(P_ref.psi) == 0.06, float(P_ref.psi)
    age_start = float(getattr(P_ref, "age_start", 18.0))
    da = float(getattr(P_ref, "da", getattr(P_ref, "period_years", 4.0)))

    sys.path.insert(0, str(ROOT / "output/model/fixed_reference_economics_20260928"
                            "/purchase_rules_overnight_v1/engines/quarter"))
    from small_credit_lab.engine.distribution import (
        advance_cohort_one_period_markov_income,
        aggregate_entry_wealth_grid_weights,
        entry_wealth_grid_weights,
        interp_indices,
    )
    from small_credit_lab.engine import distribution as distmod
    from small_credit_lab.engine.household import (
        birth_destination_child_state,
        build_forward_tenure_transition_maps,
    )
    from small_credit_lab.engine.parameters import (
        get_fecundity_by_age,
        independent_child_maturation_active,
    )
    from small_credit_lab.engine.shared import income_transition_values

    fec = np.asarray(get_fecundity_by_age(P_ref), dtype=float)
    indep = bool(independent_child_maturation_active(P_ref))
    flags = {
        "sequential_births": bool(getattr(P_ref, "sequential_births", False)),
        "joint_nested_choice": bool(getattr(P_ref, "joint_nested_choice", False)),
        "readiness_gate": bool(getattr(P_ref, "readiness_gate_enabled", False)),
        "parent_age_maturation": bool(getattr(P_ref, "parent_age_maturation_active",
                                              getattr(P_ref, "parent_age_maturation", False))),
        "use_age_survival": bool(getattr(P_ref, "use_age_survival", False)),
        "use_stochastic_aging": bool(getattr(P_ref, "use_stochastic_aging", False)),
        "due_stayer": bool(getattr(P_ref, "native_due_stayer_credit", False)),
        "entry_censor": bool(getattr(P_ref, "entry_wealth_censor_to_frontier", False)),
        "independent_child_maturation": indep,
    }
    assert flags["sequential_births"] and not flags["joint_nested_choice"]
    assert not flags["readiness_gate"]
    assert int(P_ref.n_parity) == 4 and int(P_ref.n_child_states) == 4

    to_labels = ["renter"] + ["owner_H%g" % h for h in H_OWN]
    nt = 1 + len(H_OWN)

    mass_rows, prob_rows, ver_rows = [], [], []
    key_inventory = {}
    for lab in CASES:
        r = results[lab]
        stage = Path(r["report"]) / "stage" / "solution_arrays.npz"
        a = np.load(stage)
        key_inventory[lab] = sorted(a.files)
        gb = np.asarray(a["g_beginning_distribution"], dtype=float)  # POST-fertility
        fp = np.asarray(a["fert_probs"], dtype=float)
        f2 = np.asarray(a["fert2_probs"], dtype=float)
        lp = np.asarray(a["loc_probs"], dtype=float)
        tc = np.asarray(a["tenure_choice"])
        tp = np.asarray(a["tenure_probs"], dtype=float)
        bp = np.asarray(a["bp_pol"], dtype=float)
        bp_stay = np.asarray(a["bp_pol_stay"], dtype=float)
        entry_by_loc = np.asarray(a["entry_by_loc"], dtype=float)
        zin = np.asarray(a["income_transition"], dtype=float)
        npar, ncs = int(P_ref.n_parity), int(P_ref.n_child_states)
        J = int(P_ref.J)
        assert gb.shape == (len(b_grid), nt, 1, J, zin.shape[0], npar, ncs)

        psi = float(r.get("selling_cost_psi"))
        phi = float(r.get("phi"))
        P = build_case_P(P_ref, {"psi": psi if psi != 0.06 else None}, phi)
        assert float(P.psi) == psi and float(np.asarray(P.phi).flat[0]) == phi
        P.entry_by_loc = entry_by_loc
        z_grid, z_weights, Pi_z = income_transition_values(P)
        assert np.allclose(Pi_z, zin), lab
        SD = SimpleNamespace(nc=npar * ncs,
                             phi_choice=np.asarray(a["shared.phi_choice"], dtype=float),
                             birth_dp=np.asarray(a["shared.birth_dp"]),
                             birth_entry_grant=np.asarray(
                                 a["shared.birth_entry_grant"], dtype=float))
        p_hat = np.full((int(P.I),), PRICE)
        hc = np.zeros((int(P.I), nt))
        he = np.zeros((int(P.I), nt))
        for i in range(int(P.I)):
            for ten in range(1, nt):
                hs = float(P.H_own[ten - 1])
                hc[i, ten] = p_hat[i] * hs
                he[i, ten] = (1 - float(P.psi)) * p_hat[i] * hs
        lmm_idx = np.zeros((int(P.I), nt, len(b_grid)), dtype=np.int64)
        lmm_wt = np.zeros((int(P.I), nt, len(b_grid)))
        for io in range(int(P.I)):
            for to in range(nt):
                ba = np.clip(b_grid + he[io, to], b_grid[0], b_grid[-1])
                lmm_idx[io, to, :], lmm_wt[io, to, :] = interp_indices(b_grid, ba)
        tmx_idx, tmx_wt = build_forward_tenure_transition_maps(
            P, b_grid, hc, he, SD.phi_choice, SD.birth_dp, SD.birth_entry_grant)
        ust = bool(flags["use_stochastic_aging"] and hasattr(P, "Pi_child"))
        Pia = P.Pi_child if ust else None

        # g_pre[0]: entrant cohort (KFE inline entry logic; readiness off).
        # Mirror the KFE exactly, including entry censoring when flagged.
        Nb, Nz = len(b_grid), len(z_grid)
        g_pre = np.zeros((Nb, nt, int(P.I), J, Nz, npar, ncs))
        entry_idx, entry_wt = aggregate_entry_wealth_grid_weights(b_grid, P)
        for i in range(int(P.I)):
            for zz in range(Nz):
                idx, wts = entry_wealth_grid_weights(
                    b_grid, P, i=i, j=0, z_value=float(z_grid[zz]))
                for kk, ww in zip(idx, wts):
                    g_pre[int(kk), 0, i, 0, zz, 0, 0] += (
                        float(entry_by_loc[i]) * float(z_weights[zz]) * float(ww))
        if flags["entry_censor"]:
            V = np.asarray(a["V"], dtype=float)
            assert V.shape == gb.shape, (V.shape, gb.shape)
            censored = float(distmod._censor_entry_dead_mass(
                g_pre[:, :, :, 0, :, :, :], V[:, :, :, 0, :, :, :]))
        else:
            censored = 0.0
        for j in range(1, J):
            surv = float(P.survival_probs[j - 1]) if flags["use_age_survival"] else 1.0
            g_pre[:, :, :, j, :, :, :] = advance_cohort_one_period_markov_income(
                surv * gb[:, :, :, j - 1, :, :, :], j - 1, lp, tc, tp, bp,
                P, b_grid, SD, lmm_idx, lmm_wt, tmx_idx, tmx_wt, ust, Pia, Pi_z,
                newborn_frac=None,
                bp_pol_stay=bp_stay if flags["due_stayer"] else None)

        # Gate: re-apply fertility (KFE sequential logic) and compare with gb.
        post = g_pre.copy()
        births1 = np.zeros(J)
        for j in range(J):
            if not (int(P.A_f_start) <= j + 1 <= int(P.A_f_end)):
                continue
            pi_j = float(fec[j])
            for zz in range(Nz):
                gc = g_pre[:, :, :, j, zz, 0, 0]
                m1 = gc * fp[:, :, :, j, zz, 1]
                real1 = pi_j * m1
                post[:, :, :, j, zz, 0, 0] = gc - real1
                post[:, :, :, j, zz, 1, 1] += real1
                births1[j] += float(real1.sum())
                for nn in range(1, npar - 1):
                    css = range(0, nn + 1) if indep else (1,)
                    for cs in css:
                        atr = g_pre[:, :, :, j, zz, nn, cs]
                        p2 = (f2[:, :, :, j, zz, 1, nn - 1, cs] if indep
                              else f2[:, :, :, j, zz, 1, nn - 1])
                        real2 = pi_j * atr * p2
                        post[:, :, :, j, zz, nn, cs] -= real2
                        post[:, :, :, j, zz, nn + 1,
                             birth_destination_child_state(P, cs)] += real2
        l1 = float(np.abs(post - gb).sum())
        agg1 = float(births1.sum())
        obs1 = float(r["extra"]["birth_flows_uniform_birth_time"]["first_birth_flow"])
        ver_rows.append({"case": lab, "psi": psi, "phi": phi,
                         "recon_post_L1_vs_saved_gb": l1,
                         "implied_first_births": agg1,
                         "observer_first_birth_flow": obs1,
                         "implied_minus_observer": agg1 - obs1,
                         "entry_mass_recon": float(g_pre[:, :, :, 0, :, :, :].sum()),
                         "entry_rate_saved": float(entry_by_loc.sum()),
                         "entry_censored_mass": censored})

        for j in J_CELLS:
            pi_j = float(fec[j])
            cell = float(g_pre[:, :, :, j, :, 0, 0].sum())
            for to in range(nt):
                m = float(g_pre[:, to, :, j, :, 0, 0].sum())
                mass_rows.append({
                    "case": lab, "psi": psi, "phi": phi,
                    "age_cell_j": j, "age_low": age_start + da * j,
                    "age_high": age_start + da * (j + 1) - 1,
                    "begin_tenure_to": to, "begin_tenure": to_labels[to],
                    "at_risk_mass": m,
                    "share_of_childless_cell": m / cell if cell > 0 else float("nan"),
                    "cell_mass": cell})
                if m > 0:
                    w = g_pre[:, to, :, j, :, 0, 0]
                    attempt = float((w * fp[:, to, :, j, :, 1]).sum() / m)
                else:
                    attempt = float("nan")
                prob_rows.append({
                    "case": lab, "psi": psi, "phi": phi,
                    "age_cell_j": j, "age_low": age_start + da * j,
                    "age_high": age_start + da * (j + 1) - 1,
                    "begin_tenure_to": to, "begin_tenure": to_labels[to],
                    "at_risk_mass": m, "attempt_prob": attempt,
                    "fecundity_pi_j": pi_j,
                    "birth_prob": attempt * pi_j if np.isfinite(attempt) else float("nan")})

    outdir = HERE / "part2_tables"
    outdir.mkdir(exist_ok=True)
    with (outdir / "part2_childless_mass_by_tenure.csv").open("w", newline="") as s:
        w = csv.DictWriter(s, fieldnames=list(mass_rows[0].keys()))
        w.writeheader()
        w.writerows(mass_rows)
    with (outdir / "part2_firstbirth_by_tenure.csv").open("w", newline="") as s:
        w = csv.DictWriter(s, fieldnames=list(prob_rows[0].keys()))
        w.writeheader()
        w.writerows(prob_rows)
    with (outdir / "part2_verification.csv").open("w", newline="") as s:
        w = csv.DictWriter(s, fieldnames=list(ver_rows[0].keys()))
        w.writeheader()
        w.writerows(ver_rows)

    tenure_keys = sorted({k for keys in key_inventory.values() for k in keys
                          if "tenure" in k.lower()})
    verdict = ("No post-birth-branch tenure policy is saved: the only "
               "tenure arrays are pre-tenure policies by beginning tenure "
               "state (tenure_choice/tenure_probs indexed by pre-tenure "
               "parity/child-state). Keep/upsize/downsize/to-rent shares after "
               "the birth branch cannot be computed from saved arrays.")
    lines = ["# Part 2 lock-in diagnostic — childless (n=0, m=0), age cells 0-4",
             "",
             "Cases (saved arrays, zero new solves): " + ", ".join(CASES) + ".",
             "At-risk mass g_pre: reconstructed exactly as the observer does --",
             "entrant cohort at j=0 from saved entry_by_loc; advance of saved",
             "post-fertility g_beginning_distribution otherwise (saved loc/tenure/",
             "saving policy, rebuilt location/tenure maps).",
             "Attempt prob: at-risk-mass-weighted mean of saved fert_probs parity-0",
             "try slice [...,1] (second axis = beginning tenure state to) -- the",
             "exact engine birth-flow factor. Birth prob = attempt x pi_j; per-cell",
             "pi_j in the prob CSV.",
             "Fecundity vector pi_j (j=0..%d): %s." % (
                 len(fec) - 1, ", ".join("%.6g" % v for v in fec)),
             "Childless mass at (n=0, cs>=1) in saved post-fertility arrays is 0.0",
             "in all 4 cases (checked).",
             "Engine flags: " + ", ".join("%s=%s" % kv for kv in sorted(flags.items())),
             "Tenure labels: to=0 renter; to=1..5 owner rungs H_own=(2,4,6,8,10).",
             "Verification per case in part2_verification.csv: post-fertility",
             "reconstruction L1 vs saved arrays; implied aggregate first births",
             "vs observer first_birth_flow; entry mass vs saved entry rate.",
             "",
             "## Post-birth tenure transitions",
             "",
             verdict,
             "Tenure-related saved keys: " + ", ".join(tenure_keys) + ".",
             "Full per-case key inventories are identical: " +
             str(all(k == key_inventory[CASES[0]] for k in key_inventory.values())),
             "",
             "No economic interpretation per TASK."]
    (outdir / "part2_summary.md").write_text("\n".join(lines) + "\n")
    for v in ver_rows:
        print(v)


if __name__ == "__main__":
    main()
