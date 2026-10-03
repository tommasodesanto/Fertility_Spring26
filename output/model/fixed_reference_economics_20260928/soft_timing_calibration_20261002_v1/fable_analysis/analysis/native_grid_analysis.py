"""Read-only native-grid analysis of the two verified soft-timing winners.

No model solve. Loads the exact-repeat saved solution arrays of the chosen
48-chain winners (original chain 15, alternative chain 13), reproduces the
early-fertility observer from the arrays, decomposes the age-25 gap, and
slices fertility/tenure policies and the debt floor on the native grid.

Run from anywhere:
  MPLBACKEND=Agg code/model/.venv/bin/python <this file>
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FA = HERE.parent
PACKET = FA.parent
ROOT = PACKET.parents[2]
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)

ARMS = {
    "original": PACKET / "collection/production_original_chain_15/run/native_postcheck/selected_postcheck/phase_b_ge",
    "alternative": PACKET / "collection/production_alternative_chain_13/run/native_postcheck/selected_postcheck/phase_b_ge",
}
# Fixed objects from CALIBRATION_STATUS.md / parameters.csv (identical in both arms).
R = 1.08243216
PHI = 0.8
PSI_SELL = 0.06
HOLD = 0.05545379079326218 + 0.042393443095490375  # period depreciation + property tax
C_MIN = 0.04
H_OWN = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
H_P = {"original": 2.6, "alternative": 2.593759507364224}
# Four-year income by age cell (run_model.py EXTERNAL_INPUTS; same profile in both arms)
INCOME_AGE = np.array([
    2.650830656801071, 2.650830656801071, 3.4664708588937074, 3.4664708588937074,
    4.078201010463186, 4.078201010463186, 4.078201010463186, 4.017027995306238,
    4.017027995306238, 4.017027995306238, 3.8131179447830785, 3.8131179447830785,
    0.917784047463731, 0.917784047463731, 0.917784047463731, 0.917784047463731,
    0.917784047463731])
AGES = 18 + 4 * np.arange(17)

# CPS age-25 empirical object (early_fertility_target.json, authoritative candidate 25)
DATA25 = dict(mean_capped3=0.8095276384290021, any_birth=0.45725426281169845,
              share_above3=0.029163597533520928, se=0.02803450449153517)
DATA_BY_AGE = {24: 0.6701713286931946, 25: 0.8095276384290021, 26: 0.9226758406059775}


def wsum(a, axes):
    return np.sum(a, axis=axes)


def parity_by_age(dist):
    # axes: b, tenure, loc, age, income, n, m
    return wsum(dist, (0, 1, 2, 4, 6))


def analyze(arm: str) -> dict:
    root = ARMS[arm]
    z = np.load(root / "selected_repeat/stage/solution_arrays.npz", allow_pickle=True)
    obs = json.load(open(root / "selected_root/observers.json"))["fertility"]["uniform_birth_time"]
    acc = obs["accounting"]
    b = z["b_grid"]
    p = float(z["p_eq"][0])
    g = z["g"]                                   # realized post-tenure cross-section
    g_beg = z["g_beginning_distribution"]        # post-fertility, pre-tenure
    g_cross = z["g_cross_sectional_wealth_distribution"]
    g_stay = z["g_stay_distribution"]
    fert = z["fert_probs"]                       # [..., 1] = first-birth attempt prob
    tprob = z["tenure_probs"].astype(float)      # [..., inherited tenure, ..., chosen 6]
    tchoice = z["tenure_choice"]
    bp = z["bp_pol"]; bp_stay = z["bp_pol_stay"]
    zval = z["type_values"]
    out: dict = {"arm": arm, "price": p}

    # ---- 1. Which saved distribution is pre-fertility? Identity check vs observers.json
    pre_obs = np.array(acc["pre_parity_mass_by_age"]); post_obs = np.array(acc["post_parity_mass_by_age"])
    cands = {"g": g, "g_beginning_distribution": g_beg, "g_cross_sectional_wealth_distribution": g_cross}
    ident = {}
    for name, d in cands.items():
        pa = parity_by_age(d)
        ident[name] = {"max_abs_diff_vs_pre": float(np.max(np.abs(pa - pre_obs))),
                       "max_abs_diff_vs_post": float(np.max(np.abs(pa - post_obs)))}
    out["distribution_identity"] = ident
    # Finding: every saved distribution reproduces the POST-fertility stock; no
    # pre-fertility array is saved. Use the verified observer accounting for the
    # pre/post stocks and the saved post-fertility, pre-tenure array for weights.
    out["pre_fertility_array"] = "not saved; observers.json pre_parity_mass_by_age"
    out["post_fertility_array"] = "g_beginning_distribution (post-fertility, pre-tenure)"
    pre = pre_obs; post = post_obs

    # ---- 2. Reproduce the early-fertility observer and decompose
    proj25 = 0.125 * pre[1] + 0.875 * post[1]
    shares25 = proj25 / proj25.sum()
    mean25 = float(np.dot([0, 1, 2, 3], shares25))
    out["early_fertility_reproduced"] = mean25
    out["early_fertility_observer"] = obs["moments"]["mean_children_ever_born_capped3_age25"]
    any25 = float(shares25[1:].sum())
    out["age25"] = {
        "model_shares_0_1_2_3": shares25.tolist(), "model_any_birth": any25,
        "model_children_per_mother": mean25 / any25,
        "data_any_birth": DATA25["any_birth"],
        "data_children_per_mother_capped3": DATA25["mean_capped3"] / DATA25["any_birth"],
        "data_share_two_plus_bounds": [(DATA25["mean_capped3"] - DATA25["any_birth"]) / 2.0,
                                       DATA25["mean_capped3"] - DATA25["any_birth"] - DATA25["share_above3"]],
        "model_share_two_plus": float(shares25[2:].sum()),
    }
    # Hazards and flows
    flows = np.array(acc["parity_birth_flows_by_age"])  # J x 3: first, second, third+ entry
    mass = pre.sum(axis=1)
    h1 = flows[:, 0] / np.where(pre[:, 0] > 0, pre[:, 0], np.nan)
    h2 = flows[:, 1] / np.where(pre[:, 1] > 0, pre[:, 1], np.nan)
    h3 = flows[:, 2] / np.where(pre[:, 2] > 0, pre[:, 2], np.nan)
    out["first_birth_hazard_by_cell"] = np.nan_to_num(h1[:8]).tolist()
    out["second_birth_hazard_by_cell"] = np.nan_to_num(h2[:8]).tolist()
    out["third_birth_hazard_by_cell"] = np.nan_to_num(h3[:8]).tolist()
    fb_share = flows[:, 0] / flows[:, 0].sum()
    out["first_birth_cell_shares"] = fb_share[:8].tolist()
    out["first_birth_mean_age_midpoint"] = float(np.dot(AGES + 2, fb_share))
    # Ceiling: every 18-21 mother has a second birth at 22-25; first-birth timing unchanged
    P0 = flows[0, 0] / mass[0]; P1 = flows[1, 0] / mass[1]
    out["age25_ceiling_given_first_birth_timing"] = float(P0 * (1 + 0.875) + 0.875 * P1)
    # Share of women needing a first birth at 18-21 to reach the target when
    # every such mother has a second birth at 22-25 and the 22-25 first-birth
    # hazard of the remaining childless is held at its model value.
    h1_cell1 = flows[1, 0] / pre[1, 0]
    out["age25_first_birth_18_21_share_needed_for_target"] = float(
        (DATA25["mean_capped3"] - 0.875 * h1_cell1) / (1.875 - 0.875 * h1_cell1))
    out["age25_first_birth_18_21_share_model"] = float(P0)
    # Alternative projections (robustness of the 0.875 convention)
    out["age25_if_post_cell_stock"] = float(np.dot([0, 1, 2, 3], post[1]) / mass[1])
    out["age25_if_pre_cell_stock"] = float(np.dot([0, 1, 2, 3], pre[1]) / mass[1])
    # Children ever born by age cell (pre and post), for the lifecycle figure
    out["ceb_by_cell_pre"] = (pre @ np.array([0, 1, 2, 3]) / mass).tolist()
    out["ceb_by_cell_post"] = (post @ np.array([0, 1, 2, 3]) / mass).tolist()

    # ---- 3. Young childless renters on the native grid.
    # Weights: post-fertility, pre-tenure distribution (inherited-tenure axis).
    # Childless mass here excludes households whose birth occurred this period,
    # so weighted mean attempt probabilities are selection-understated; the
    # by-income-state and by-wealth policy values themselves are exact.
    w = g_beg
    young = {}
    pre_identity = {}
    for j in (0, 1, 2):
        # Invert the within-period first-birth selection: post childless mass
        # equals pre mass times (1 - pi_j * attempt), with the engine's
        # fecundity pi_j = 1 - omega1 exp(omega2 (age-18)) (parameters.py).
        pi_j = 1.0 - 0.02 * np.exp(0.134 * (AGES[j] - 18))
        all_post = w[:, :, 0, j, :, 0, 0]
        all_att = fert[:, :, 0, j, :, 1]
        all_pre = all_post / np.clip(1.0 - pi_j * all_att, 1e-9, None)
        pre_identity[int(AGES[j])] = {"reconstructed_pre_childless_mass": float(all_pre.sum()),
                                      "observer_pre_childless_mass": float(pre[j, 0]),
                                      "post_childless_mass": float(all_post.sum()),
                                      "fecundity_pi_j": float(pi_j)}
        wj = all_pre[:, 0]                       # inherited renters, childless, PRE-fertility: b x z
        fj = fert[:, 0, 0, j, :, 1]              # first-birth attempt prob: b x z
        mj = wj.sum()
        cum = np.cumsum(wj.sum(axis=1)) / mj
        lo = b[np.searchsorted(cum, 0.005)]; hi = b[min(np.searchsorted(cum, 0.995), b.size - 1)]
        # attempt prob by income state
        by_z = np.where(wj.sum(0) > 0, (wj * fj).sum(0) / np.maximum(wj.sum(0), 1e-300), np.nan)
        # by wealth tercile (pooled)
        wb = wj.sum(1); cb = np.cumsum(wb) / mj
        terc = np.digitize(cb, [1 / 3, 2 / 3])
        by_terc = [float((wj[terc == t] * fj[terc == t]).sum() / wj[terc == t].sum()) for t in range(3)]
        # purchase feasibility of the 4-room rung (smallest rung above the parent floor)
        Q4 = 4.0 * p
        y = INCOME_AGE[j] * zval                 # four-year income by state
        B = b[:, None]; Y = y[None, :]
        if arm == "original":
            feas_soft = B + Y / R >= (1 - PHI / R + HOLD / R) * Q4 + C_MIN / R  # end floor + c_min
        else:
            feas_soft = B + (Y - C_MIN - HOLD * Q4) / R >= -(PHI - 1) * Q4 / R
        feas_screen = R * B + Y >= (1 - PHI) * Q4   # executed soft closing screen
        feas_stock = B >= (1 - PHI) * Q4            # DUE/SSV-style cash-in-hand test (not implemented)
        own_prob = tprob[:, 0, 0, j, :, 0, 0, 1:].sum(-1)
        def share(mask): return float((wj * mask).sum() / mj)
        def fert_given(mask):
            m = (wj * mask).sum(); return float((wj * mask * fj).sum() / m) if m > 0 else float("nan")
        young[int(AGES[j])] = {
            "mass": float(mj), "occupied_b_range_99pct": [float(lo), float(hi)],
            "share_b_le_0": float(wj[b <= 1e-12].sum() / mj),
            "attempt_prob_mean": float((wj * fj).sum() / mj),
            "attempt_prob_by_income_state": np.nan_to_num(by_z).tolist(),
            "attempt_prob_by_wealth_tercile": by_terc,
            "income_state_mass_shares": (wj.sum(0) / mj).tolist(),
            "Q_4room": Q4,
            "share_feasible_end_floor_rule": share(feas_soft),
            "share_passing_closing_screen": share(feas_screen),
            "share_passing_stock_20pct_down_test": share(feas_stock),
            "attempt_prob_if_stock_feasible": fert_given(feas_stock),
            "attempt_prob_if_stock_infeasible": fert_given(~feas_stock),
            "attempt_prob_if_endfloor_feasible": fert_given(feas_soft),
            "attempt_prob_if_endfloor_infeasible": fert_given(~feas_soft),
            "ownership_choice_prob_mean": float((wj * own_prob).sum() / mj),
            "ownership_choice_prob_by_income_state": np.nan_to_num(
                (wj * own_prob).sum(0) / np.maximum(wj.sum(0), 1e-300)).tolist(),
        }
        # Within-income wealth gradient of the attempt probability at the modal income states
        grad = {}
        for zz in (3, 4, 5):
            col = wj[:, zz]; cz = np.cumsum(col) / col.sum()
            tz = np.digitize(cz, [1 / 3, 2 / 3])
            grad[int(zz + 1)] = [float((col[tz == t] * fj[tz == t, zz]).sum() / col[tz == t].sum())
                                 if col[tz == t].sum() > 0 else None for t in range(3)]
        young[int(AGES[j])]["attempt_prob_wealth_terciles_within_income_state"] = grad
    out["young_childless_renters"] = young
    out["pre_childless_identity_check"] = pre_identity

    # ---- 4. Debt floor: owners at b' = -phi * p * h (realized post-tenure, by age)
    floor = -PHI * p * H_OWN
    bind = []
    for j in range(17):
        tot = 0.0; at = 0.0; near = 0.0
        for k in range(1, 6):
            gs = g_stay[:, k, 0, j]; gm = g[:, k, 0, j] - gs
            gm = np.clip(gm, 0, None)
            for dist, pol in ((gs, bp_stay[:, k, 0, j]), (gm, bp[:, k, 0, j])):
                tot += dist.sum()
                at += (dist * (pol <= floor[k - 1] + 1e-6)).sum()
                near += (dist * (pol <= floor[k - 1] + 0.05 * p * H_OWN[k - 1])).sum()
        bind.append({"age": int(AGES[j]), "owner_mass": float(tot),
                     "share_at_floor": float(at / tot) if tot > 0 else None,
                     "share_within_5pct_of_Q_of_floor": float(near / tot) if tot > 0 else None})
    out["owner_debt_floor_by_age"] = bind
    # Occupied region overall
    pooled = wsum(w, (1, 2, 3, 4, 5, 6)); cum = np.cumsum(pooled) / pooled.sum()
    out["occupied_region_pooled_99pct"] = [float(b[np.searchsorted(cum, 0.005)]),
                                           float(b[min(np.searchsorted(cum, 0.995), b.size - 1)])]
    out["endpoint_mass"] = [float(pooled[0] / pooled.sum()), float(pooled[-1] / pooled.sum())]
    # Births by realized tenure among new mothers at 22-25 and 26-29 (post-tenure g, n=1,m=1 minus inherited)
    out["own_rate_new_mothers_proxy"] = {}
    for j in (1, 2):
        new_m = g[:, :, 0, j, :, 1, 1].sum((0, 2))
        out["own_rate_new_mothers_proxy"][int(AGES[j])] = float(new_m[1:].sum() / new_m.sum())
    out["own_rate_childless_by_cell"] = [float(g[:, 1:, 0, j, :, 0, 0].sum() / max(g[:, :, 0, j, :, 0, 0].sum(), 1e-300)) for j in range(4)]
    out["mean_b_childless_by_cell"] = [float((w[:, :, 0, j, :, 0, 0].sum((1, 2)) * b).sum() / max(w[:, :, 0, j, :, 0, 0].sum(), 1e-300)) for j in range(4)]
    return out


def main():
    res = {arm: analyze(arm) for arm in ARMS}
    json.dump(res, open(OUT / "native_grid_analysis.json", "w"), indent=1)
    # ---- figures
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    o, a = res["original"], res["alternative"]
    # F1 age-25 decomposition
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.8))
    x = np.arange(4); wd = 0.27
    ax[0].bar(x - wd, o["age25"]["model_shares_0_1_2_3"], wd, label="model, original timing")
    ax[0].bar(x, a["age25"]["model_shares_0_1_2_3"], wd, label="model, alternative timing")
    lo2, hi2 = DATA25_two = o["age25"]["data_share_two_plus_bounds"]
    ax[0].bar([0 + wd], [1 - DATA25["any_birth"]], wd, color="k", alpha=.6, label="CPS 2004/06, age 25")
    ax[0].errorbar([2 + wd], [(lo2 + hi2) / 2], yerr=[[(hi2 - lo2) / 2]], fmt="ko", capsize=4)
    ax[0].text(2 + wd, hi2 + 0.02, "CPS 2+ (bounds)", ha="center", fontsize=8)
    ax[0].set_xticks(x); ax[0].set_xticklabels(["0", "1", "2", "3+"]); ax[0].set_xlabel("children ever born at age 25")
    ax[0].set_ylabel("share of women"); ax[0].legend(fontsize=7); ax[0].set_title("Age-25 children-ever-born shares", fontsize=10)
    labels = ["CPS", "orig.", "alt.", "bound"]
    vals = [DATA25["mean_capped3"], o["early_fertility_reproduced"], a["early_fertility_reproduced"], o["age25_ceiling_given_first_birth_timing"]]
    ax[1].bar(labels, vals, color=["k", "C0", "C1", "C7"])
    ax[1].axhline(DATA25["mean_capped3"] - 2 * DATA25["se"], ls=":", c="k", label="CPS ±2 working bootstrap SE")
    ax[1].axhline(DATA25["mean_capped3"] + 2 * DATA25["se"], ls=":", c="k")
    ax[1].tick_params(axis="x", labelrotation=20)
    ax[1].legend(fontsize=7, loc="upper right")
    ax[1].set_ylabel("mean children ever born, capped at 3"); ax[1].set_title("Target, model, and bound at held first-birth flows", fontsize=10)
    fig.suptitle("Supplemental: early-fertility decomposition (no model solve)", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / "F1_age25_decomposition.png", dpi=150); plt.close(fig)
    # F2 children ever born by age cell vs CPS points
    fig, ax = plt.subplots(figsize=(6, 4))
    for r, c in ((o, "C0"), (a, "C1")):
        ax.step(AGES[:8], r["ceb_by_cell_pre"][:8], where="post", color=c, alpha=.5)
        ax.plot(AGES[:8] + 4, r["ceb_by_cell_post"][:8], "o--", color=c, ms=4, label=f"model {r['arm']} (end of cell)")
    ax.plot(list(DATA_BY_AGE), list(DATA_BY_AGE.values()), "ks", label="CPS 2004/06 completed interview age")
    ax.plot([25.5], [o["early_fertility_reproduced"]], "C0*", ms=12, label="completed age 25 [25,26); model 0.875 interpolation")
    ax.set_xlabel("age"); ax.set_ylabel("mean children ever born (capped 3)"); ax.legend(fontsize=7)
    ax.set_title("Supplemental: children ever born by age, model cells vs CPS", fontsize=10)
    fig.tight_layout(); fig.savefig(OUT / "F2_ceb_by_age.png", dpi=150); plt.close(fig)
    # F3 first-birth attempt probability by income state and wealth tercile, ages 18-21 and 22-25
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for i, age in enumerate((18, 22)):
        for r, c in ((o, "C0"), (a, "C1")):
            y = r["young_childless_renters"][age]
            ax[i].plot(range(1, 10), y["attempt_prob_by_income_state"], "o-", color=c, label=f"{r['arm']}: by income state")
            ax[i].bar(range(1, 10), y["income_state_mass_shares"], color=c, alpha=.15)
        ax[i].set_title(f"Childless inherited renters, ages {age}-{age+3}", fontsize=10); ax[i].set_xlabel("income state (1-9); shaded = mass share")
    ax[0].set_ylabel("first-birth attempt probability (pre-tenure weights)"); ax[0].legend(fontsize=7)
    fig.suptitle("Supplemental: first-birth attempt probability on the native grid", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / "F3_first_birth_by_income.png", dpi=150); plt.close(fig)
    # F4 debt floor by age
    fig, ax = plt.subplots(figsize=(6, 4))
    for r, c in ((o, "C0"), (a, "C1")):
        d = r["owner_debt_floor_by_age"]
        ax.plot([e["age"] for e in d], [e["share_at_floor"] or 0 for e in d], "o-", color=c, label=f"{r['arm']}: at floor")
        ax.plot([e["age"] for e in d], [e["share_within_5pct_of_Q_of_floor"] or 0 for e in d], "s--", color=c, alpha=.5, label=f"{r['arm']}: within 5% of Q")
    ax.set_xlabel("age cell start"); ax.set_ylabel("share of realized owners"); ax.legend(fontsize=7)
    ax.set_title("Supplemental: owner ending-debt floor b' = -phi Q, binding mass by age", fontsize=10)
    fig.tight_layout(); fig.savefig(OUT / "F4_debt_floor_by_age.png", dpi=150); plt.close(fig)
    print(json.dumps({k: {kk: v[kk] for kk in ("price", "pre_fertility_array", "post_fertility_array", "early_fertility_reproduced", "early_fertility_observer", "age25", "age25_ceiling_given_first_birth_timing", "age25_first_birth_18_21_share_needed_for_target", "age25_first_birth_18_21_share_model", "age25_if_post_cell_stock", "first_birth_hazard_by_cell", "second_birth_hazard_by_cell", "first_birth_cell_shares", "occupied_region_pooled_99pct", "endpoint_mass", "own_rate_new_mothers_proxy", "own_rate_childless_by_cell", "mean_b_childless_by_cell")} for k, v in res.items()}, indent=1))
    for arm in res:
        print("====", arm)
        for age, y in res[arm]["young_childless_renters"].items():
            print(age, json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in y.items() if k not in ("attempt_prob_wealth_terciles_within_income_state",)}))
            print("   within-income terciles:", y["attempt_prob_wealth_terciles_within_income_state"])
        print("floor by age:", [(e["age"], round(e["share_at_floor"] or 0, 4)) for e in res[arm]["owner_debt_floor_by_age"][:8]])
        print("distribution identity:", res[arm]["distribution_identity"])
        print("pre-childless identity:", json.dumps(res[arm]["pre_childless_identity_check"]))


if __name__ == "__main__":
    main()
