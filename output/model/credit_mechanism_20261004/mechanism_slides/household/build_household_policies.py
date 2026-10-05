"""Plot saved, tenure-inclusive first-birth policies at occupied common states.

This is a zero-solve extraction.  It uses the engine's own forward tenure map
and the checked pre-birth exposure inversion.  Choices and controls are from the
October 3 post-interest chain-13, no-Estate-A, phi=.8 cached fixed-price case.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
DIAG = ROOT / "output/model/credit_mechanism_20261004/diagnostics"
sys.dont_write_bytecode = True
sys.path[:0] = [str(DIAG), str(ROOT / "code/model")]
import extract_common_states as ecs  # noqa: E402
from production.engine.household import build_forward_tenure_transition_maps  # noqa: E402

SOURCE = ecs.DEFAULT / "phi_080"
AGE_INDEX = 1                # age 22 (four-year age cell 22--25)
INCOME_INDEX = 4             # fifth of nine Markov earnings states
EXPOSURE_QUANTILE = 0.995
WAIT = "Waiting"
SUCCESS = "Successful first birth"
BLUE = "#1f4e79"
ORANGE = "#bf6b2b"


def main() -> None:
    p, a, pre, _, ages, fec, inversion = ecs.load(SOURCE)
    if (p["I"], p["native_purchase_income"], p["native_due_stayer_credit"],
        p["child_state_mode"], p.get("child_maturation_mode", "constant")) != (
            1, True, True, "independent_count", "constant"):
        raise ValueError("Saved state/timing contract has changed")
    if bool(p["use_pti_constraint"]) or bool(p["birth_dp_grant"]):
        raise ValueError("Unexpected purchase override")

    needed = ("tenure_probs", "hR_pol", "bp_pol", "bp_pol_stay", "c_pol",
              "c_pol_stay", "shared.phi_choice", "shared.birth_dp",
              "shared.birth_entry_grant", "shared.cb_flat", "shared.hb_flat")
    with np.load(SOURCE / "solution_arrays.npz", allow_pickle=False) as z:
        v = {key: z[key].copy() for key in needed}
    grid = a["b_grid"]
    price = float(a["p_eq"][0])
    house = np.asarray(p["H_own"], dtype=float)
    R = float(p["R_gross"])
    hc = np.r_[0., price * house][None, :]
    he = np.r_[0., (1. - float(p["psi"])) * price * house][None, :]
    pp = SimpleNamespace(native_purchase_income=True,
                         propagate_birth_entry_grant=p["propagate_birth_entry_grant"],
                         R_gross=R)
    idx, wt = build_forward_tenure_transition_maps(
        pp, grid, hc, he, v["shared.phi_choice"],
        v["shared.birth_dp"], v["shared.birth_entry_grant"])

    j, z = AGE_INDEX, INCOME_INDEX
    exposure = pre[:, 0, 0, j, z, 0, 0]
    assert exposure.sum() > 0 and fec[j] > 0
    age_income_mass = pre[:, :, 0, j, :, 0, 0].sum(axis=(0, 1))
    median_income_index = int(np.searchsorted(np.cumsum(age_income_mass), age_income_mass.sum() / 2))
    if z != median_income_index:
        raise ValueError("Illustrative earnings state is no longer the age-22 all-tenure childless-exposure median")
    positive = np.flatnonzero(exposure > 0)
    cutoff_pos = int(np.searchsorted(np.cumsum(exposure), EXPOSURE_QUANTILE * exposure.sum()))
    displayed = positive[positive <= cutoff_pos]
    if len(displayed) < 10:
        raise ValueError("Too little occupied support for policy figure")
    # Drop no positive-mass state: every retained node must have both menus.
    rows = []
    max_budget_error = 0.
    max_tenure_sum_error = 0.
    min_consumption_surplus = np.inf
    min_housing_surplus = np.inf
    min_saving_floor_slack = np.inf
    max_renter_stayer_control_difference = 0.
    excluded_exposure = 0.
    y = float(np.asarray(p["income"])[0, j] * a["type_values"][z]
               + p["property_tax_lump_sum_transfer"])
    annual_gross_earnings = (y - float(p["property_tax_lump_sum_transfer"])) / (
        float(p["period_years"]) * (1. - float(p["tau_pay"])))
    if annual_gross_earnings <= 0:
        raise ValueError("Annual gross earnings denominator must be positive")
    rent = (R - 1 + float(p["delta"]) + float(p["tau_H"])) * price
    for bi in positive:
        b = float(grid[bi])
        node = dict(beginning_wealth=b, age_cell_start=float(ages[j]),
                    beginning_wealth_over_annual_earnings=b / annual_gross_earnings,
                    income_state=z + 1, income_type_multiplier=float(a["type_values"][z]),
                    four_year_aftertax_income=y, annual_gross_earnings=annual_gross_earnings,
                    prebirth_exposure=float(exposure[bi]),
                    first_birth_probability=float(fec[j] * a["fert_probs"][bi, 0, 0, j, z, 1]),
                    displayed=int(bi in displayed))
        complete = True
        for label, n, m in ((WAIT, 0, 0), (SUCCESS, 1, 1)):
            raw = np.asarray(v["tenure_probs"][bi, 0, 0, j, z, n, m], dtype=float)
            total = float(raw.sum())
            if not np.isfinite(total) or total <= 0:
                complete = False
                break
            pr = raw / total
            max_tenure_sum_error = max(max_tenure_sum_error, abs(float(pr.sum()) - 1))
            expected_h = expected_bp = expected_c = 0.
            for tn, q in enumerate(pr):
                if q <= 0:
                    continue
                lo = int(idx[0, 0, tn, n, m, bi])
                w = float(wt[0, 0, tn, n, m, bi])
                parts = ((lo, 1. - w), (lo + 1, w))
                carray = v["c_pol_stay"] if tn == 0 else v["c_pol"]
                barray = v["bp_pol_stay"] if tn == 0 else v["bp_pol"]
                c = sum(ww * float(carray[ii, tn, 0, j, z, n, m]) for ii, ww in parts if ww > 0)
                bp = sum(ww * float(barray[ii, tn, 0, j, z, n, m]) for ii, ww in parts if ww > 0)
                if tn == 0:
                    h = sum(ww * float(v["hR_pol"][ii, tn, 0, j, z, n, m])
                            for ii, ww in parts if ww > 0)
                    cost = rent * h
                    floor = 0.
                    max_renter_stayer_control_difference = max(
                        max_renter_stayer_control_difference,
                        max(abs(v["c_pol_stay"][ii, tn, 0, j, z, n, m] -
                                v["c_pol"][ii, tn, 0, j, z, n, m]) for ii, ww in parts if ww > 0),
                        max(abs(v["bp_pol_stay"][ii, tn, 0, j, z, n, m] -
                                v["bp_pol"][ii, tn, 0, j, z, n, m]) for ii, ww in parts if ww > 0))
                else:
                    h = float(house[tn - 1])
                    cost = (float(p["delta"]) + float(p["tau_H"])) * price * h
                    floor = -float(v["shared.phi_choice"][0, tn, n, m]) * price * h
                cb = float(v["shared.cb_flat"][0, n + int(p["n_parity"]) * m])
                hb = float(v["shared.hb_flat"][0, n + int(p["n_parity"]) * m])
                if tn > 0:
                    hb *= float(p["owner_h_bar_scale"])
                min_consumption_surplus = min(min_consumption_surplus, c - cb)
                min_housing_surplus = min(min_housing_surplus, h - hb)
                min_saving_floor_slack = min(min_saving_floor_slack, bp - floor)
                if c - cb <= -1e-8 or h - hb <= -1e-8 or bp - floor < -1e-7:
                    raise ValueError(f"Positive-probability infeasible control at b={b}, branch={label}, tn={tn}")
                budget = R * b + y + he[0, 0] - hc[0, tn] - cost - c - bp
                max_budget_error = max(max_budget_error, abs(float(budget)))
                if abs(budget) > 1e-7:
                    raise ValueError(f"Positive-probability budget fails at b={b}, branch={label}, tn={tn}: {budget}")
                expected_h += q * h
                expected_bp += q * bp
                expected_c += q * c
            suffix = "wait" if n == 0 else "success"
            node[f"housing_{suffix}"] = expected_h
            node[f"ending_assets_{suffix}"] = expected_bp
            node[f"ending_assets_{suffix}_over_annual_earnings"] = expected_bp / annual_gross_earnings
            node[f"consumption_{suffix}"] = expected_c
            node[f"owner_probability_{suffix}"] = float(pr[1:].sum())
            node[f"renter_probability_{suffix}"] = float(pr[0])
        if complete:
            rows.append(node)
        else:
            excluded_exposure += float(exposure[bi])
    if excluded_exposure > 1e-12:
        raise ValueError("Positive prebirth exposure has a missing tenure menu")
    if len(rows) != len(positive):
        raise ValueError("Unexpected missing occupied state")
    with (HERE / "household_policies.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)

    shown = [row for row in rows if row["displayed"]]
    plt.rcParams.update({"font.size": 14, "axes.labelsize": 15, "axes.titlesize": 15,
                         "xtick.labelsize": 12, "ytick.labelsize": 12,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(11, 4), layout="constrained")
    x = [row["beginning_wealth_over_annual_earnings"] for row in shown]
    for ax, stem, ylabel in zip(
            axes, ("housing", "ending_assets", "owner_probability"),
            ("Physical rooms", "Ending wealth / annual earnings", "Probability")):
        for suffix, name, color in (("wait", WAIT, BLUE), ("success", SUCCESS, ORANGE)):
            field = f"{stem}_{suffix}_over_annual_earnings" if stem == "ending_assets" else f"{stem}_{suffix}"
            ax.plot(x, [row[field] for row in shown],
                    color=color, lw=2.6, label=name)
        ax.set_ylabel(ylabel)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="0.88", lw=.7)
        ax.set_xlim(float(x[0]), float(x[-1]))
    for ax, title in zip(axes, ("Housing", "Saving", "Ownership")):
        ax.set_title(title)
    axes[0].legend(loc="best", frameon=False, fontsize=11)
    fig.supxlabel("Beginning wealth / annual earnings", fontsize=14)
    fig.savefig(HERE / "household_policies.pdf")
    fig.savefig(HERE / "household_policies.png", dpi=220)
    plt.close(fig)

    # The native all-age matched-state comparison already exists in the checked
    # T1 tabulation.  Plot its age cells and verify its overall row from the
    # seven age rows before using it as the broad-exposure counterpart.
    t1_file = ROOT / "output/model/credit_mechanism_20261004/birth_tenure_tabulations/T1_i_matched_state_ownership_rooms.csv"
    with t1_file.open(newline="") as f:
        t1 = list(csv.DictReader(f))
    age_rows = [r for r in t1 if r["group_type"] == "age"]
    overall = next(r for r in t1 if r["group_type"] == "overall")
    if len(age_rows) != 7:
        raise ValueError("T1 age-cell table has changed")
    attempt = np.array([float(r["attempt_mass"]) for r in age_rows])
    if abs(attempt.sum() - float(overall["attempt_mass"])) > 1e-12:
        raise ValueError("T1 attempt mass does not sum")
    for field in ("rooms_mean_success", "rooms_mean_nobirth", "P_own_success", "P_own_nobirth"):
        calculated = np.average([float(r[field]) for r in age_rows], weights=attempt)
        if abs(calculated - float(overall[field])) > 1e-12:
            raise ValueError(f"T1 age rows do not recover {field}")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    xx = np.arange(len(age_rows))
    labels = [r["group"].removeprefix("age ") for r in age_rows]
    for ax, first, second, title, ylabel, factor in (
            (axes[0], "rooms_mean_nobirth", "rooms_mean_success", "Physical rooms", "Rooms", 1),
            (axes[1], "P_own_nobirth", "P_own_success", "Ownership", "Probability", 1)):
        ax.plot(xx, [factor * float(r[first]) for r in age_rows], "o-", color=BLUE, lw=2.4, label=WAIT)
        ax.plot(xx, [factor * float(r[second]) for r in age_rows], "o-", color=ORANGE, lw=2.4, label=SUCCESS)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xticks(xx, labels, rotation=35, ha="right")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="0.88", lw=.7)
    axes[0].legend(loc="best", frameon=False, fontsize=11)
    fig.supxlabel("Age at beginning of four-year period", fontsize=14)
    fig.savefig(HERE / "household_by_age.pdf")
    fig.savefig(HERE / "household_by_age.png", dpi=220)
    plt.close(fig)

    receipt = {
        "reference": "Post-interest soft chain 13, no Estate A, phi=0.8, fixed-price diagnostic",
        "zero_solves": True,
        "source_case": str(SOURCE.resolve()),
        "arrays_sha256": ecs.sha(SOURCE / "solution_arrays.npz"),
        "parameters_sha256": ecs.sha(SOURCE / "executed_P.json"),
        "age_cell_start": float(ages[j]),
        "age_cell_end": float(ages[j] + p["da"] - 1),
        "income_state_one_based": z + 1,
        "age22_all_tenure_childless_exposure_weighted_median_income_state_one_based": median_income_index + 1,
        "income_type_multiplier": float(a["type_values"][z]),
        "four_year_aftertax_income": y,
        "annual_gross_earnings": annual_gross_earnings,
        "wealth_figure_normalization": "Both beginning b and ending b-prime divide by this household's current annual gross earnings: (four-year aftertax income minus property-tax transfer)/[period_years*(1-tau_pay)]. Raw model-unit values remain in the CSV.",
        "earnings_definition": "Executed period after-tax P.income[location, age] * type_values[z] plus property-tax lump-sum transfer",
        "beginning_state": "Childless renter, location 0, same beginning b, age and Markov earnings state",
        "branches": {"waiting": "postbirth (n,m)=(0,0)", "successful_first_birth": "postbirth (n,m)=(1,1); constant child-maturation mode"},
        "controls": "Normalized tenure probabilities and the engine's two-node post-transaction wealth map; renter stay uses saved stayer controls; owner purchases use ordinary controls; renter rooms are included only when renter chosen",
        "first_birth_probability": "fecundity(age) * saved attempt probability, before tenure choice",
        "display_quantile": EXPOSURE_QUANTILE,
        "displayed_wealth_min": float(x[0]),
        "displayed_wealth_max": float(x[-1]),
        "displayed_wealth_max_model_units": float(shown[-1]["beginning_wealth"]),
        "occupied_nodes": len(rows),
        "displayed_nodes": len(shown),
        "positive_prebirth_exposure": float(exposure.sum()),
        "displayed_prebirth_exposure": float(sum(row["prebirth_exposure"] for row in shown)),
        "displayed_share_of_selected_age_income_exposure": float(sum(row["prebirth_exposure"] for row in shown) / exposure.sum()),
        "selected_age_income_share_of_all_fertile_childless_renter_exposure": float(exposure.sum() / pre[:, 0, 0, :int(p["A_f_end"]), :, 0, 0].sum()),
        "excluded_positive_exposure": excluded_exposure,
        "checks": {"inversion": inversion,
                   "max_tenure_probability_sum_error": max_tenure_sum_error,
                   "max_positive_choice_budget_error": max_budget_error,
                   "min_positive_choice_consumption_surplus": float(min_consumption_surplus),
                   "min_positive_choice_housing_surplus": float(min_housing_surplus),
                   "min_positive_choice_saving_floor_slack": float(min_saving_floor_slack),
                   "max_renter_stayer_control_difference": float(max_renter_stayer_control_difference)},
        "source_hashes": {str(path.resolve()): ecs.sha(path) for path in (
            Path(__file__), DIAG / "extract_common_states.py",
            ROOT / "code/model/production/engine/household.py",
            ROOT / "code/model/production/engine/distribution.py", t1_file)},
        "aggregate_t1": {
            "weight": "Baseline prebirth childless attempt mass across all fertile ages, all beginning wealth, tenure and earnings states",
            "attempt_mass": float(overall["attempt_mass"]),
            "physical_rooms_success": float(overall["rooms_mean_success"]),
            "physical_rooms_wait": float(overall["rooms_mean_nobirth"]),
            "ownership_probability_success": float(overall["P_own_success"]),
            "ownership_probability_wait": float(overall["P_own_nobirth"]),
            "source": str(t1_file.resolve()),
        },
    }
    (HERE / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"displayed_wealth_max": receipt["displayed_wealth_max"],
                      "coverage": receipt["displayed_share_of_selected_age_income_exposure"],
                      "checks": receipt["checks"]}, indent=2))


if __name__ == "__main__":
    main()
