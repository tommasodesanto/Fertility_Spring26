"""Saved h_P-only effect on age-22 median-earner renters; no model solve."""
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
sys.dont_write_bytecode = True
sys.path[:0] = [str(ROOT / "output/model/credit_mechanism_20261004/diagnostics"), str(ROOT / "code/model")]
import extract_common_states as ecs  # noqa: E402
from production.engine.household import build_forward_tenure_transition_maps  # noqa: E402

CASES = {"baseline": HERE.parent.parent / "price_decomposition/cases/baseline_resolve",
         "higher_floor": HERE.parent.parent / "price_decomposition/cases/arm_floor_only"}
J = 1  # ages 22-25
Z = 4  # fifth of nine Markov earnings states


def extract(case: Path, state_indices: np.ndarray) -> tuple[list[dict], dict]:
    p, a, pre, rate, ages, _, checks = ecs.load(case)
    if (p["I"], p["native_purchase_income"], p["native_due_stayer_credit"],
        p["child_state_mode"], p.get("child_maturation_mode", "constant")) != (
        1, True, True, "independent_count", "constant"):
        raise ValueError("Unexpected saved timing/state contract")
    if p["use_pti_constraint"] or p["birth_dp_grant"]:
        raise ValueError("Unexpected purchase override")
    keys = ("tenure_probs", "hR_pol", "bp_pol", "bp_pol_stay", "c_pol", "c_pol_stay",
            "shared.phi_choice", "shared.birth_dp", "shared.birth_entry_grant", "shared.cb_flat", "shared.hb_flat")
    with np.load(case / "solution_arrays.npz", allow_pickle=False) as z:
        v = {key: z[key].copy() for key in keys}
    grid = a["b_grid"]
    price = float(a["p_eq"][0])
    house = np.asarray(p["H_own"], dtype=float)
    R = float(p["R_gross"])
    hc = np.r_[0., price * house][None, :]
    he = np.r_[0., (1 - float(p["psi"])) * price * house][None, :]
    pp = SimpleNamespace(native_purchase_income=True,
                         propagate_birth_entry_grant=p["propagate_birth_entry_grant"], R_gross=R)
    idx, wt = build_forward_tenure_transition_maps(
        pp, grid, hc, he, v["shared.phi_choice"], v["shared.birth_dp"], v["shared.birth_entry_grant"])
    y = float(np.asarray(p["income"])[0, J] * a["type_values"][Z]
              + p["property_tax_lump_sum_transfer"])
    rent = (R - 1 + float(p["delta"]) + float(p["tau_H"])) * price
    rows = []
    max_budget_error = 0.
    min_consumption_surplus = np.inf
    min_housing_surplus = np.inf
    min_saving_floor_slack = np.inf
    for bi in state_indices:
        b = float(grid[bi])
        n, m = 1, 1  # successful first-birth branch
        raw = np.asarray(v["tenure_probs"][bi, 0, 0, J, Z, n, m], dtype=float)
        total = float(raw.sum())
        if total <= 0 or not np.isfinite(total):
            raise ValueError(f"Missing successful-birth tenure menu at b={b}")
        pr = raw / total
        expected_h = expected_bp = expected_c = 0.
        for tn, q in enumerate(pr):
            if q <= 0:
                continue
            lo = int(idx[0, 0, tn, n, m, bi])
            w = float(wt[0, 0, tn, n, m, bi])
            parts = ((lo, 1 - w), (lo + 1, w))
            carray = v["c_pol_stay"] if tn == 0 else v["c_pol"]
            barray = v["bp_pol_stay"] if tn == 0 else v["bp_pol"]
            c = sum(ww * float(carray[ii, tn, 0, J, Z, n, m]) for ii, ww in parts if ww > 0)
            bp = sum(ww * float(barray[ii, tn, 0, J, Z, n, m]) for ii, ww in parts if ww > 0)
            if tn == 0:
                h = sum(ww * float(v["hR_pol"][ii, tn, 0, J, Z, n, m])
                        for ii, ww in parts if ww > 0)
                cost = rent * h
                floor = 0.
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
                raise ValueError(f"Infeasible positive-probability choice at b={b}, tn={tn}")
            budget = R * b + y - hc[0, tn] - cost - c - bp
            max_budget_error = max(max_budget_error, abs(float(budget)))
            if abs(budget) > 1e-7:
                raise ValueError(f"Budget fails at b={b}, tn={tn}: {budget}")
            expected_h += q * h
            expected_bp += q * bp
            expected_c += q * c
        rows.append(dict(beginning_wealth=b, prebirth_exposure=float(pre[bi, 0, 0, J, Z, 0, 0]),
                         first_birth_probability=float(rate[bi, 0, 0, J, Z, 0, 0]),
                         successful_birth_physical_rooms=expected_h,
                         successful_birth_ending_assets=expected_bp,
                         successful_birth_consumption=expected_c,
                         successful_birth_owner_probability=float(pr[1:].sum())))
    receipt = dict(case=str(case.resolve()), arrays_sha256=ecs.sha(case / "solution_arrays.npz"),
                   parameters_sha256=ecs.sha(case / "executed_P.json"), price=price,
                   h_P=float(p["hbar_first_child_jump"]), hR_max=float(p["hR_max"]),
                   max_budget_error=max_budget_error,
                   min_consumption_surplus=min_consumption_surplus,
                   min_housing_surplus=min_housing_surplus,
                   min_saving_floor_slack=min_saving_floor_slack,
                   inversion_checks=checks)
    return rows, receipt


def main() -> None:
    p0, a0, pre0, _, _, _, _ = ecs.load(CASES["baseline"])
    p1, a1, _, _, _, _, _ = ecs.load(CASES["higher_floor"])
    if not np.array_equal(a0["b_grid"], a1["b_grid"]) or not np.array_equal(a0["p_eq"], a1["p_eq"]):
        raise ValueError("Unmatched grid or price")
    diff = [k for k in p0 if not k.startswith("_") and p0[k] != p1[k]]
    if diff != ["hbar_first_child_jump"]:
        raise ValueError(f"Non-floor differences: {diff}")
    four_year_aftertax_income = float(
        np.asarray(p0["income"])[0, J] * a0["type_values"][Z]
        + p0["property_tax_lump_sum_transfer"]
    )
    annual_gross_earnings = four_year_aftertax_income / (
        float(p0["period_years"]) * (1 - float(p0["tau_pay"]))
    )
    if annual_gross_earnings <= 0:
        raise ValueError("Invalid earnings denominator")
    exposure = pre0[:, 0, 0, J, Z, 0, 0]
    positive = np.flatnonzero(exposure > 0)
    cutoff = int(np.searchsorted(np.cumsum(exposure), .995 * exposure.sum()))
    displayed = positive[(positive <= cutoff) & (a0["b_grid"][positive] >= 0)]
    if len(displayed) < 10 or float(exposure[displayed].sum() / exposure.sum()) < .5:
        raise ValueError("Insufficient occupied plotted support")
    values, receipts = {}, {}
    for label, case in CASES.items():
        values[label], receipts[label] = extract(case, positive)
        if [r["beginning_wealth"] for r in values[label]] != [float(a0["b_grid"][i]) for i in positive]:
            raise ValueError("Misaligned wealth states")
    with (HERE / "main_floor_state.csv").open("w", newline="") as f:
        fields = ["beginning_wealth", "beginning_wealth_to_annual_earnings",
                  "prebirth_exposure", "displayed"]
        for label in CASES:
            fields += [f"{label}_{k}" for k in values[label][0] if k not in ("beginning_wealth", "prebirth_exposure")]
            fields += [f"{label}_successful_birth_ending_assets_to_annual_earnings"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for bi, row0, row1 in zip(positive, values["baseline"], values["higher_floor"]):
            out = dict(beginning_wealth=row0["beginning_wealth"],
                       beginning_wealth_to_annual_earnings=row0["beginning_wealth"] / annual_gross_earnings,
                       prebirth_exposure=row0["prebirth_exposure"],
                       displayed=int(bi in displayed))
            for label, row in (("baseline", row0), ("higher_floor", row1)):
                out.update({f"{label}_{k}": v for k, v in row.items()
                            if k not in ("beginning_wealth", "prebirth_exposure")})
                out[f"{label}_successful_birth_ending_assets_to_annual_earnings"] = (
                    row["successful_birth_ending_assets"] / annual_gross_earnings
                )
            writer.writerow(out)

    plt.rcParams.update({"font.size": 14, "axes.labelsize": 15, "xtick.labelsize": 13,
                         "ytick.labelsize": 13, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(11, 4), layout="constrained")
    names = ("first_birth_probability", "successful_birth_physical_rooms", "successful_birth_ending_assets")
    panel_titles = ("First birth", "Housing after birth", "Saving after birth")
    ylabels = ("Probability (%)", "Physical rooms", "Ending wealth /\nannual earnings")
    raw_x = [float(a0["b_grid"][bi]) for bi in displayed]
    x = [v / annual_gross_earnings for v in raw_x]
    for ax, key, title, ylabel in zip(axes, names, panel_titles, ylabels):
        for label, color, line_label in (("baseline", "#1f4e79", "Baseline"),
                                         ("higher_floor", "#bf6b2b", "Higher space requirement")):
            lookup = {r["beginning_wealth"]: r for r in values[label]}
            y = [lookup[b][key] * (100 if key == "first_birth_probability" else
                                   1 / annual_gross_earnings if key == "successful_birth_ending_assets" else 1)
                 for b in raw_x]
            ax.plot(x, y, color=color, lw=2.6, label=line_label)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="0.88", lw=.7)
        ax.set_xlim(float(x[0]), float(x[-1]))
    axes[0].legend(frameon=False, fontsize=11, loc="best")
    fig.supxlabel("Beginning wealth / annual earnings")
    fig.savefig(HERE / "main_floor_state.pdf")
    fig.savefig(HERE / "main_floor_state.png", dpi=220)
    plt.close(fig)

    receipt = dict(reference="Post-interest chain 13, no Estate-A; fixed-price h_P-only diagnostic",
                   age_cell_start=22, earnings_state_index_zero_based=Z,
                   initial_tenure="renter", initial_children_ever_born=0,
                   initial_children_at_home=0, success_branch=(1, 1),
                   floor_change=[p0["hbar_first_child_jump"], p1["hbar_first_child_jump"]],
                   four_year_aftertax_income=four_year_aftertax_income,
                   annual_gross_earnings_denominator=annual_gross_earnings,
                   earnings_denominator_definition="four-year aftertax income / [period_years * (1 - tau_pay)]; baseline denominator shared across arms",
                   exact_public_input_diff=diff, plotted_positive_states=len(displayed),
                   all_positive_states=len(positive),
                   plotted_baseline_prebirth_exposure_share=float(exposure[displayed].sum() / exposure.sum()),
                   all_positive_exposure_covered=float(exposure[positive].sum() / exposure.sum()),
                   case_receipts=receipts,
                   interpretation="All curves evaluate the same baseline-occupied beginning states. Physical rooms and b-prime integrate each arm's tenure probabilities and post-transaction controls conditional on successful first birth. No GE price adjustment or recalibration.")
    (HERE / "main_floor_state_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k: receipt[k] for k in ("floor_change", "plotted_positive_states", "plotted_baseline_prebirth_exposure_share", "all_positive_exposure_covered")}, indent=2))
    print(json.dumps({k: {x: v for x, v in r.items() if x.startswith("max_") or x.startswith("min_")}
                      for k, r in receipts.items()}, indent=2))


if __name__ == "__main__":
    main()
