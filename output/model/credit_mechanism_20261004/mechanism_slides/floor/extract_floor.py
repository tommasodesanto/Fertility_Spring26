"""Read the saved chain-13 fixed-price h_P-only diagnostic; no model solve."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
DECOMP = ROOT / "output/model/credit_mechanism_20261004/price_decomposition"
BASE = DECOMP / "cases/baseline_resolve"
FLOOR = DECOMP / "cases/arm_floor_only"
sys.path.insert(0, str(ROOT / "output/model/credit_mechanism_20261004/diagnostics"))
from extract_common_states import load  # noqa: E402 - exact pre-birth inversion


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    p0, a0, pre0, rate0, ages, _, check0 = load(BASE)
    p1, a1, pre1, rate1, _, _, check1 = load(FLOOR)
    changes = [k for k in p0 if not k.startswith("_") and p0[k] != p1[k]]
    if changes != ["hbar_first_child_jump"]:
        raise ValueError(f"Non-floor input differences: {changes}")
    if not np.array_equal(a0["b_grid"], a1["b_grid"]) or not np.array_equal(a0["p_eq"], a1["p_eq"]):
        raise ValueError("Price or wealth grid differs")
    audit = json.loads((FLOOR / "arm_spec_and_audit.json").read_text())["audit"]
    if audit["problems"] or list(audit["changed_fields"]) != changes:
        raise ValueError("Executed-input audit failed")

    # Identical boundaries for both arms, using the baseline's childless
    # pre-birth exposure throughout fertile ages, as in extract_common_states.
    earnings = np.asarray(p0["income"])[0, :, None] * a0["type_values"][None, :] / (
        p0["period_years"] * (1 - p0["tau_pay"])
    )
    age_income_mass = pre0[..., 0, 0].sum(axis=(0, 1, 2))
    age_income_mass[p0["A_f_end"] :] = 0
    xx, ww = earnings.ravel(), age_income_mass.ravel()
    order = np.argsort(xx, kind="stable")
    cumul = np.cumsum(ww[order])
    cuts = [float(xx[order[np.searchsorted(cumul, q * cumul[-1])]]) for q in (1 / 3, 2 / 3)]
    income_bin = np.searchsorted(cuts, earnings, side="left")
    b = a0["b_grid"]
    wealth_labels = ["≤ −1", "(−1, 0)", "0", "(0, .3]", "(.3, .6]", "(.6, 1]", "> 1"]
    wealth_bin = np.select(
        [b <= -1, (b > -1) & (b < 0), b == 0, (b > 0) & (b <= .3),
         (b > .3) & (b <= .6), (b > .6) & (b <= 1), b > 1],
        list(range(len(wealth_labels))), default=-1,
    )
    if np.any(wealth_bin < 0):
        raise ValueError("Unassigned wealth node")
    rows = []
    for j in (1, 3):
        for inc in range(3):
            zmask = income_bin[j] == inc
            group_mass = float(pre0[:, :, :, j, zmask, 0, 0].sum())
            for wb, label in enumerate(wealth_labels):
                bmask = wealth_bin == wb
                w = pre0[bmask, :, :, j, :, 0, 0][:, :, :, zmask]
                r0 = rate0[bmask, :, :, j, :, 0, 0][:, :, :, zmask]
                r1 = rate1[bmask, :, :, j, :, 0, 0][:, :, :, zmask]
                w1 = pre1[bmask, :, :, j, :, 0, 0][:, :, :, zmask]
                mass = float(w.sum())
                rows.append(dict(
                    age=int(ages[j]), earnings_group=inc + 1, wealth_group=label,
                    baseline_exposure=mass, baseline_group_share=mass / group_mass if group_mass else np.nan,
                    baseline_first_birth_probability=float(np.sum(w * r0) / mass) if mass else np.nan,
                    floor_first_birth_probability_at_baseline_states=float(np.sum(w * r1) / mass) if mass else np.nan,
                    policy_response_pp=float(100 * np.sum(w * (r1 - r0)) / mass) if mass else np.nan,
                    floor_exposure=float(w1.sum()),
                    price=float(a0["p_eq"][0]),
                ))
    HERE.mkdir(parents=True, exist_ok=True)
    with (HERE / "first_birth_by_age_earnings_wealth.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 6.5), sharey="row", layout="constrained")
    for j, age in enumerate((22, 30)):
        for inc in range(3):
            ax = axes[j, inc]
            rr = [r for r in rows if r["age"] == age and r["earnings_group"] == inc + 1]
            # Bins below 0.5% of this group's baseline exposure are omitted
            # from the picture; the complete weighted table retains them.
            x = [i for i, r in enumerate(rr) if r["baseline_group_share"] >= .005]
            y0 = [100 * rr[i]["baseline_first_birth_probability"] for i in x]
            y1 = [100 * rr[i]["floor_first_birth_probability_at_baseline_states"] for i in x]
            ax.plot(x, y0, "o-", color="#2c5b86", lw=1.8, label="Chain 13")
            ax.plot(x, y1, "s-", color="#b6503c", lw=1.8, label=r"$h_P$ +7.236%")
            ax.set_xticks(range(len(wealth_labels)), wealth_labels, rotation=45, ha="right")
            ax.set_title(f"Age {age} · earnings group {inc + 1}")
            ax.grid(alpha=.2)
            if inc == 0:
                ax.set_ylabel("First-birth probability (%)")
            if j == 1:
                ax.set_xlabel("Beginning net financial wealth")
    axes[0, 2].legend(frameon=False, fontsize=9)
    fig.suptitle("Higher family housing floor lowers first-birth probability at the same states", fontsize=13)
    fig.text(.5, -.005, "Baseline pre-birth distribution weights; fixed house price 0.776; 4-year age cells. Bins below 0.5% exposure omitted from plot.", ha="center", fontsize=8)
    fig.savefig(HERE / "first_birth_floor_common_states.pdf", bbox_inches="tight")
    fig.savefig(HERE / "first_birth_floor_common_states.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    first0 = float((pre0[..., 0, 0] * rate0[..., 0, 0]).sum())
    first1 = float((pre1[..., 0, 0] * rate1[..., 0, 0]).sum())
    policy = float((pre0[..., 0, 0] * (rate1[..., 0, 0] - rate0[..., 0, 0])).sum())
    composition = float(((pre1[..., 0, 0] - pre0[..., 0, 0]) * rate1[..., 0, 0]).sum())
    receipt = dict(reference="post-interest chain 13, no Estate-A; fixed-price lifecycle diagnostic",
                   baseline=str(BASE), floor_case=str(FLOOR),
                   executed_public_input_diff={"hbar_first_child_jump": [p0["hbar_first_child_jump"], p1["hbar_first_child_jump"]]},
                   floor_multiplier=float(p1["hbar_first_child_jump"] / p0["hbar_first_child_jump"]),
                   rental_cap_physical_rooms=p0["hR_max"], price=float(a0["p_eq"][0]),
                   first_birth_flow_baseline=first0, first_birth_flow_floor=first1,
                   first_birth_flow_policy_at_baseline_exposure=policy,
                   first_birth_flow_composition=composition,
                   flow_identity_error=first1 - first0 - policy - composition,
                   inversion_checks=[check0, check1], earnings_cuts_annual_gross=cuts,
                   source_hashes={str(case): {file: sha(case / file) for file in ("executed_P.json", "solution_arrays.npz")}
                                  for case in (BASE, FLOOR)},
                   interpretation="The floor only changes one nonzero public P field. This is not a price effect or market-clearing equilibrium. First-birth probabilities are conditional on exact common beginning states; physical housing and ending assets require the tenure/transaction branch maps and are not inferred from conditional renter arrays.")
    (HERE / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k: receipt[k] for k in ("executed_public_input_diff", "first_birth_flow_baseline", "first_birth_flow_floor", "first_birth_flow_policy_at_baseline_exposure", "first_birth_flow_composition", "flow_identity_error")}, indent=2))


if __name__ == "__main__":
    main()
