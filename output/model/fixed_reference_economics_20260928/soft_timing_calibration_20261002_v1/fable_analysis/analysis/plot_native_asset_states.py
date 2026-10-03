"""Plot saved native-node policies and corresponding occupied wealth mass; no solve.

At ages 22--25, condition on inherited renter tenure, location 0, children
ever born n=0 and children at home m=0, and income states 4 and 6. Fertility
is a pre-tenure first-birth attempt probability. Consumption and next-period
financial wealth are conditional on choosing the renter branch. Reconstruct
pre-fertility state mass from saved post-fertility, pre-tenure childless mass
using the executed first-birth selection probability; normalize it within
each arm and income state. Plot only native b nodes with mass > 1e-14.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
PACKET = BASE.parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
RECEIPT = json.loads((BASE / "revision1/saved_array_source_receipt.json").read_text())
ARMS = {
    "original": PACKET / "collection/production_original_chain_15/run/native_postcheck/selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz",
    "alternative": PACKET / "collection/production_alternative_chain_13/run/native_postcheck/selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz",
}

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def main() -> None:
    age_index = 1  # model cell 22--25
    fecundity = 1.0 - 0.02 * np.exp(0.134 * (22 - 18))
    data = {}
    for arm, path in ARMS.items():
        expected = RECEIPT["arrays"][arm]
        if path.resolve() != Path(expected["path"]).resolve() or sha256(path) != expected["sha256"]:
            raise RuntimeError(f"saved-array identity mismatch for {arm}")
        with np.load(path, allow_pickle=True) as z:
            b = z["b_grid"]
            arm_data = {}
            for income_state in (4, 6):
                iz = income_state - 1
                fert = z["fert_probs"][:, 0, 0, age_index, iz, 1]
                post_mass = z["g_beginning_distribution"][:, 0, 0, age_index, iz, 0, 0]
                pre_mass = post_mass / np.clip(1.0 - fecundity * fert, 1e-9, None)
                keep = pre_mass > 1e-14
                mass_conditional = pre_mass / pre_mass.sum()
                arm_data[str(income_state)] = {
                    "b": b[keep].tolist(),
                    "first_birth_attempt": fert[keep].tolist(),
                    "renter_consumption": z["c_pol"][:, 0, 0, age_index, iz, 0, 0][keep].tolist(),
                    "renter_next_financial_wealth": z["bp_pol"][:, 0, 0, age_index, iz, 0, 0][keep].tolist(),
                    "pre_fertility_mass_conditional_income": mass_conditional[keep].tolist(),
                    "total_pre_fertility_mass": float(pre_mass.sum()),
                    "positive_mass_cutoff": 1e-14,
                }
            data[arm] = arm_data
    meta = {
        "conditioning": "ages 22--25; location 0; inherited renter; n=0 children ever born; m=0 at home; income state 4 or 6",
        "fertility_policy": "pre-tenure first-birth attempt probability",
        "consumption_saving_policy": "conditional on choosing the renter branch; b prime is next-period net financial wealth",
        "mass": "reconstructed pre-fertility childless inherited-renter mass at native b nodes, divided by total mass within that arm/income state",
        "source_receipt": str((BASE / "revision1/saved_array_source_receipt.json").resolve()),
        "data": data,
    }
    (OUT / "native_asset_policy_mass.json").write_text(json.dumps(meta, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8), sharex=True)
    fields = ["first_birth_attempt", "renter_consumption", "renter_next_financial_wealth", "pre_fertility_mass_conditional_income"]
    titles = ["First-birth attempt, pre-tenure", "Consumption if renting", "Next financial wealth if renting", "Conditional wealth-node mass"]
    colors = {"original": "C0", "alternative": "C1"}
    markers = {4: "o", 6: "s"}
    max_b = 0.0
    min_b = 0.0
    for arm, arm_data in data.items():
        for state in (4, 6):
            row = arm_data[str(state)]
            x = np.asarray(row["b"])
            max_b = max(max_b, float(x.max()))
            min_b = min(min_b, float(x.min()))
            for ax, field in zip(axes.flat, fields):
                ax.scatter(x, row[field], s=20, color=colors[arm], marker=markers[state],
                           label=f"{arm}, income {state}", alpha=.8)
    for ax, title in zip(axes.flat, titles):
        ax.set_title(title)
        ax.grid(alpha=.2)
        ax.set_xlim(min_b - .25, max_b + .25)
    for ax in axes[1]:
        ax.set_xlabel("beginning liquid financial wealth b (native nodes)")
    axes[0, 0].set_ylabel("probability")
    axes[0, 1].set_ylabel("consumption, model units")
    axes[1, 0].set_ylabel("b prime, model units")
    axes[1, 1].set_ylabel("fraction of state mass (log scale)")
    axes[1, 1].set_yscale("log")
    axes[0, 0].legend(fontsize=8, ncol=2, loc="best")
    fig.suptitle("Ages 22–25, childless inherited renters, income states 4 and 6\n"
                 "Native occupied asset nodes only; mass conditional on arm and income state", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "F5_native_asset_policy_mass.png", dpi=160)
    plt.close(fig)
    print(json.dumps({arm: {state: {"nodes": len(row["b"]), "state_mass": row["total_pre_fertility_mass"]}
                            for state, row in arm_data.items()} for arm, arm_data in data.items()}, indent=2))

if __name__ == "__main__":
    main()
