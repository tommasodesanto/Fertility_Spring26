"""Matched native-node fixed-price policy overlay after both solutions pass."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ARMS = {"phi_080": ("80% financed", "#2563a6"),
        "phi_100": ("100% financed", "#d27023")}
ROOMS = np.array([2., 4., 6., 8., 10.])


def main():
    if not (HERE / "completed.json").is_file():
        raise RuntimeError("two-solve completion receipt absent")
    arrays = {}
    for arm in ARMS:
        with np.load(HERE / arm / "solution_arrays.npz", allow_pickle=False) as saved:
            arrays[arm] = {key: saved[key].copy() for key in
                           ("b_grid", "fert_probs", "fert2_probs", "tenure_probs", "hR_pol",
                            "g_beginning_distribution", "g", "shared.h_bar", "shared.phi_choice")}
    b = arrays["phi_080"]["b_grid"]
    if not np.array_equal(b, arrays["phi_100"]["b_grid"]):
        raise RuntimeError("asset grid differs across financing arms")
    pi = 1.0 - 0.02 * np.exp(0.134 * (22 - 18))
    results = {"design": "one fixed-price lifecycle solve per arm at p=0.7760569760205563; only uniform phi changes 0.8 to 1.0",
               "conditioning": "age22-25, location0, inherited renter, income state4/6; conditional policy schedules at identical beginning b",
               "weight": "baseline phi=.8 pre-fertility childless inherited-renter native-node mass, normalized separately by income state; identical weights across both arms",
               "population": "fixed-price credit counterfactual; each arm has its own solved stationary distribution. Conditional family-state schedules are not forced-birth transitions; no GE or dated transition",
               "income_states": {}}
    rung_decomposition = {}
    fig, axes = plt.subplots(4, 2, figsize=(12.8, 12), sharex="col")
    for col, state in enumerate((4, 6)):
        iz = state - 1; j = 1
        reference = arrays["phi_080"]
        attempt_ref = reference["fert_probs"][:, 0, 0, j, iz, 1]
        post_mass = reference["g_beginning_distribution"][:, 0, 0, j, iz, 0, 0]
        pre_mass = post_mass / (1.0 - pi * attempt_ref)
        occupied = pre_mass > 1e-14
        weights = pre_mass[occupied] / pre_mass[occupied].sum()
        x = b[occupied]
        keep = x <= 5.5
        state_data = {"occupied_node_range": [float(x[0]), float(x[-1])],
                      "baseline_prefertility_childless_mass": float(pre_mass.sum()),
                      "arm_means_fixed_baseline_weights": {}}
        rung_decomposition[str(state)] = {"childless": {}, "one_child": {}}
        for arm, (label, color) in ARMS.items():
            z = arrays[arm]
            if not np.allclose(z["shared.phi_choice"], .8 if arm == "phi_080" else 1.0, atol=0, rtol=0):
                raise RuntimeError("financed share did not reach policy solver")
            attempt = z["fert_probs"][:, 0, 0, j, iz, 1][occupied]
            axes[0, col].plot(x[keep], attempt[keep], "o-", color=color, lw=1.2, ms=3, label=label)
            summaries = {"first_birth_attempt": float(weights @ attempt),
                         "solved_g_mass": float(z["g"].sum())}
            for n, m, fam in ((0, 0, "childless"), (1, 1, "one_child")):
                tp = z["tenure_probs"][:, 0, 0, j, iz, n, m, :][occupied].astype(float)
                if not np.allclose(tp.sum(1), 1, atol=1e-5):
                    raise RuntimeError("tenure probabilities fail to sum to one")
                rent = z["hR_pol"][:, 0, 0, j, iz, n, m][occupied]
                own = tp[:, 1:].sum(1)
                expected = tp[:, 0] * rent + tp[:, 1:] @ ROOMS
                style = "-" if fam == "childless" else "--"
                text_label = f"{label}, {fam.replace('_', ' ')}"
                for ax, y in ((axes[1, col], own), (axes[2, col], rent), (axes[3, col], expected)):
                    ax.plot(x[keep], y[keep], linestyle=style, marker="o", color=color,
                            lw=1.15, ms=2.8, label=text_label)
                summaries[f"{fam}_ownership"] = float(weights @ own)
                summaries[f"{fam}_rental_rooms"] = float(weights @ rent)
                summaries[f"{fam}_expected_rooms"] = float(weights @ expected)
                rung_decomposition[str(state)][fam][arm] = {
                    str(room): float(weights @ tp[:, k])
                    for k, room in enumerate((0, 2, 4, 6, 8, 10))
                }
            state_data["arm_means_fixed_baseline_weights"][arm] = summaries
        results["income_states"][str(state)] = state_data
        for ax in axes[:, col]:
            ax.set_xlim(-.08, 5.55); ax.grid(alpha=.2)
        axes[0, col].set_title(f"Income state {state}; inherited renters, ages 22–25")
        axes[3, col].set_xlabel("beginning net financial wealth b (model units; native nodes)")
    for ax in axes[0]: ax.set_ylim(-.03, 1.03)
    for ax in axes[1]: ax.set_ylim(-.03, 1.03)
    for ax in axes[2]: ax.set_ylim(2.0, 6.15)
    for ax in axes[3]: ax.set_ylim(2.0, 8.65)
    axes[0, 0].set_ylabel("first-birth attempt probability")
    axes[1, 0].set_ylabel("probability of owning")
    axes[2, 0].set_ylabel("rental rooms if renting")
    axes[3, 0].set_ylabel("expected rooms across tenure choices")
    axes[1, 0].legend(fontsize=7, loc="upper left", ncol=2)
    axes[0, 1].legend(fontsize=8, loc="lower right")
    fig.suptitle("Fixed-price credit relaxation at the same native beginning-wealth nodes\n"
                 "Solid: childless; dashed: one child at home. Policy schedules are conditional, not transitions.", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, .955))
    fig.savefig(HERE / "matched_native_policy_overlay.png", dpi=170)
    plt.close(fig)
    (HERE / "overlay_data.json").write_text(json.dumps(results, indent=2) + "\n")
    (HERE / "owner_room_decomposition.json").write_text(json.dumps(rung_decomposition, indent=2) + "\n")
    print(json.dumps(results["income_states"], indent=2))


if __name__ == "__main__":
    main()
