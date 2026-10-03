"""Supplemental native-grid mechanism views from pinned chain-13 saved arrays.

No solve. All policy comparisons hold beginning net financial wealth, inherited
renter tenure, location, age and income state fixed. Family states are conditional
policy schedules, not forced-birth counterfactual outcomes.
"""
from __future__ import annotations

import hashlib
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PACKET = HERE.parents[1]
SOURCE = PACKET / "collection/production_alternative_chain_13/run/native_postcheck/selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz"
PIN = json.loads((PACKET / "fable_analysis/revision1/saved_array_source_receipt.json").read_text())["arrays"]["alternative"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    if SOURCE.resolve() != Path(PIN["path"]).resolve() or sha256(SOURCE) != PIN["sha256"]:
        raise RuntimeError("chain-13 saved-array identity differs from pinned receipt")
    rows = {}
    with np.load(SOURCE, allow_pickle=True) as z:
        b = np.asarray(z["b_grid"])
        p = float(z["p_eq"][0])
        rooms = np.array([2., 4., 6., 8., 10.])
        type_values = np.asarray(z["type_values"])
        j = 1  # ages 22--25
        pi = 1.0 - 0.02 * np.exp(0.134 * (22 - 18))
        floor = np.asarray(z["shared.h_bar"])
        active_parent = np.array([floor[n, m] for n in range(1, 4) for m in range(1, n + 1)])
        if not (np.allclose(floor[0], 0) and np.allclose(active_parent, floor[1, 1])):
            raise RuntimeError("unexpected saved child housing floor")
        for state in (4, 6):
            iz = state - 1
            attempt = np.asarray(z["fert_probs"][:, 0, 0, j, iz, 1])
            post = np.asarray(z["g_beginning_distribution"][:, 0, 0, j, iz, 0, 0])
            pre = post / (1.0 - pi * attempt)
            occupied = pre > 1e-14
            if not occupied.any():
                raise RuntimeError(f"income state {state} has no occupied childless native nodes")
            cumulative = np.cumsum(pre) / pre.sum()
            b999 = float(b[np.searchsorted(cumulative, .999)])
            policies = {}
            for n, m, label in ((0, 0, "childless"), (1, 1, "one_child"), (2, 2, "two_children")):
                prob = np.asarray(z["tenure_probs"][:, 0, 0, j, iz, n, m, :], dtype=float)
                if not np.allclose(prob.sum(axis=1)[occupied], 1.0, atol=1e-5):
                    raise RuntimeError("tenure probabilities fail to sum to one")
                rent = np.asarray(z["hR_pol"][:, 0, 0, j, iz, n, m])
                expected_rooms = prob[:, 0] * rent + prob[:, 1:] @ rooms
                policies[label] = {
                    "renter_rooms_if_renting": rent[occupied].tolist(),
                    "ownership_probability": prob[occupied, 1:].sum(axis=1).tolist(),
                    "expected_rooms_over_tenure_choices": expected_rooms[occupied].tolist(),
                }
            rows[str(state)] = {
                "beginning_b": b[occupied].tolist(),
                "first_birth_attempt": attempt[occupied].tolist(),
                "pre_fertility_childless_mass_conditional_on_income_state": (pre[occupied] / pre.sum()).tolist(),
                "total_pre_fertility_childless_mass": float(pre.sum()),
                "post_fertility_childless_mass": float(post.sum()),
                "occupied_range": [float(b[occupied][0]), float(b[occupied][-1])],
                "conditional_mass_99_9pct_upper_b": b999,
                "policies": policies,
            }
        assert np.isclose(p, 0.7760569760205564, rtol=1e-5)

    meta = {
        "source": str(SOURCE.resolve()), "source_sha256": PIN["sha256"],
        "reference": "working revised-timing soft-constraint chain 13, not certified paper baseline",
        "age_cell": "22--25", "location": 0, "inherited_tenure": "renter",
        "income_states": [4, 6], "price_per_owner_room": p,
        "income_state_multipliers": {str(s): float(type_values[s - 1]) for s in (4, 6)},
        "owner_room_rungs": rooms.tolist(), "rental_room_cap": 6.0,
        "executed_parent_room_floor": float(floor[1, 1]),
        "floor_definition": "0 if m=0; 2.593759507... rooms if m>=1, for feasible n,m",
        "fecundity_probability": pi,
        "mass_definition": "g_beginning_distribution is post-fertility/pre-tenure. At n=m=0, pre-fertility mass = saved post mass / (1 - fecundity * first-birth attempt). Normalize separately within income state.",
        "interpretation": "All housing and tenure comparisons are conditional policies at the same beginning b, age, income and inherited renter tenure; no causal birth experiment.",
        "data": rows,
    }
    (HERE / "native_grid_data.json").write_text(json.dumps(meta, indent=2) + "\n")
    with (HERE / "selected_native_nodes.csv").open("w", newline="") as stream:
        fields = ["income_state", "beginning_b", "first_birth_attempt", "pre_fertility_mass_share",
                  "renter_rooms_childless", "renter_rooms_one_child", "renter_rooms_two_children",
                  "ownership_childless", "ownership_one_child", "ownership_two_children",
                  "expected_rooms_childless", "expected_rooms_one_child", "expected_rooms_two_children"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for state in (4, 6):
            row = rows[str(state)]
            for index, bval in enumerate(row["beginning_b"]):
                if not (np.isclose(bval, 0) or np.isclose(bval, 1)):
                    continue
                out = {"income_state": state, "beginning_b": bval,
                       "first_birth_attempt": row["first_birth_attempt"][index],
                       "pre_fertility_mass_share": row["pre_fertility_childless_mass_conditional_on_income_state"][index]}
                for family in ("childless", "one_child", "two_children"):
                    out[f"renter_rooms_{family}"] = row["policies"][family]["renter_rooms_if_renting"][index]
                    out[f"ownership_{family}"] = row["policies"][family]["ownership_probability"][index]
                    out[f"expected_rooms_{family}"] = row["policies"][family]["expected_rooms_over_tenure_choices"][index]
                writer.writerow(out)

    colors = {"childless": "#2563a6", "one_child": "#e07826", "two_children": "#67823a"}
    labels = {"childless": "no children", "one_child": "one child at home", "two_children": "two children at home"}
    fig, axes = plt.subplots(4, 2, figsize=(12.5, 12), sharex="col")
    for col, state in enumerate((4, 6)):
        row = rows[str(state)]; x = np.asarray(row["beginning_b"])
        keep = x <= 5.5
        axes[0, col].plot(x[keep], np.asarray(row["first_birth_attempt"])[keep], "o-", color="black", ms=3, lw=1)
        for family, data in row["policies"].items():
            for row_ix, field in ((1, "renter_rooms_if_renting"), (2, "ownership_probability"), (3, "expected_rooms_over_tenure_choices")):
                axes[row_ix, col].plot(x[keep], np.asarray(data[field])[keep], "o-", color=colors[family], ms=3, lw=1, label=labels[family])
        for ax in axes[:, col]:
            ax.set_xlim(-.08, 5.55); ax.grid(alpha=.22)
            ax.axvline(row["conditional_mass_99_9pct_upper_b"], ls=":", color="gray", lw=1)
        axes[0, col].set_title(f"Income state {state} · ages 22–25 · inherited renters")
        axes[3, col].set_xlabel("beginning net financial wealth b (model units; native nodes)")
    axes[0, 0].set_ylabel("first-birth attempt probability")
    axes[1, 0].set_ylabel("rental rooms if renting")
    axes[2, 0].set_ylabel("probability of owning")
    axes[3, 0].set_ylabel("expected rooms across tenure choices")
    for ax in axes[1]: ax.axhline(6, color="gray", ls="--", lw=.8)
    for ax in axes[0]: ax.set_ylim(-.03, 1.03)
    for ax in axes[1]: ax.set_ylim(2.0, 6.15)
    for ax in axes[2]: ax.set_ylim(-.03, 1.03)
    for ax in axes[3]: ax.set_ylim(2.0, 8.65)
    axes[1, 0].legend(fontsize=8, loc="upper left")
    fig.suptitle("Conditional housing and tenure schedules at identical beginning wealth\n"
                 "Dashed line: six-room rental cap; dotted vertical: 99.9% childless pre-fertility mass cutoff", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, .955))
    fig.savefig(HERE / "housing_tenure_fertility_native.png", dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), sharex="col", sharey=True)
    for col, state in enumerate((4, 6)):
        row = rows[str(state)]; x = np.asarray(row["beginning_b"]); keep = x <= 5.5
        axes[col].stem(x[keep], np.asarray(row["pre_fertility_childless_mass_conditional_on_income_state"])[keep], basefmt=" ", linefmt="C0-", markerfmt="C0o")
        axes[col].set_yscale("log")
        axes[col].set_xlim(-.08, 5.55); axes[col].grid(alpha=.22)
        axes[col].axvline(row["conditional_mass_99_9pct_upper_b"], ls=":", color="gray", lw=1)
        axes[col].set_title(f"Income state {state} · inherited renters, ages 22–25")
        axes[col].set_xlabel("beginning net financial wealth b (model units; native nodes)")
    axes[0].set_ylabel("pre-fertility childless mass share\nwithin income state (log)")
    fig.suptitle("Occupied native wealth nodes; dotted line marks 99.9% mass cutoff", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, .93))
    fig.savefig(HERE / "childless_mass_native.png", dpi=170)
    plt.close(fig)
    print(json.dumps({s: {"nodes": len(v["beginning_b"]), "range": v["occupied_range"],
                           "b99.9": v["conditional_mass_99_9pct_upper_b"]} for s, v in rows.items()}, indent=2))


if __name__ == "__main__":
    main()
