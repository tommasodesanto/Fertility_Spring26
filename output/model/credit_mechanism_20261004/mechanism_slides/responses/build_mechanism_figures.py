"""Draw talk-size mechanism figures from the retained common-state CSVs.

No model solve is performed. Run from the repository root with code/model/.venv/bin/python.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[5]
DIAGNOSTICS = ROOT / "output/model/credit_mechanism_20261004/diagnostics"
OUT = Path(__file__).resolve().parent
BLUE = "#1f4e79"
ORANGE = "#bf6b2b"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def source_checks() -> dict[str, object]:
    price_receipt = json.loads((DIAGNOSTICS / "price110_decomposition/receipt.json").read_text())
    credit_receipt = json.loads((DIAGNOSTICS / "phi095_decomposition/receipt.json").read_text())
    assert price_receipt["kind"] == "price110", "Price CSV is not the matched price case"
    assert credit_receipt["pair_changes"] == ["phi"], "Credit CSV is not the phi-only case"
    assert abs(price_receipt["price1"] / price_receipt["price0"] - 1.10) < 1e-12
    assert abs(credit_receipt["price"] - price_receipt["price0"]) < 1e-12
    assert price_receipt["sources"][0]["arrays_sha256"] == credit_receipt["sources"][0]["arrays_sha256"]
    assert price_receipt["sources"][0]["parameters_sha256"] == credit_receipt["sources"][0]["parameters_sha256"]
    return {
        "baseline_arrays_sha256": price_receipt["sources"][0]["arrays_sha256"],
        "price_arrays_sha256": price_receipt["sources"][1]["arrays_sha256"],
        "credit_arrays_sha256": credit_receipt["sources"][1]["arrays_sha256"],
        "price0": price_receipt["price0"],
        "price1": price_receipt["price1"],
        "wealth_group0": credit_receipt["wealth_groups"][0],
        "income_cuts": credit_receipt["income_cuts"],
    }


def price_rows() -> tuple[list[dict[str, object]], dict[str, float]]:
    source = read_csv(DIAGNOSTICS / "price110_decomposition/common_state_groups.csv")
    grouped = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])
    for row in source:
        if int(row["children_ever_born"]) != 0 or int(row["children_at_home"]) != 0:
            continue
        age, tenure = int(row["age"]), int(row["beginning_tenure"]) != 0
        exposure = float(row["exposure0"])
        cell = grouped[age, tenure]
        cell[0] += exposure
        cell[1] += exposure * float(row["birth_rate0"])
        cell[2] += exposure * float(row["birth_rate1_common"])
        cell[3] += float(row["policy"])
    rows = []
    for (age, tenure), (exposure, birth0, birth1, policy) in sorted(grouped.items()):
        assert exposure > 0
        assert abs(policy - (birth1 - birth0)) < 1e-11
        rows.append({
            "age_start": age,
            "age_end": age + 3,
            "beginning_tenure": "owner" if tenure else "renter",
            "baseline_prebirth_exposure": exposure,
            "baseline_first_birth_probability_pct": 100 * birth0 / exposure,
            "price110_first_birth_probability_at_baseline_states_pct": 100 * birth1 / exposure,
            "response_pp": 100 * (birth1 - birth0) / exposure,
        })
    source_policy = sum(float(r["policy"]) for r in source if int(r["children_ever_born"]) == 0 and int(r["children_at_home"]) == 0)
    plotted_policy = sum(float(r["baseline_prebirth_exposure"]) * float(r["response_pp"]) / 100 for r in rows)
    assert abs(source_policy - plotted_policy) < 1e-11
    return rows, {"first_birth_policy_flow_change": plotted_policy}


def credit_rows() -> list[dict[str, object]]:
    source = read_csv(DIAGNOSTICS / "phi095_decomposition/first_birth_income_wealth_summary.csv")
    result = []
    for row in source:
        if int(row["wealth_group"]) != 0:
            continue
        wait = float(row["wait_credit_gain"])
        success = float(row["success_credit_gain"])
        gap = float(row["success_gap_delta"])
        assert abs((success - wait) - gap) < 1e-10
        result.append({
            "income_group": ("Low", "Middle", "High")[int(row["income_group"])],
            "beginning_assets": "b<=0",
            "baseline_prebirth_exposure": float(row["exposure0"]),
            "baseline_first_birth_probability_pct": 100 * float(row["birth_rate0"]),
            "logit_supported_exposure_share": float(row["gap_supported_share"]),
            "wait_gain_lifetime_utility": wait,
            "successful_first_birth_gain_lifetime_utility": success,
            "success_minus_wait_gain_lifetime_utility": gap,
        })
    assert [r["income_group"] for r in result] == ["Low", "Middle", "High"]
    return result


def cap_rows(price0: float) -> list[dict[str, object]]:
    cap_dir = ROOT / "output/model/credit_mechanism_20261004/price_decomposition/lead_followup/cap"
    source = {row["case"]: row for row in read_csv(cap_dir / "results.csv")}
    result = []
    baseline = source["cached_baseline"]
    for name, cap in [("cap4", 4), ("cached_baseline", 6), ("cap20", 20)]:
        row = source[name]
        assert abs(float(row["price"]) - price0) < 1e-12
        if name != "cached_baseline":
            status = json.loads((cap_dir / "cases" / name / "status.json").read_text())
            audit = json.loads((cap_dir / "cases" / name / "arm_spec_and_audit.json").read_text())
            assert status["status"] == "completed_fixed_price_diagnostic"
            assert audit["spec"]["hR_max"] == cap
            assert list(audit["audit"]["changed_fields"]) == ["hR_max"]
            assert not audit["audit"]["problems"]
        births_change = 100 * (float(row["explicit_births"]) / float(baseline["explicit_births"]) - 1)
        owner_change = 100 * (float(row["ownership_18_29"]) - float(baseline["ownership_18_29"]))
        assert abs(births_change - float(row["pct_change_vs_baseline"])) < 1e-10
        assert abs(owner_change - float(row["own_pp_change"])) < 1e-10
        result.append({
            "renter_room_cap": cap,
            "explicit_births": float(row["explicit_births"]),
            "explicit_births_pct_change_vs_cap6": births_change,
            "ownership_18_29": float(row["ownership_18_29"]),
            "ownership_18_29_pp_change_vs_cap6": owner_change,
        })
    return result


def style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 14,
        "axes.labelsize": 15,
        "axes.titlesize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 13,
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white")
    fig.savefig(OUT / f"{stem}.png", dpi=220, facecolor="white")
    plt.close(fig)


def plot_price(rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, tenure, label in zip(axes, ["renter", "owner"], ["Beginning renters", "Beginning owners"]):
        cell = [r for r in rows if r["beginning_tenure"] == tenure]
        x = np.array([int(r["age_start"]) for r in cell])
        base = [float(r["baseline_first_birth_probability_pct"]) for r in cell]
        shock = [float(r["price110_first_birth_probability_at_baseline_states_pct"]) for r in cell]
        ax.plot(x, base, color=BLUE, marker="o", lw=2.7, ms=6, label="Baseline")
        ax.plot(x, shock, color=ORANGE, marker="s", lw=2.7, ms=6, label="Price +10%")
        ax.set_title(label, loc="left", pad=12)
        ax.set_xticks([18, 22, 26, 30, 34, 38, 42])
        ax.set_xticklabels(["18–21", "22–25", "26–29", "30–33", "34–37", "38–41", "42–45"], rotation=35, ha="right")
        ax.grid(axis="y", color="#d9dfe5", lw=0.8)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("First birth probability (%)")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.5, 0.99), frameon=False, ncol=2)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.22, top=0.80, wspace=0.16)
    save(fig, "first_birth_price110_by_beginning_tenure")


def plot_credit(rows: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(11, 4))
    x = np.arange(len(rows))
    width = 0.30
    ax.bar(x - width / 2, [1000 * float(r["wait_gain_lifetime_utility"]) for r in rows], width, color=BLUE, label="Wait")
    ax.bar(x + width / 2, [1000 * float(r["successful_first_birth_gain_lifetime_utility"]) for r in rows], width, color=ORANGE, label="Successful first birth")
    ax.set_xticks(x, [str(r["income_group"]) for r in rows])
    ax.set_xlabel("Current earnings group, beginning assets b ≤ 0")
    ax.set_ylabel("Lifetime utility gain (×10⁻³)")
    ax.set_ylim(0, 22)
    ax.grid(axis="y", color="#d9dfe5", lw=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper right")
    fig.subplots_adjust(left=0.13, right=0.99, bottom=0.21, top=0.96)
    save(fig, "credit_phi095_wait_success_gains_debt_states")


def plot_cap(rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    x = np.arange(3)
    labels = [str(row["renter_room_cap"]) for row in rows]
    colors = [ORANGE, BLUE, ORANGE]
    series = [
        ("explicit_births_pct_change_vs_cap6", "Explicit births", "Change from cap 6 (%)", (-3.2, 1.4)),
        ("ownership_18_29_pp_change_vs_cap6", "Ownership, ages 18–29", "Change from cap 6 (pp)", (-14, 19)),
    ]
    for ax, (key, title, ylabel, ylim) in zip(axes, series):
        values = [float(row[key]) for row in rows]
        ax.bar(x, values, width=0.56, color=colors)
        ax.plot(1, 0, marker="o", color=BLUE, ms=8, zorder=3)
        ax.set_xticks(x, labels)
        ax.set_xlabel("Maximum renter rooms")
        ax.set_title(title, loc="left", pad=10)
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        ax.axhline(0, color="#424b54", lw=1)
        ax.grid(axis="y", color="#d9dfe5", lw=0.8)
        ax.set_axisbelow(True)
        for i, value in enumerate(values):
            offset = 0.12 if key.startswith("explicit") else 0.7
            display = "0" if i == 1 else (f"{value:+.2f}" if key.startswith("explicit") else f"{value:+.1f}")
            ax.text(i, value + (offset if value >= 0 else -offset), display,
                    ha="center", va="bottom" if value >= 0 else "top", fontsize=12)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.22, top=0.86, wspace=0.32)
    save(fig, "renter_room_cap_births_ownership")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checks = source_checks()
    price, reconciliation = price_rows()
    credit = credit_rows()
    cap = cap_rows(checks["price0"])
    write_csv(OUT / "price110_first_birth_by_age_tenure.csv", price)
    write_csv(OUT / "phi095_branch_value_gains_b_nonpositive.csv", credit)
    write_csv(OUT / "renter_room_cap_responses.csv", cap)
    style()
    plot_price(price)
    plot_credit(credit)
    plot_cap(cap)
    (OUT / "plot_receipt.json").write_text(json.dumps({**checks, **reconciliation, "rows_price": len(price), "rows_credit": len(credit), "rows_cap": len(cap)}, indent=2) + "\n")


if __name__ == "__main__":
    main()
