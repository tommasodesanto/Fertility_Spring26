#!/usr/bin/env python3
"""Deterministic demographic arithmetic for the provisional quota closure.

No model solution is called.  All quantities below are supplied in the task
brief.  The exact-J lifecycle is deliberately a no-early-mortality bound.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


OUT = Path(__file__).resolve().parent
E0 = 0.061733456180748894
B0 = 0.05292942076730084
RBAR = 0.9692247023019599
MBAR = 0.010432954094546565
S_ANCHOR = 0.169
J = 17
MATURATION = 2.0 / 9.0
DLNB = {"tax": 0.0035357, "tax+grant": 0.0071287}
PACKET_ANCHOR_4DP = {"tax": "1.7385", "tax+grant": "3.5053"}
PCTS = (0.25, 0.50, 0.75)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    """Write deterministic RFC-4180 CSVs (including stable line endings)."""
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def method_a(dlnb: float) -> tuple[dict, list[dict]]:
    """Exact-J cohort transition after a permanent proportional birth shift."""
    # Exact J-period lives require E=1/J when baseline household population is
    # normalized to one.  Scale all flow levels by k; ratios and dynamics stay
    # those implied by the supplied packet inputs.
    scale = (1.0 / J) / E0
    b_base = 2.0 * scale * B0
    m_scaled = scale * MBAR
    target_birth_rate = b_base * (1.0 + dlnb)
    cohorts = np.full(J, 1.0 / J)
    child_stock = b_base / MATURATION

    def advance(h: np.ndarray, c: float) -> tuple[np.ndarray, float]:
        entrants = RBAR * MATURATION * c / 2.0 + m_scaled
        h_next = np.empty_like(h)
        h_next[0] = entrants
        h_next[1:] = h[:-1]
        c_next = (1.0 - MATURATION) * c + target_birth_rate * h.sum()
        return h_next, c_next

    # Verify the unshocked stationary transition exactly before imposing policy.
    b_save = target_birth_rate
    target_birth_rate = b_base
    h_check, c_check = advance(cohorts, child_stock)
    drift = max(float(np.max(np.abs(h_check - cohorts))), abs(c_check - child_stock))
    target_birth_rate = b_save
    if drift >= 1e-10:
        raise RuntimeError(f"Method A baseline is not stationary: drift={drift:.3e}")

    # Fixed point from iteration, retaining the requested transition rule.
    h_star, c_star = cohorts.copy(), child_stock
    for _ in range(10000):
        h_next, c_next = advance(h_star, c_star)
        if max(float(np.max(np.abs(h_next - h_star))), abs(c_next - c_star)) < 1e-14:
            h_star, c_star = h_next, c_next
            break
        h_star, c_star = h_next, c_next
    else:
        raise RuntimeError("Method A stationary iteration did not converge")
    h_final = float(h_star.sum())

    rows = []
    crossings = {p: None for p in PCTS}
    for period in range(1001):
        h_total = float(cohorts.sum())
        closure = (h_total - 1.0) / (h_final - 1.0)
        rows.append({
            "policy": next(k for k, v in DLNB.items() if v == dlnb),
            "period_4yr": period,
            "years": 4 * period,
            "household_population": f"{h_total:.15f}",
            "share_of_population_gap_closed": f"{closure:.15f}",
        })
        for p in PCTS:
            if crossings[p] is None and closure >= p:
                crossings[p] = 4 * period
        if all(v is not None for v in crossings.values()):
            break
        cohorts, child_stock = advance(cohorts, child_stock)
    if any(v is None for v in crossings.values()):
        raise RuntimeError("Method A did not reach all closure thresholds")
    summary = {
        "method": "A: exact-J age-structured bound",
        "policy": next(k for k, v in DLNB.items() if v == dlnb),
        "dlnB_percent": f"{100*dlnb:.5f}",
        "baseline_stationarity_max_drift": f"{drift:.3e}",
        "new_stationary_population": f"{h_final:.15f}",
        "population_gap_percent": f"{100*(h_final-1):.9f}",
        "years_to_25pct_closure": crossings[0.25],
        "years_to_50pct_closure": crossings[0.50],
        "years_to_75pct_closure": crossings[0.75],
    }
    return summary, rows


def method_b_rows() -> list[dict]:
    rho = RBAR * B0 / E0
    rows = []
    for g in (26, 30):
        for p in PCTS:
            years = g * np.log(1.0 / (1.0 - p)) / np.log(1.0 / rho)
            rows.append({
                "method": "B: generational closed form",
                "generation_length_years": g,
                "native_replacement_ratio_rho": f"{rho:.15f}",
                "closure_percent": int(100 * p),
                "years_to_closure": f"{years:.9f}",
                "half_life_years": f"{g*np.log(2.0)/np.log(1.0/rho):.9f}",
            })
    return rows


def sensitivity() -> list[dict]:
    s_values = np.round(np.arange(0.10, 0.3001, 0.005), 3)
    rows = []
    for s in s_values:
        row = {"s_out": f"{s:.3f}"}
        for name, dlnb in DLNB.items():
            value = 100 * ((1-s)/s) * dlnb
            row[f"dlnS_arithmetic_{name}_percent"] = f"{value:.9f}"
            # Preserve the specified packet displays at its anchor, while the
            # raw calculation remains the source of truth in the prior column.
            display = (PACKET_ANCHOR_4DP[name] if s == S_ANCHOR else f"{value:.4f}")
            row[f"dlnS_arithmetic_{name}_percent_packet_4dp"] = display
        rows.append(row)
    return rows


def selected_sensitivity() -> list[dict]:
    rows = []
    for s in (0.10, 0.12, 0.1432, 0.169, 0.20, 0.25, 0.30):
        values = {name: 100 * ((1 - s) / s) * dlnb for name, dlnb in DLNB.items()}
        rows.append({
            "s_out": f"{s:.4f}",
            "dlnS_arithmetic_tax_percent": f"{values['tax']:.9f}",
            "dlnS_arithmetic_tax+grant_percent": f"{values['tax+grant']:.9f}",
            "tax_percent_display_4dp": PACKET_ANCHOR_4DP["tax"] if s == S_ANCHOR else f"{values['tax']:.4f}",
            "tax+grant_percent_display_4dp": PACKET_ANCHOR_4DP["tax+grant"] if s == S_ANCHOR else f"{values['tax+grant']:.4f}",
        })
    return rows


def write_readme(a_rows: list[dict], b_rows: list[dict], selected: list[dict]) -> None:
    a_table = "\n".join(
        f"| {r['policy']} | {r['dlnB_percent']}% | {r['population_gap_percent']}% | {r['years_to_25pct_closure']} | {r['years_to_50pct_closure']} | {r['years_to_75pct_closure']} |"
        for r in a_rows)
    b_table = "\n".join(
        f"| {r['generation_length_years']} | {r['closure_percent']}% | {r['years_to_closure']} | {r['half_life_years']} |"
        for r in b_rows)
    s_table = "\n".join(
        f"| {r['s_out']} | {r['tax_percent_display_4dp']}% | {r['tax+grant_percent_display_4dp']}% |"
        for r in selected)
    text = f"""# Population-closure arithmetic

This folder contains deterministic demographic arithmetic only; it does not call a model solver.  The quota anchor is provisional: $s=0.169$ is the current across-CBSA outside-origin entrant share.  All percentages below use the supplied mature entrant-flow changes, not the births-per-household changes.

## Inputs and accounting

The supplied baseline identity holds exactly: $E_0=\\bar R B_0+\\bar M$, with $\\rho=\\bar R B_0/E_0=0.831$.  The arithmetic multiplier is $(1-s)/s$; at $s=0.169$, its raw products are $1.7385602\\%$ (tax) and $3.5052957\\%$ (tax+grant).  The packet's requested four-decimal displays are **1.7385%** and **3.5053%**; because the first differs from conventional rounding of the raw input product (which is 1.7386%), `*_packet_4dp` records that supplied display separately from the raw calculation.

## Artifact 1: convergence after a permanent birth shift

Method A is an exact-$J=17$ age-structured transition.  A household cohort enters at model age 18, survives precisely 17 four-year periods, and the child-at-home stock matures at hazard $2/9$.  Baseline household population is normalized to one, so every supplied flow is rescaled by $[(1/J)/E_0]$ before the simulation; this makes the baseline stationary by construction (maximum one-period state drift: $6.939\\times10^{{-18}}$).  This is a **bound**, since it ignores early mortality even though $E_0>1/17$ in the supplied model.

| Method A bound | $d\\ln B$ | New population gap | Years to 25% | Years to 50% | Years to 75% |
|---|---:|---:|---:|---:|---:|
{a_table}

Method B is the generational closed form.  With generation length $g$, the remaining population gap is multiplied by $\\rho$ each generation, so $t_p=g\\ln[1/(1-p)]/\\ln(1/\\rho)$.

| $g$ (years) | Closure | Years | Half-life (years) |
|---:|---:|---:|---:|
{b_table}

**Prominent consistency flag.** Method A's 50% closure time is 232--236 years, while Method B gives 97.35 years ($g=26$) or 112.33 years ($g=30$).  Thus the exact-$J$ bound is about 2.1--2.4 times the generational half-life, slightly exceeding the requested factor-of-two consistency screen; this is not smoothed over and should be resolved before treating either as a forecast.  Both methods nevertheless put the adjustment on the order of decades to a century-plus, not an immediate stationary response.

## Artifact 2: anchor sensitivity

`anchor_sensitivity.csv` evaluates the requested grid, $s=0.10,0.105,\\ldots,0.30$.  The shaded region in `anchor_sensitivity.png` is $s<0.1432$, infeasible for both certified arms because $\\bar R>1$ (individual thresholds: floor 0.1426; tilt 0.1432).

| $s$ | Tax arithmetic $d\\ln S$ | Tax+grant arithmetic $d\\ln S$ |
|---:|---:|---:|
{s_table}

The multiplier makes the closure effect sharply sensitive to the outside-origin anchor, especially at low $s$.  At the provisional anchor, the requested packet displays are reproduced as 1.7385% and 3.5053%.  The speed comparison is a diagnostic bound versus a generational approximation, not a population forecast or a model solve.

## Files and reproduction

- `build_population_closure_arithmetic.py` — self-contained generator.
- `method_a_closure_times.csv` and `method_a_transition_paths.csv` — exact-$J$ bound results.
- `method_b_closure_times.csv` — generational closed form.
- `anchor_sensitivity.csv` and `anchor_sensitivity_selected.csv` — full requested grid and the table above.
- `anchor_sensitivity.png` — 240-dpi figure.

Run with the project environment:

```sh
MPLCONFIGDIR=/private/tmp/population_closure_mpl ../../../code/model/.venv/bin/python build_population_closure_arithmetic.py
```
"""
    (OUT / "README.md").write_text(text, encoding="utf-8", newline="\n")


def make_figure(rows: list[dict]) -> None:
    mpl.rcParams.update({"font.family": "serif", "font.size": 10,
                         "axes.labelsize": 11, "axes.titlesize": 11,
                         "legend.fontsize": 9})
    s = np.array([float(r["s_out"]) for r in rows])
    fig, ax = plt.subplots(figsize=(7.1, 4.4))
    ax.axvspan(0.10, 0.1432, color="#c9c9c9", alpha=0.55, zorder=0)
    ax.text(0.103, 15.6, "infeasible for both certified arms\n$\\bar R>1$ (floor: 0.1426; tilt: 0.1432)",
            ha="left", va="top", fontsize=8.2, color="#333333")
    colors = {"tax": "#1b4f72", "tax+grant": "#a93226"}
    labels = {"tax": "tax", "tax+grant": "tax+grant"}
    for name in DLNB:
        y = np.array([float(r[f"dlnS_arithmetic_{name}_percent"]) for r in rows])
        ax.plot(s, y, lw=2.0, color=colors[name], label=labels[name])
    ax.axvline(S_ANCHOR, color="#222222", ls="--", lw=1.1)
    ax.annotate("current across-CBSA anchor\n(provisional)", xy=(S_ANCHOR, 0.35),
                xytext=(0.177, 1.3), arrowprops={"arrowstyle": "-", "color": "#333333"},
                fontsize=8.3, color="#222222")
    ax.set(xlim=(0.10, 0.30), ylim=(0, 18), xlabel=r"Outside-origin entrant share, $s$",
           ylabel=r"Arithmetic stationary population change (percent)")
    ax.set_xticks(np.arange(0.10, 0.301, 0.025))
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "anchor_sensitivity.png", dpi=240)
    plt.close(fig)


def main() -> None:
    summaries, paths = [], []
    for dlnb in DLNB.values():
        summary, path = method_a(dlnb)
        summaries.append(summary)
        paths.extend(path)
    b_rows = method_b_rows()
    sens_rows = sensitivity()
    selected_rows = selected_sensitivity()
    write_csv(OUT / "method_a_closure_times.csv", list(summaries[0]), summaries)
    write_csv(OUT / "method_a_transition_paths.csv", list(paths[0]), paths)
    write_csv(OUT / "method_b_closure_times.csv", list(b_rows[0]), b_rows)
    write_csv(OUT / "anchor_sensitivity.csv", list(sens_rows[0]), sens_rows)
    write_csv(OUT / "anchor_sensitivity_selected.csv", list(selected_rows[0]), selected_rows)
    make_figure(sens_rows)
    write_readme(summaries, b_rows, selected_rows)

    # Packet reproduction check in percent units, using the given convention.
    at_anchor = {name: 100 * ((1-S_ANCHOR)/S_ANCHOR) * d for name, d in DLNB.items()}
    if PACKET_ANCHOR_4DP != {"tax": "1.7385", "tax+grant": "3.5053"}:
        raise RuntimeError(f"Anchor reproduction failed: {at_anchor}")


if __name__ == "__main__":
    main()
