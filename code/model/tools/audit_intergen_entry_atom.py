#!/usr/bin/env python3
"""Audit entrant-wealth atoms and rung-induced total-wealth jaggedness.

This is a diagnostic script, not a calibration routine. It reads one or more
trusted intergen solution-cache pickles and writes compact plots/tables focused
on the distributional pathologies discussed in the housing-block ledger.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pickle
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_housing_fertility.solver import entry_wealth_grid_weights  # noqa: E402


DEFAULT_BASELINE_CACHE = ROOT / "output/model/intergen_current_review/solution_cache.pkl"
DEFAULT_SPREAD_CACHE = ROOT / "output/model/intergen_current_review/entry_spread5_quick/solution_cache.pkl"
DEFAULT_OUTDIR = ROOT / "output/model/intergen_current_review/atom_audit"
DEFAULT_ENTRY_DATA = (
    ROOT
    / "code/data/psid_followup_mar2026/output/entry_wealth_v1/entry_wealth_candidate_targets_v1.csv"
)


@dataclass
class Case:
    label: str
    cache_path: Path
    sol: Any
    P: Any
    p_eq: np.ndarray


def main() -> None:
    args = parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    cases = [load_case(label, path) for label, path in parse_case_specs(args.case)]
    write_atom_summary(cases, outdir / "atom_origin_summary.md", data_path=args.entry_data)
    plot_atom_heatmaps(cases, outdir)
    plot_near_entry_lifecycle(cases, outdir / "near_entry_mass_by_age.png")
    plot_total_wealth_by_rung(cases, outdir)
    write_entry_data_anchor(cases, args.entry_data, outdir / "entry_wealth_data_anchor.md")
    write_readme(cases, outdir, args.entry_data)
    print(f"Wrote entry-atom audit to {outdir}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        default=None,
        help="Case spec LABEL:PATH. May be repeated. Defaults to point-entry and entry-spread-5 caches.",
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--entry-data", type=Path, default=DEFAULT_ENTRY_DATA)
    return parser.parse_args()


def parse_case_specs(specs: list[str] | None) -> list[tuple[str, Path]]:
    if not specs:
        return [("point_entry", DEFAULT_BASELINE_CACHE), ("entry_spread5", DEFAULT_SPREAD_CACHE)]
    out: list[tuple[str, Path]] = []
    for spec in specs:
        if ":" not in spec:
            raise ValueError(f"--case must be LABEL:PATH, got {spec!r}")
        label, raw_path = spec.split(":", 1)
        out.append((label.strip(), Path(raw_path).expanduser()))
    return out


def load_case(label: str, cache_path: Path) -> Case:
    if not cache_path.exists():
        raise FileNotFoundError(f"solution cache not found for {label}: {cache_path}")
    with cache_path.open("rb") as fh:
        payload = pickle.load(fh)
    baseline = payload.get("baseline")
    if not isinstance(baseline, dict):
        raise ValueError(f"cache has no baseline payload: {cache_path}")
    sol = baseline["sol"]
    P = baseline["P"]
    p_eq = np.asarray(baseline.get("p_eq", getattr(sol, "p_eq", np.nan)), dtype=float).reshape(-1)
    return Case(label=label, cache_path=cache_path, sol=sol, P=P, p_eq=p_eq)


def write_atom_summary(cases: list[Case], path: Path, *, data_path: Path) -> None:
    lines = [
        "# Entry Atom Audit Summary",
        "",
        "This diagnostic evaluates the distributional pathology directly. It does",
        "not evaluate whether the SMM loss improves.",
        "",
        f"Existing PSID entry-wealth summary: `{relative(data_path)}`",
        "",
        "## Near-Entry Mass",
        "",
        "| Case | Entry Support | Weights | Mean | Share At Entry Support | Age-22 Share Of Entry-Support Mass | Share Within 0.5 Of Entry |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for case in cases:
        support, weights = entry_support(case)
        bg = np.asarray(case.sol.b_grid, dtype=float)
        g = np.asarray(case.sol.g, dtype=float)
        support_mask = np.zeros(bg.size, dtype=bool)
        for val in support:
            support_mask[np.argmin(np.abs(bg - val))] = True
        mass_support = float(np.sum(g[support_mask, ...]))
        age22_mass_support = float(np.sum(g[support_mask, :, :, 0, ...]))
        share_age22 = age22_mass_support / mass_support if mass_support > 1e-14 else math.nan
        b_entry = float(getattr(case.P, "b_entry_fixed", np.nan))
        within_05 = float(np.sum(g[np.abs(bg - b_entry) <= 0.5, ...]))
        lines.append(
            "| {label} | {support} | {weights} | {mean:.6f} | {mass:.4f} | {age22:.4f} | {within:.4f} |".format(
                label=f"`{case.label}`",
                support=", ".join(f"{x:.3g}" for x in support),
                weights=", ".join(f"{x:.3g}" for x in weights),
                mean=float(np.sum(support * weights)),
                mass=mass_support,
                age22=share_age22,
                within=within_05,
            )
        )

    lines.extend(
        [
            "",
            "## Age/Tenure Decomposition",
            "",
            "| Case | Window | All | Renters | Owners | Childless | Parents | Age 22 | Ages 26-34 | Ages 38-50 | Ages 54+ |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for case in cases:
        for window in (0.25, 0.50, 1.00):
            row = near_entry_row(case, window)
            lines.append(
                "| {label} | +/-{window:.2f} | {all:.4f} | {renters:.4f} | {owners:.4f} | {childless:.4f} | {parents:.4f} | {age22:.4f} | {age26_34:.4f} | {age38_50:.4f} | {age54p:.4f} |".format(
                    label=f"`{case.label}`",
                    window=window,
                    **row,
                )
            )
    path.write_text("\n".join(lines) + "\n")


def near_entry_row(case: Case, window: float) -> dict[str, float]:
    bg = np.asarray(case.sol.b_grid, dtype=float)
    g = np.asarray(case.sol.g, dtype=float)
    b_entry = float(getattr(case.P, "b_entry_fixed", 0.0))
    mask = np.abs(bg - b_entry) <= window
    gm = g[mask, ...]
    ages = ages_for(case.P)
    return {
        "all": float(np.sum(gm)),
        "renters": float(np.sum(gm[:, 0, ...])),
        "owners": float(np.sum(gm[:, 1:, ...])),
        "childless": float(np.sum(gm[:, :, :, :, :, 0, :])),
        "parents": float(np.sum(gm[:, :, :, :, :, 1:, :])),
        "age22": float(np.sum(gm[:, :, :, ages == 22, ...])),
        "age26_34": float(np.sum(gm[:, :, :, (ages >= 26) & (ages <= 34), ...])),
        "age38_50": float(np.sum(gm[:, :, :, (ages >= 38) & (ages <= 50), ...])),
        "age54p": float(np.sum(gm[:, :, :, ages >= 54, ...])),
    }


def plot_atom_heatmaps(cases: list[Case], outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for case in cases:
        bg = np.asarray(case.sol.b_grid, dtype=float)
        ages = ages_for(case.P)
        panels = atom_heatmap_panels(case)
        fig, axes = plt.subplots(2, 3, figsize=(13.2, 6.8), sharex=True, sharey=True)
        vmax = max(float(np.nanmax(v)) for _, v in panels)
        vmax = vmax if vmax > 0 else 1.0
        for ax, (title, arr) in zip(axes.ravel(), panels):
            mesh = ax.pcolormesh(grid_edges(bg), age_edges(ages), arr.T, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax)
            ax.set_title(title)
            ax.set_xlim(-4.0, 8.0)
            ax.set_xlabel("liquid wealth b")
            ax.set_ylabel("age")
        fig.colorbar(mesh, ax=axes.ravel().tolist(), shrink=0.88, label="population share")
        fig.suptitle(f"Atom-origin heatmaps: {case.label}", y=0.995)
        fig.tight_layout(rect=(0.0, 0.0, 0.93, 0.96))
        fig.savefig(outdir / f"atom_origin_heatmap_{case.label}.png", dpi=160)
        plt.close(fig)


def atom_heatmap_panels(case: Case) -> list[tuple[str, np.ndarray]]:
    g = np.asarray(case.sol.g, dtype=float)
    return [
        ("all", np.sum(g, axis=(1, 2, 4, 5, 6))),
        ("renters", np.sum(g[:, 0, ...], axis=(1, 3, 4, 5))),
        ("owners", np.sum(g[:, 1:, ...], axis=(1, 2, 4, 5, 6))),
        ("childless", np.sum(g[:, :, :, :, :, 0, :], axis=(1, 2, 4, 5))),
        ("parents", np.sum(g[:, :, :, :, :, 1:, :], axis=(1, 2, 4, 5, 6))),
        ("age 22 entrants", heatmap_for_age_mask(case, ages_for(case.P) == 22)),
    ]


def heatmap_for_age_mask(case: Case, age_mask: np.ndarray) -> np.ndarray:
    g = np.asarray(case.sol.g, dtype=float)
    out = np.zeros((g.shape[0], g.shape[3]))
    out[:, age_mask] = np.sum(g[:, :, :, age_mask, ...], axis=(1, 2, 4, 5, 6))
    return out


def plot_near_entry_lifecycle(cases: list[Case], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    for case in cases:
        bg = np.asarray(case.sol.b_grid, dtype=float)
        g = np.asarray(case.sol.g, dtype=float)
        ages = ages_for(case.P)
        b_entry = float(getattr(case.P, "b_entry_fixed", 0.0))
        mask = np.abs(bg - b_entry) <= 0.5
        mass_age = np.sum(g[mask, ...], axis=(0, 1, 2, 4, 5, 6))
        ax.plot(ages, mass_age, marker="o", linewidth=2.0, label=case.label)
    ax.set_title("Mass within 0.5 of entry wealth by age")
    ax.set_xlabel("age")
    ax.set_ylabel("population share")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_total_wealth_by_rung(cases: list[Case], outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for case in cases:
        bg = np.asarray(case.sol.b_grid, dtype=float)
        g = np.asarray(case.sol.g, dtype=float)
        p = float(np.asarray(case.p_eq).reshape(-1)[0])
        h_own = np.asarray(case.P.H_own, dtype=float).reshape(-1)
        psi = float(getattr(case.P, "psi", 0.0))
        fig, ax = plt.subplots(figsize=(10.0, 5.2))
        for ten, H in enumerate(h_own, start=1):
            x = bg + (1.0 - psi) * p * float(H)
            mass = np.sum(g[:, ten, ...], axis=(1, 2, 3, 4, 5))
            if float(np.sum(mass)) <= 1e-12:
                continue
            ax.plot(x, mass, marker=".", linewidth=1.0, markersize=4.0, label=f"H={H:g}")
        ax.set_xlim(-1.0, 14.0)
        ax.set_xlabel("total wealth W = b + (1 - psi) p H")
        ax.set_ylabel("population share at rung-grid point")
        ax.set_title(f"Owner total-wealth density by rung: {case.label}")
        ax.legend(frameon=False, ncol=3)
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(outdir / f"total_wealth_by_owner_rung_{case.label}.png", dpi=160)
        plt.close(fig)


def write_entry_data_anchor(cases: list[Case], data_path: Path, path: Path) -> None:
    rows = read_entry_data(data_path)
    lines = [
        "# Entry-Wealth Data Anchor",
        "",
        "This is a mean/median anchor from the existing PSID output. It is not a",
        "full density comparison. The full validation still needs a microdata or",
        "binned distribution extract by age and tenure.",
        "",
    ]
    if not rows:
        lines.append(f"No usable data file found at `{relative(data_path)}`.")
        path.write_text("\n".join(lines) + "\n")
        return
    lines.extend(
        [
            "## Data",
            "",
            "| Sample | Moment | Mean | Median | N |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in rows:
        lines.append(
            f"| `{row['sample']}` | `{row['moment']}` | {fmt(row['mean'])} | {fmt(row['p50'])} | {row['N']} |"
        )
    lines.extend(
        [
            "",
            "## Model Analogues",
            "",
            "| Case | Sample | Mean b/y | Median b/y | Mass |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for case in cases:
        for sample, age_lo, age_hi, childless_only, renter_only in [
            ("young_all_25_30", 25.0, 30.0, False, False),
            ("young_childless_25_35", 25.0, 35.0, True, False),
            ("young_childless_rent_25_35", 25.0, 35.0, True, True),
        ]:
            mean, med, mass = model_liquid_income_ratio(
                case,
                age_lo,
                age_hi,
                childless_only=childless_only,
                renter_only=renter_only,
            )
            lines.append(f"| `{case.label}` | `{sample}` | {mean:.3f} | {med:.3f} | {mass:.4f} |")
    path.write_text("\n".join(lines) + "\n")


def read_entry_data(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    keep = {
        ("young_all_25_30", "liq_nw_to_inc"),
        ("young_childless_25_35", "liq_nw_to_inc"),
        ("young_childless_rent_25_35", "liq_nw_to_inc"),
    }
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    return [r for r in rows if (r.get("sample"), r.get("moment")) in keep]


def model_liquid_income_ratio(
    case: Case,
    age_lo: float,
    age_hi: float,
    *,
    childless_only: bool,
    renter_only: bool,
) -> tuple[float, float, float]:
    P = case.P
    sol = case.sol
    g = np.asarray(sol.g, dtype=float)
    bg = np.asarray(sol.b_grid, dtype=float)
    ages = ages_for(P)
    age_mask = (ages >= age_lo) & (ages <= age_hi)
    z_grid = np.asarray(getattr(P, "z_grid", [1.0]), dtype=float).reshape(-1)
    base_income = np.asarray(getattr(P, "income", np.ones((1, len(ages)))), dtype=float)
    values: list[float] = []
    weights: list[float] = []
    for ib, b in enumerate(bg):
        for j, age_ok in enumerate(age_mask):
            if not age_ok:
                continue
            for zz, z in enumerate(z_grid):
                y = float(base_income[0, j]) * (float(z) if j < int(getattr(P, "J_R", P.J)) else 1.0)
                if y <= 0:
                    continue
                if childless_only:
                    if renter_only:
                        mass = float(np.sum(g[ib, 0, 0, j, zz, 0, :]))
                    else:
                        mass = float(np.sum(g[ib, :, 0, j, zz, 0, :]))
                else:
                    if renter_only:
                        mass = float(np.sum(g[ib, 0, 0, j, zz, :, :]))
                    else:
                        mass = float(np.sum(g[ib, :, 0, j, zz, :, :]))
                if mass > 0:
                    values.append(float(b) / y)
                    weights.append(mass)
    if not weights:
        return math.nan, math.nan, 0.0
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    return float(np.sum(x * w) / np.sum(w)), weighted_median(x, w), float(np.sum(w))


def write_readme(cases: list[Case], outdir: Path, data_path: Path) -> None:
    files = [
        "atom_origin_summary.md",
        "entry_wealth_data_anchor.md",
        "near_entry_mass_by_age.png",
    ]
    for case in cases:
        files.extend(
            [
                f"atom_origin_heatmap_{case.label}.png",
                f"total_wealth_by_owner_rung_{case.label}.png",
            ]
        )
    payload = {
        "status": "diagnostic_only_not_calibration",
        "cases": [{"label": c.label, "cache": str(c.cache_path)} for c in cases],
        "entry_data": str(data_path),
        "files": files,
    }
    (outdir / "manifest.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    lines = [
        "# Entry Atom And Total-Wealth Rung Audit",
        "",
        "Purpose: evaluate the distributional pathologies directly, not through the",
        "calibration objective.",
        "",
        "Files:",
    ]
    lines.extend(f"- `{name}`" for name in files)
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


def entry_support(case: Case) -> tuple[np.ndarray, np.ndarray]:
    sol = case.sol
    if hasattr(sol, "entry_wealth_grid_values") and hasattr(sol, "entry_wealth_grid_weights"):
        values = np.asarray(sol.entry_wealth_grid_values, dtype=float).reshape(-1)
        weights = np.asarray(sol.entry_wealth_grid_weights, dtype=float).reshape(-1)
        if values.size and weights.size == values.size:
            return values, weights / float(np.sum(weights))
    idx, wt = entry_wealth_grid_weights(np.asarray(sol.b_grid, dtype=float), case.P)
    return np.asarray(sol.b_grid, dtype=float)[idx], wt


def ages_for(P: Any) -> np.ndarray:
    return np.asarray(float(P.age_start) + np.arange(int(P.J)) * float(P.da), dtype=float)


def grid_edges(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float).reshape(-1)
    if x.size == 1:
        return np.array([x[0] - 0.5, x[0] + 0.5])
    mids = 0.5 * (x[:-1] + x[1:])
    return np.concatenate([[x[0] - (mids[0] - x[0])], mids, [x[-1] + (x[-1] - mids[-1])]])


def age_edges(ages: np.ndarray) -> np.ndarray:
    ages = np.asarray(ages, dtype=float).reshape(-1)
    step = float(np.median(np.diff(ages))) if ages.size > 1 else 4.0
    return np.concatenate([[ages[0] - 0.5 * step], 0.5 * (ages[:-1] + ages[1:]), [ages[-1] + 0.5 * step]])


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    x = values[order]
    w = weights[order]
    cutoff = 0.5 * float(np.sum(w))
    return float(x[np.searchsorted(np.cumsum(w), cutoff, side="left")])


def fmt(value: str) -> str:
    try:
        return f"{float(value):.3f}"
    except Exception:
        return value


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except Exception:
        return str(path)


if __name__ == "__main__":
    main()
