#!/usr/bin/env python3
"""Read-only calibration autopsy for the certified E5b winner."""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
SOURCE = ROOT / "output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json"
CERTIFIED = SOURCE.parent / "target_fit_full.csv"
NCHS = ROOT / "code/data/nchs_natality_timing/first_birth_counts_year_age.csv"
OUTDIR = ROOT / "output/model/eqscale_seq_e5b_autopsy_20260726"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty output {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[keep], weights[keep]
    if values.size == 0:
        return math.nan
    order = np.argsort(values); values, weights = values[order], weights[order]
    return float(np.interp(q * weights.sum(), np.cumsum(weights), values))


def load_chain() -> Any:
    os.environ["E3_L4"] = "1"; os.environ["E5"] = "1"
    os.environ["E3_TFR_TOP_BIN_WEIGHT"] = "3.602359422009"
    if str(MODEL_ROOT) not in sys.path:
        sys.path.insert(0, str(MODEL_ROOT))
    from intergen_eqscale_seq_optimized import run_e1_chain
    chain = importlib.reload(run_e1_chain)
    chain.load_runtime()
    return chain


def theta() -> dict[str, float]:
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    raw = data["winners"]["E1"]["theta"]
    return {key: float(value) for key, value in raw.items()}


def solve(chain: Any, base_theta: dict[str, float], h0_scale: float = 1.0) -> tuple[Any, Any, np.ndarray, dict[str, Any]]:
    overrides = chain.common_overrides(argparse.Namespace(J=17, Nb=120, max_iter_eq=40, tol_eq=2.5e-5))
    overrides.update(base_theta)
    overrides["H0"] = float(overrides["H0"]) * h0_scale
    sol, P, p = chain.run_model_cp_dt(overrides, verbose=False)
    return sol, P, np.asarray(p, dtype=float).reshape(-1), overrides


def verify(chain: Any, sol: Any, P: Any) -> dict[str, float]:
    moments = chain.extract_moments(sol, P)
    with CERTIFIED.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 12:
        raise ValueError(f"Expected 12 certified rows, found {len(rows)}")
    failed = []
    for row in rows:
        name, expected = row["moment"], float(row["model"])
        actual = float(moments[name]); diff = abs(actual - expected)
        print(f"VERIFY {name}: solved={actual:.12f} certified={expected:.12f} abs_diff={diff:.3e}", flush=True)
        if diff > 1e-6:
            failed.append((name, actual, expected, diff))
    if failed:
        raise RuntimeError(f"Mandatory certified-moment gate failed: {failed}")
    return moments


def ages(P: Any) -> np.ndarray:
    return float(P.age_start) + float(P.da) * np.arange(int(P.J))


def state_arrays(sol: Any) -> tuple[np.ndarray, np.ndarray]:
    g, fp = np.asarray(sol.g, dtype=float), np.asarray(sol.fert_probs, dtype=float)
    # The Markov-income solution must expose the same income-state axis in
    # the stationary distribution and entry policy; fail loudly if it does not.
    if g.ndim != 7 or fp.ndim != 6 or g.shape[2] != 1 or fp.shape[2] != 1:
        raise RuntimeError(f"Unexpected E5b tensor shapes: g={g.shape}, fert_probs={fp.shape}")
    # Solver layout: (b, tenure, location, age, z, parity, child-stage).
    # Autopsy layout: (b, tenure, z, age, parity, child-stage).
    return g[:, :, 0, :, :, :, :].transpose(0, 1, 3, 2, 4, 5), fp[:, :, 0, :, :, :].transpose(0, 1, 3, 2, 4)


def fertility_anatomy(sol: Any, P: Any, outdir: Path) -> tuple[float, float]:
    from intergen_eqscale_seq_optimized.parameters import get_fecundity_by_age
    g, fp = state_arrays(sol); aa = ages(P); fec = np.asarray(get_fecundity_by_age(P), dtype=float)
    rows, flows = [], np.zeros(int(P.J))
    for j in range(int(P.A_f_start) - 1, int(P.A_f_end)):
        mass = g[:, :, :, j, 0, 0]
        attempt = 1.0 - fp[:, :, :, j, 0]
        attempt_share = float(np.sum(mass * attempt) / max(np.sum(mass), 1e-15))
        flows[j] = float(np.sum(mass) * fec[j] * attempt_share)
        rows.append({"age": float(aa[j]), "attempt_share_of_childless": attempt_share,
                     "realized_first_births_mass": flows[j], "share_of_all_first_births": math.nan})
    total = float(flows.sum())
    for row in rows:
        row["share_of_all_first_births"] = row["realized_first_births_mass"] / total
    write_csv(outdir / "first_births_by_age.csv", rows)

    counts: dict[int, float] = {}
    with NCHS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            year, age = int(row["year"]), int(row["age"])
            if 1979 <= year - age <= 1984:
                counts[age] = counts.get(age, 0.0) + float(row["n_first_births"])
    denom = sum(counts.values())
    nchs_rows = [{"granularity": "single_year", "age": age, "age_bin_start": age,
                  "age_bin_end": age, "count": count, "share": count / denom}
                 for age, count in sorted(counts.items())]
    for start in range(18, 47, 4):
        count = sum(v for age, v in counts.items() if start <= age < start + 4)
        nchs_rows.append({"granularity": "model_4_year_bin", "age": math.nan, "age_bin_start": start,
                          "age_bin_end": start + 3, "count": count, "share": count / denom})
    write_csv(outdir / "nchs_first_births_by_age.csv", nchs_rows)

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    model_by_bin = {int(row["age"]): row["share_of_all_first_births"] / float(P.da) for row in rows}
    ax.step(list(model_by_bin), list(model_by_bin.values()), where="mid", lw=2, label="Model (4-year-bin per-year density)")
    x, y = zip(*[(a, c / denom) for a, c in sorted(counts.items()) if 14 <= a <= 46])
    ax.plot(x, y, lw=1.8, label="NCHS cohort shares")
    ax.set(xlim=(14, 46), xlabel="Age", ylabel="Share / year", title="First-birth age distribution")
    ax.grid(alpha=.25); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(outdir / "first_birth_age_distribution.png", dpi=180); plt.close(fig)
    return float(sum(v for a, v in model_by_bin.items() if a >= 30) * float(P.da)), float(sum(v for a, v in counts.items() if a >= 30) / denom)


def childless_tables(sol: Any, P: Any, p: np.ndarray, outdir: Path) -> None:
    g, _ = state_arrays(sol); aa = ages(P); b = np.asarray(sol.b_grid, dtype=float)
    from intergen_eqscale_seq_optimized.parameters import get_fecundity_by_age
    fec = np.asarray(get_fecundity_by_age(P), dtype=float)
    eligible = np.where((fec == 0) & (np.arange(P.J) >= int(P.A_f_end)) & (np.arange(P.J) < int(P.J_R)))[0]
    if eligible.size == 0: raise RuntimeError("No post-fertile working ages")
    block = np.take(g, eligible, axis=3)
    total_wealth = np.empty_like(block)
    for ten in range(g.shape[1]):
        equity = 0.0 if ten == 0 else (1.0 - float(P.psi)) * float(p[0]) * float(P.H_own[ten - 1])
        total_wealth[:, ten, :, :, :, :] = b[:, None, None, None, None] + equity
    mass = block.reshape(-1); wealth = total_wealth.reshape(-1)
    cuts = [weighted_quantile(wealth, mass, q) for q in (.2, .4, .6, .8)]
    rows = []
    for z in range(g.shape[2]):
        m = float(np.sum(block[:, :, z, :, :, :])); c = float(np.sum(block[:, :, z, :, 0, :]))
        rows.append({"group_type": "income_state", "group": z + 1, "z_value": float(np.asarray(sol.type_values)[z]),
                     "wealth_quintile": math.nan, "childless_share": c / m, "mass": m})
    for q in range(5):
        lo, hi = ([-math.inf] + cuts)[q], (cuts + [math.inf])[q]
        mask = (total_wealth >= lo) & (total_wealth <= hi if q == 4 else total_wealth < hi)
        m = float(np.sum(block * mask)); c = float(np.sum(block[:, :, :, :, 0, :] * mask[:, :, :, :, 0, :]))
        rows.append({"group_type": "wealth_quintile", "group": q + 1, "z_value": math.nan,
                     "wealth_quintile": q + 1, "childless_share": c / m, "mass": m})
    write_csv(outdir / "childless_by_income_wealth.csv", rows)


def entry_surface(sol: Any, P: Any, outdir: Path, make_plot: bool = True) -> list[dict[str, Any]]:
    g, fp = state_arrays(sol); aa, b = ages(P), np.asarray(sol.b_grid, dtype=float); rows = []
    for j in range(int(P.A_f_start) - 1, int(P.A_f_end)):
        for z in range(g.shape[2]):
            mass = g[:, 0, z, j, 0, 0]
            for k, wealth in enumerate(b):
                rows.append({"age": float(aa[j]), "z": z + 1, "wealth_grid_point": float(wealth),
                             "prob_try": float(1.0 - fp[k, 0, z, j, 0]), "stationary_mass": float(mass[k])})
    if make_plot:
        fig, ax = plt.subplots(figsize=(7.5, 4.5)); picks = [0, g.shape[2] // 2, g.shape[2] - 1]
        for z in picks:
            series = []
            for j in range(int(P.A_f_start) - 1, int(P.A_f_end)):
                m = g[:, 0, z, j, 0, 0]; med = weighted_quantile(b, m, .5); k = int(np.argmin(abs(b - med)))
                series.append((aa[j], 1.0 - fp[k, 0, z, j, 0]))
            ax.plot(*zip(*series), marker="o", label=f"z state {z + 1}")
        ax.set(xlabel="Age", ylabel="P(try)", title="Entry probability at mass-weighted median renter wealth")
        ax.grid(alpha=.25); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(outdir / "entry_prob_by_age.png", dpi=180); plt.close(fig)
    write_csv(outdir / "entry_prob_by_age_wealth_income.csv", rows)
    return rows


def child_cost(sol: Any, P: Any, p: np.ndarray, theta0: dict[str, float], outdir: Path) -> None:
    from intergen_eqscale_seq_optimized.solver import income_at_state
    g, _ = state_arrays(sol); aa, b = ages(P), np.asarray(sol.b_grid, dtype=float)
    ix = np.where((aa >= 22) & (aa <= 30))[0]
    budgets, masses = [], []
    for z, z_value in enumerate(np.asarray(sol.type_values, dtype=float)):
        for j in ix:
            budgets.extend(float(P.R_gross) * b + income_at_state(P, 0, int(j), float(z_value)))
            masses.extend(g[:, 0, z, j, 0, 0])
    budgets, masses = np.asarray(budgets), np.asarray(masses)
    xs = [weighted_quantile(budgets, masses, q) for q in (.25, .5, .75)]
    if min(xs) <= 0:
        raise RuntimeError(f"Child-cost budget quantiles must be positive, got {xs}")
    alpha0, alpha1, sigma = .733, .733 - theta0["delta_alpha_jump"] - theta0["delta_alpha"], 2.0
    e1, psi = (2.7 / 2.0) ** .7, theta0["psi_child"]
    rent = float(P.user_cost_rate * p[0])
    def util(x: float, alpha: float, e: float, child_psi: float, r: float) -> float:
        c, h = alpha * x, (1 - alpha) * x / r
        composite = c ** alpha * h ** (1 - alpha)
        return e * composite ** (1 - sigma) / (1 - sigma) + child_psi
    rows=[]
    for x in xs:
        for mult in (1.0, 1.10):
            r = rent * mult; vals = [util(x, alpha0, 1, 0, r), util(x, alpha1, 1, 0, r), util(x, alpha1, e1, 0, r), util(x, alpha1, e1, psi, r)]
            rows.append({"budget_x": x, "rent_multiplier": mult, "rent": r, "u_childless": vals[0], "u_tilt_only": vals[1], "u_scale": vals[2], "u_parent": vals[3],
                         "increment_tilt": vals[1]-vals[0], "increment_scale": vals[2]-vals[1], "increment_psi": vals[3]-vals[2]})
    for ix in range(0, len(rows), 2):
        base, high = rows[ix], rows[ix+1]
        for key in ("increment_tilt", "increment_scale", "increment_psi"):
            base[f"change_{key}_rent_plus10"] = high[key] - base[key]
            high[f"change_{key}_rent_plus10"] = high[key] - base[key]
    write_csv(outdir / "child_cost_decomposition.csv", rows)


def passthrough(chain: Any, t: dict[str, float], baseline: tuple[Any, Any, np.ndarray], outdir: Path) -> None:
    rows=[]
    for label, scale, record in [("baseline", 1., baseline), ("H0_minus_7pct", .93, None), ("H0_plus_7pct", 1.07, None)]:
        print(f"PASSTHROUGH_START {label}", flush=True)
        sol, P, p = record if record is not None else solve(chain, t, scale)[:3]
        print(f"PASSTHROUGH_SOLVED {label}", flush=True)
        moments = chain.extract_moments(sol, P); surface = entry_surface(sol, P, outdir, make_plot=False)
        middle = (np.asarray(sol.type_values).size - 1) // 2
        profile=[]
        for age in sorted({r["age"] for r in surface}):
            cells=[r for r in surface if r["age"] == age and r["z"] == middle + 1]; med=weighted_quantile(np.array([r["wealth_grid_point"] for r in cells]), np.array([r["stationary_mass"] for r in cells]), .5)
            profile.append(float(min(cells, key=lambda r: abs(r["wealth_grid_point"]-med))["prob_try"]))
        rows.append({"case":label, "H0_scale":scale, "equilibrium_price_vector":json.dumps(p.tolist()), "rent":float(P.user_cost_rate*p[0]), "tfr":float(moments["tfr"]), "childless_rate":float(moments["childless_rate"]), "mean_age_first_birth":float(moments["mean_age_first_birth"]), "own_rate":float(moments["own_rate"]), "try_profile_middle_z_median_wealth":json.dumps(profile)})
        print(f"PASSTHROUGH_WRITTEN {label}", flush=True)
    write_csv(outdir / "price_passthrough.csv", rows)


def wealth_anatomy(sol: Any, P: Any, p: np.ndarray, outdir: Path) -> None:
    g, _ = state_arrays(sol); aa,b=ages(P),np.asarray(sol.b_grid,dtype=float)
    annual=[]; weights=[]
    for z in range(g.shape[2]):
        for j in range(int(P.J_R)):
            annual.append(float(P.income[0,j] * float(sol.type_values[z]) / float(P.period_years) / (1-float(P.tau_pay)))); weights.append(float(g[:,:,z,j,:,:].sum()))
    ratio=weighted_quantile(np.array(annual),np.array(weights),.9)/weighted_quantile(np.array(annual),np.array(weights),.5)
    print("INCOME_P90_P50", f"{ratio:.12f}", "z_grid", np.asarray(sol.type_values), "stationary_weights", np.array(weights).reshape(g.shape[2],int(P.J_R)).sum(axis=1)/sum(weights), flush=True)
    rows=[]
    for name, lo, hi in [("30_45",30,45),("46_65",46,65),("66_85",66,85),("76_84",76,84)]:
        js=np.where((aa>=lo)&(aa<=hi))[0]; vals=[]; w=[]
        for ten in range(g.shape[1]):
            eq=0 if ten==0 else (1-float(P.psi))*float(p[0])*float(P.H_own[ten-1])
            vals.extend(np.repeat(b+eq, 1)); w.extend(np.take(g[:,ten,:,:,:,:], js, axis=2).sum(axis=(1,2,3,4)))
        v,w=np.asarray(vals),np.asarray(w)
        rows.append({"age_band":name,"p50":weighted_quantile(v,w,.5),"p90":weighted_quantile(v,w,.9),"p99":weighted_quantile(v,w,.99),"mass":float(w.sum())})
    write_csv(outdir / "wealth_percentiles_by_age.csv", rows)


def loss_forensics(outdir: Path) -> None:
    with CERTIFIED.open(newline="",encoding="utf-8") as f: source=list(csv.DictReader(f))
    schemes={"certified":{},"share30plus_se_x2":{"share_first_births_age30plus":.016},"timing_5pct_target":{"mean_age_first_birth":.05*25.310560799362,"share_first_births_age30plus":.05*.270062376851342}}
    rows=[]
    for scheme, changes in schemes.items():
        ranked=[]
        for r in source:
            se=changes.get(r["moment"]); weight=1/se**2 if se else float(r["weight"]); loss=weight*float(r["gap"])**2
            ranked.append((loss,r,weight))
        for rank,(loss,r,weight) in enumerate(sorted(ranked,reverse=True,key=lambda x:x[0]),1):
            rows.append({"scheme":scheme,"rank":rank,"moment":r["moment"],"model":float(r["model"]),"target":float(r["target"]),"gap":float(r["gap"]),"weight":weight,"loss_contribution":loss})
    write_csv(outdir / "loss_forensics.csv", rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-passthrough", action="store_true", help="Reuse an existing pass-through CSV during interrupted-run recovery.")
    args = parser.parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    chain=load_chain(); t=theta(); sol,P,p,_=solve(chain,t); moments=verify(chain,sol,P)
    model30,nchs30=fertility_anatomy(sol,P,OUTDIR)
    print(f"CHILDLESS_SPLIT chosen={sol.childless_chosen_45:.12f} clock={sol.childless_clock_45:.12f} total={moments['childless_rate']:.12f}",flush=True)
    print(f"FIRST_BIRTH_30PLUS model={model30:.12f} nchs={nchs30:.12f}",flush=True)
    childless_tables(sol,P,p,OUTDIR); entry_surface(sol,P,OUTDIR); child_cost(sol,P,p,t,OUTDIR)
    if not args.skip_passthrough:
        passthrough(chain,t,(sol,P,p),OUTDIR)
    print("WEALTH_ANATOMY_START", flush=True); wealth_anatomy(sol,P,p,OUTDIR)
    print("LOSS_FORENSICS_START", flush=True); loss_forensics(OUTDIR)
    print("AUTOPSY_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
