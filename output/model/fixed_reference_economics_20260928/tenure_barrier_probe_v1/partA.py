"""Part A: zero-solve decomposition on saved need1 arrays (no model solves).

Population: pre-tenure renter parents (origin tenure to=0, children at home
m>=1) whose conditional renter policy is at the 6-room cap, weighted by saved
pre-tenure mass g_beginning_distribution (verified: gb x tenure_probs[...,0]
summed over origins reproduces realized renter mass to ~3e-10; see RECEIPT).

Per (m, age cell j): rung feasibility under the kernel purchase screen,
tenure choice probabilities, mean b / income state / grid-bound proximity.
Writes partA_tables/ CSVs + MD. Read-only w.r.t. engine + saved arrays.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
PREV = ROOT / "output/model/fixed_reference_economics_20260928/per_child_need_probe_v1"

PRICE = 0.6744838540900874
H_OWN = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
CAP = 6.0
ARMS = (("need1_phi10", 1.0), ("need1_phi08", 0.8))


def main():
    sys.path.insert(0, str(PREV))
    import run_probe as rp
    native = rp.setup_imports(HERE / "partA_setup_tmp")
    out_pre = HERE / "partA_setup_tmp"
    out_pre.mkdir(parents=True, exist_ok=True)
    _, _, _, _, P_ref, grid, _ = rp.build_base(native, out_pre)
    assert np.array_equal(np.asarray(P_ref.H_own, dtype=float), H_OWN), "H_own drift"
    assert int(P_ref.n_house) == 5, "n_house drift"
    assert bool(getattr(P_ref, "native_purchase_income", False)) is True, "purchase timing drift"
    b_grid = np.asarray(grid, dtype=float)
    age_start = float(getattr(P_ref, "age_start", 18.0))
    da = float(getattr(P_ref, "da", getattr(P_ref, "period_years", 4.0)))
    z_grid = np.asarray(P_ref.z_grid, dtype=float)
    Rg = float(P_ref.R_gross)
    sys.path.insert(0, str(ROOT / "output/model/fixed_reference_economics_20260928"
                            "/purchase_rules_overnight_v1/engines/quarter"))
    from small_credit_lab.engine.shared import income_at_state
    income_jz = np.array([[income_at_state(P_ref, 0, j, float(zv)) for zv in z_grid]
                          for j in range(int(P_ref.J))])  # (J, Nz)

    tables = {}
    for lab, phi in ARMS:
        rep = PREV / lab / "phase_b_ge" / (lab + "_selected")
        stage = rep / "stage" / "solution_arrays.npz"
        a = np.load(stage)
        gb = np.asarray(a["g_beginning_distribution"], dtype=float)  # pre-tenure
        tp = np.asarray(a["tenure_probs"], dtype=float)  # (b,to,loc,j,z,n,cs,choice)
        hR = np.asarray(a["hR_pol"], dtype=float)  # (b,to,loc,j,z,n,cs)
        assert gb.shape == hR.shape == tp.shape[:-1], (gb.shape, hR.shape, tp.shape)
        assert bool(a["shared.birth_dp"].any()) is False
        assert float(np.abs(a["shared.birth_entry_grant"]).max()) == 0.0
        assert set(np.unique(a["shared.phi_choice"]).tolist()) == {phi}, \
            np.unique(a["shared.phi_choice"])

        hc = PRICE * H_OWN  # purchase price per rung (household.py:181)
        dp0 = (1.0 - phi) * hc  # base down-payment screen (household.py:189)
        bm0 = -phi * hc  # base borrowing floor after purchase (household.py:190)
        # Operative screen counts current income (native_purchase_income=True):
        # dp_choice = dp - y/Rg, floor = max(bmo - y/Rg, b_grid[0])
        # (household.py:907-912); renter-origin test kernels.py:403-405.

        ncs = gb.shape[-1]
        rows = []
        for m in (1, 2, 3):
            for j in range(gb.shape[3]):
                w = gb[:, 0, :, j, :, :, m].copy()  # origin renters, child state m
                # cap indicator on conditional renter policy (to=0)
                at_cap = hR[:, 0, :, j, :, :, m] >= CAP - 1e-9
                wcap = w * at_cap
                mass = float(wcap.sum())
                row = {"arm": lab, "phi": phi, "m": m, "age_cell_j": j,
                       "age_low": age_start + da * j, "age_high": age_start + da * (j + 1) - 1,
                       "cell_mass": mass}
                if mass > 0:
                    # 1. rung feasibility under the operative income-augmented
                    # purchase screen + borrowing floors.
                    def feas_mask(t):
                        zz = np.arange(gb.shape[4])
                        y = income_jz[j][zz]  # (Nz,)
                        dpc = dp0[t - 1] - y / Rg
                        bmc = np.maximum(bm0[t - 1] - y / Rg, b_grid[0])
                        b = b_grid[:, None, None, None]  # (b,loc,z,n)->broadcast
                        f = (b >= dpc[None, None, :, None]) & \
                            ((b - hc[t - 1]) >= bmc[None, None, :, None])
                        return np.broadcast_to(f, wcap.shape)
                    for t, H in enumerate(H_OWN, start=1):
                        row["feas_rung_H%g" % H] = float((wcap * feas_mask(t)).sum() / mass)
                    # cross-check: positive choice prob on infeasible cells
                    worst = 0.0
                    for t, H in enumerate(H_OWN, start=1):
                        f = feas_mask(t)
                        bad = float(((wcap * tp[:, 0, :, j, :, :, m][..., t])[~f]).sum())
                        worst = max(worst, bad)
                    row["infeas_positive_prob_mass"] = worst
                    # 2. tenure choice probs, mass-weighted mean
                    pr = tp[:, 0, :, j, :, :, m, :]  # (b,loc,z,choice)
                    for t in range(tp.shape[-1]):
                        row["prob_choice_t%d" % t] = float((wcap * pr[..., t]).sum() / mass)
                    # 3. financial position / income / bounds
                    row["mean_b"] = float((wcap * b_grid[:, None, None, None]).sum() / mass)
                    row["share_b_at_bottom_node"] = float(wcap[0].sum() / mass)
                    row["share_b_at_top_node"] = float(wcap[-1].sum() / mass)
                    zidx = np.arange(gb.shape[4])
                    zm = wcap.sum(axis=(0, 1, 3))
                    row["mean_z_index"] = float((zm * zidx).sum() / mass)
                    row["mean_income"] = float((wcap * income_jz[j][None, None, :, None]).sum() / mass)
                    for z in zidx:
                        row["zshare_%d" % z] = float(zm[z] / mass)
                rows.append(row)
        tables[lab] = rows

    outdir = HERE / "partA_tables"
    outdir.mkdir(exist_ok=True)
    for lab, rows in tables.items():
        keys = list(rows[0].keys())
        with (outdir / ("partA_" + lab + ".csv")).open("w", newline="") as s:
            w = csv.DictWriter(s, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
    # pooled rows
    with (outdir / "partA_pooled.csv").open("w", newline="") as s:
        w = csv.writer(s)
        hdr = ["arm", "phi", "m", "mass", "mean_b", "mean_income"] + \
            ["feas_rung_H%g" % H for H in H_OWN] + \
            ["prob_choice_t%d" % t for t in range(6)] + \
            ["mean_z_index", "share_b_bottom", "share_b_top"]
        w.writerow(hdr)
        for lab, rows in tables.items():
            for m in (1, 2, 3):
                sub = [r for r in rows if r["m"] == m and r["cell_mass"] > 0]
                mass = sum(r["cell_mass"] for r in sub)
                if mass == 0:
                    continue
                def avg(k):
                    return sum(r[k] * r["cell_mass"] for r in sub) / mass
                w.writerow([lab, dict(ARMS)[lab], m, mass, avg("mean_b"), avg("mean_income")] +
                           [avg("feas_rung_H%g" % H) for H in H_OWN] +
                           [avg("prob_choice_t%d" % t) for t in range(6)] +
                           [avg("mean_z_index"), avg("share_b_at_bottom_node"),
                            avg("share_b_at_top_node")])
    # MD summary
    lines = ["# Part A tables — capped renter parents (pre-tenure weights)",
             "",
             "Population: origin-tenure renters (to=0) with children at home m>=1 whose",
             "conditional renter policy hR_pol is at the 6-room cap, weighted by saved",
             "pre-tenure mass g_beginning_distribution. Per-cell CSVs: partA_<arm>.csv;",
             "pooled: partA_pooled.csv. Choice-specific tenure values are NOT saved in",
             "stage/solution_arrays.npz (kernels return VH/tcj/probs only), so no",
             "rent-vs-best-rung value gap is reported.",
             "",
             "Operative screen (native_purchase_income=True): renter buying rung tn is",
             "feasible iff b >= (1-phi)*p*H - y/Rg AND b - p*H >= max(-phi*p*H - y/Rg, b_grid[0]).",
             "",
             "z_grid: " + ", ".join("%.4g" % v for v in z_grid),
             "Rg: %.6g" % Rg]
    (outdir / "partA_summary.md").write_text("\n".join(lines) + "\n")
    print("wrote", outdir, "rows:", {k: len(v) for k, v in tables.items()})


if __name__ == "__main__":
    main()
