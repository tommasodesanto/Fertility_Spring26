#!/usr/bin/env python3
"""Diagnostic: children currently at home (m) by parent age, model vs. ACS.

    PYTHONPATH=. sandbox/.venv-relative python code/model/sandbox/diagnostics_dependents_by_age.py

Reuses run_ss.py's own solve machinery (same retained theta, same
build_overrides/apply_spec path, same package solver) for the model side --
no reimplementation of the Bellman/root solve. Reuses the ACS 2005-2006
extract and household-head filters already established in
output/model/e5f_matched_pf_20260909a/design_research/housing/inspect_early_housing.py
for the data side.

Writes exactly three files to --out
(default output/model/sandbox/dependents_by_parent_age/):
    dependents_by_age.png, dependents_by_age.csv, README.md

Does not touch code/model/intergen_eqscale_seq_optimized/ or code/model/tools/,
uses no git commands, and launches no cluster jobs.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

SANDBOX_ROOT = Path(__file__).resolve().parent
MODEL_ROOT = SANDBOX_ROOT.parent
REPO_ROOT = MODEL_ROOT.parents[1]
TOOLS_ROOT = MODEL_ROOT / "tools"
DEFAULT_OUT = REPO_ROOT / "output/model/sandbox/dependents_by_parent_age"
ACS_SOURCE = REPO_ROOT / "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
ACS_YEARS = (2005, 2006)

sys.path[:0] = [str(MODEL_ROOT), str(SANDBOX_ROOT), str(TOOLS_ROOT)]

import run_ss  # noqa: E402  (sandbox/run_ss.py -- reused, not modified)

SPEC_NAME = "baseline_psi_fixed"
AGE_BIN_WIDTH = 4.0  # matches P.da


# --------------------------------------------------------------------------
# Model side
# --------------------------------------------------------------------------

def solve_model() -> tuple[Any, Any]:
    """Reproduces exactly the solve path run_ss.py --spec baseline_psi_fixed
    takes (same theta, same overrides, same fix_psi/root logic), warm-started
    from the existing output/model/sandbox/baseline_psi_fixed run if present.
    Does not write anything under output/model/sandbox/baseline_psi_fixed/.
    """
    spec = run_ss.load_spec(SANDBOX_ROOT / "specs" / f"{SPEC_NAME}.yaml")

    import audit_closed_reproductive_closure as closure
    import run_e5f_transition_calibration as calib

    chain = closure.load_chain(profile=run_ss.PROFILE)
    theta, retained_psi, _candidate = run_ss.load_retained_theta()

    nb = run_ss.FULL_NB
    overrides = run_ss.build_overrides(closure, chain, nb, theta)
    overrides, _switches = run_ss.apply_spec(overrides, spec)

    psi_mode = str(spec.get("psi_mode", "root")).lower()
    fix_psi = psi_mode == "fixed" or bool(spec.get("fix_psi", False))
    initial_psi = float(spec.get("psi_child", retained_psi)) if fix_psi else retained_psi

    baseline_dir = run_ss.DEFAULT_OUT_ROOT / SPEC_NAME
    if baseline_dir.exists():
        run_ss.warm_start_price(overrides, baseline_dir)

    t0 = time.perf_counter()
    sol, P, _price, diagnostics, evaluations = run_ss.solve_stationary_state(
        chain, calib, overrides, initial_psi=initial_psi, fix_psi=fix_psi,
    )
    elapsed = time.perf_counter() - t0
    print(f"[model] solved spec={SPEC_NAME} evaluations={evaluations} elapsed={elapsed:.1f}s "
          f"status={diagnostics.get('status')}")
    return sol, P


def m_lookup_table(n_parity: int, n_child_states: int) -> np.ndarray:
    """current children at home m(n,cs) under child_state_mode=independent_count
    (solver.py:current_child_bin_dt, independent_count branch):
    m = cs if cs <= n else 0.
    """
    lut = np.zeros((n_parity, n_child_states), dtype=float)
    for nn in range(n_parity):
        for cs in range(n_child_states):
            lut[nn, cs] = cs if cs <= nn else 0.0
    return lut


def model_series(sol: Any, P: Any) -> pd.DataFrame:
    """g's axis order is (Nb, n_tenure, I, J, [Nz,] n_parity, n_child_states):
    solver.py builds it as `np.zeros((Nb, nt, I, J, npar, ncs))` (6D, no
    permanent-income sub-grid) or `np.zeros((Nb, nt, I, J, Nz, npar, ncs))`
    (7D, under the e5f 15-state income-entry profile this spec's overrides
    use -- confirmed here: measured g.shape = (120, 6, 1, 17, 15, 4, 4)).
    Age is always axis 3; the child axes are always the last two, regardless
    of how many intermediate axes (I, Nz) sit in between, so this locates
    axes by position rather than assuming a fixed ndim.
    """
    g = np.asarray(sol.g, dtype=float)
    J = int(P.J)
    n_parity = int(P.n_parity)
    n_child_states = int(P.n_child_states)
    assert g.shape[3] == J, f"axis 3 of sol.g (size {g.shape[3]}) is not the age axis (P.J={J})"
    assert g.shape[-2] == n_parity, f"second-to-last axis of sol.g (size {g.shape[-2]}) != P.n_parity ({n_parity})"
    assert g.shape[-1] == n_child_states, f"last axis of sol.g (size {g.shape[-1]}) != P.n_child_states ({n_child_states})"
    lut = m_lookup_table(n_parity, n_child_states)
    ages = float(P.age_start) + np.arange(J, dtype=float) * float(P.da)

    rows = []
    for j in range(J):
        age_slice = np.take(g, j, axis=3)  # drops the age axis; last two axes remain (n_parity, n_child_states)
        collapse_axes = tuple(range(age_slice.ndim - 2))
        mass_nm = age_slice.sum(axis=collapse_axes)  # (n_parity, n_child_states)
        total = float(mass_nm.sum())
        if total <= 0:
            rows.append(dict(age=ages[j], mean_m=np.nan, share_m_pos=np.nan,
                              mean_m_given_n_pos=np.nan, share_m0=np.nan, share_m1=np.nan,
                              share_m2=np.nan, share_m3plus=np.nan, total_mass=0.0))
            continue
        mean_m = float((mass_nm * lut).sum() / total)
        share_m_pos = float(mass_nm[lut > 0].sum() / total)
        mass_n_pos = float(mass_nm[1:, :].sum())
        mean_m_given_n_pos = (
            float((mass_nm[1:, :] * lut[1:, :]).sum() / mass_n_pos) if mass_n_pos > 0 else np.nan
        )
        share_m0 = float(mass_nm[lut == 0].sum() / total)
        share_m1 = float(mass_nm[lut == 1].sum() / total)
        share_m2 = float(mass_nm[lut == 2].sum() / total)
        share_m3plus = float(mass_nm[lut >= 3].sum() / total)
        rows.append(dict(age=ages[j], mean_m=mean_m, share_m_pos=share_m_pos,
                          mean_m_given_n_pos=mean_m_given_n_pos, share_m0=share_m0,
                          share_m1=share_m1, share_m2=share_m2, share_m3plus=share_m3plus,
                          total_mass=total))
    df = pd.DataFrame(rows)
    df.insert(0, "source", "model")
    return df


# --------------------------------------------------------------------------
# ACS side
# --------------------------------------------------------------------------

class HeaderReader(pd.io.stata.StataReader):
    """Identical minimal header-only reader used by
    design_research/housing/inspect_early_housing.py, to avoid pandas 1.5.3's
    full in-memory copy of a 9.2 GiB .dta file."""

    def __init__(self, path: Path) -> None:
        pd.io.stata.StataParser.__init__(self)
        self.col_sizes = []
        self._convert_dates = False
        self._convert_categoricals = False
        self._index_col = None
        self._convert_missing = False
        self._preserve_dtypes = True
        self._columns = None
        self._order_categoricals = True
        self._encoding = ""
        self._chunksize = 1
        self._using_iterator = False
        self._has_string_data = False
        self._missing_values = False
        self._can_read_value_labels = False
        self._column_selector_set = False
        self._value_labels_read = False
        self._data_read = False
        self._dtype = None
        self._lines_read = 0
        self._native_byteorder = "<" if sys.byteorder == "little" else ">"
        self.path_or_buf = path.open("rb")
        self._read_header()
        self._setup_dtype()


def acs_series(age_edges: np.ndarray) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Household-head sample, ACS 2005-2006 pooled, same head filter as
    inspect_early_housing.py (gq in {1,2}, pernum==1, relate==1, hhwt>0,
    age 18-85) but WITHOUT its ownershp/rooms>0 restriction, since dependents
    are wanted for the whole head population, not just the housing-row
    subsample. NCHILD = own children of any age currently in the household
    (IPUMS ACS definition: "children at home", any age); YNGCH/ELDCH give the
    age of the youngest/eldest own child but not a per-child age list, so
    "own children under 18" is a lower bound: NCHILD counted only for
    households where ELDCH<18 (eldest child is a minor, so every child NCHILD
    counts is under 18); mixed adult+minor-child households (ELDCH>=18 and
    YNGCH<18) are flagged separately and excluded from that lower bound.
    """
    r = HeaderReader(ACS_SOURCE)
    dtype = r._dtype
    fields = dict(zip(r.varlist, dtype.names))
    data = np.memmap(ACS_SOURCE, dtype=dtype, mode="r", offset=r.data_location, shape=(r.nobs,))

    def lower(y: int) -> int:
        lo, hi = 0, r.nobs
        while lo < hi:
            mid = (lo + hi) // 2
            if int(data[fields["year"]][mid]) < y:
                lo = mid + 1
            else:
                hi = mid
        return lo

    required = ("year", "sample", "gq", "pernum", "relate", "hhwt", "age", "nchild", "yngch", "eldch")
    assert all(x in fields for x in required)

    frames = []
    meta = []
    for y in ACS_YEARS:
        lo, hi = lower(y), lower(y + 1)
        assert hi - lo <= 4_000_000
        block = data[lo:hi]
        a = {x: np.asarray(block[fields[x]]) for x in required}
        assert np.all(a["year"] == y)
        keep = (
            (a["sample"] == y * 100 + 1)
            & np.isin(a["gq"], (1, 2))
            & (a["pernum"] == 1)
            & (a["relate"] == 1)
            & (a["hhwt"] > 0)
            & (a["age"] >= 18)
            & (a["age"] <= 85)
        )
        h = pd.DataFrame({x: a[x][keep].astype(float) for x in required})
        frames.append(h)
        meta.append({"year": y, "raw_records": int(hi - lo), "head_records": int(keep.sum()),
                      "head_weight": float(h.hhwt.sum())})
    heads = pd.concat(frames, ignore_index=True)

    heads["nchild_any_age"] = heads["nchild"]
    heads["has_minor_child"] = (heads["nchild"] > 0) & (heads["yngch"] < 18)
    # Lower-bound "own children under 18": full NCHILD only when the eldest
    # child is also under 18 (so no adult child is being double counted as a
    # minor); ambiguous mixed-age families get NaN and are excluded from the
    # under-18 mean, not zeroed.
    all_minor = (heads["nchild"] > 0) & (heads["eldch"] < 18)
    heads["nchild_under18_lb"] = np.where(
        heads["nchild"] == 0, 0.0, np.where(all_minor, heads["nchild"], np.nan)
    )
    mixed_family = (heads["nchild"] > 0) & (heads["yngch"] < 18) & (heads["eldch"] >= 18)

    rows = []
    for lo_edge in age_edges:
        hi_edge = lo_edge + AGE_BIN_WIDTH
        bucket = heads[(heads.age >= lo_edge) & (heads.age < hi_edge)]
        w = bucket.hhwt.to_numpy()
        wsum = w.sum()
        if wsum <= 0:
            rows.append(dict(age=lo_edge, mean_m=np.nan, share_m_pos=np.nan,
                              mean_m_given_n_pos=np.nan, share_m0=np.nan, share_m1=np.nan,
                              share_m2=np.nan, share_m3plus=np.nan, total_mass=0.0,
                              mean_m_under18_lb=np.nan, mixed_family_share=np.nan))
            continue
        nchild = bucket.nchild_any_age.to_numpy()
        mean_m = float((w * nchild).sum() / wsum)
        share_m_pos = float(w[nchild > 0].sum() / wsum)
        pos_mask = nchild > 0
        mean_m_given_n_pos = float((w[pos_mask] * nchild[pos_mask]).sum() / w[pos_mask].sum()) if pos_mask.any() else np.nan
        share_m0 = float(w[nchild == 0].sum() / wsum)
        share_m1 = float(w[nchild == 1].sum() / wsum)
        share_m2 = float(w[nchild == 2].sum() / wsum)
        share_m3plus = float(w[nchild >= 3].sum() / wsum)
        under18 = bucket.nchild_under18_lb.to_numpy()
        u_mask = ~np.isnan(under18)
        mean_m_u18 = float((w[u_mask] * under18[u_mask]).sum() / w[u_mask].sum()) if u_mask.any() else np.nan
        mixed_mask = mixed_family.loc[bucket.index].to_numpy()
        mixed_share = float(w[mixed_mask].sum() / wsum) if len(bucket) else np.nan
        rows.append(dict(age=lo_edge, mean_m=mean_m, share_m_pos=share_m_pos,
                          mean_m_given_n_pos=mean_m_given_n_pos, share_m0=share_m0,
                          share_m1=share_m1, share_m2=share_m2, share_m3plus=share_m3plus,
                          total_mass=float(wsum), mean_m_under18_lb=mean_m_u18,
                          mixed_family_share=mixed_share))
    df = pd.DataFrame(rows)
    df.insert(0, "source", "acs_2005_2006")
    receipt = {"source": str(ACS_SOURCE), "years": ACS_YEARS, "year_meta": meta,
               "n_heads_total": int(len(heads)), "weight_total": float(heads.hhwt.sum())}
    return df, receipt


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

def write_png(model_df: pd.DataFrame, acs_df: pd.DataFrame, out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    ax.plot(model_df.age, model_df.mean_m, "o-", label="Model E[m | age]")
    ax.plot(acs_df.age, acs_df.mean_m, "s-", label="ACS mean NCHILD (any age)")
    ax.plot(acs_df.age, acs_df.mean_m_under18_lb, "s--", alpha=0.6, label="ACS mean children <18 (lower bound)")
    ax.set_xlabel("age")
    ax.set_ylabel("mean children at home")
    ax.set_title("Mean children at home by age")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.plot(model_df.age, model_df.share_m_pos, "o-", label="Model share(m>0)")
    ax.plot(acs_df.age, acs_df.share_m_pos, "s-", label="ACS share(NCHILD>0)")
    ax.set_xlabel("age")
    ax.set_ylabel("share with any child at home")
    ax.set_title("Share with any child at home, by age")
    ax.legend(fontsize=8)

    fig.suptitle("Children at home by parent age: model (baseline_psi_fixed) vs. ACS 2005-2006")
    fig.tight_layout()
    fig.savefig(out_dir / "dependents_by_age.png", dpi=150)
    plt.close(fig)


def write_readme(out_dir: Path, mu: float, A_m: float, period_years: float,
                  model_df: pd.DataFrame, acs_df: pd.DataFrame, acs_receipt: dict[str, Any]) -> None:
    key_ages = (26, 30, 42, 58, 66, 74)

    def row_at(df: pd.DataFrame, age: float) -> pd.Series:
        idx = (df.age - age).abs().idxmin()
        return df.loc[idx]

    lines = []
    lines.append("# Children at home by parent age: model vs. ACS")
    lines.append("")
    lines.append("**Model**: spec `baseline_psi_fixed` (retained E5F calibration, psi_child held fixed, "
                 "one GE stationary evaluation, Nb=120/J=17). `m` = children currently at home, read off "
                 "`sol.g[b,tenure,i,j,n,cs]` under `child_state_mode=independent_count` "
                 "(solver.py `current_child_bin_dt`): m(n,cs) = cs if cs<=n else 0, where n = children "
                 "ever born (parity axis) and cs = child-state axis. Model ages are "
                 "18 + j*4 for j=0..16 (P.age_start, P.da).")
    lines.append("")
    lines.append("**Data**: ACS 2005-2006 pooled (extract27.dta), household heads only "
                 "(gq in {1,2}, pernum==1, relate==1, hhwt>0, age 18-85), same head filter as "
                 "output/model/e5f_matched_pf_20260909a/design_research/housing/inspect_early_housing.py "
                 "minus its ownershp/rooms restriction (not relevant to dependents). "
                 "\"Children of any age at home\" = IPUMS NCHILD directly. "
                 "\"Children under 18 (lower bound)\" = NCHILD only for households where ELDCH<18 "
                 "(so every counted child is confirmed a minor); mixed adult+minor-child households "
                 "(YNGCH<18 and ELDCH>=18) are excluded from that lower bound rather than zeroed, "
                 "since ACS gives no per-child age list to split them. 4-year age bins matching model ages.")
    lines.append("")
    lines.append(f"**mu (implied per-period exit probability)**: A_m = {A_m:g} years "
                 f"(expected duration a child stays at home), period length = {period_years:g} years, "
                 f"so under the model's memoryless per-period exit process "
                 f"mu = period_years / A_m = {mu:.4f} per 4-year period "
                 "(a per-child-period hazard, not an age-dependent schedule).")
    lines.append("")
    lines.append("## Key numbers (mean children at home)")
    lines.append("| Age | Model E[m] | ACS mean NCHILD (any age) | ACS mean <18 (lower bound) | "
                 "Model share(m>0) | ACS share(NCHILD>0) |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for age in key_ages:
        mr = row_at(model_df, age)
        ar = row_at(acs_df, age)
        lines.append(f"| {age} | {mr.mean_m:.3f} | {ar.mean_m:.3f} | {ar.mean_m_under18_lb:.3f} | "
                     f"{mr.share_m_pos:.3f} | {ar.share_m_pos:.3f} |")
    lines.append("")
    lines.append("## Reading")
    y_m, y_a = row_at(model_df, 30).mean_m, row_at(acs_df, 30).mean_m_under18_lb
    m_m, m_a = row_at(model_df, 42).mean_m, row_at(acs_df, 42).mean_m_under18_lb
    o_m, o_a = row_at(model_df, 66).mean_m, row_at(acs_df, 66).mean_m_under18_lb
    lines.append(f"- Young (22-34): model {y_m:.2f} vs. ACS-under-18 {y_a:.2f} at age 30 -- "
                 f"{'model overstates' if y_m > y_a else 'model understates'} dependents at this margin.")
    lines.append(f"- Middle (38-54): model {m_m:.2f} vs. ACS-under-18 {m_a:.2f} at age 42 -- "
                 f"{'model overstates' if m_m > m_a else 'model understates'} dependents at this margin.")
    lines.append(f"- Old (58+): model {o_m:.2f} vs. ACS-under-18 {o_a:.2f} at age 66 -- "
                 f"the model's memoryless mu-exit process ({mu:.3f}/period) implies a geometric tail of "
                 f"children still coded 'at home' at old parent ages that ACS households do not show; "
                 f"{'model overstates' if o_m > o_a else 'model understates'} dependents at this margin.")
    lines.append("")
    lines.append(f"ACS receipt: {json.dumps(acs_receipt)}")
    (out_dir / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    args = parser.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    sol, P = solve_model()
    model_df = model_series(sol, P)

    age_edges = model_df.age.to_numpy()
    acs_df, acs_receipt = acs_series(age_edges)

    combined = pd.concat([model_df, acs_df], ignore_index=True)
    combined.to_csv(out_dir / "dependents_by_age.csv", index=False)

    write_png(model_df, acs_df, out_dir)

    A_m = float(getattr(P, "A_m", 18.0))
    period_years = float(P.period_years)
    mu = period_years / A_m
    write_readme(out_dir, mu, A_m, period_years, model_df, acs_df, acs_receipt)

    print(f"Wrote {out_dir}/{{dependents_by_age.png,dependents_by_age.csv,README.md}}")


if __name__ == "__main__":
    main()
