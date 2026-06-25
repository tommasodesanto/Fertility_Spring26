#!/usr/bin/env python3
"""Reproduce + quantify the Markov room-moment collapse bug (read-only).

Root cause hypothesis: compute_markov_statistics() collapses the income (z)
dimension BEFORE compute_statistics() applies nonlinear (threshold / median)
operators to the renter housing policy. For a continuous renter policy that
varies across z within a (b,j,n,cs) cell, 1{E_z[h]>=t} != E_z[1{h>=t}] and
median(E_z[h]) != weighted-median over z. Means are linear and survive.

This script computes each affected renter moment three ways:
  stored      : sol.<attr> from the cache (the buggy reported value)
  buggy_repro : recompute the collapsed-then-threshold path (should == stored)
  corrected   : compute the threshold/median on the FULL z-resolved g and hR
"""
from __future__ import annotations
import pickle, sys
from pathlib import Path
import numpy as np

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
sys.path.insert(0, str(REPO / "code/model"))
from intergen_housing_fertility.solver import (  # noqa: E402
    collapse_markov_policy, current_child_bin_dt, age_to_index,
    weighted_median_from_cells, income_transition_values,
)

d = pickle.load(open(REPO / "output/model/intergen_model_run_current/solution_cache.pkl", "rb"))
sol = d["baseline"]["sol"]; P = d["baseline"]["P"]
g = np.asarray(sol.g, float)          # (Nb, nt, I, J, Nz, npar, ncs)
hR = np.asarray(sol.hR_pol, float)
Nb, nt, I, J, Nz, npar, ncs = g.shape
z_grid, z_weights, _ = income_transition_values(P)
g_tot = np.sum(g, axis=4)
hR_tot = collapse_markov_policy(hR, g, z_weights)
dep_last = P.n_child_stages
a25s, a45e = age_to_index(P, 25), age_to_index(P, 45)
a30s, a55e = age_to_index(P, 30), age_to_index(P, 55)
hRmax = float(P.hR_max)
print(f"age idx: 25->{a25s} 45->{a45e} 30->{a30s} 55->{a55e}; Nz={Nz}; hR_max={hRmax}")

def share_ge_collapsed(jlo, jhi, thr, child_bins):
    num = den = 0.0
    for j in range(jlo, jhi + 1):
        for i in range(I):
            for nn in range(npar):
                for cs in range(ncs):
                    if current_child_bin_dt(nn, cs, dep_last) not in child_bins: continue
                    gr = g_tot[:, 0, i, j, nn, cs]; hr = hR_tot[:, 0, i, j, nn, cs]
                    kr = (gr > 0) & np.isfinite(hr) & (hr > 0)
                    den += float(np.sum(gr[kr])); num += float(np.sum(gr[kr][hr[kr] >= thr - 1e-8]))
    return num / max(den, 1e-12)

def share_ge_correct(jlo, jhi, thr, child_bins):
    num = den = 0.0
    for j in range(jlo, jhi + 1):
        for i in range(I):
            for zz in range(Nz):
                for nn in range(npar):
                    for cs in range(ncs):
                        if current_child_bin_dt(nn, cs, dep_last) not in child_bins: continue
                        gr = g[:, 0, i, j, zz, nn, cs]; hr = hR[:, 0, i, j, zz, nn, cs]
                        kr = (gr > 0) & np.isfinite(hr) & (hr > 0)
                        den += float(np.sum(gr[kr])); num += float(np.sum(gr[kr][hr[kr] >= thr - 1e-8]))
    return num / max(den, 1e-12)

def median_collapsed(jlo, jhi):
    vals, wts = [], []
    for j in range(jlo, jhi + 1):
        for i in range(I):
            for nn in range(npar):
                for cs in range(ncs):
                    if current_child_bin_dt(nn, cs, dep_last) != 2: continue
                    gr = g_tot[:, 0, i, j, nn, cs]; hr = hR_tot[:, 0, i, j, nn, cs]
                    kr = (gr > 0) & np.isfinite(hr) & (hr > 0)
                    if np.any(kr): vals.append(hr[kr]); wts.append(gr[kr])
    return weighted_median_from_cells(vals, wts)

def median_correct(jlo, jhi):
    vals, wts = [], []
    for j in range(jlo, jhi + 1):
        for i in range(I):
            for zz in range(Nz):
                for nn in range(npar):
                    for cs in range(ncs):
                        if current_child_bin_dt(nn, cs, dep_last) != 2: continue
                        gr = g[:, 0, i, j, zz, nn, cs]; hr = hR[:, 0, i, j, zz, nn, cs]
                        kr = (gr > 0) & np.isfinite(hr) & (hr > 0)
                        if np.any(kr): vals.append(hr[kr]); wts.append(gr[kr])
    return weighted_median_from_cells(vals, wts)

def mean_rooms_correct(jlo, jhi, child_bins):
    num = den = 0.0
    for j in range(jlo, jhi + 1):
        for i in range(I):
            for zz in range(Nz):
                for nn in range(npar):
                    for cs in range(ncs):
                        if current_child_bin_dt(nn, cs, dep_last) not in child_bins: continue
                        gr = g[:, 0, i, j, zz, nn, cs]; hr = hR[:, 0, i, j, zz, nn, cs]
                        kr = (gr > 0) & np.isfinite(hr) & (hr > 0)
                        den += float(np.sum(gr[kr])); num += float(np.sum(gr[kr] * hr[kr]))
    return num / max(den, 1e-12)

rows = [
    ("prime30_55_childless_renter_share_rooms_ge6", "TARGET (some sets, w=25)",
     getattr(sol, "prime30_55_childless_renter_share_rooms_ge6", np.nan),
     share_ge_collapsed(a30s, a55e, 6.0, {2}), share_ge_correct(a30s, a55e, 6.0, {2})),
    ("prime_childless_renter_median_rooms", "TARGET (some sets, w=10)",
     getattr(sol, "prime_childless_renter_median_rooms", np.nan),
     median_collapsed(a25s, a45e), median_correct(a25s, a45e)),
    ("renter25_45_all_cap_share", "diagnostic",
     getattr(sol, "renter25_45_all_cap_share", np.nan),
     share_ge_collapsed(a25s, a45e, hRmax, {2, 3, 4}), share_ge_correct(a25s, a45e, hRmax, {2, 3, 4})),
    ("renter25_45_current0_cap_share", "diagnostic",
     getattr(sol, "renter25_45_current0_cap_share", np.nan),
     share_ge_collapsed(a25s, a45e, hRmax, {2}), share_ge_correct(a25s, a45e, hRmax, {2})),
    ("renter25_45_current1_cap_share", "diagnostic",
     getattr(sol, "renter25_45_current1_cap_share", np.nan),
     share_ge_collapsed(a25s, a45e, hRmax, {3}), share_ge_correct(a25s, a45e, hRmax, {3})),
]

print(f"\n{'moment':<46}{'role':<26}{'stored':>9}{'buggy_repro':>12}{'CORRECTED':>11}")
for name, role, stored, buggy, corr in rows:
    print(f"{name:<46}{role:<26}{stored:>9.4f}{buggy:>12.4f}{corr:>11.4f}")

# means must be invariant (linear)
print("\n--- linearity check (means should be ~equal stored vs corrected) ---")
rm_stored = getattr(sol, "prime30_55_childless_renter_mean_rooms", np.nan)
rm_corr = mean_rooms_correct(a30s, a55e, {2})
gap_stored = getattr(sol, "prime30_55_childless_owner_minus_renter_mean_rooms", np.nan)
print(f"renter mean rooms 30-55 childless: stored={rm_stored:.4f}  corrected={rm_corr:.4f}  (ACTIVE-TARGET gap stored={gap_stored:.4f})")
print(f"  => mean gap (active target, weight 12) uses this renter mean; if invariant, the active loss is NOT corrupted.")
