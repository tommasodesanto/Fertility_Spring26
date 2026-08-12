#!/usr/bin/env python3
"""Solve the quota-closure equilibrium from the traced loci (diagnostic).

Demographic locus S_d(r) = Mbar / (E - Rbar * B(r)) with the production
floor-arm quota objects; clearing locus S_c(r) = H^S(r)/D(r) from the trace.
Intersect on the trace grid, then vary Mbar to quantify Proposition 5's
perturbation: how long-run population and rent respond to the outside inflow.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

SCRATCH = Path(__file__).resolve().parent
RBAR = 0.9692247023019599
MBAR0 = 0.010432954094546565

rows = [
    row
    for row in csv.DictReader((SCRATCH / "trace_points.csv").open())
    if row["status"] == "ok"
]
rows.sort(key=lambda row: float(row["rent"]))
rent = np.array([float(row["rent"]) for row in rows])
B = np.array([float(row["mature_entrants_B"]) for row in rows])
E = float(rows[0]["entry_rate_E"])
Sc = np.array([float(row["clearing_scale_S"]) for row in rows])

lgrid = np.linspace(np.log(rent[1]), np.log(rent[-1]), 4000)
Bg = np.interp(lgrid, np.log(rent), B)
Scg = np.exp(np.interp(lgrid, np.log(rent), np.log(Sc)))

print("min over grid of E - Rbar*B(r):", np.min(E - RBAR * Bg))
base_idx = None
results = []
for frac in [1.0, 0.9, 0.75, 0.5, 0.25, 0.1]:
    Sd = frac * MBAR0 / (E - RBAR * Bg)
    diff = np.log(Scg) - np.log(Sd)
    k = int(np.argmin(np.abs(diff)))
    r_star, s_star = float(np.exp(lgrid[k])), float(Scg[k])
    if frac == 1.0:
        base_idx = (r_star, s_star)
    results.append((frac, r_star, s_star))
    print(
        f"Mbar x{frac:4.2f}: rent*={r_star:.5f}  S*={s_star:.5f}  "
        f"(rel to Mbar x1: rent {r_star/base_idx[0]-1:+.3%}, S {s_star/base_idx[1]-1:+.3%})"
    )

# Local elasticities at the baseline intersection.
k0 = int(np.argmin(np.abs(np.log(Scg) - np.log(MBAR0 / (E - RBAR * Bg)))))
h = 60
eta_c = (np.log(Scg[k0 + h]) - np.log(Scg[k0 - h])) / (lgrid[k0 + h] - lgrid[k0 - h])
phi = (
    np.log(E - RBAR * Bg[k0 + h]) - np.log(E - RBAR * Bg[k0 - h])
) / (lgrid[k0 + h] - lgrid[k0 - h])
print(f"\nlocal clearing elasticity eta_c = {eta_c:.3f}")
print(f"local demographic-feedback elasticity phi = {phi:.3f}")
print(f"dlnS*/dlnM = eta_c/(eta_c+phi) = {eta_c/(eta_c+phi):.3f}")
print(f"dln rent*/dlnM = 1/(eta_c+phi) = {1/(eta_c+phi):.3f}")
print(f"feedback offset share phi/(eta_c+phi) = {phi/(eta_c+phi):.3f}")
