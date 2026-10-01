"""Fail-closed saved-state audit for the isolated renter-only credit diagnostic."""
from __future__ import annotations

import numpy as np


def audit_renter_support(P, evaluation, solution, grid, trace):
    bg = np.asarray(grid, dtype=float)
    mass = np.asarray(evaluation.g_current, dtype=float)
    saving = np.asarray(solution.bp_pol, dtype=float)
    beginning = np.asarray(solution.g_beginning_distribution, dtype=float)
    if (mass.shape != saving.shape or mass.shape != beginning.shape
            or mass.ndim != 7 or mass.shape[0] != len(bg)
            or mass.shape[1] != 1 + int(P.n_house)
            or mass.shape[2] != int(P.I) or mass.shape[3] != int(P.J)
            or mass.shape[4] != int(P.Nz)
            or not all(np.isfinite(x).all() for x in (bg, mass, saving, beginning))
            or np.any(mass < 0) or np.any(beginning < 0)):
        raise ValueError("Renter support audit requires aligned finite native arrays")
    if len(trace) != int(P.J) * int(P.Nz):
        raise ValueError("Renter support trace does not cover exactly one complete Bellman solve")
    if bool(getattr(P, "use_age_survival", False)):
        survival = np.asarray(P.survival_probs, dtype=float)
        if survival.shape != (int(P.J) - 1,):
            raise ValueError("Renter support audit survival shape changed")
    else:
        survival = np.ones(int(P.J) - 1)
    death = np.r_[1.0 - survival, 1.0]
    expected = {(j, z) for j in range(int(P.J)) for z in range(int(P.Nz))}
    seen = set()
    below_floor_mass = negative_estate_mass = lower_grid_mass = boundary_support_mass = 0.0
    for j, z_value, floor_flat, dead_flat in trace:
        z = int(np.argmin(np.abs(np.asarray(P.z_grid, dtype=float) - z_value)))
        if abs(float(P.z_grid[z]) - z_value) > 1e-12 or (j, z) in seen:
            raise ValueError("Renter support trace has duplicate or unknown age/income state")
        seen.add((j, z))
        floor = np.asarray(floor_flat, dtype=float)
        dead = np.asarray(dead_flat, dtype=bool)
        if (floor.shape != (int(P.n_parity) * int(P.n_child_states),)
                or dead.shape != floor.shape or not np.isfinite(floor).all()):
            raise ValueError("Renter support floor has wrong family-state shape")
        floor = floor.reshape((int(P.n_parity), int(P.n_child_states)), order="F")
        dead = dead.reshape(floor.shape, order="F")
        g = mass[:, 0, 0, j, z]
        bp = saving[:, 0, 0, j, z]
        if np.any(g[:, dead] > 0):
            raise RuntimeError("Positive renter mass occupies an all-dead continuation state")
        below_floor_mass += float(g[bp < floor[None, :, :] - 1e-9].sum())
        if death[j] > 0:
            negative_estate_mass += float(g[bp < -1e-10].sum())
        lower_grid_mass += float(g[bp <= bg[0] + 1e-10].sum())
        family_mass = beginning[:, :, 0, j, z].sum(axis=(0, 1))
        boundary_support_mass += float(family_mass[floor <= bg[0] + 1e-10].sum())
    if seen != expected:
        raise ValueError("Renter support trace omitted age/income states")
    result = dict(status="occupied_support_pass_unoccupied_alternatives_unverified",
                  trace_branches=len(trace), grid_lower=float(bg[0]),
                  below_natural_floor_mass=below_floor_mass,
                  negative_renter_death_estate_mass=negative_estate_mass,
                  occupied_lower_grid_mass=lower_grid_mass,
                  beginning_family_mass_at_lower_support=boundary_support_mass,
                  finite_grid_only=True, value_cutoff_approximate=True,
                  natural_support_certified=False)
    if max(below_floor_mass, negative_estate_mass, lower_grid_mass, boundary_support_mass) > 0:
        raise RuntimeError("Renter support or estate gate failed: " + str(result))
    return result
