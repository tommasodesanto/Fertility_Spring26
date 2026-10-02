"""Tiny owner-choice contract check; no equilibrium or model solve."""
from types import SimpleNamespace

import numpy as np

from refactor_lab.engine.household import native_due_death_floor
from refactor_lab.engine.kernels import full_owner_block_kernel


def owner_policy(phi, death_floor):
    grid = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    n = grid.size
    z = np.zeros(1)
    value, saving, consumption = full_owner_block_kernel(
        np.full(n, 5.0), np.full(n, 5.0), np.zeros((n, 1)),
        np.zeros((n, 1)), 0, grid, z, z, z, z,
        np.array([0.5]), np.ones(1), np.array([-phi]),
        0.0, 1.0, 0.0, 1.0, 0.01, 0.5, -1.0, 0.95,
        0.0, 0.0, 0.382, 0.618, 1e-5,
        buyer_death_floor=death_floor,
    )
    assert np.isfinite(value).all() and np.isfinite(consumption).all()
    assert np.min(saving) >= max(-phi, death_floor) - 1e-10
    return saving


def main():
    P = SimpleNamespace(J=2, psi=0.1, use_age_survival=False)
    death = native_due_death_floor(P, 1, 1.0, 1.0)
    assert death == -0.9
    hard = owner_policy(1.0, death)
    conventional = owner_policy(0.8, death)
    assert np.min(hard) >= -0.9 - 1e-10
    assert np.min(conventional) >= -0.8 - 1e-10
    # At phi=.8 and psi=.1, -phi*pH > -(1-psi)*pH.
    assert -0.8 > death
    print("buyer floor fixture PASS; hard100 estate floor and phi=.8 redundancy")


if __name__ == "__main__":
    main()
