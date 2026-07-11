"""Decisive identity test for the July 9 timing repair (audit task D).

Builds a TINY one-market markov-income model, solves it once at a fixed price
under both use_postdecision_current_distribution arms (policies are identical
across arms because the flag only changes the KFE observation step), and then
checks, from first principles:

  T1  Age-by-age (and age x n x cs) mass conservation of the post-decision
      cross-section vs the beginning-of-period distribution.
  T2  Contemporaneous-tenure identity: the ownership rate by age computed from
      the post-decision distribution equals the ownership rate implied by
      applying the age-j tenure policy to the age-j beginning-of-period mass
      (I=1 so the location step is trivial). This is the exact statement that
      the moment integrand picks the CURRENT choice, not the inherited state.
  T3  Inherited vs contemporaneous tenure genuinely differ somewhere (the test
      has power): max_j |own_by_age_post - own_by_age_pre| > 0.
  T4  Birth-triggered same-period purchases: among the newborn child state
      (n>=1, cs=1) at fertile ages, owner mass under the post-decision
      distribution vs the inherited-tenure distribution.
  T5  Market clearing and moments use the same object: sol.housing_demand
      (the clearing integrand) recomputed from the RETURNED distribution and
      hR_pol matches exactly.
  T6  Window audit: which model ages enter own_rate_2534 and
      old_age_own_rate_6575.

Read-only with respect to production code; writes nothing outside stdout.
"""

from __future__ import annotations

import numpy as np

from intergen_housing_fertility.parameters import setup_parameters, apply_overrides
from intergen_housing_fertility.solver import (
    age_to_index,
    solve_markov_income_at_prices,
)
from intergen_housing_fertility.utils import make_grid


def tiny_overrides(flag: bool) -> dict:
    return {
        "J": 17,
        "Nb": 28,
        "H_own": np.array([2.0, 6.0]),
        "hR_max": 6.0,
        "use_income_types": True,
        "income_type_transition": "markov",
        "z_grid": np.array([0.8, 1.2]),
        "z_weights": np.array([0.5, 0.5]),
        "income_shock_persistence": 0.85,
        "max_iter_eq": 1,
        "tol_eq": 1e-4,
        "interp_method": "linear",
        "use_postdecision_current_distribution": flag,
        "propagate_birth_entry_grant": True,
        "housing_event_horizon": 0,
        "tenure_choice_kappa": 0.0,  # deterministic argmax -> exact identity
    }


def solve(flag: bool):
    P = apply_overrides(setup_parameters(), tiny_overrides(flag))
    bg = make_grid(P)
    p = np.array([2.0])
    sol = solve_markov_income_at_prices(p, P, bg, verbose=False, fast_stats=False)
    return P, bg, sol


def main() -> None:
    P_post, bg, sol_post = solve(True)
    P_pre, _, sol_pre = solve(False)

    g_post = np.asarray(sol_post.g)   # returned distribution, flag True  -> realized
    g_pre = np.asarray(sol_pre.g)     # returned distribution, flag False -> beginning-of-period
    tc = np.asarray(sol_pre.tenure_choice)  # identical across arms
    assert np.array_equal(tc, np.asarray(sol_post.tenure_choice)), "policies differ across arms!"
    np.testing.assert_allclose(
        np.asarray(sol_pre.hR_pol), np.asarray(sol_post.hR_pol), rtol=0, atol=0,
    )

    J = P_post.J
    nt = 1 + P_post.n_house
    npar = P_post.n_parity
    ncs = P_post.n_child_states
    Nz = len(P_post.z_grid)

    print("=" * 72)
    print("T1  mass conservation, age by age (post-decision vs inherited)")
    max_age_gap = 0.0
    max_cell_gap = 0.0
    for j in range(J):
        m_post = float(np.sum(g_post[:, :, :, j]))
        m_pre = float(np.sum(g_pre[:, :, :, j]))
        max_age_gap = max(max_age_gap, abs(m_post - m_pre))
        for nn in range(npar):
            for cs in range(ncs):
                a = float(np.sum(g_post[:, :, :, j, :, nn, cs]))
                b = float(np.sum(g_pre[:, :, :, j, :, nn, cs]))
                max_cell_gap = max(max_cell_gap, abs(a - b))
    print(f"    max |mass_post - mass_pre| by age          = {max_age_gap:.3e}")
    print(f"    max |mass_post - mass_pre| by (age,n,cs)   = {max_cell_gap:.3e}")
    t1 = max_age_gap < 1e-12 and max_cell_gap < 1e-12
    print(f"    T1 {'PASS' if t1 else 'FAIL'}")

    print("=" * 72)
    print("T2  contemporaneous tenure identity (I=1, deterministic tenure)")
    # Manually apply the age-j tenure policy to the age-j pre-decision mass:
    # owner share after this period's choice.
    max_dev = 0.0
    rows = []
    for j in range(J):
        num = 0.0
        den = 0.0
        for zz in range(Nz):
            for to in range(nt):
                for nn in range(npar):
                    for cs in range(ncs):
                        gs = g_pre[:, to, 0, j, zz, nn, cs]
                        if gs.sum() <= 0.0:
                            continue
                        chosen = tc[:, to, 0, j, zz, nn, cs]
                        num += float(np.sum(gs * (chosen > 0)))
                        den += float(np.sum(gs))
        manual_own = num / max(den, 1e-300)
        post_own = float(
            np.sum(g_post[:, 1:, :, j]) / max(np.sum(g_post[:, :, :, j]), 1e-300)
        )
        pre_own = float(
            np.sum(g_pre[:, 1:, :, j]) / max(np.sum(g_pre[:, :, :, j]), 1e-300)
        )
        dev = abs(manual_own - post_own)
        max_dev = max(max_dev, dev)
        rows.append((j, pre_own, manual_own, post_own, dev))
    print("    j   own_pre(inherited)  own_manual(choice@j)  own_post(reported)  |diff|")
    for j, a, b, c, d in rows:
        print(f"    {j:2d}  {a:18.6f}  {b:20.6f}  {c:18.6f}  {d:.2e}")
    t2 = max_dev < 1e-10
    print(f"    max deviation manual-vs-post = {max_dev:.3e}  -> T2 {'PASS' if t2 else 'FAIL'}")

    print("=" * 72)
    print("T3  power check: inherited vs contemporaneous tenure differ somewhere")
    gap = max(abs(r[1] - r[3]) for r in rows)
    print(f"    max_j |own_pre - own_post| = {gap:.6f}")
    t3 = gap > 1e-4
    print(f"    T3 {'PASS (test has power)' if t3 else 'FAIL (no state where they differ)'}")

    print("=" * 72)
    print("T4  birth-triggered same-period purchases visible at the birth age")
    for j in range(0, min(P_post.A_f_end, J - 1)):
        nb_post = float(np.sum(g_post[:, 1:, :, j, :, 1:, 1]))
        nb_pre = float(np.sum(g_pre[:, 1:, :, j, :, 1:, 1]))
        tot = float(np.sum(g_post[:, :, :, j, :, 1:, 1]))
        if tot > 1e-14:
            print(
                f"    j={j:2d}: newborn-state mass={tot:.5f}  owner share post={nb_post / tot:.4f}"
                f"  inherited={nb_pre / tot:.4f}"
            )
    print("    (post >= inherited at birth ages indicates same-period purchases counted)")

    print("=" * 72)
    print("T5  market clearing integrand == moment integrand (same distribution)")
    hR = np.asarray(sol_post.hR_pol)
    renter_demand = float(np.sum(g_post[:, 0, 0] * hR[:, 0, 0]))
    owner_demand = 0.0
    for ten in range(1, nt):
        owner_demand += float(np.sum(g_post[:, ten, 0])) * float(P_post.H_own[ten - 1])
    manual_hd = renter_demand + owner_demand
    reported_hd = float(np.asarray(sol_post.housing_demand).reshape(-1)[0])
    # housing_demand_normalizer: check both raw and normalized
    from intergen_housing_fertility.solver import housing_demand_normalizer

    norm = housing_demand_normalizer(P_post)
    dev5 = abs(manual_hd / norm - reported_hd)
    print(f"    recomputed demand from returned g = {manual_hd / norm:.10f}")
    print(f"    sol.housing_demand (clearing)     = {reported_hd:.10f}")
    t5 = dev5 < 1e-10
    print(f"    |diff| = {dev5:.3e} -> T5 {'PASS' if t5 else 'FAIL'}")
    agg_rooms = float(sol_post.aggregate_housing_demand) / float(sol_post.total_mass)
    print(f"    aggregate_mean_occupied_rooms(18-85) moment = {agg_rooms:.6f} (same object)")

    print("=" * 72)
    print("T6  age windows actually used by the ownership moments")
    a25 = age_to_index(P_post, 25)
    a34 = age_to_index(P_post, 34)
    a65 = age_to_index(P_post, 65)
    a75 = age_to_index(P_post, 75)
    def ages(j0, j1):
        return [f"[{18 + 4 * j},{18 + 4 * j + 4})" for j in range(j0, j1 + 1)]
    print(f"    own_rate_2534 uses j in [{a25},{a34}] -> age bins {ages(a25, a34)}")
    print(f"    old_age_own_rate_6575 uses j in [{a65},{a75}] -> age bins {ages(a65, a75)}")
    print("    (data windows are 25-34 and 65-75)")

    print("=" * 72)
    ok = t1 and t2 and t3 and t5
    print(f"OVERALL: {'ALL DECISIVE CHECKS PASS' if ok else 'AT LEAST ONE CHECK FAILED'}")


if __name__ == "__main__":
    main()
