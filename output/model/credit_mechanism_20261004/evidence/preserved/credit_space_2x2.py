"""Fixed-price credit x family-space experiment at the chain-13 inputs.

Household solve + stationary distribution at the chain-13 GE price; no price
root, no market clearing, no recalibration, no writes to the repository.
Credit: financed share phi 0.80 vs 0.95. Housing setups: baseline (renter cap
6, owner rungs 2-10), renter cap 4, no 2-room owner rung, both.
"""
import os, sys, json, time
for v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
OUT = os.path.dirname(os.path.abspath(__file__))
os.environ["NUMBA_CACHE_DIR"] = os.path.join(OUT, "numba_cache")
sys.dont_write_bytecode = True
import numpy as np

ROOT = "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
sys.path.insert(0, ROOT + "/code/model")
from production.inputs import load_inputs, DEFAULT_PRICE
from production.equilibrium import solve_at_price
from production.engine.parameters import get_fecundity_by_age

OBS = json.load(open(ROOT + "/output/model/local_solution/latest/native/phase_b_ge/selected_root/observers.json"))
OBS_F = OBS["fertility"]["uniform_birth_time"]


def run(phi, hR, drop2):
    P, grid = load_inputs(external_inputs={"phi": [phi] * 4, "hR_max": float(hR)})
    if drop2:
        P.H_own = np.array([4.0, 6.0, 8.0, 10.0]); P.n_house = 4
    t0 = time.time()
    out = solve_at_price(P, grid, DEFAULT_PRICE)
    sol, Q = out["solution"], out["P"]
    secs = time.time() - t0
    J = int(Q.J)
    fec = np.asarray(get_fecundity_by_age(Q), float)
    post = np.asarray(sol.g_beginning_distribution)          # post-fertility, inherited tenure
    gch = np.asarray(sol.g)                                   # post-tenure choice
    F1 = np.asarray(Q._first_births_by_age, float)
    F2 = np.asarray(Q._second_births_by_age, float)
    F3 = np.asarray(Q._third_births_by_age, float)
    # parity stocks by age, post and pre fertility
    pp_post = post.sum(axis=(0, 1, 2, 4, 6))                  # (J, npar)
    pp_pre = pp_post.copy()
    pp_pre[:, 0] += F1; pp_pre[:, 1] += -F1 + F2; pp_pre[:, 2] += -F2 + F3; pp_pre[:, 3] += -F3
    mass = pp_post.sum(1)
    ceb_post = pp_post @ np.arange(4) / mass
    s_pre = pp_pre[1] / pp_pre[1].sum(); s_post = pp_post[1] / pp_post[1].sum()
    sh25 = 0.125 * s_pre + 0.875 * s_post
    age25 = float(sh25 @ np.arange(4))
    # first-birth attempt by income state, pre-fertility childless weights
    fp1 = np.asarray(sol.fert_probs)[..., 1]                  # (Nb, nt, I, J, Nz)
    post0 = post[..., 0, 0]
    pre0 = post0 / np.clip(1.0 - fec[None, None, None, :, None] * fp1, 1e-12, None)
    flow_rec = (fec[None, None, None, :, None] * pre0 * fp1).sum(axis=(0, 1, 2, 4))
    rec_err = float(np.max(np.abs(flow_rec - F1)))
    w0 = pre0.sum(axis=(0, 1, 2))                             # (J, Nz)
    att_z = (pre0 * fp1).sum(axis=(0, 1, 2)) / np.maximum(w0, 1e-300)
    # ownership after the tenure choice
    m_jz = gch.sum(axis=(0, 1, 2, 5, 6)); o_jz = gch[:, 1:].sum(axis=(0, 1, 2, 5, 6))
    own_j = o_jz.sum(1) / m_jz.sum(1)
    own_z = o_jz / np.maximum(m_jz, 1e-300)
    # owner share among parents vs childless, ages 22-33
    mom = gch[..., 1:, :].sum(axis=(0, 2, 4, 5, 6)); cl = gch[..., 0, :].sum(axis=(0, 2, 4, 5))
    fb_share = F1 / F1.sum()
    mid = 18 + 4 * np.arange(J) + 2
    res = dict(
        phi=phi, hR_max=hR, drop2=drop2, seconds=round(secs, 1), flow_reconstruction_err=rec_err,
        age25_ceb=age25, age25_any=float(1 - sh25[0]), age25_two_plus=float(sh25[2:].sum()),
        mean_age_first_birth=float(fb_share @ mid),
        first_births=float(F1.sum()), second_births=float(F2.sum()), third_births=float(F3.sum()),
        explicit_births=float(F1.sum() + F2.sum() + F3.sum()),
        ceb_by_age=ceb_post[:8].tolist(), childless_end=float(pp_post[7, 0] / mass[7]),
        first_birth_flow_by_age=F1[:8].tolist(),
        own_by_age=own_j[:10].tolist(),
        attempt_22_by_z=att_z[1].tolist(), attempt_18_by_z=att_z[0].tolist(),
        own_22_by_z=own_z[1].tolist(), own_26_by_z=own_z[2].tolist(),
        z_mass_22=(w0[1] / w0[1].sum()).tolist(),
        own_parents_by_age=(mom[1:].sum(0) / np.maximum(mom.sum(0), 1e-300))[:8].tolist(),
        own_childless_by_age=(cl[1:].sum(0) / np.maximum(cl.sum(0), 1e-300))[:8].tolist(),
        renewal_residual=float(getattr(sol, "adult_entry_stationary_residual", np.nan)),
    )
    return res


if __name__ == "__main__":
    cases = [(phi, hR, d2) for (hR, d2) in [(6.0, False), (4.0, False), (6.0, True), (4.0, True)] for phi in (0.80, 0.95)]
    results = []
    for phi, hR, d2 in cases:
        r = run(phi, hR, d2)
        results.append(r)
        print(json.dumps({k: r[k] for k in ("phi", "hR_max", "drop2", "seconds", "flow_reconstruction_err", "age25_ceb", "explicit_births")}), flush=True)
        if len(results) == 1:  # validate baseline against the saved chain-13 observer
            m = OBS_F["moments"]
            print("VALIDATION age25 model", r["age25_ceb"], "saved", m["mean_children_ever_born_capped3_age25"],
                  "| mean age", r["mean_age_first_birth"], "saved", m["period_mean_age_first_birth"], flush=True)
        json.dump(results, open(os.path.join(OUT, "credit_space_2x2.json"), "w"), indent=1)
