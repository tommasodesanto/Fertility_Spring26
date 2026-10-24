"""Section 6.3 numbers from the free-entry eta 0.63 runs (Oct 24 2026): who uses the credit. Runs on Torch next to the saved solutions
(pol_063fi/out/runs/P063fi_<ARM>/solution_best.npz: full distributions for dates 0-2 = 2023, 2027, 2031; read only).
Distribution axes (b, tenure, loc, age, z, n, cs); tenure 0 = rent, 1..5 = own H_own {2,4,6,8,10}; age cell j = ages 18+4j at the date
(as stats.own_18_29 = cells 0:3 in driver_policy_v10_rc.py and cohort_capture.py). b_grid = the newacct engine grid (120 nodes),
passed as b_grid_newacct.npy. Prices from the run's result.json rows.
  (i)   2023 after the tenure choice (g_current_000), ages 22-29 (cells 1:3): share renting and owning each size.
  (ii)  start of 2027 (g_pre_001), ages 22-25 (cell 1): ownership and two-room ownership.
  (iii) start of 2031 (g_pre_002), ages 22-37 (cells 1:5): ownership, owners' mean rooms, owners' mortgage debt / home value
        (debt = max(-b, 0) at the start of the period, value = P_t * h).
usage: python s6_owners_fi.py RUNS_DIR B_GRID.npy OUT.json"""
import json, sys, zipfile
import numpy as np
RUNS, BG, OUTF = sys.argv[1], sys.argv[2], sys.argv[3]
b = np.load(BG); H = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
def arr(z, name):
    with z.open(name + '.npy') as f: return np.lib.format.read_array(f)
res = {}
for arm in ('BASE', 'PHI100', 'SUPMC', 'PPHI100'):
    try:
        z = zipfile.ZipFile(f'{RUNS}/P063fi_{arm}/solution_best.npz')
    except FileNotFoundError:
        continue
    rows = json.load(open(f'{RUNS}/P063fi_{arm}/result.json'))['last']['rows']; P = [r['asset_price'] for r in rows]
    gc0 = arr(z, 'g_current_000'); gp1 = arr(z, 'g_pre_001'); gp2 = arr(z, 'g_pre_002')
    assert gc0.shape[0] == b.size and gc0.shape[1] == 6, gc0.shape
    by_ten = lambda g, js: g[:, :, :, js].sum(axis=tuple(i for i in range(g.ndim) if i != 1))
    t0 = by_ten(gc0, slice(1, 3)); t0 = t0 / t0.sum()
    t1 = by_ten(gp1, slice(1, 2)); t1 = t1 / t1.sum()
    g2 = gp2[:, :, :, 1:5]; m_bt = g2.sum(axis=(2, 3, 4, 5, 6))          # (b, tenure)
    own = m_bt[:, 1:]; own_mass = own.sum()
    rooms = float((own.sum(axis=0) * H).sum() / own_mass)
    debt = float((own.sum(axis=1) * np.maximum(-b, 0.0)).sum()); value = float(P[2] * (own.sum(axis=0) * H).sum())
    res[arm] = dict(tenure_22_29_2023=t0.tolist(), own_22_29_2023=float(t0[1:].sum()),
                    own_22_25_2027=float(t1[1:].sum()), own2room_22_25_2027=float(t1[1]),
                    own_22_37_2031=float(own_mass / m_bt.sum()), owner_rooms_22_37_2031=rooms, debt_value_22_37_2031=debt / value,
                    price=P[:3])
json.dump(res, open(OUTF, 'w'), indent=1); print(json.dumps(res, indent=1))
