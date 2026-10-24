"""Long-run credit test for Section 6.2 (Oct 24 2026): does removing the down payment raise births at the long-run price, and how
much less at a higher price or with a larger rental cap? Stationary household solve at a fixed price, newacct point, rebate balanced
(../20261021/scripts/mf_common.py), with the long-run preference of the free-entry history: psi_inf = 0.148414 (last date of
data_fi_torch/P063fi_BASE.json); long-run price 0.407567 (its stationary endpoint).
  python lr_credit_cell.py PRICE PHI HRMAX   -> out/p{PRICE}_phi{PHI}_cap{HRMAX}.json
Reports completed fertility 2.1*(1 + renewal residual) and the renewal residual (births relative to replacement)."""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '20261021', 'scripts'))
import mf_common as M
PSI_INF = 0.148414
price, phi, cap = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
f = os.path.join(HERE, 'out', f'p{price:.6f}_phi{phi:.2f}_cap{cap:g}.json')
if os.path.exists(f): print(open(f).read()); sys.exit(0)
t0 = time.time()
P, grid = M.build_P(psi=PSI_INF, phi=phi)
P.hR_max = cap
r = M.solve_balanced(P, grid, price, start=float(P.property_tax_lump_sum_transfer))
st = M.renewal_and_demand(r['P'], r['solution'])
rec = dict(price=price, phi=phi, hR_max=cap, psi=PSI_INF, cf=2.1 * (1 + st['renewal_residual']), transfer=float(r['fiscal']['transfer']),
           seconds=time.time() - t0, **st)
open(f, 'w').write(json.dumps(rec, indent=1) + '\n'); print(json.dumps(rec), flush=True)
