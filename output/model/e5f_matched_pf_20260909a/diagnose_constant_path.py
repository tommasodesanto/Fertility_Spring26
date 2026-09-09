"""Two rental-arithmetic replays of the saved original-grid stationary policy."""
import argparse, copy, gzip, json, os, pickle, sys, time
from pathlib import Path
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
p = argparse.ArgumentParser()
p.add_argument('--arm', choices=['sequential','nested'], required=True)
p.add_argument('--contract', type=Path, required=True)
p.add_argument('--contract-sha256', required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
root = Path(json.loads(a.contract.read_text())['source_root'])
sys.path[:0] = [str(root/'code/model'), str(root/'code/model/tools')]
import numpy as np
import run_e5f_matched_pf_smoke as smoke
c = smoke.load_contract(a.contract, a.contract_sha256)
smoke.transition.configure_sequential_model()
smoke.calendar.model = smoke.model
a.output.mkdir(parents=True, exist_ok=False)
with gzip.open(c['checkpoint'],'rb') as f:
    packet = pickle.load(f)
P = copy.deepcopy(packet['parameters'])
P.joint_nested_choice = P.fertility_nest_choice = a.arm == 'nested'
P.two_shock_choice = False
P.exhaustive_saving_control = True
P.property_tax_lump_sum_transfer = 0.
saved = np.load(root/'output/model/e5f_matched_pf_20260909a/pilot_03'/a.arm/'stationary_arrays.npz')
grid = saved['wealth_grid']; V = saved['V']; occupied = saved['g_pre'] > 1e-12
price = float(packet['evaluation'].policy.price[0])
rent_pf = float(smoke.pf.rents_from_asset_prices([price],price,P)[0])
rent_ss = float(P.user_cost_rate * price)
result = {'arm':a.arm, 'rents':{'pf':rent_pf, 'stationary':rent_ss, 'difference':rent_pf-rent_ss},
          'scope':'fixed-price Bellman diagnostic; no changed model, grid, or gates', 'replays':[]}
for label, rent in [('pf_arithmetic',rent_pf), ('stationary_arithmetic',rent_ss)]:
    start=time.monotonic();params=copy.deepcopy(P)
    shared=smoke.model.precompute_shared(params,grid)
    policy=smoke.pf.solve_date_policy(price=price,rent=rent,P=params,b_grid=grid,
                                    shared=shared,continuation_V=V)
    fields={}
    for name,actual in smoke.policy_arrays(policy).items():
        expected=saved[name]; difference=np.abs(actual-expected)
        row={'max_absolute':float(difference.max()),'above_2e10':int((difference>2e-10).sum()),
             'exact_equal':bool(np.array_equal(actual,expected))}
        if name=='V':
            index=tuple(int(i) for i in np.unravel_index(difference.argmax(),difference.shape))
            row.update(worst_index=index,expected=float(expected[index]),actual=float(actual[index]),
                       occupied_max=float(difference[occupied].max()),
                       occupied_above_2e10=int(((difference>2e-10)&occupied).sum()))
        fields[name]=row
    result['replays'].append({'label':label,'seconds':time.monotonic()-start,'fields':fields})
    smoke.pf.write_json(a.output/'diagnosis.json',result)
    print(label,json.dumps(fields['V']),flush=True)
