"""Local, read-only explorer of authenticated saved household solutions.

Raw conditional policy curves use saved asset-grid coordinates by default.
The alternate view averages the saved tenure lottery after applying its
transaction map. Owner-stayer consumption/saving remain separate from purchase
policies. This server never solves the model or changes a saved solution.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import mimetypes
import time
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import parse_qs, urlparse

import numpy as np

try:
    from .model_policy_tools import aggregate_solution
except ImportError:
    from model_policy_tools import aggregate_solution


class SavedCase:
    def __init__(self, spec, common):
        started = time.perf_counter()
        self.spec, self.common = spec, common
        path = Path(spec['arrays'])
        if hashlib.sha256(path.read_bytes()).hexdigest() != spec['sha256']:
            raise ValueError('Saved array fingerprint differs: ' + spec['id'])
        names = ('b_grid', 'V', 'c_pol', 'hR_pol', 'bp_pol', 'c_pol_stay',
                 'bp_pol_stay', 'tenure_probs', 'fert_probs', 'fert2_probs',
                 'g_beginning_distribution', 'g', 'g_stay_distribution', 'type_values')
        with np.load(path, allow_pickle=False) as arrays:
            self.a = {k: arrays[k] for k in names}
        self.b = self.a['b_grid']
        self.shape = self.a['V'].shape
        assert len(self.shape) == 7 and self.shape[2] == 1
        assert self.a['tenure_probs'].shape == self.shape + (self.shape[1],)
        assert self.shape[0] == len(self.b) and np.all(np.diff(self.b) > 0)
        assert self.shape[1] == 1 + len(common['houses'])
        for k in ('c_pol', 'bp_pol', 'hR_pol', 'g_beginning_distribution'):
            assert self.a[k].shape == self.shape
        mass = self.a['g_beginning_distribution']
        assert np.isfinite(mass).all() and mass.min() >= 0
        self.total = float(mass.sum())
        pooled = mass.sum(axis=tuple(range(1, 7)))
        cdf = np.cumsum(pooled) / self.total
        self.central = (self.b[max(0, np.searchsorted(cdf, .0001)-1)],
                        self.b[min(len(self.b)-1, np.searchsorted(cdf, .995)+1)])
        self.load_seconds = time.perf_counter() - started

    @lru_cache(maxsize=256)
    def slice(self, age, income, tenure, children, at_home, view='central',
              policy_mode='raw', owner_policy='buying', policy_tenure=None,
              asset_display='raw'):
        c, a, b = self.common, self.a, self.b
        age_index = (age-c['age_start']) / c['period_years']
        if not age_index.is_integer():
            raise ValueError('Age is not a model decision date')
        j, z, old, n, m = int(age_index), income, tenure, children, at_home
        policy_old = old if policy_tenure is None else int(policy_tenure)
        for value, size in ((j,self.shape[3]),(z,self.shape[4]),(old,self.shape[1]),
                            (policy_old,self.shape[1]),
                            (n,self.shape[5]),(m,self.shape[6])):
            if not 0 <= value < size:
                raise ValueError('Selector out of range')
        if view not in ('central', 'all'):
            raise ValueError('Unknown wealth range')
        if policy_mode not in ('raw', 'average'):
            raise ValueError('Unknown policy mode')
        if owner_policy not in ('buying', 'staying'):
            raise ValueError('Unknown owner policy branch')
        if asset_display not in ('raw', 'normalized', 'node_index'):
            raise ValueError('Unknown asset display')
        ix = (slice(None),old,0,j,z,n,m)
        probs = np.asarray(a['tenure_probs'][ix], float)
        valid = (a['V'][ix] > -1e9) & (probs.sum(axis=1) > .99)
        if m > n:
            valid[:] = False
        mass = a['g_beginning_distribution'][ix]
        expected_c, expected_b, expected_h = [np.zeros(len(b)) for _ in range(3)]
        costs = np.r_[0., float(self.spec['price'])*np.asarray(c['houses'])]
        sales = (1-c['selling_cost'])*costs
        divisor = c['R_gross'] if self.spec['timing'] == 'inherited_only' else 1.
        if policy_mode == 'average':
            for new in range(self.shape[1]):
                x = b if old == new else b+(sales[old]-costs[new])/divisor
                dest = (slice(None),new,0,j,z,n,m)
                staying = old == new and old > 0
                cc = a['c_pol_stay' if staying else 'c_pol'][dest]
                bp = a['bp_pol_stay' if staying else 'bp_pol'][dest]
                hh = a['hR_pol'][dest] if new == 0 else np.full(len(b),c['houses'][new-1])
                expected_c += probs[:,new]*np.interp(x,b,cc)
                expected_b += probs[:,new]*np.interp(x,b,bp)
                expected_h += probs[:,new]*np.interp(x,b,hh)
        raw_valid = np.ones(len(b), dtype=bool)
        if policy_mode == 'raw':
            branch = 'c_pol_stay' if policy_old > 0 and owner_policy == 'staying' else 'c_pol'
            saving_branch = 'bp_pol_stay' if policy_old > 0 and owner_policy == 'staying' else 'bp_pol'
            branch_ix = (slice(None),policy_old,0,j,z,n,m)
            policy_c = a[branch][branch_ix]
            policy_b = a[saving_branch][branch_ix]
            policy_h = (np.full(len(b), c['houses'][policy_old-1]) if policy_old > 0
                        else a['hR_pol'][branch_ix])
            if asset_display == 'node_index':
                policy_x = np.arange(1, len(b)+1, dtype=float)
                policy_xlabel = 'Asset-grid node index (1-based)'
            else:
                policy_x = b.copy()
                policy_xlabel = ('Assets b (model units)' if asset_display == 'raw'
                                 else 'Assets b (already normalized by mean annual earnings)')
            policy_name = ('renter' if policy_old == 0 else
                           f'owner with {c["houses"][policy_old-1]:g} rooms')
            branch_name = ('staying' if policy_old > 0 and owner_policy == 'staying' else 'buying')
            policy_condition = (f'Raw saved policy for {policy_name}' +
                (f', {branch_name} branch' if policy_old > 0 else '') +
                f'; age {age:g}, income state {z+1}, n={n} children ever born, m={m} currently at home.')
            display_note = ('The normalized earnings display reuses the same saved b values; it does not rescale them again.'
                            if asset_display == 'normalized' else '')
            policy_note = (f'The first three charts show raw conditional policy arrays for {policy_name} at saved asset-grid '
                'coordinates. The inherited tenure selector still controls the bottom three charts. Raw transaction branches '
                'can use a coordinate that is not transformed physical wealth. The arrays include initialized entries for '
                'infeasible states; zero entries are not certified choices. For owners, the housing curve repeats the selected '
                'fixed house size and is not a solver choice at infeasible nodes. The bottom charts show inherited-state '
                'objects: ownership and birth-attempt probabilities are before the tenure choice, and population mass is after '
                f'fertility and before tenure. {display_note}')
            raw_markers = True
        else:
            policy_c, policy_b, policy_h = expected_c, expected_b, expected_h
            policy_name = ('renter' if old == 0 else f'owner with {c["houses"][old-1]:g} rooms')
            policy_condition = (f'Tenure-choice lottery average conditional on inherited {policy_name}; age {age:g}, '
                f'income state {z+1}, n={n} children ever born, m={m} currently at home.')
            policy_x = (np.arange(1, len(b)+1, dtype=float) if asset_display == 'node_index'
                        else b.copy())
            policy_xlabel = ('Asset-grid node index (1-based)' if asset_display == 'node_index' else
                'Assets b (model units)' if asset_display == 'raw' else
                'Assets b (already normalized by mean annual earnings)')
            policy_note = ('Consumption, saving, housing and ownership condition on the family state after fertility and average over tenure choices. '
                'The fertility curve is the attempt probability for that family state before the birth decision. '
                'Mass is after fertility, before tenure; zero-mass states are hypothetical. The top child count is capped at 3.')
            raw_markers = False
        if n == 0:
            attempt = a['fert_probs'][:,old,0,j,z,1]
        elif n < self.shape[5]-1 and m <= n:
            attempt = a['fert2_probs'][:,old,0,j,z,1,n-1,m]
        else:
            attempt = np.zeros(len(b))
        keep = np.ones(len(b),bool) if view == 'all' else ((b >= self.central[0]) & (b <= self.central[1]))

        def values(x):
            return [float(v) if np.isfinite(v) else None for v in np.asarray(x)[keep]]

        inherited_x = np.arange(1, len(b)+1, dtype=float) if asset_display == 'node_index' else b.copy()
        inherited_xlabel = ('Asset-grid node index (1-based)' if asset_display == 'node_index' else
            'Assets b (model units)' if asset_display == 'raw' else
            'Assets b (already normalized by mean annual earnings)')
        inherited_name = ('renter' if old == 0 else f'owner with {c["houses"][old-1]:g} rooms')
        inherited_condition = (f'Inherited tenure: {inherited_name}; age {age:g}, income state {z+1}, '
            f'n={n} children ever born, m={m} currently at home.')

        def chart(title, unit, y, mask=valid, markers=False, x=None, xlabel=None, conditioning=None):
            return dict(title=title,ylabel=unit,conditioning=conditioning or inherited_condition,
                valid=np.asarray(mask)[keep].tolist(),markers=markers,
                x=values(inherited_x if x is None else x), xlabel=xlabel or inherited_xlabel,
                series=[dict(label=title,values=values(y))])

        return dict(label=f"{self.spec['label']} · age {age:g} · income state {z+1} · inherited {inherited_name} · {n} children ever born, {m} currently at home",
            note=policy_note,
            slice_mass=float(mass.sum()/self.total),
            xlabel=inherited_xlabel,
            x=values(inherited_x), valid=valid[keep].tolist(), charts=[
                chart('Raw conditional consumption' if policy_mode == 'raw' else 'Consumption', 'Four-year consumption / mean annual earnings', policy_c, raw_valid if policy_mode == 'raw' else valid, raw_markers, policy_x, policy_xlabel, policy_condition),
                chart('Raw conditional next wealth' if policy_mode == 'raw' else 'Next-period financial wealth', 'Wealth / mean annual earnings', policy_b, raw_valid if policy_mode == 'raw' else valid, raw_markers, policy_x, policy_xlabel, policy_condition),
                chart('Raw conditional housing' if policy_mode == 'raw' else 'Housing services', 'Rooms', policy_h, raw_valid if policy_mode == 'raw' else valid, raw_markers, policy_x, policy_xlabel, policy_condition),
                chart('Probability of owning', 'Probability', probs[:,1:].sum(axis=1), conditioning=inherited_condition+' Ownership lottery before the tenure choice.'),
                chart('Birth attempt probability', 'Probability', attempt, conditioning=inherited_condition+' Attempt probability before the birth decision.'),
                chart('Household wealth distribution', 'Mass at node (% of total population)', 100*mass/self.total, conditioning=inherited_condition+' Population mass after fertility and before tenure.')])

    @lru_cache(maxsize=8)
    def aggregates(self):
        sol = SimpleNamespace(**self.a)
        return aggregate_solution(sol, houses=self.common['houses'],
            age_start=self.common['age_start'], period_years=self.common['period_years'])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', required=True, type=Path)
    ap.add_argument('--port', type=int, default=8765)
    args = ap.parse_args()
    config_path = args.config.resolve(strict=True)
    config_bytes = config_path.read_bytes()
    config = json.loads(config_bytes)
    cases = {x['id']: SavedCase(x,config['common']) for x in config['cases']}
    first = next(iter(cases.values()))
    c = config['common']
    meta = dict(title='Housing and fertility — model explorer',
        config_path=str(config_path),
        config_sha256=hashlib.sha256(config_bytes).hexdigest(),
        note='Saved solutions. Changing selectors does not solve or recalibrate the model.',
        links=([dict(label='All moments, parameters and diagnostics',url='/report/readout.html')]
               if 'report_root' in config and (Path(config['report_root'])/'readout.html').is_file() else []),
        cases=[dict(id=k,label=v.spec['label']) for k,v in cases.items()],
        ages=[c['age_start']+c['period_years']*j for j in range(first.shape[3])],
        incomes=[dict(index=i,label=f'{i+1}: productivity {z:.3g}') for i,z in enumerate(first.a['type_values'])],
        tenures=[dict(index=0,label='Renter')]+[dict(index=i+1,label=f'Owner: {h:g} rooms') for i,h in enumerate(c['houses'])],
        children=list(range(first.shape[5])),at_home=list(range(first.shape[6])),
        default=dict(age=30,income=first.shape[4]//2,tenure=0,policy_tenure=0,children=0,at_home=0))

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *_):
            pass

        def do_GET(self):
            url = urlparse(self.path)
            try:
                if url.path == '/':
                    body = Path(__file__).with_suffix('.html').read_bytes()
                    kind = 'text/html; charset=utf-8'
                elif url.path == '/api/meta':
                    body = json.dumps(meta,allow_nan=False).encode(); kind = 'application/json'
                elif url.path == '/api/slice':
                    q = {k:v[0] for k,v in parse_qs(url.query).items()}
                    result = cases[q['case']].slice(float(q['age']),int(q['income']),int(q['tenure']),
                        int(q['children']),int(q['at_home']),q.get('view','central'),
                        q.get('policy_mode','raw'),q.get('owner_policy','buying'),
                        int(q['policy_tenure']) if 'policy_tenure' in q else None,
                        q.get('asset_display','raw'))
                    body = json.dumps(result,allow_nan=False).encode(); kind = 'application/json'
                elif url.path == '/api/aggregates':
                    q = {k:v[0] for k,v in parse_qs(url.query).items()}
                    result = cases[q['case']].aggregates()
                    body = json.dumps(result,allow_nan=False).encode(); kind = 'application/json'
                elif url.path.startswith('/report/') and 'report_root' in config:
                    root = Path(config['report_root']).resolve()
                    path = (root / url.path.removeprefix('/report/')).resolve()
                    if not path.is_relative_to(root) or not path.is_file() or path.suffix not in ('.html','.png','.csv','.json'):
                        self.send_error(404); return
                    body = path.read_bytes()
                    kind = mimetypes.guess_type(path.name)[0] or 'application/octet-stream'
                else:
                    self.send_error(404); return
                self.send_response(200)
                self.send_header('Content-Type',kind)
                self.send_header('Content-Length',str(len(body)))
                self.send_header('Cache-Control','no-store')
                self.end_headers(); self.wfile.write(body)
            except (ValueError,KeyError,IndexError) as exc:
                self.send_error(400,str(exc))
            except (BrokenPipeError,ConnectionResetError):
                pass

    server = ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    print(json.dumps(dict(url=f'http://127.0.0.1:{args.port}',
        config_path=meta['config_path'],config_sha256=meta['config_sha256'],
        loaded_seconds={k:v.load_seconds for k,v in cases.items()},native_solves=0)),flush=True)
    server.serve_forever()


if __name__ == '__main__':
    main()
