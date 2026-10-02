"""Isolated continuation for prematurely stopped 0--47 search chains.

Only a candidate that exhausts its *own* 32-solve price search while ample
global time remains is censored. The optimizer receives +inf, but no model
loss is saved and the candidate can never become the selected fit.
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import os
import sys
import time
from pathlib import Path

PACKET = Path(os.environ.get('RESTART_PACKET_ROOT', Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(PACKET))


def read(path):
    return json.loads(Path(path).read_text())


def require_launcher_scaffold(out):
    """Accept only the empty directories/log created by launch_torch.sh."""
    out = Path(out)
    assert out.is_dir() and not out.is_symlink(), 'Missing launcher result scaffold'
    assert {p.name for p in out.iterdir()} == {'numba_cache', 'matplotlib', 'search.log'}, 'Existing scientific results or unexpected launcher scaffold'
    for name in ('numba_cache', 'matplotlib'):
        folder = out / name
        assert folder.is_dir() and not folder.is_symlink() and not any(folder.iterdir()), 'Nonempty or invalid launcher cache'
    log = out / 'search.log'
    assert log.is_file() and not log.is_symlink(), 'Missing launcher search log'


def original_contract(parent, chain):
    parent = Path(parent)
    launch = read(parent / 'launcher_start.json')
    search = read(parent / 'search/search_completed.json')
    assert launch['chain'] == chain
    assert launch['deadline_epoch'] == launch['start_epoch'] + 14400
    assert search['search_stop_reason'] == 'native_evaluation_budget_exhausted'
    assert 0 < search['objective_calls'] < 250
    assert search['completed_full_ge'] <= search['objective_calls']
    assert (parent / 'search/input_contract.json').is_file()
    assert (parent / 'search/cases.json').is_file()
    cases = read(parent / 'search/cases.json')
    assert len(cases) == search['completed_full_ge']
    assert cases[-1]['status'] == 'budget_exhausted'
    assert cases[-1]['reason'] == 'uncomputed_bounded_budget'
    selected = search.get('selected')
    if selected is not None:
        assert selected['status'] == 'passed' and selected['label'] in {c['label'] for c in cases}
        assert selected['loss'] >= 0
        postcheck = read(parent / 'postcheck/completed.json')
        assert postcheck['status'] == 'selected_numerically_verified'
        assert postcheck['selected']['label'] == selected['label']
        assert postcheck['selected']['loss'] == selected['loss']
    return launch, search


def install_censored_optimizer(module, remaining_calls):
    source = inspect.getsource(module.optimize)
    old = """        if key in cache:
            write(out/'latest.json',dict(status='exact_vector_cache_hit',objective_calls=calls,completed_full_ge=len(cases),loss=cache[key]['objective'],parameters=point))
            return cache[key]['objective']"""
    new = """        if key in cache:
            previous=cache[key]
            write(out/'latest.json',dict(status='exact_vector_cache_hit',objective_calls=calls,completed_full_ge=len(cases),loss=previous['objective'],parameters=point,optimizer_only_censor=previous.get('optimizer_only_censor',False)))
            return float('inf') if previous.get('optimizer_only_censor') else previous['objective']"""
    assert source.count(old) == 1
    source = source.replace(old, new)
    old = """        elif result['status']=='budget_exhausted':
            cases.append(row);write(out/'cases.json',cases)
            raise BudgetStop('native_evaluation_budget_exhausted')"""
    new = """        elif result['status']=='budget_exhausted':
            # Native can return the same status for a candidate cap or a global
            # clock stop. Censor only the exact candidate-cap signature, with
            # the original 900s final reserve plus native 700s gate intact.
            candidate_cap=(row.get('reason')=='uncomputed_bounded_budget'
                           and row.get('lifecycle_solves',0)>=31
                           and time.time()+RESERVE+700 < deadline)
            if candidate_cap:
                row.update(objective=None,loss=None,computed_valid_loss=False,
                           optimizer_only_censor=True,censor_reason='per_candidate_price_root_cap_32',
                           model_moment_status='uncomputed')
                new_best=False
                write(out/label/'candidate_censored.json',row)
            else:
                row.update(objective=None,loss=None,computed_valid_loss=False,
                           optimizer_only_censor=False)
                cases.append(row);write(out/'cases.json',cases)
                raise BudgetStop('native_or_global_evaluation_budget_exhausted')"""
    assert source.count(old) == 1
    source = source.replace(old, new)
    old = """        return row['objective']"""
    new = """        return float('inf') if row.get('optimizer_only_censor') else row['objective']"""
    assert source.count(old) == 1
    source = source.replace(old, new)
    namespace = module.__dict__
    exec(compile(source, str(Path(__file__)), 'exec'), namespace)
    patched = namespace['optimize']
    def bounded(*args, **kwargs):
        kwargs['maxeval'] = remaining_calls
        return patched(*args, **kwargs)
    module.optimize = bounded
    return hashlib.sha256(source.encode()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--chain', type=int, choices=range(48), required=True)
    ap.add_argument('--parent', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    launch, search = original_contract(args.parent, args.chain)
    remaining = 250 - search['objective_calls']
    if time.time() >= launch['deadline_epoch'] - 900:
        raise SystemExit('Original four-hour deadline leaves no search time; no restart')
    contract = dict(status='prepared',chain=args.chain,parent=str(args.parent),
                    parent_search_receipt_sha256=hashlib.sha256((args.parent/'search/search_completed.json').read_bytes()).hexdigest(),
                    parent_launcher_receipt_sha256=hashlib.sha256((args.parent/'launcher_start.json').read_bytes()).hexdigest(),
                    parent_postcheck_receipt_sha256=hashlib.sha256((args.parent/'postcheck/completed.json').read_bytes()).hexdigest() if search.get('selected') else None,
                    original_start_epoch=launch['start_epoch'],original_deadline_epoch=launch['deadline_epoch'],
                    original_objective_calls=search['objective_calls'],remaining_objective_calls=remaining,
                    original_selected_loss=search.get('selected',{}).get('loss') if search.get('selected') else None,
                    original_target_contract_sha256=read(args.parent/'search/input_contract.json')['parameter_contract']['target_contract_sha256'],
                    no_clock_reset=True,no_call_count_reset=True)
    if args.dry_run:
        print(json.dumps(contract,indent=2));return
    require_launcher_scaffold(args.out)
    (args.out/'restart_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    # Import the immutable original controller with --chain present so its
    # isolated engine is loaded before the frozen integration.
    sys.argv=['run_psi.py','--chain',str(args.chain),'--out',str(args.out/'search'),
              '--deadline-epoch',str(launch['deadline_epoch']),'--fast-objective']
    import run_psi as original
    selected=search.get('selected')
    if selected is not None:
        original.CONFIG['nearby_starts'][args.chain]['parameters']=selected['parameters']
    contract['optimizer_source_sha256']=install_censored_optimizer(original,remaining)
    (args.out/'restart_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    original.main()
    new=read(args.out/'search/search_completed.json')
    assert search['objective_calls']+new['objective_calls']<=250
    assert new['selected'] is None or new['selected']['status']=='passed'
    winner='restart' if new['selected'] is not None and (selected is None or new['selected']['loss']<selected['loss']) else 'parent'
    receipt=dict(status='restart_search_finished',chain=args.chain,winner=winner,
                 cumulative_objective_calls=search['objective_calls']+new['objective_calls'],
                 original_deadline_epoch=launch['deadline_epoch'],parent_selected_loss=selected['loss'] if selected else None,
                 restart_selected_loss=new['selected']['loss'] if new['selected'] else None,
                 original_search_stop=search['search_stop_reason'],restart_search_stop=new['search_stop_reason'],
                 selected_requires_fresh_postcheck=winner=='restart')
    (args.out/'restart_summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))


if __name__=='__main__': main()
