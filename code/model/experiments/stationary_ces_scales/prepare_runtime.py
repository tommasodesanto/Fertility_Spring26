"""Prepare isolated reviewable CES runtime; never submit or solve."""
from __future__ import annotations
import argparse, hashlib, json, shutil
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
PACKET=ROOT/'output/model/fixed_reference_economics_20260928/utility_ces_scales_v1'
OLD=ROOT/'output/model/fixed_reference_economics_20260928/utility_floor_round2_v1'
NM=ROOT/'output/model/fixed_reference_economics_20260928/utility_floor_nm_v1'
PSI=ROOT/'output/model/fixed_reference_economics_20260928/utility_floor_psi_v1'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def replace(source,old,new):
    if source.count(old)!=1:raise RuntimeError('Ambiguous source boundary '+old[:90])
    return source.replace(old,new,1)
def write(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--refresh',action='store_true');a=ap.parse_args()
    runtime=PACKET/'runtime'
    if runtime.exists():
        if not a.refresh:raise RuntimeError('Runtime already exists; explicit refresh required')
        shutil.rmtree(runtime)
    runtime.mkdir(parents=True);tree=runtime/'source';adapter=runtime/'utility_adapter';adapter.mkdir()
    inherited=json.loads((PSI/'source_pins.json').read_text());copied={}
    for rel,digest in inherited.items():
        src=ROOT/rel
        if sha(src)!=digest:raise RuntimeError('Inherited dependency drift '+rel)
        if src.suffix in ('.py','.json','.csv','.sh') and src.stat().st_size<2_000_000:
            dest=tree/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest);copied[rel]=digest
    matched=tree/'output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed'
    base=tree/'output/model/publication_refactor_20260929/grid_resolution_v1/credit053_v2/runner'
    for src in (OLD/'runner.py',OLD/'inputs.py',OLD/'phase_b_pilot.py',OLD/'plan.json'):
        shutil.copy2(src,adapter/src.name)
    src=(adapter/'inputs.py').read_text()
    # Retain exact pinned bundle/data paths, replace only experimental binding.
    src=replace(src,'ROOT=HERE.parents[3]',f'ROOT=Path({str(ROOT)!r})')
    src=replace(src,"Q.child_room_floor=arm=='floor'","Q.ces_enabled=True;Q.exhaustive_saving_control=True;Q.ces_eta=.487;Q.lambda_housing=float(point['lambda_housing']);Q.alpha_cons=.733\n    Q.child_room_floor=False")
    (adapter/'inputs.py').write_text(src)
    src=(adapter/'runner.py').read_text();src=replace(src,'ROOT=HERE.parents[3]',f'ROOT=Path({str(ROOT)!r})')
    src=replace(src,"BASE=ROOT/'output/model/publication_refactor_20260929/grid_resolution_v1/credit053_v2/runner'",f'BASE=Path({str(base)!r})')
    src=replace(src,"    x.update(point,wealth_grid_nodes", "    point={k:v for k,v in point.items() if k!='lambda_housing'}\n    x.update(point,wealth_grid_nodes")
    src=replace(src,"h_P=float(point['h_P']) if arm=='floor' else 0.","h_P=0.")
    src=replace(src,"delta_alpha_jump=float(point['delta_alpha_jump']) if arm=='no_A' else 0.","delta_alpha_jump=0.")
    start=src.index('def verify_sources():');end=src.index('\n\ndef compare_repeated',start)
    src=src[:start]+'''def verify_sources():
    packet=HERE.parents[1]
    for rel,digest in json.loads((packet/'source_pins.json').read_text()).items():
        require(inputs.sha(packet/rel)==digest,'Private dependency drift: '+rel)
'''+src[end:]
    src=replace(src,"    require(float(live_sd.h_bar[1,1])==expected['h_P'] and float(live_sd.h_bar[1,0])==0.,'Executed physical floor differs')",'''    require(np.all(live_sd.h_bar==0.),'CES must have no physical floor')
    for key in ('ces_ec_flat','ces_eh_flat'):
        np.testing.assert_array_equal(getattr(live_sd,key),getattr(checked_sd,key))
    require(P.ces_enabled and P.ces_eta==.487 and P.alpha_cons==.733 and P.lambda_housing==seed['lambda_housing'],'CES binding drift')
    write(out/'ces_binding.json',dict(eta=P.ces_eta,alpha=P.alpha_cons,lambda_housing=P.lambda_housing,consumption_scales=live_sd.ces_ec_flat.tolist(),housing_scales=live_sd.ces_eh_flat.tolist(),floor_zero=True,lifecycle_solves=0))''')
    src=replace(src, 'def native_evaluator(out,lane,P,grid,deadline,price_start=None):', '_CES_AUTH_CONTEXT = None\n\ndef native_evaluator(out,lane,P,grid,deadline,price_start=None):\n    global _CES_AUTH_CONTEXT')
    src=replace(src, '    base.authored.authenticate_frozen(ctx)', '''    if _CES_AUTH_CONTEXT is None:
        base.authored.authenticate_frozen(ctx)
        _CES_AUTH_CONTEXT = dict(ctx)
    else:
        ctx.update({k:_CES_AUTH_CONTEXT[k] for k in ('fp','prepared','manifest','objective','runtime','reference')})
        write(out/'reused_authenticated_context.json',dict(authenticated_once=True,reused_solutions=False,reused_price_results=False,lifecycle_solves=0))''')

    src=replace(src, "    actual=ctx['fp'].actual_parameters(ctx['prepared'],P,grid)", '''    actual=ctx['fp'].actual_parameters(ctx['prepared'],P,grid)
    actual_free={k:(float(P.lambda_housing) if k=='lambda_housing' else float(actual[k])) for k in coordinates}
    intended_free={k:(float(P.beta)**(1./float(P.period_years)) if k=='beta_annual' else float(seed[k])) for k in coordinates}
    require(actual_free==intended_free,'All ten effective coordinates must match bound input')
    require(target_identity(ctx['manifest']['full_target_table'])==PLAN['target_contract'],'Authenticated base target weights differ')
    write(out/'effective_contract_identity.json',dict(effective_free=actual_free,effective_free_sha256=inputs.canonical(actual_free),target_weight_sha256=inputs.canonical(PLAN['target_contract']),CES_eta=float(P.ces_eta),CES_enabled=bool(P.ces_enabled),alpha=float(P.alpha_cons),lifecycle_solves=0))''')
    src=replace(src, "    write(out/'native_initializer_verification.json'", '''    executed_dir=Path(executed_kernels.__file__).parent
    checked_dir=Path(checked_kernels.__file__).parent
    engine_hashes={p.name:inputs.sha(p) for p in checked_dir.glob('*.py')}
    require(engine_hashes=={p.name:inputs.sha(p) for p in executed_dir.glob('*.py')},'All executed engine files must match CES clone')
    write(out/'executed_engine_identity.json',dict(executed_engine_path=str(executed_dir),checked_engine_path=str(checked_dir),engine_hashes=engine_hashes,full_indexed_kernel_required=True,exhaustive_saving_control=bool(P.exhaustive_saving_control),ces_enabled=bool(P.ces_enabled),lifecycle_solves=0))
    write(out/'native_initializer_verification.json' ''')
    # Metadata for lambda belongs in a separate CES receipt; native 31-row identity remains exact.
    src=replace(src,"        actual=candidate['fp'].actual_parameters(candidate['prepared'],candidate['P'],grid)",'''        require(candidate['P'].ces_enabled and candidate['P'].ces_eta==.487 and candidate['P'].lambda_housing==point['lambda_housing'],'CES candidate binding drift')
        write(directory/'ces_parameters.json',dict(eta=.487,alpha=.733,lambda_housing=point['lambda_housing'],bounds=bounds['lambda_housing'],experimental_not_adopted=True))
        actual=candidate['fp'].actual_parameters(candidate['prepared'],candidate['P'],grid)''')
    (adapter/'runner.py').write_text(src)
    # Base only selects the private indexed package; observer authentication stays frozen.
    src=(base/'run_comparison.py').read_text();src=replace(src,"ROOT=PACKET.parents[4]",f'ROOT=Path({str(ROOT)!r})')
    src=replace(src,"MATCHED=ROOT/'output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed'",f'MATCHED=Path({str(matched)!r})')
    src=replace(src,"str(ROOT/'code/model')",f'str(Path({str(tree / "code/model")!r}))')
    (base/'run_comparison.py').write_text(src)
    # Private engines must be freshly frozen from the experiment clone, never active sources.
    engine=HERE/'refactor_lab'
    if not (engine/'engine/kernels.py').is_file():raise RuntimeError('Parent CES clone not ready')
    private=tree/'code/model/refactor_lab';shutil.copytree(engine,private,dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','*.nbc','*.nbi'))
    small=matched/'source/small_credit_lab';shutil.copytree(private/'engine',small/'engine',dirs_exist_ok=True)
    for p in (private/'engine').glob('*.py'):
        if sha(p)!=sha(small/'engine'/p.name):raise RuntimeError('Private engine pair differs '+p.name)
    fast=(NM/'fast_objective.py').read_text();fast=replace(fast,'NATIVE = HERE.parent / "utility_calibration_round1_v1"','NATIVE = HERE / "utility_adapter"');(runtime/'fast_objective.py').write_text(fast)
    # Native expected-parameter reporting never receives an unreported lambda key.
    oldplan=json.loads((adapter/'plan.json').read_text());lane=oldplan['lanes']['floor_s0']
    saved=json.loads((PSI/'chain_11/results/0026_nm/proposed_parameters.json').read_text());seed=dict(saved['free']);seed.pop('h_P');seed['lambda_housing']=.2
    bounds={k:v for k,v in lane['bounds'].items() if k!='h_P'};bounds['psi_child']=[.01,.5];bounds['lambda_housing']=[0.,1.]
    oldplan['lanes']['floor_s0'].update(seed=seed,bounds=bounds,free_coordinates=list(seed),arm='ces')
    write(adapter/'plan.json',oldplan)
    contract=oldplan['target_contract'];fingerprint=hashlib.sha256(json.dumps(contract,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    write(PACKET/'plan.json',dict(seed=seed,bounds=bounds,seed_source=str(PSI/'chain_11/results/0026_nm/proposed_parameters.json'),seed_source_sha256=sha(PSI/'chain_11/results/0026_nm/proposed_parameters.json'),target_contract=contract,target_weight_fingerprint=fingerprint,price_start=float(json.loads((PSI/'chain_11/results/0026_nm/phase_b_ge/selected_root/closure.json').read_text())['price']),maximum_full_GE=20,maximum_search_proposals=17,budget_seconds=3600,reserve_seconds=900,closure='Price clears birth renewal; population clears absolute physical housing supply',dimensions=[120,9],credit=0.,fixed_eta=.487,fixed_alpha=.733,economic_changes=['Experimental CES replaces Cobb-Douglas','Experimental consumption/housing child scales; lambda_housing estimated','Experimental floors zero; compensation off; h_P removed'],experimental_not_adopted=True))
    write(PACKET/'inherited_source_pins.json',copied)
    pins={str(p.relative_to(PACKET)):sha(p) for p in runtime.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    shutil.copy2(HERE/'run_ces.py',PACKET/'run_ces.py')
    pins['run_ces.py']=sha(PACKET/'run_ces.py')
    pins['plan.json']=sha(PACKET/'plan.json');pins['inherited_source_pins.json']=sha(PACKET/'inherited_source_pins.json')
    write(PACKET/'source_pins.json',pins)
    write(PACKET/'deployment/preparation.json',dict(status='prepared_no_submission_no_solves',source_files=len(pins),lifecycle_solves=0,private_engine_pair_exact=True,authenticated_inputs_read_only=True))
    print(json.dumps(dict(source_files=len(pins),runtime=str(runtime),lifecycle_solves=0)))
if __name__=='__main__':main()
