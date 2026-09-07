#!/usr/bin/env python3
"""Bounded saved-state checks of reporting or saving, with explicit source pins."""
import argparse, copy, hashlib, json, sys, time
from pathlib import Path
import numpy as np

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model-root',type=Path,required=True)
    ap.add_argument('--bundle-sha256',required=True)
    ap.add_argument('--checkpoint',type=Path,required=True)
    ap.add_argument('--checkpoint-sha256',required=True)
    ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--reference',type=Path)
    ap.add_argument('--saving-diagnosis',action='store_true',help='Compare local and exhaustive saving at the checkpoint price, without another policy change.')
    a=ap.parse_args();root=a.model_root.resolve();out=a.outdir.resolve()
    if out.exists():raise RuntimeError('Refusing an existing diagnostic directory')
    out.mkdir(parents=True)
    sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
    import run_e5f_independent_numerical_audit as audit
    import run_e5f_transition_calibration as calibration
    from intergen_eqscale_seq_optimized import solver
    if calibration.code_fingerprint_contract(solver)['bundle_sha256']!=a.bundle_sha256:
        raise RuntimeError('Wrong scientific source')
    if audit.digest(a.checkpoint)!=a.checkpoint_sha256:raise RuntimeError('Changed inherited checkpoint')
    packet=audit.load_checkpoint(a.checkpoint);P=copy.deepcopy(packet['parameters']);bg=packet['b_grid']
    audit.transition.configure_sequential_model()
    audit.calendar.apply_fertility=audit.transition.apply_sequential_fertility
    audit.calendar.advance_calendar_distribution=audit.transition.advance_sequential_calendar_distribution
    if a.saving_diagnosis:
        return saving_diagnosis(packet,audit,out,a)
    audit.policy.apply_policy(P,audit.policy.POLICIES['supply-plus-20'])
    rule=audit.policy.policy_supply_rule(packet['supply_rule'],audit.policy.POLICIES['supply-plus-20'])
    shared=solver.precompute_shared(P,bg)
    price=np.array([.52918721636948]) # exact failed-policy price, already cleared
    start=time.time()
    evaluation=audit.calendar.evaluate_period(price,packet['evaluation'].inherited_g_pre,P,bg,shared,audit.calendar.SolveCounter(),rule)
    current=dict(parameters=P,b_grid=bg,evaluation=evaluation,shared=shared,supply_rule=rule)
    budget=audit.budget_audit(current,out);arrays=audit.policy_array_audit(current,out)
    audit.standard_diagnostics(current,out,validate_production_young=False)
    p=evaluation.policy
    fields={name:getattr(p,name) for name in ('V','hR_pol','bp_pol','tenure_choice','tenure_probs','loc_probs','fert_probs','fert_value','fert2_probs')}
    fields.update({name:getattr(evaluation,name) for name in ('g_pre','g_post_fertility','g_current')})
    fields.update(joint_probabilities=p.joint_choice.probabilities,joint_products=p.joint_choice.products,joint_wait=p.joint_choice.wait_probabilities)
    hashes={name:hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest() for name,value in fields.items()}
    quantities=dict(births=evaluation.births,demand=evaluation.demand_by_loc.tolist(),supply=evaluation.supply_by_loc.tolist(),price=price.tolist())
    result=dict(status='reference_reporting_defect_reproduced',bundle=a.bundle_sha256,
        elapsed_seconds=time.time()-start,array_hashes=hashes,quantities=quantities,budget=budget,policy_arrays=arrays,
        consumption_sha256=hashlib.sha256(np.ascontiguousarray(p.c_pol).tobytes()).hexdigest(),
        graph_hashes={f.name:audit.digest(f) for f in sorted((out/'standard_diagnostics').glob('*.png'))})
    if a.reference:
        ref=json.loads(a.reference.read_text())
        if hashes!=ref['array_hashes'] or quantities!=ref['quantities']:
            raise RuntimeError('Reporting repair changed a value, saving, housing, choice, distribution or market quantity')
        if result['graph_hashes']!=ref['graph_hashes']:raise RuntimeError('Standard graph packet unexpectedly changed')
        if budget['budget_excess_mass']>2e-10 or arrays['occupied_negative_steps']!=0:
            raise RuntimeError('Repaired replay still fails unchanged numerical gates')
        result.update(status='pass',exact_unchanged_arrays=len(hashes),exact_standard_graphs=len(result['graph_hashes']),reference_sha256=audit.digest(a.reference))
    elif abs(budget['budget_excess_mass']-1.9390478108254123e-7)>1e-18:
        raise RuntimeError('The original failed policy budget was not reproduced')
    audit.save_json(out/'reporting_check.json',result)
    print(json.dumps(dict(status=result['status'],elapsed_seconds=result['elapsed_seconds'],budget_excess_mass=budget['budget_excess_mass'])),flush=True)


def saving_diagnosis(packet,audit,out,args):
    """Two full Bellmans; the second substitutes previously audited saving kernels."""
    import gzip,pickle
    import run_e5f_global_saving_quantification as saving
    P=packet['parameters'];bg=packet['b_grid'];old=packet['evaluation']
    inherited=old.inherited_g_pre.copy();price=old.policy.price.copy()
    rows=[];original_v=None;start=time.time()
    for method in ('local','global'):
        target=out/method;target.mkdir()
        audit.save_json(out/'heartbeat.json',dict(phase=method,epoch=time.time(),elapsed_seconds=time.time()-start))
        saving.set_method(method)
        shared=audit.model.precompute_shared(P,bg)
        result=audit.calendar.evaluate_period(price,inherited.copy(),P,bg,shared,audit.calendar.SolveCounter(),packet['supply_rule'])
        current=dict(packet,evaluation=result,shared=shared)
        arrays=audit.policy_array_audit(current,target);budget=audit.budget_audit(current,target)
        audit.standard_diagnostics(current,target,validate_production_young=False)
        with gzip.open(target/'dated_state.pkl.gz','wb',compresslevel=1) as stream:pickle.dump(current,stream,protocol=5)
        if method=='local':
            for name in ('V','c_pol','hR_pol','bp_pol','tenure_choice','tenure_probs','loc_probs','fert_probs','fert_value','fert2_probs'):
                if not np.array_equal(getattr(result.policy,name),getattr(old.policy,name)):raise RuntimeError('Local policy replay changed '+name)
            for name in ('g_pre','g_post_fertility','g_current'):
                if not np.array_equal(getattr(result,name),getattr(old,name)):raise RuntimeError('Local population replay changed '+name)
            original_v=result.policy.V.copy()
            dominance=None
        else:
            dominance=float((result.policy.V-original_v)[old.g_pre>1e-12].min())
            if dominance < -1e-7:raise RuntimeError('Exhaustive saving lowers an occupied value')
        row=dict(method=method,quantities=saving.quantities(result,P),policy_arrays=arrays,budget=budget,
                 min_occupied_value_gain=dominance,elapsed_seconds=time.time()-start)
        rows.append(row);audit.save_json(out/'latest_completed_case.json',row)
        audit.save_json(out/'best_so_far.json',dict(scope='diagnostic, not calibration selection',latest_method=method))
        print(json.dumps(dict(method=method,occupied_drops=arrays['occupied_negative_steps'],elapsed_seconds=time.time()-start)),flush=True)
    saving.set_method('local')
    left,right=[row['quantities'] for row in rows]
    summary=dict(status='comparison_complete',checkpoint_sha256=args.checkpoint_sha256,scientific_bundle=args.bundle_sha256,
        driver_sha256=audit.digest(__file__),global_helper_sha256=audit.digest(saving.__file__),oracle_helper_sha256=audit.digest(audit.__file__),
        original_population_preserved=True,local_thirteen_arrays_exact=True,rows=rows,
        births_global_minus_local_percent=100*(right['adjusted_births']/left['adjusted_births']-1),
        ownership_global_minus_local_pp=100*(right['ownership']-left['ownership']),
        rooms_global_minus_local_percent=100*(right['rooms_per_household']/left['rooms_per_household']-1),
        scope='Same inherited population and prices; full lifecycle continuation is solved under each saving method. Global prices are not recleared. No repaired history or policy path is certified.',production_changed=False)
    audit.save_json(out/'saving_diagnosis.json',summary)
    print(json.dumps(dict(status=summary['status'],births_difference_percent=summary['births_global_minus_local_percent'])),flush=True)

if __name__=='__main__':main()
