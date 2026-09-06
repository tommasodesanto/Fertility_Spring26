#!/usr/bin/env python3
"""One fixed-price replay checks that owner reporting repair changes no choice."""
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

if __name__=='__main__':main()
