"""One isolated Estate-A stationary GE, using current fixed parameters."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
PROJECT = Path(__file__).resolve().parents[4]
OUTPUT = PROJECT/'output/model/experiments/birth_count_choice/estate_a_v1'


def selected_inputs(birth_cap):
    from run import selected_inputs as canonical_inputs
    from model.estate_contract import apply_experiment_flags, experiment_flags
    config,P,grid=canonical_inputs()
    apply_experiment_flags(P,experiment_flags(birth_cap))
    return config,P,grid


def preflight(birth_cap):
    import tempfile,time
    from types import SimpleNamespace
    from run import selected_inputs as canonical_inputs
    from model.calibration import effective_input_fingerprint,_contract
    from model.estate_contract import contract,experiment_flags,apply_experiment_flags,rescore_rows,OBSERVER_CONTRACT
    from model.reporting import build_context,production_model_facade
    config,P,grid=canonical_inputs()
    fields=dict(vars(P)); before=effective_input_fingerprint(P,grid)
    flags=experiment_flags(birth_cap); apply_experiment_flags(P,flags)
    restored=dict(vars(P))
    for key in flags:
        if key in fields: restored[key]=fields[key]
        else: restored.pop(key)
    if effective_input_fingerprint(SimpleNamespace(**restored),grid)!=before: raise RuntimeError('Nonflag input drift')
    base,old_target,old_weight,_=_contract()
    for key,val in [('target_fingerprint',old_target),('weight_fingerprint',old_weight)]:
        if config['provenance'][key]!=val: raise RuntimeError('Canonical parameter provenance drift')
    new,target_pin,weight_pin,bounds=contract()
    with tempfile.TemporaryDirectory(prefix='estate_a_preflight_') as directory:
        context=build_context(P,grid,directory,price_start=config['price_guess'],deadline=time.time()+120,max_lifecycle=32,closure='fixed_h0')
        facade=production_model_facade()
        from model.engine import solver
        if facade.solve_markov_income_at_prices is not solver.solve_markov_income_at_prices: raise RuntimeError('Wrong engine')
        if len(context['manifest']['full_parameter_table'])!=31: raise RuntimeError('Parameter table changed')
        adapters=context.get('birth_count_observer_adapters',[])
    return dict(status='passed_zero_solves',lifecycle_solves=0,birth_cap=birth_cap,experiment_flags=flags,
        canonical_input_fields=len(fields),changed_input_fields=sorted(flags),canonical_fields_preserved=True,
        canonical_input_fingerprint=before,target_rows=len(new),parameter_rows=31,
        target_fingerprint=target_pin,weight_fingerprint=weight_pin,native_original_target_fingerprint=old_target,
        native_original_weight_fingerprint=old_weight,observer_contract=OBSERVER_CONTRACT,
        budget_seconds=900,max_lifecycle=32,closure='fixed_h0',fixed_H0=P.H0.tolist(),
        output_root=str(OUTPUT/('single' if birth_cap==1 else 'multiple')),
        observer_functions_adapted=[r['function'] for r in adapters],search_beta_bounds=bounds['beta_annual'])


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--birth-cap',type=int,choices=(1,3),required=True)
    parser.add_argument('--preflight',action='store_true')
    args=parser.parse_args(argv)
    sys.path.insert(0,str(Path(__file__).resolve().parent))
    if args.preflight:
        print(json.dumps(preflight(args.birth_cap),indent=2)); return
    from model.workflow import run_stationary
    from model.estate_contract import experiment_flags
    config,_,_=selected_inputs(args.birth_cap)
    result,case=run_stationary(config['parameters'],config['external_inputs'],config['native_overrides'],
        price_guess=config['price_guess'],budget_seconds=900,max_lifecycle=32,closure='fixed_h0',
        output_root=OUTPUT/('single' if args.birth_cap==1 else 'multiple'),experiment_flags=experiment_flags(args.birth_cap),
        parameter_file_metadata={key:config[key] for key in ('config_source','config_sha256','config_text','provenance')})
    receipt=json.loads((case/'estate_a_rescore_receipt.json').read_text())
    print(f'Converged Estate-A cap={args.birth_cap}: {case}; price={result.price:.12g}; new-contract loss={receipt["loss"]:.12g}')
    return case

if __name__=='__main__': main()
