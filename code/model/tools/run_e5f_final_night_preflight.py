"""Pinned saved-policy demographic checks; no new equilibrium is claimed."""
from pathlib import Path
from dataclasses import asdict, replace
import argparse, gzip, hashlib, json, pickle, sys, time
import numpy as np

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a=ap.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    def save(name,data):
        p=a.output/name;q=p.with_suffix('.tmp')
        q.write_text(json.dumps(data,indent=2,default=lambda v:v.tolist() if hasattr(v,'tolist') else str(v))+'\n');q.replace(p)
    p=json.loads(a.plan.read_text());root=Path(p['source_root'])
    sys.path.extend([str(root/'code/model/tools'),str(root/'code/model')])
    try:
        import run_e5f_candidate_terminal as td
        c=p['terminal_template'];td.validate_contract(c);td.verify_sources(c)
        td.verify(c['initial_checkpoint']['path'],c['initial_checkpoint']['sha256'])
        with gzip.open(c['initial_checkpoint']['path'],'rb') as f:packet=pickle.load(f)
        import e5f_balanced_history as balanced
        from e5f_overnight_demography import dependent_units_by_location, _dependent_deaths_by_location
        _,primitive,_,_=balanced._runtime()
        import run_e5f_open_population_transition as transition
        P=packet['parameters'];e=packet['evaluation'];grid=packet['b_grid'];shared=packet['shared']
        save('heartbeat.json',dict(phase='actual_saved_policy_operator',elapsed=time.monotonic()-started))
        advanced,mature,exits,mass_error=transition.advance_sequential_calendar_distribution(e,np.zeros(P.I),P,grid,shared)
        children=float(dependent_units_by_location(e.g_post_fertility).sum())
        deaths=float(_dependent_deaths_by_location(e.g_post_fertility,P).sum())
        next_children=float(dependent_units_by_location(advanced).sum())
        child_error=next_children-(children-deaths-float(mature.sum()))
        if max(abs(mass_error),abs(child_error))>2e-9:
            raise RuntimeError('Actual saved-policy demographic identity failed')
        demographics=packet['demographic_seed']
        zero=replace(demographics,net_migration={year:np.zeros_like(value) for year,value in demographics.net_migration.items()})
        if any(np.any(x!=0) for x in zero.net_migration.values()):raise RuntimeError('Nonzero closed migration')
        # No arbitrary formation conversion is used. Report the restriction
        # needed for this particular state to replace household exits.
        result=dict(status='passed_saved_policy_demographic_preflight',
            post_birth_dependents=children,dependent_deaths=deaths,surviving_maturation=float(mature.sum()),
            household_exits=float(exits),household_identity_error=float(mass_error),child_identity_error=child_error,
            implied_entries_per_model_child_for_constant_heads=float(exits/mature.sum()),
            implied_conversion_is_diagnostic_only=True,zero_migration_cells_verified=True,
            original_fiscal_rule='unrebated seed, used only for demographic operator check',
            B_production_formation_contract='outstanding',B_plus_migrant_allocation='outstanding',
            equilibrium_solved=False,production_eligible=False,elapsed=time.monotonic()-started)
        for name in ('latest_completed.json','best_so_far.json','summary.json'):save(name,result)
    except Exception as exc:
        save('failure.json',dict(status='failed',error=str(exc),type=type(exc).__name__,elapsed=time.monotonic()-started));raise

if __name__=='__main__':main()
