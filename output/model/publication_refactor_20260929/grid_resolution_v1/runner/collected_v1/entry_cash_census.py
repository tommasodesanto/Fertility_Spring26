"""Input-only necessary entrant cash bound; no Bellman/KFE/model solve."""
import json,sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import run_comparison as runner
import phase_a
prep=runner.verify_sources()
context=runner.authored.context_from_bundle(SimpleNamespace(bundle=runner.ROOT/'output/model/publication_refactor_20260929/local_export_v1/inputs',reference_root=runner.ROOT,out=HERE))
results={}
for name in ['control_160x15','proposal_120x9']:
    if name.startswith('proposal'):runner.apply_proposal(context,prep)
    bound=phase_a.entrant_budget_bound(context,price_factors=(1.,))
    for key in ('selected_d_bar','operational_mesh','native_upper_search_buffer','margin_above_requirement','robustness_price_factors'):bound.pop(key,None)
    P=context['P'];C=P.fixed_reference_entry_conditional;grid=context['b_grid']
    from small_credit_lab.engine.shared import income_at_state
    from small_credit_lab.engine import solver
    sd=solver.precompute_shared(P,grid)
    assert float(sd.gb_flat[0,0])==0.,'Childless transfer floor invalidates cash formula'
    cb=float(sd.cb_flat[0,0]);hb=float(sd.hb_flat[0,0]);rent=P.user_cost_rate*context['q_ref']
    infeasible=[];mass=0.
    for bi,zi in zip(*np.nonzero(C>0)):
        income=float(income_at_state(P,0,0,float(P.z_grid[zi])))
        slack=float(P.R_gross*grid[bi]+income+runner.D-cb-rent*hb)
        if slack<=1e-6:
            weight=float(C[bi,zi]*P.z_weights[zi]);mass+=weight
            infeasible.append(dict(b_index=int(bi),income_index=int(zi),wealth=float(grid[bi]),income_multiplier=float(P.z_grid[zi]),income=income,entrant_joint_mass=weight,necessary_renter_surplus_at_credit_floor=slack))
    results[name]=dict(necessary_budget_bound=bound,necessary_renter_surplus_failure_mass=mass,failing_cells=infeasible,wealth_endpoints=grid[[0,-1]].tolist(),childless_transfer_floor=float(sd.gb_flat[0,0]),income_function_includes_property_rebate_and_active_estate_transfer=True,entry_location_weights_sum=float(np.sum(P.entry_by_loc)),lifecycle_solves=0)
receipt=dict(credit_magnitude_selection=False,scope='Necessary current entrant renter-cash condition only; no continuation/Bellman frontier or relocation census exists. Cannot establish which states the failed KFE relocated or how much.',d_fixed=runner.D,reference_price=context['q_ref'],arms=results,lifecycle_solves=0)
(HERE/'entry_cash_census.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:dict(max_necessary_d=v['necessary_budget_bound']['unrounded_requirement'],renter_cash_failure_mass=v['necessary_renter_surplus_failure_mass'],failing_cells=len(v['failing_cells'])) for k,v in results.items()},indent=2))
