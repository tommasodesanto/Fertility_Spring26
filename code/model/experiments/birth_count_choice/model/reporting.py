"""Authenticated frozen observers, bound to the canonical production engine.

The empirical observers and certificates remain a read-only compatibility
stack. Household and distribution callbacks are explicitly production callbacks.
"""
from __future__ import annotations
import copy, hashlib, importlib.util, inspect, json, sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from .inputs import DEFAULT_PARAMETERS
ROOT = Path(__file__).resolve().parents[5]
PACKETS = ROOT/'output/model/fixed_reference_economics_20260928'
AUTHORED = ROOT/'output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed'
_AUTHENTICATED_TEMPLATE = None

def _load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return module

def _install_read_only_overlay():
    # Same two hash-checked frozen source redirects used by the authenticated
    # local timing oracle. No writes or source mutation is performed.
    overlay=PACKETS/'purchase_rules_overnight_v1/local_runtime/bootstrap.py'
    marker='_fertility_model_playground_overlay'
    if getattr(sys,marker,None):
        if getattr(sys,marker)!=str(overlay): raise RuntimeError('conflicting frozen overlay')
        return
    source=overlay.read_text().split("if '--preflight-context' in sys.argv:",1)[0]
    if 'DIGESTS=' not in source or 'Frozen overlay write forbidden' not in source: raise RuntimeError('frozen overlay source shape drift')
    exec(compile(source,str(overlay),'exec'),{'__file__':str(overlay),'__name__':'production_read_only_overlay'})
    setattr(sys,marker,str(overlay))

def production_model_facade():
    """Flatten the same stage definitions expected by native calendar observers."""
    from .engine import solver, parameters, utils, adult_entry
    values={}
    for module in (parameters,utils,adult_entry,solver):
        values.update({k:v for k,v in vars(module).items() if not k.startswith('__')})
    values['__file__']=solver.__file__
    return SimpleNamespace(**values)

def _install_timing_audit(context,out):
    rt=context['prepared'].rt
    if getattr(rt['accounting'],'_transaction_timing_installed',False): return
    old=rt['accounting']._inherited_audit();before=inspect.getsource(old)
    replacements={
        'x = grid if old == new else grid + sale[old] - costs[new]':'x = grid if old == new else grid + (sale[old] - costs[new]) / float(P.R_gross)',
        'invalid = x + y / float(P.R_gross) < floor - 1e-10':'invalid = float(P.R_gross) * x + y < floor - 1e-10',
    }
    after=before
    for a,b in replacements.items():
        if after.count(a)!=1: raise RuntimeError('Purchase audit source drift: '+a)
        after=after.replace(a,b)
    namespace=dict(old.__globals__);path=Path(out)/'timing_purchase_audit.py';path.write_text(after)
    exec(compile(after,str(path),'exec'),namespace)
    rt['accounting']._INHERITED=namespace['audit_purchase_accounting']
    rt['accounting']._transaction_timing_installed=True
    context['fp'].write(Path(out)/'timing_observer_receipt.json',{
        'original_audit_sha256':hashlib.sha256(before.encode()).hexdigest(),
        'timing_audit_sha256':hashlib.sha256(after.encode()).hexdigest(),
        'forward_map_file':inspect.getsourcefile(rt['model'].build_forward_tenure_transition_maps),
        'acceptance_tolerances_unchanged':True,'audit_budget':'R*b + net_sale - purchase + income - consumption - costs'})

def _install_recent_observer(context,facade,out):
    # The authenticated wrapper retains its reviewed dead-tail options. Its
    # enclosed recent observer has a local import used in an identity check.
    # Replace only that import with the same canonical facade as calendar.model.
    wrapped=context['prepared'].rt['observe_recent_parent_flow']
    recent=next((cell.cell_contents for cell in (wrapped.__closure__ or ())
                 if hasattr(cell.cell_contents,'observe_recent_parent_flow')),None)
    if recent is None: raise RuntimeError('Authenticated recent-observer closure drift')
    original=getattr(recent,'_production_original_observe',recent.observe_recent_parent_flow)
    recent._production_original_observe=original
    before=inspect.getsource(original)
    old='    from intergen_eqscale_seq_optimized import solver as model'
    if before.count(old)!=1: raise RuntimeError('Recent-observer model import drift')
    after=before.replace(old,'    model = _production_model_facade')
    namespace=dict(original.__globals__);namespace['_production_model_facade']=facade
    path=Path(out)/'production_recent_parent_observer.py';path.write_text(after)
    exec(compile(after,str(path),'exec'),namespace)
    recent.observe_recent_parent_flow=namespace['observe_recent_parent_flow']
    context['fp'].write(Path(out)/'recent_observer_adapter_receipt.json',dict(
        original_source_path=inspect.getsourcefile(original),clone_source_path=str(path),
        original_sha256=hashlib.sha256(before.encode()).hexdigest(),
        clone_sha256=hashlib.sha256(after.encode()).hexdigest(),
        edit='Only local solver import becomes the bound production facade; identity guard and all numerical statements unchanged',
        acceptance_tolerances_unchanged=True))

def build_context(P,grid,out,*,price_start,deadline,max_lifecycle,closure):
    """Authenticate reporting once; never build or replace production primitives."""
    global _AUTHENTICATED_TEMPLATE
    _install_read_only_overlay()
    from .frozen_sources import install_recovered_sources
    install_recovered_sources()
    if _AUTHENTICATED_TEMPLATE is None:
        sys.path[:0]=[str(AUTHORED),str(AUTHORED/'source'),str(ROOT/'code/model')]
        authored=_load(AUTHORED/'driver.py','production_authored_reporting')
        initial=authored.context_from_bundle(SimpleNamespace(bundle=ROOT/'output/model/publication_refactor_20260929/local_export_v1/inputs',reference_root=ROOT,out=Path(out)))
        authored.authenticate_frozen(initial)
        initial['_native_actual_parameters']=initial['fp'].actual_parameters
        _AUTHENTICATED_TEMPLATE=initial
    context=dict(_AUTHENTICATED_TEMPLATE)
    context['out']=Path(out)
    # The utility adapter changes only the recorded physical-room-floor field.
    fp=context['fp'];native_actual=context['_native_actual_parameters']
    def actual(prepared,Q,b):
        values=native_actual(prepared,Q,b)
        values['h_P']=float(Q.hbar_first_child_jump) if Q.child_room_floor else 0.
        if float(Q.delta)!=float(PRECISE_BASE_DELTA): values['annual_depreciation']=1-(1-float(Q.delta))**(1/float(Q.period_years))
        if float(Q.tau_H)!=float(PRECISE_BASE_TAX): values['annual_property_tax']=float(Q.tau_H)/float(Q.period_years)
        return values
    fp.actual_parameters=actual
    rt=context['prepared'].rt;facade=production_model_facade()
    rt['model']=facade
    cal=rt['primitive'].pf.calendar;cal.model=facade
    _install_timing_audit(context,out)
    _install_recent_observer(context,facade,out)
    from .observer_adapters import install_birth_count_observers
    install_birth_count_observers(context, facade, out)
    from .estate_audit_adapter import install_estate_audit
    install_estate_audit(context, P, out)
    from .estate_observer_adapter import install_estate_observer
    install_estate_observer(context, P, facade, out)
    context.update(P=copy.deepcopy(P),b_grid=np.asarray(grid).copy(),selected_d_bar=float(P.unsecured_credit_limit),
        reference_psi=float(P.psi_child),expected_dimensions={'wealth_grid_nodes':int(P.Nb),'income_states':int(P.Nz)},
        free_coordinates=list(DEFAULT_PARAMETERS),price_start=float(price_start),phase_b_max_new_lifecycle=int(max_lifecycle),
        deadline_epoch=float(deadline),closure_mode=closure)
    parameter_rows=json.loads((ROOT/'code/model/production/reference_inputs/parameter_table.json').read_text())
    for row in parameter_rows:
        row['estimate']=row['reference_estimate']
    context['manifest']=copy.deepcopy(context['manifest'])
    context['manifest']['full_parameter_table']=parameter_rows
    context['expected_parameters']=fp.actual_parameters(context['prepared'],P,grid)
    # Preserve historical restrictions as provenance, without applying them to
    # ordinary solves. Current effective values are checked, never substituted.
    rows=copy.deepcopy(json.loads((ROOT/'code/model/production/reference_inputs/bundle.json').read_text()))
    context['reference_metadata']={k:rows[k] for k in ('source','source_sha256','target_fingerprint','weight_fingerprint')}
    return context

PRECISE_BASE_DELTA=0.05545379079326218
PRECISE_BASE_TAX=0.042393443095490375

def historical_reporting_dependency():
    return {'observer':str(PACKETS/'sources/fixed_price_v1/run_fixed_price.py'),
            'household_engine':'code/model/experiments/birth_count_choice/model/engine',
            'status':'authenticated read-only empirical/accounting compatibility'}
