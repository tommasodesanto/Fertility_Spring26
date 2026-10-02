"""Isolated N0=1 observer adapter; unchanged household solves and renewal root."""
from __future__ import annotations
import copy, difflib, hashlib, inspect, json, sys
from pathlib import Path
from types import ModuleType
import numpy as np
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'utility_floor_nm_v1'
sys.path.insert(0,str(OLD))
import fast_objective as fast

class H0BoundError(RuntimeError): pass

def normalize_observation(context, live, P, ev, cal, price, final):
    """Derive H0 at fixed policies, then update every actual accounting object."""
    demand=np.asarray(ev.demand_by_loc,dtype=float)
    if demand.shape!=(1,) or not np.isfinite(demand).all() or demand[0]<=0:
        raise RuntimeError('Invalid one-market normalized demand')
    if float(P.N_target)!=1. or abs(float(np.asarray(ev.g_current).sum())-1.)>5e-9 or abs(float(live['sol'].total_mass)-1.)>5e-9:
        raise RuntimeError('Actual normalized household mass differs from N0=1')
    factor=(P.user_cost_rate*price/P.r_bar)**P.xi_supply
    if not np.isfinite(factor).all() or (factor<=0).any():raise RuntimeError('Nonpositive normalized supply factor')
    H0=demand/factor
    if final and not .2<=float(H0[0])<=80.:
        raise H0BoundError('Derived normalized H0 outside original [.2,80] bounds')
    P=copy.copy(P);P.H0=H0.copy()
    supply_values=P.H0*factor
    supply=cal.HousingSupplyRule('static-elastic',float(price[0]),float(supply_values[0]),float(P.xi_supply[0]))
    ev.supply_by_loc=supply_values.copy()
    excess=demand-supply_values
    metric=float(np.max(np.abs(excess)/np.maximum(np.abs(supply_values),1e-12)))
    ev.relative_market_residual=metric
    sol=live['sol']
    # Packed fixed-price solution has old supply fields; policies and mass remain exact.
    sol.housing_supply=supply_values.copy()
    sol.aggregate_housing_supply=float(supply_values.sum())
    sol.aggregate_housing_demand=float(demand.sum())
    sol.aggregate_housing_excess=float(excess.sum())
    sol.market_housing_demand=demand.copy()
    sol.market_aggregate_housing_demand=float(demand.sum())
    sol.market_housing_excess=excess.copy()
    sol.market_aggregate_housing_excess=float(excess.sum())
    sol.best_max_abs_rel_excess=sol.best_market_metric=metric
    sol.converged=bool(metric<=getattr(P,'tol_eq',1e-4))
    live['P']=P
    ctx=dict(context);ctx['P']=P
    ctx['expected_parameters']=dict(context['expected_parameters'],H0=float(H0[0]))
    # Retained native selected/repeat stage arrays must contain normalized quantities.
    stage=live.get('stage_dir')
    if stage is not None:
        stage=Path(stage);arrays=stage/'solution_arrays.npz'
        if arrays.exists():
            with np.load(arrays,allow_pickle=False) as saved: values={k:saved[k].copy() for k in saved.files}
            for k,v in vars(sol).items():
                if isinstance(v,np.ndarray) and v.dtype!=object: values[k]=v
            for name in ('g_pre','g_post_fertility','g_current','g_stay_distribution'):
                value=getattr(ev,name,None)
                if isinstance(value,np.ndarray) and value.dtype!=object:
                    values['distribution.'+name]=value.copy()
            values['parameters.H0']=P.H0;values['normalization.population']=np.array([1.])
            np.savez_compressed(arrays,**values)
        context['fp'].write(stage/'normalization.json',dict(population=1.,H0=float(H0[0]),demand=float(demand[0]),supply=float(supply_values[0]),residual=float(excess[0]),derived_H0_bounds=[.2,80.]))
    return ctx,P,ev,supply

def observer_variant(source):
    replace=fast._replace_once
    source=replace(source,'    _require(np.array_equal(P.H0, context["P"].H0) and np.array_equal(P.xi_supply, context["P"].xi_supply),\n             "Physical housing supply curve drift")',
        '    _require(np.array_equal(P.xi_supply, context["P"].xi_supply), "Housing supply elasticity drift")')
    source=replace(source,'    renter_floor = audit_realized_renter_floor(P, grid, policy, ev)',
        '    context, P, ev, supply = _normalize_observation(context, live, P, ev, cal, price, final)\n    renter_floor = audit_realized_renter_floor(P, grid, policy, ev)')
    source=replace(source,'    population = physical_supply / demand','    population = 1.0')
    source=replace(source,'                  occupied_renter_floor=renter_floor)',
        '                  occupied_renter_floor=renter_floor, normalized_population=1.0,\n                  H0_derived=float(P.H0[0]), H0_bounds=[.2,80.],\n                  housing_supply_coefficient_role="derived calibrated coefficient at N0=1")')
    source=replace(source,'                row["status"] = "fixed population scale; not identified by per-household targets"',
        '                row["status"] = "derived calibrated housing supply coefficient at N0=1"\n                row["lower"], row["upper"] = "0.2", "80.0"\n                row["near_bound"] = str(min(float(P.H0[0])-.2,80.-float(P.H0[0])) <= .01*79.8)')
    # Supply and demand already share N0=1 units; preserve standard diagnostic graph set.
    source=replace(source,'        result["standard_plot_supply_units"] = "physical supply divided by endogenous household population"',
        '        result["standard_plot_supply_units"] = "physical housing supply at normalized N0=1"')
    return source

def make_evaluator(out,lane,P,grid,deadline,price_start=None,*,native_runner,exploratory=False):
    out=Path(out);native_runner.verify_sources()
    # Same initialization/import route as the reviewed fast evaluator.
    sys.path.insert(0,str(native_runner.BASE))
    fast._load(native_runner.BASE/'run_comparison.py','normalized_path_initializer')
    sys.path.insert(0,str(fast.NATIVE))
    original=fast._load(fast.NATIVE/'phase_b_pilot.py','normalized_original_ge')
    originals=[inspect.getsource(original.observe_price),inspect.getsource(original.run_phase_b),inspect.getsource(native_runner.native_evaluator)]
    if exploratory:
        sources=list(fast.variants(*originals))
        # The fast adapter expects the old plotting-units text when removing plots.
        sources[0]=observer_variant(sources[0]+ '\n') if 'standard_plot_supply_units' in sources[0] else observer_variant_no_plot(sources[0])
    else: sources=list(originals);sources[0]=observer_variant(sources[0])
    sources[2]=fast._replace_once(sources[2],'    import phase_b_pilot as ge\n','    ge = _get_normalized_ge(out)\n') if not exploratory else fast._replace_once(sources[2],'    ge = _get_fast_ge(out)\n','    ge = _get_normalized_ge(out)\n')
    sources[2]=fast._replace_once(sources[2],'        except TimeoutError as exc:',
        '        except H0BoundError as exc:\n'
        '            return dict(status=\"inadmissible_numerical\",reason=str(exc),lifecycle_solves=budget.used_lifecycle,derived_H0_bound_rejection=True,rejection_kind=\"derived_H0_constraint\")\n'
        '        except TimeoutError as exc:')
    review=out/'normalization_source_review';review.mkdir(parents=True,exist_ok=True)
    receipt=dict(normalized_population=1.,derived_H0_bounds=[.2,80.],exploration_unverified=exploratory,household_solver_unchanged=True,renewal_root_unchanged=True,source_sha256={})
    for name,before,after in zip(('observe_price','run_phase_b','native_evaluator'),originals,sources):
        (review/(name+'.py')).write_text(after)
        (review/(name+'.diff')).write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='original/'+name,tofile='normalized/'+name)))
        receipt['source_sha256'][name]={'original':hashlib.sha256(before.encode()).hexdigest(),'normalized':hashlib.sha256(after.encode()).hexdigest()}
    native_runner.write(review/'receipt.json',receipt)
    ge=ModuleType('normalized_ge');ge.__dict__.update(original.__dict__);ge._normalize_observation=normalize_observation
    exec(compile(sources[0],str(review/'observe_price.py'),'exec'),ge.__dict__)
    exec(compile(sources[1],str(review/'run_phase_b.py'),'exec'),ge.__dict__)
    exec(compile(inspect.getsource(original._observe_with_deadline),'normalized_deadline_wrapper','exec'),ge.__dict__)
    namespace=dict(native_runner.__dict__);namespace['_get_normalized_ge']=lambda _out:ge;namespace['H0BoundError']=H0BoundError
    exec(compile(sources[2],str(review/'native_evaluator.py'),'exec'),namespace)
    evaluate=namespace['native_evaluator'](out,lane,P,grid,deadline,price_start)
    def checked(label,point,end):
        result=evaluate(label,point,end)
        if result['status']=='passed':
            closure=json.loads((Path(result['report'])/'closure.json').read_text())
            if result['population']!=1. or closure['population_scale']!=1.:raise RuntimeError('N0=1 accounting drift')
            result['H0_derived']=closure['H0_derived'];result['normalization']='N0=1'
        return result
    return checked

def observer_variant_no_plot(source):
    # Fast variants remove this exact line; use the full transform then remove it.
    line='        result["standard_plot_supply_units"] = "physical supply divided by endogenous household population"\n'
    marker='        fp.write(out / "closure.json", result)\n'
    pos=source.rfind(marker)
    if pos<0:raise RuntimeError('Missing fast observer final receipt')
    full=source[:pos]+line+source[pos:]
    transformed=observer_variant(full)
    return transformed.replace('        result["standard_plot_supply_units"] = "physical housing supply at normalized N0=1"\n','',1)
