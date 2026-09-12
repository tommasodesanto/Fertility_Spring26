"""Experimental successive permanent surprises, with a full PF solve per vintage.

No target loading, calibration, launch, or default shock sequence. This module
runs against the approved utility/PAYGO runtime in tmp/e5f_matched_pf. Its caller
owns pinned inputs, endpoint solves, process watchdogs and checkpoint storage.
Historical demographic conditioning is preserved, not re-estimated.
"""
from __future__ import annotations
import copy
from dataclasses import dataclass, replace
from types import SimpleNamespace
import time
import numpy as np


@dataclass
class InheritedState:
    year: int
    households: object


@dataclass
class SurpriseResult:
    path: object
    root_receipt: dict
    next_state: InheritedState | None
    realized_row: dict | None


def local_conditioning(original, year):
    """Reindex observed bridges, preserving original 2007 mass and closure."""
    if year not in (2007, 2011, 2015, 2019, 2023):
        raise ValueError('Surprises currently supported on 2007--2023 four-year dates')
    return replace(original, start_year=year, observer=None,
        next_age_targets={(y-year)//4:y for y in range(year+4,2024,4)})


def _runtime():
    import e5f_balanced_history as balanced
    joined, primitive, checks, rent_domain = balanced._runtime()
    return balanced, joined, primitive, checks, rent_domain


def evaluate_forecast(*, inherited, old_state, demographics, prices, pensions,
                      psi, terminal, observer=None):
    """A constant-preference forecast, with inherited states and the 2023 bridge.

    Future *realized* surprise levels are deliberately absent from this API.
    """
    _, joined, primitive, _, _ = _runtime()
    pf, person = joined.pf, joined.person_pf
    P, grid = old_state.parameters, old_state.b_grid
    start=inherited.year
    if start not in (2007,2011,2015,2019,2023):
        raise ValueError('No restart beyond the fixed historical observation window')
    p=np.asarray(prices,dtype=float);benefits=np.asarray(pensions,dtype=float)
    h=(2023-start)//4
    if p.ndim!=1 or len(p)<h+2 or benefits.shape!=p.shape:
        raise ValueError('Forecast must include history and at least two tail dates')
    if not np.isfinite(psi) or not np.isfinite(p).all() or np.any(p<=0):
        raise ValueError('Finite preference and positive finite prices required')
    psi_path=np.full(len(p),float(psi));transfers=np.zeros(len(p));taxes=np.full(len(p),.179)
    years=start+4*np.arange(len(p))
    tail_values,backward=pf.backward_value_path(prices=p[h:],
        rents=pf.rents_from_asset_prices(p[h:],float(terminal.policy.price[0]),P),
        psi_path=psi_path[h:],terminal_V=terminal.policy.V,base_parameters=P,b_grid=grid,
        transfer_path=transfers[h:],pension_path=benefits[h:],payroll_tax_path=taxes[h:])
    if h:
        conditioning=replace(local_conditioning(old_state.historical_conditioning,start),observer=observer)
        history=pf.evaluate_path_at_prices(prices=p[:h],psi_path=psi_path[:h],
            transfer_path=transfers[:h],terminal_price=float(p[h]),terminal_V=tail_values[0],
            base_parameters=P,b_grid=grid,initial_state=inherited.households,
            supply_rule=old_state.supply_rule,birth_to_entry_conversion=1/2.1,
            historical_conditioning=conditioning,pension_path=benefits[:h],payroll_tax_path=taxes[:h])
        g=history.terminal_state.g_pre
    else:
        g=inherited.households.g_pre
        history=SimpleNamespace(rows=[],values=[],bellman_solves=0,
            maximum_mass_accounting_error=0.,maximum_policy_reproduction_error=0.,
            maximum_feasibility_projection_mass=0.)
    people=demographics.initial_person_state.validated()
    if people.year!=2023:
        raise ValueError('Retained demographic anchor must remain 2023')
    if not h and hasattr(inherited.households,'persons'):
        if inherited.households.persons.year!=2023 or not np.array_equal(inherited.households.persons.persons,people.persons) or not np.array_equal(inherited.households.persons.heads,people.heads):
            raise ValueError('Inherited 2023 person anchor differs; no silent reset allowed')
        people=inherited.households.persons
    heads=person.aggregate_heads_to_model_age_cells(people,age_start=int(P.age_start),
        cell_width=int(P.da),number_of_cells=int(P.J))
    gap=float(np.max(np.abs(g.sum(axis=(0,1,2,4,5,6))-heads)))
    if not np.isfinite(gap) or gap>2e-9:
        raise RuntimeError('Inherited 2023 head-age bridge failed')
    def tail_observer(i,*args):
        if observer is not None:observer(h+i,*args)
    tail=person.evaluate_path_at_prices_person_demography(prices=p[h:],psi_path=psi_path[h:],
        transfer_path=transfers[h:],terminal_price=float(terminal.policy.price[0]),
        terminal_V=terminal.policy.V,base_parameters=P,b_grid=grid,
        initial_state=person.PersonPFState(g_pre=g.copy(),persons=people),
        demographic_primitives=demographics,supply_rule=old_state.supply_rule,
        precomputed_value_path=tail_values,observer=tail_observer if observer else None,
        pension_path=benefits[h:],payroll_tax_path=taxes[h:])
    result=joined.ConditionalHistoryEvaluation(history=history,person_tail=tail,
        rows=[dict(r) for r in history.rows]+[dict(r,period=int(r['period'])+h) for r in tail.rows],
        values=(history.values[:-1] if h else [])+tail.values,
        bellman_solves=history.bellman_solves+backward+tail.bellman_solves,
        initial_2023_age_head_gap=gap)
    joined.check_smoke_gates(result,expected_years=years.tolist())
    return result


def first_period_state(*, inherited, old_state, demographics, path, prices, pensions, psi):
    """Replay only the implemented period under its ORIGINAL expected prices/V.

    Never use the final stationary boundary for this one-period replay; never
    use next vintage's revised house price to reconstruct the preceding rent.
    """
    _,joined,_,_,_=_runtime();pf,person=joined.pf,joined.person_pf
    P=old_state.parameters
    common=dict(prices=[float(prices[0])],psi_path=[float(psi)],transfer_path=[0.],
        terminal_price=float(prices[1]),terminal_V=path.values[1],base_parameters=P,
        b_grid=old_state.b_grid,initial_state=inherited.households,supply_rule=old_state.supply_rule,
        pension_path=[float(pensions[0])],payroll_tax_path=[.179])
    if inherited.year<2023:
        replay=pf.evaluate_path_at_prices(**common,birth_to_entry_conversion=1/2.1,
            historical_conditioning=local_conditioning(old_state.historical_conditioning,inherited.year))
        if replay.maximum_mass_accounting_error>2e-8 or replay.maximum_policy_reproduction_error>2e-10 or replay.maximum_feasibility_projection_mass>1e-6:
            raise RuntimeError('First-period historical replay gate failed')
        state=replay.terminal_state
        if inherited.year==2019:
            state=person.PersonPFState(g_pre=state.g_pre,persons=copy.deepcopy(demographics.initial_person_state))
    else:
        if not hasattr(inherited.households,'persons'):
            raise ValueError('2023 surprise requires the inherited joint household/person state')
        replay=person.evaluate_path_at_prices_person_demography(**common,
            demographic_primitives=demographics,precomputed_value_path=path.values[:2])
        state=replay.terminal_state
        for name,limit in [('maximum_person_identity_error',2e-9),('maximum_head_identity_error',2e-9),
                           ('maximum_household_person_head_gap',2e-9),('maximum_age_head_gap',2e-9),
                           ('maximum_policy_reproduction_error',2e-10),('maximum_feasibility_projection_mass',1e-6)]:
            if not np.isfinite(getattr(replay,name)) or getattr(replay,name)>limit:
                raise RuntimeError('First-period person replay gate failed: '+name)
    if not np.allclose(replay.values[0],path.values[0],rtol=0,atol=2e-10):
        raise RuntimeError('First-period continuation value differs from accepted forecast')
    for key in ('asset_price','renter_price','housing_demand','owner_rate','birth_children_topcode_adjusted',
                'pension_period_units','payroll_tax_revenue','pension_outlays'):
        if not np.isclose(float(replay.rows[0][key]),float(path.rows[0][key]),rtol=0,atol=2e-10):
            raise RuntimeError('First-period replay differs: '+key)
    return InheritedState(inherited.year+4,state)


def solve_surprise(*, inherited, psi, old_state, terminal, terminal_root_receipt,
                   demographic_primitives, terminal_demographic_primitives,
                   count, initial_prices, initial_pensions, audit_controls,
                   root_controls, deadline_monotonic, pension_tail_tolerance,
                   callback=None, observer=None):
    """Solve one full forecast and return its first-period state only if accepted.

    Terminal must be solved at THIS psi, not the ultimate realized historical
    psi. Source/hash/target contracts remain caller-owned and unchanged.
    """
    balanced,joined,primitive,checks,rent_domain=_runtime()
    from e5f_social_security import fiscal_accounts
    from e5f_social_security_root import solve_social_security_path
    from e5f_balanced_terminal import _household_checks,TerminalAuditControls
    if not np.isfinite(psi) or not np.isclose(float(terminal.parameters.psi_child),psi,rtol=0,atol=1e-14):
        raise ValueError('Each surprise requires its own constant-preference terminal')
    # Reuse all existing reference primitive, endpoint, grid, supply and initial
    # normalization checks. This synthetic line is validation-only, never solved.
    reference=copy.copy(old_state)
    reference.psi_path=np.linspace(float(old_state.parameters.psi_child),psi,5)
    reference.diagnostics=dict(old_state.diagnostics,preference_change_2023=psi-float(old_state.parameters.psi_child))
    balanced._validate(reference,terminal,terminal_root_receipt,
        demographic_primitives,terminal_demographic_primitives,count)
    if not isinstance(audit_controls,TerminalAuditControls):raise ValueError('Explicit household audits required')
    for name,limit in dict(reconstruction_tolerance=5e-9,feasibility_projection_tolerance=1e-6,
        probability_tolerance=1e-12,occupied_mass_tolerance=1e-12,value_drop_tolerance=1e-7).items():
        value=getattr(audit_controls,name)
        if not np.isfinite(value) or not 0<=value<=limit:raise ValueError('Invalid household audit '+name)
    rc=dict(root_controls)
    if type(rc['max_evaluations']) is not int or not 2<=rc['max_evaluations']<=8:
        raise ValueError('Retain the bounded2--8 mapping budget; use explicit continuation for more')
    if not 0<rc['market_tolerance']<=2e-4 or not 0<rc['fiscal_tolerance']<=1e-6 or not 0<=rc['final_reproduction_tolerance']<=2e-10:
        raise ValueError('Do not relax market/fiscal/replay tolerances')
    if not np.isfinite(pension_tail_tolerance) or not 0<pension_tail_tolerance<=.01:
        raise ValueError('Explicit terminal pension tolerance no larger than1% required')
    if not np.isfinite(deadline_monotonic) or deadline_monotonic<=time.monotonic():raise ValueError('Expired deadline')
    p0=np.asarray(initial_prices,dtype=float);b0=np.asarray(initial_pensions,dtype=float)
    if p0.shape!=(count,) or b0.shape!=(count,):raise ValueError('Complete dated price/pension guesses required')
    lo,hi=rc.pop('price_bounds');pb=rc.pop('pension_bounds')
    if not np.isfinite([lo,hi]).all() or not 0<lo<hi or not np.isfinite(p0).all() or np.any(p0<lo) or np.any(p0>hi):
        raise ValueError('Finite dated guesses inside explicit positive price bounds required')
    endpoint=SimpleNamespace(parameters=terminal.parameters,asset_price=float(terminal.policy.price[0]),
        renter_price=float(terminal.parameters.user_cost_rate)*float(terminal.policy.price[0]),
        equal_transfer=0.,psi_child=psi,state=joined.person_pf.PersonPFState(terminal.fixed_point.g_pre,terminal.fixed_point.persons))
    final_path=None;coordinates=None;trial=0
    def project(p):
        result=rent_domain.project_price_path_to_positive_rents(np.clip(p,lo,hi),terminal=endpoint,minimum_rent_share=1e-6)[0]
        if np.any(result>hi):raise ValueError('Positive-rent projection exceeds price bounds')
        return result
    def evaluate(p,b):
        nonlocal final_path,coordinates,trial
        final_path=None;coordinates=None;trial+=1;accounts=[];audits=[]
        rents=joined.pf.rents_from_asset_prices(p,endpoint.asset_price,old_state.parameters)
        def observe(i,e,P,grid,shared):
            if i!=len(accounts) or P.pension!=float(b[i]) or P.tau_pay!=.179:raise RuntimeError('Wrong dated fiscal path')
            d,gates=_household_checks(e,P,shared,grid,float(rents[i]),primitive,audit_controls)
            if not all(gates.values()):raise RuntimeError('Dated household audit failed')
            audits.append(dict(year=inherited.year+4*i,diagnostics=d,gates=gates));accounts.append(fiscal_accounts(e.g_current,P))
            if observer:observer(i,e,P,grid,shared)
        result=evaluate_forecast(inherited=inherited,old_state=old_state,demographics=demographic_primitives,
            prices=p,pensions=b,psi=psi,terminal=terminal,observer=observe)
        if len(accounts)!=count:raise RuntimeError('Missing dated fiscal ledger')
        market=[];fiscal=[]
        for row,a in zip(result.rows,accounts):
            if not np.isfinite(row['housing_supply']) or row['housing_supply']<=0:raise RuntimeError('Invalid supply')
            market.append((row['housing_demand']-row['housing_supply'])/row['housing_supply'])
            rev,out=a['payroll_tax_revenue'],a['pension_outlays']
            fiscal.append((rev-out)/max(abs(rev),abs(out),1e-12))
        distance=checks.terminal_convergence_diagnostics(result.person_tail,terminal=endpoint,psi_path=np.full(len(result.person_tail.rows),psi))
        gap=abs(float(b[-1])-float(terminal.parameters.pension))/float(terminal.parameters.pension)
        distance['last_pension_relative_gap']=gap
        distance.setdefault('metrics',{})['pension_relative_gap']=gap
        distance.setdefault('tolerances',{})['pension_relative_gap']=pension_tail_tolerance
        distance['checks']['pension_relative_gap']=gap<=pension_tail_tolerance
        distance['all_checks_pass']=all(distance['checks'].values())
        distance['status']='passed' if distance['all_checks_pass'] else 'not_converged'
        final_path=result;coordinates=(p.copy(),b.copy(),trial)
        return dict(market_residual=market,fiscal_residual=fiscal,mapping_valid=True,
            payload=dict(trial=trial,terminal_distance=distance,dated_household_audits=audits,fiscal_accounts=accounts))
    receipt=solve_social_security_path(closure='fixed_tax',initial_prices=p0,initial_fiscal_values=b0,
        evaluate=evaluate,project_prices=project,fiscal_bounds=pb,deadline_monotonic=deadline_monotonic,
        callback=callback,**rc)
    final=receipt.get('final')
    matched=bool(final is not None and coordinates is not None and final_path is not None
        and final['payload']['trial']==coordinates[2] and np.array_equal(final['prices'],coordinates[0])
        and np.array_equal(final['fiscal_values'],coordinates[1]))
    finite=bool(receipt.get('converged') and matched)
    horizon=bool(finite and final['payload']['terminal_distance']['all_checks_pass'])
    receipt.update(schema='e5f_permanent_surprise_v1',start_year=inherited.year,psi=psi,
        expected_psi_path=[float(psi)]*count,finite_horizon_market_fiscal_converged=finite,
        terminal_distance_passed=horizon,horizon_verified=False,production_eligible=False,
        information='Current preference persists; subsequent preference surprises not anticipated',
        demographic_information='Existing observed historical head-age bridge and frozen2023 person anchor retained')
    if not horizon:return SurpriseResult(final_path,receipt,None,None)
    next_state=first_period_state(inherited=inherited,old_state=old_state,demographics=demographic_primitives,
        path=final_path,prices=final['prices'],pensions=final['fiscal_values'],psi=psi)
    realized=dict(final_path.rows[0],forecast_vintage_year=inherited.year,
        expected_next_asset_price=float(final['prices'][1]),expected_constant_psi=float(psi))
    return SurpriseResult(final_path,receipt,next_state,realized)


def run_sequence(*, initial_state, dated_preferences, solve_episode, persist):
    """Sequential execution with explicit per-stage solve and persistence hooks.

    A solve callback sees only its inherited state and current preference, never
    later surprises. The persisted result contains the full forecast vintage.
    No fitted targets or optimizer are selected here.
    """
    events=list(dated_preferences)
    if not events or any(type(y) is not int or y!=initial_state.year+4*i or y>2023
                         or not np.isfinite(psi) for i,(y,psi) in enumerate(events)):
        raise ValueError('Explicit consecutive four-year surprise dates through2023 required')
    state=copy.deepcopy(initial_state);realized=[]
    for year,psi in events:
        try:
            result=solve_episode(inherited=copy.deepcopy(state),psi=float(psi))
        except Exception as exc:
            persist(year,SurpriseResult(None,dict(start_year=year,psi=float(psi),
                finite_horizon_market_fiscal_converged=False,terminal_distance_passed=False,
                error_type=type(exc).__name__,error=str(exc)),None,None))
            raise
        persist(year,result)
        if not result.root_receipt.get('finite_horizon_market_fiscal_converged') or not result.root_receipt.get('terminal_distance_passed') or result.next_state is None:
            raise RuntimeError(f'Surprise{year} is incomplete; no subsequent realized state accepted')
        if result.next_state.year!=year+4:raise RuntimeError('Must advance exactly one period')
        if result.realized_row is None or int(result.realized_row['calendar_year'])!=year:
            raise RuntimeError('Missing current-vintage observation')
        realized.append(copy.deepcopy(result.realized_row));state=copy.deepcopy(result.next_state)
    return dict(realized_rows=realized,final_state=state,historical_path_fitted=False,production_eligible=False)


def persist_episode(directory, year, result, *, provenance):
    """Save forecast vintage and full next state; no overwrite or automatic resume.

    Caller supplies its pinned scientific/input/target contract as provenance.
    No failed forecast creates an accepted next-state checkpoint.
    """
    import csv
    import gzip
    import hashlib
    import json
    import pickle
    from pathlib import Path
    if not provenance:raise ValueError('Pinned run provenance required')
    out=Path(directory)/str(year);out.mkdir(parents=True,exist_ok=False)
    def serial(value):
        if isinstance(value,np.ndarray):return value.tolist()
        if isinstance(value,np.generic):return value.item()
        raise TypeError(type(value).__name__)
    def clean(value):
        if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
        if isinstance(value,(list,tuple)):return [clean(v) for v in value]
        if isinstance(value,np.ndarray):return clean(value.tolist())
        if isinstance(value,np.generic):return clean(value.item())
        if isinstance(value,float) and not np.isfinite(value):return str(value)
        return value
    def save(name,value):
        (out/name).write_text(json.dumps(clean(value),default=serial,indent=2,allow_nan=False)+'\n')
    save('provenance.json',provenance);save('root_receipt.json',result.root_receipt)
    if result.path is not None and result.path.rows:
        keys=list(dict.fromkeys(k for row in result.path.rows for k in row))
        with (out/'expected_transition.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(result.path.rows)
    if result.next_state is not None:
        save('realized_period.json',result.realized_row)
        checkpoint=out/'next_state.pkl.gz'
        with gzip.open(checkpoint,'wb',compresslevel=1) as f:pickle.dump(result.next_state,f,protocol=5)
        with gzip.open(checkpoint,'rb') as f:reloaded=pickle.load(f)
        if reloaded.year!=year+4:raise RuntimeError('Checkpoint date verification failed')
        np.testing.assert_array_equal(reloaded.households.g_pre,result.next_state.households.g_pre)
        if hasattr(reloaded.households,'persons'):
            for name in ('persons','heads'):
                np.testing.assert_array_equal(getattr(reloaded.households.persons,name),getattr(result.next_state.households.persons,name))
        else:
            for name in ('scheduled_entries','scheduled_raw_entries'):
                np.testing.assert_array_equal(getattr(reloaded.households,name),getattr(result.next_state.households,name))
        with checkpoint.open('rb') as f:checksum=hashlib.file_digest(f,'sha256').hexdigest()
        save('checkpoint_verification.json',dict(year=reloaded.year,reloaded_exactly=True,sha256=checksum))
    return out
