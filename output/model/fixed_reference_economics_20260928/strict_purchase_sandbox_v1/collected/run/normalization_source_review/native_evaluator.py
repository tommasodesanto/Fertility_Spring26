def native_evaluator(out,lane,P,grid,deadline,price_start=None):
    coordinates=inputs.parameters(lane)
    arm=inputs.LANES[lane]['arm'];dims=inputs.LANES[lane]['dimensions']
    # Import the same indexed full-GE integration, never the default CLI.
    sys.path.insert(0,str(BASE))
    spec=importlib.util.spec_from_file_location('pilot_matched_workflow',BASE/'run_comparison.py')
    base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
    sys.path.insert(0,str(HERE))
    ge = _get_normalized_ge(out)
    ctx=base.authored.context_from_bundle(SimpleNamespace(bundle=ROOT/'output/model/publication_refactor_20260929/local_export_v1/inputs',reference_root=ROOT,out=out))
    install_reporter_on_authored(base.authored)
    base.authored.authenticate_frozen(ctx)
    ctx.update(P=P,b_grid=grid,selected_d_bar=inputs.ARMS[arm],reference_psi=float(P.psi_child),expected_dimensions={'wealth_grid_nodes':dims[0],'income_states':dims[1]},free_coordinates=list(coordinates),price_start=float(price_start) if price_start is not None else float(ctx['q_ref']),phase_b_max_new_lifecycle=32)
    from small_credit_lab.engine.shared import annual_gross_income_at_state
    from small_credit_lab.engine import solver
    sd=solver.precompute_shared(P,grid)
    require(float(sd.cb_flat[0,0])==float(sd.hb_flat[0,0])==float(sd.gb_flat[0,0])==0.,'Necessary cash preflight missing childless floors')
    actual_income=np.asarray([annual_gross_income_at_state(P,0,0,float(z)) for z in P.z_grid])
    np.testing.assert_array_equal(actual_income,P.income[0,0]*P.z_grid/P.period_years/(1-P.tau_pay))
    cal=ctx['prepared'].rt['primitive'].pf.calendar
    cohort=cal.entrant_cohort(np.asarray([1.]),P,grid)
    np.testing.assert_allclose(cohort.sum(axis=(1,2,4,5)),P.fixed_reference_entry_conditional*P.z_weights[None,:],rtol=0,atol=2e-16)
    require(abs(cohort.sum()-1)<2e-12,'Entrant mass changed')
    write(out/'native_entry_verification.json',dict(calendar_joint_maximum_error=float(abs(cohort.sum(axis=(1,2,4,5))-P.fixed_reference_entry_conditional*P.z_weights[None,:]).max()),income_formula_exact=True,childless_cash_floors_zero=True,lifecycle_solves=0))
    seed,bounds,_=inputs.seed_and_bounds(lane)
    for row in ctx['manifest']['full_parameter_table']:
        if row['parameter'] in bounds:
            row['lower'],row['upper']=map(str,bounds[row['parameter']])
    install_observer_metadata(ge,arm,bounds)
    from small_credit_lab import credit
    credit.bind_engine_credit(P,'corrected',0.)
    actual=ctx['fp'].actual_parameters(ctx['prepared'],P,grid)
    expected=expected_parameters(seed,dims,arm)
    ge.validate_parameter_estimates(dict(expected_parameters=expected),PLAN['reference_parameter_table'],actual)
    live_sd=solver.precompute_shared(P,grid)
    from refactor_lab.engine import solver as checked_solver
    checked_sd=checked_solver.precompute_shared(P,grid)
    for key in ('h_bar','c_bar','g_bar','alpha_flat','psi_v','escale_flat'):
        np.testing.assert_array_equal(getattr(live_sd,key),getattr(checked_sd,key))
    require(float(live_sd.h_bar[1,1])==expected['h_P'] and float(live_sd.h_bar[1,0])==0.,'Executed physical floor differs')
    import inspect
    from small_credit_lab.engine import shared as executed_shared,child_preferences as executed_child,kernels as executed_kernels
    from refactor_lab.engine import shared as checked_shared,child_preferences as checked_child,kernels as checked_kernels
    for checked,executed in [(checked_shared,executed_shared),(checked_child,executed_child),(checked_kernels,executed_kernels)]:
        require(inputs.sha(Path(checked.__file__))==inputs.sha(Path(executed.__file__)),'Executed preference source differs')
    write(out/'native_initializer_verification.json',dict(status='passed_zero_solve',actual_parameters=actual,expected_parameters=expected,authenticated_once=True,lifecycle_solves=0))
    def evaluate(label,point,evaluation_deadline):
        directory=out/label;directory.mkdir()
        candidate=dict(ctx);candidate.update(P=inputs.bind(P,point,bounds,arm),out=directory,expected_parameters=expected_parameters(point,dims,arm),deadline_epoch=evaluation_deadline)
        # No direct field may silently fail to propagate through the imported engine.
        from small_credit_lab import credit
        credit.bind_engine_credit(candidate['P'],'corrected',inputs.ARMS[arm])
        actual=candidate['fp'].actual_parameters(candidate['prepared'],candidate['P'],grid)
        ge.validate_parameter_estimates(candidate,PLAN['reference_parameter_table'],actual)
        budget=base.ArmBudget(directory,evaluation_deadline);budget.max_lifecycle=32
        write(directory/'proposed_parameters.json',dict(free=point,actual=actual,credit=inputs.ARMS[arm],starting_price=candidate['price_start'],target_contract_sha256=inputs.canonical(PLAN['target_contract'])))
        try:
            solved=ge.run_phase_b(candidate,dict(selected_d_bar=inputs.ARMS[arm]),budget)
            if solved['status']!='passed':
                exhausted=solved['status']=='uncomputed_bounded_budget' or solved.get('price_search',{}).get('termination_reason')=='budget_or_repeat_reserve'
                return dict(status='budget_exhausted' if exhausted else 'inadmissible_numerical',reason=solved['status'],lifecycle_solves=budget.used_lifecycle,price_search=solved.get('price_search'))
        except H0BoundError as exc:
            return dict(status="inadmissible_numerical",reason=str(exc),lifecycle_solves=budget.used_lifecycle,derived_H0_bound_rejection=True,rejection_kind="derived_H0_constraint")
        except TimeoutError as exc:
            if evaluation_deadline<deadline and time.time()>=evaluation_deadline:
                return dict(status='budget_exhausted',reason=str(exc),lifecycle_solves=budget.used_lifecycle)
            raise
        except RuntimeError as exc:
            if is_search_budget_exit(exc,budget,deadline):
                write(directory/'search_budget_exhausted.json',dict(reason=str(exc),lifecycle_solves=budget.used_lifecycle))
                return dict(status='budget_exhausted',reason=str(exc),lifecycle_solves=budget.used_lifecycle)
            # Only explicit failure to locate a root can reject an exploratory proposal.
            # Accounting, feasibility, source and parameter failures halt the arm.
            if label!='000_baseline' and str(exc).startswith('Renewal root unbracketed'):
                write(directory/'inadmissible.json',dict(reason=str(exc),lifecycle_solves=budget.used_lifecycle))
                return dict(status='inadmissible_numerical',reason=str(exc),lifecycle_solves=budget.used_lifecycle)
            raise
        report=directory/'phase_b_ge/selected_root'
        rows=readtable(report/'target_fit.csv');rr=residual(rows)
        # The native selected-price repeat verifies arrays/tables; verify actual PNGs too.
        repeat=directory/'phase_b_ge/selected_repeat_final'
        write(directory/'native_selected_repeat.json',compare_repeated(report,repeat))
        return dict(status='passed',residual=rr.tolist(),report=str(report),lifecycle_solves=budget.used_lifecycle,price=solved['selected_price'],starting_price=candidate['price_start'],population=solved['selected']['population_scale'])
    return evaluate
