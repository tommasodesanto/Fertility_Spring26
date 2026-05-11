%% run_model_cp_dt.m
%
%  DISCRETE-TIME CENTER-PERIPHERY LIFECYCLE MODEL (April 2026)
%  HIGH-PERFORMANCE VERSION
%
%  Optimizations:
%    1. Howard improvement: full grid search only on select iterations;
%       cheap policy evaluation on the rest.
%    2. All (nn,cs) batched into 3D array operations (Nb×Nb×nc).
%    3. CompEcon lookup MEX for fast interpolation (falls back to discretize).
%    4. accumarray for distribution savings advance.
%  ========================================================================

function [sol, P, p_eq] = run_model_cp_dt(P_override)

    t_start = tic;

    %% Setup
    P = setup_parameters();
    if nargin >= 1 && ~isempty(P_override)
        P = apply_overrides(P, P_override);
    end
    % beta is the primitive in DT; rho is derived for compatibility
    if ~isfield(P, 'beta') || isempty(P.beta)
        if isfield(P, 'rho') && ~isempty(P.rho)
            P.beta = 1/(1+P.rho);
        else
            P.beta = 0.96;
        end
    end
    P.rho = 1/P.beta - 1;  % derived, for any legacy code that reads rho
    P.rho_hat = P.rho;
    P.user_cost_rate = P.q + P.delta + P.tau_H;
    P.R_gross = 1+P.q;
    if isscalar(P.phi), P.phi = P.phi*ones(P.n_parity,1); end
    if isfield(P,'entry_init_override') && ~isempty(P.entry_init_override)
        e0=max(P.entry_init_override(:),0); e0=e0/sum(e0);
        P.entry_shares=e0; P.entry_by_loc=(1/P.J)*e0;
    end
    if ~isfield(P,'use_stochastic_aging'), P.use_stochastic_aging=false; end
    P = finalize_location_choice_spec(P);

    % lookup.mexmaci64 should be in this directory (local copy from CompEcon)
    P.has_lookup_mex = (exist('lookup','file') == 3);  % 3 = MEX
    if P.has_lookup_mex
        fprintf('  lookup MEX: AVAILABLE\n');
    end

    fprintf('============================================================\n');
    fprintf('  DT LIFECYCLE MODEL (Howard + vectorized)\n');
    fprintf('  J=%d, beta=%.4f, R=%.4f, Nb=%d\n', P.J, P.beta, P.R_gross, P.Nb);
    fprintf('  alpha=%.2f, sigma=%.2f, chi=%.2f, hR_max=%.1f\n', ...
        P.alpha_cons, P.sigma, P.chi, P.hR_max);
    fprintf('  kappa_loc=%.2f, kappa_fert=%.2f\n', P.kappa_loc, P.kappa_fert);
    if isfield(P,'location_choice_legacy_converted') && P.location_choice_legacy_converted
        fprintf('  location_choice_form=additive_due (converted legacy multiplicative inputs)\n');
    else
        fprintf('  location_choice_form=%s\n', P.location_choice_form);
    end
    fprintf('============================================================\n\n');

    %% Grid
    b_grid = make_grid(P);
    db_f = [diff(b_grid); b_grid(end)-b_grid(end-1)];

    %% Solve
    p_init = P.r_bar / P.user_cost_rate;
    if isfield(P,'p_init_override')&&~isempty(P.p_init_override), p_init=P.p_init_override(:); end

    solve_mode = 'ge';
    if isfield(P,'solve_mode') && ~isempty(P.solve_mode)
        solve_mode = lower(char(P.solve_mode));
    end
    do_pe = any(strcmp(solve_mode, {'pe','partial','partial_equilibrium','fixed'}));

    if do_pe
        fprintf('=== Solving DT Partial Equilibrium (fixed prices/wages/entry) ===\n\n');
        if ~isfield(P,'p_fixed') || ~isfield(P,'w_fixed') || ~isfield(P,'entry_shares_fixed')
            error('PE mode requires P_override fields: p_fixed, w_fixed, entry_shares_fixed.');
        end
        [p_eq, sol, P] = solve_partial_equilibrium_dt(P.p_fixed, P.w_fixed, ...
            P.entry_shares_fixed, P, b_grid);
    else
        fprintf('=== Solving DT Equilibrium (Howard) ===\n\n');
        [p_eq, sol, P] = solve_equilibrium(p_init, P, b_grid);
    end

    %% Display
    fprintf('\n==================== RESULTS ====================\n');
    fprintf('TFR=%.2f, Own=%.1f%%, Prices=[%.3f,%.3f]\n', ...
        2*sol.mean_parity, 100*sol.own_rate, p_eq(1), p_eq(2));
    fprintf('Pop=[%.2f,%.2f], Gradient=%.3f\n', ...
        sol.pop_share(1), sol.pop_share(2), ...
        sol.mean_parity_by_loc(1)-sol.mean_parity_by_loc(2));

    elapsed = toc(t_start);
    fprintf('Total: %.1f sec\n', elapsed);

    out_dir = fullfile(fileparts(mfilename('fullpath')), 'output');
    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    save(fullfile(out_dir,'last_run_dt.mat'), 'sol','P','p_eq', '-v7');
end


%% ========================================================================
%%   PARTIAL EQUILIBRIUM WITH FIXED PRICES/WAGES/ENTRY SHARES
%% ========================================================================
function [p_eq, sol, P] = solve_partial_equilibrium_dt(p_fixed, w_fixed, entry_shares_fixed, P, b_grid)
    P0 = P;
    p_eq = p_fixed(:);
    w = w_fixed(:);
    entry_shares = entry_shares_fixed(:);

    if length(p_eq) ~= P0.I
        error('p_fixed must have length I=%d.', P0.I);
    end
    if length(w) ~= P0.I
        error('w_fixed must have length I=%d.', P0.I);
    end
    if length(entry_shares) ~= P0.I
        error('entry_shares_fixed must have length I=%d.', P0.I);
    end
    if any(~isfinite(p_eq)) || any(p_eq <= 0)
        error('p_fixed must be finite and strictly positive.');
    end
    if any(~isfinite(w)) || any(w <= 0)
        error('w_fixed must be finite and strictly positive.');
    end

    entry_shares = max(entry_shares, 0);
    em = sum(entry_shares);
    if ~(em > 0)
        error('entry_shares_fixed must have positive mass.');
    end
    entry_shares = entry_shares / em;

    P0.w_hat = w;
    P0 = set_income_given_w_and_pension(P0);
    P0.entry_shares = entry_shares;
    P0.entry_by_loc = P0.E_total * entry_shares;
    P0.eq_iter = 1;

    r = P0.user_cost_rate * p_eq;
    SD = precompute_shared(P0, b_grid);
    [V,c_pol,hR_pol,bp_pol,tc,lp_j,fp,fv] = solve_bellman_full(r,p_eq,P0,b_grid,SD);
    [g,stats] = forward_distribution(bp_pol,hR_pol,tc,lp_j,fp,r,p_eq,P0,b_grid,SD);

    sol = pack_solution(V,c_pol,hR_pol,bp_pol,tc,lp_j,fp,fv,g,stats,P0.w_hat,p_eq);
    P = P0;
end


%% ========================================================================
%%   EQUILIBRIUM LOOP WITH HOWARD IMPROVEMENT
%%
%%   Strategy:
%%     - First few iterations: full Bellman (grid search, ~2.3 sec)
%%     - Later iterations: policy evaluation (skip grid search, ~0.5 sec)
%%       Reuse stored bp_pol; recompute V, tenure, loc, fert from new prices.
%%     - Re-optimize when convergence stalls.
%% ========================================================================
function [p_eq, sol, P] = solve_equilibrium(p_init, P, b_grid)
    P0 = P; p = p_init;
    entry_shares = P0.entry_shares(:);
    entry_shares = max(entry_shares,0);
    em=sum(entry_shares);
    if ~(em>0), entry_shares=ones(P0.I,1)/P0.I;
    else, entry_shares=entry_shares/em; end
    E_total = P0.E_total;
    P0.entry_shares=entry_shares; P0.entry_by_loc=E_total*entry_shares;

    % Damping state
    lam=P0.lambda_eq;
    lp=lam*ones(P0.I,1); le=lam*ones(P0.I,1);
    dp_prev=zeros(P0.I,1); de_prev=zeros(P0.I,1);
    pts=p; ets=entry_shares;

    best_err=inf; best_p=p; best_entry=entry_shares; best_snap=struct();
    t_full=0; t_eval=0; t_dist=0; n_full=0; n_eval=0; n_dist=0;

    % Precompute shared data that only depends on grid (not prices)
    SD = precompute_shared(P0, b_grid);

    % Howard: stored policies from last full solve
    stored_bp = [];

    % Howard schedule: full on 1-3, then every howard_freq, or when stalled
    howard_freq = 5;
    if isfield(P0,'howard_freq') && ~isempty(P0.howard_freq)
        howard_freq = max(1, round(P0.howard_freq));
    end
    force_full = isfield(P0,'force_full_bellman') && P0.force_full_bellman;
    stall_count = 0;
    prev_err = inf;

    fprintf('  Howard eq: lam=%.2f, tol=%.1e, freq=%d', lam, P0.tol_eq, howard_freq);
    if force_full, fprintf(' [full-only]'); end
    fprintf('\n');

    for iter = 1:P0.max_iter_eq
        r = P0.user_cost_rate * p;
        P0.eq_iter = iter;

        % Decide: full Bellman or policy evaluation?
        do_full = force_full || (iter <= 3) || isempty(stored_bp) || ...
                  (mod(iter, howard_freq) == 0) || (stall_count >= 2);

        if do_full
            tic_f = tic;
            [V,c_pol,hR_pol,bp_pol,tc,lp_j,fp,fv] = solve_bellman_full(r,p,P0,b_grid,SD);
            t_full = t_full+toc(tic_f); n_full=n_full+1;
            stored_bp = bp_pol;  % Store for Howard evaluation
            stall_count = 0;
            mode_str = 'F';
        else
            tic_e = tic;
            [V,c_pol,hR_pol,bp_pol,tc,lp_j,fp,fv] = solve_bellman_eval(stored_bp,r,p,P0,b_grid,SD);
            t_eval = t_eval+toc(tic_e); n_eval=n_eval+1;
            mode_str = 'E';
        end

        tic_d = tic;
        [g,stats] = forward_distribution(bp_pol,hR_pol,tc,lp_j,fp,r,p,P0,b_grid,SD,true);
        t_dist = t_dist+toc(tic_d); n_dist=n_dist+1;

        % Price targets
        Hd = stats.housing_demand;
        p_target = zeros(P0.I,1);
        Hs = zeros(P0.I,1);
        for i=1:P0.I
            hd=max(Hd(i),P0.housing_demand_floor_for_supply);
            p_target(i) = P0.r_bar(i)*(hd/P0.H0(i))^(1/P0.xi_supply(i)) / P0.user_cost_rate;
            Hs(i) = P0.H0(i)*(r(i)/P0.r_bar(i))^P0.xi_supply(i);
        end
        entry_target = stats.mature_entry_shares;

        % Smoothing
        if P0.adaptive_price_damping
            fw=P0.target_filter_weight;
            pts=fw*p_target+(1-fw)*pts; ets=fw*entry_target+(1-fw)*ets;
            pt=pts; et=ets;
        else, pt=p_target; et=entry_target; end

        err_p=max(abs(p_target-p)./max(abs(p),1e-6));
        err_e=max(abs(entry_target-entry_shares));
        qe=max(abs((Hd-Hs)./max(Hs,1e-6)));
        ov=max([err_p,err_e]);

        if ov<best_err
            best_err=ov; best_p=p; best_entry=entry_shares;
            % Store lightweight snapshot; full stats computed after convergence
            best_snap = struct('V',V,'c_pol',c_pol,'hR_pol',hR_pol,'bp_pol',bp_pol,...
                'tc',tc,'lp_j',lp_j,'fp',fp,'fv',fv,'g',g,'r',r,'p',p);
        end

        fprintf('  %s%3d: ep=%.4f ee=%.4f own=%.1f%% TFR=%.2f pop=[%.2f,%.2f]\n', ...
            mode_str, iter, err_p, err_e, 100*stats.own_rate, 2*stats.mean_parity, ...
            stats.pop_share(1), stats.pop_share(2));

        if ov<P0.tol_eq, fprintf('  *** CONVERGED %d iters ***\n',iter); break; end
        if iter>10 && best_err<=P0.tol_eq, fprintf('  *** NEAR-CONV ***\n'); break; end
        % Accept if best error is small and we've been cycling
        if iter>25 && best_err<10*P0.tol_eq
            fprintf('  *** ACCEPTED at %d iters (best_err=%.2e) ***\n',iter,best_err); break;
        end
        % Hard cap: if we've spent 50 iterations and best is reasonable, stop
        if iter>50 && best_err<0.02
            fprintf('  *** CAPPED at %d iters (best_err=%.2e) ***\n',iter,best_err); break;
        end

        % Stall detection for Howard
        if ov >= prev_err * 0.999, stall_count = stall_count+1;
        else, stall_count = 0; end
        prev_err = ov;

        % Damped update (same as before)
        if P0.adaptive_price_damping
            dp=pt-p; de=et-entry_shares;
            fp_mask=(abs(dp_prev)>1e-12)&(sign(dp)~=sign(dp_prev));
            fe_mask=(abs(de_prev)>1e-12)&(sign(de)~=sign(de_prev));
            dc=P0.adaptive_damping_decay; gr=P0.adaptive_damping_grow;
            lp(fp_mask)=max(P0.lambda_price_min,dc*lp(fp_mask));
            lp(~fp_mask)=min(P0.lambda_price_max,gr*lp(~fp_mask));
            le(fe_mask)=max(P0.lambda_entry_min,dc*le(fe_mask));
            le(~fe_mask)=min(P0.lambda_entry_max,gr*le(~fe_mask));
            if sum(fp_mask)+sum(fe_mask)>=2*P0.I
                lp=max(P0.lambda_price_min,P0.cycle_guard_factor*lp);
                le=max(P0.lambda_entry_min,P0.cycle_guard_factor*le);
            end
            if iter>15&&ov>1.05*best_err
                lp=max(P0.lambda_price_min,0.9*lp); le=max(P0.lambda_entry_min,0.9*le);
            end
            p=p+lp.*dp; dp_prev=dp;
            entry_shares=entry_shares+le.*de; de_prev=de;
        else
            p=p+lam*(pt-p); entry_shares=entry_shares+lam*(et-entry_shares);
        end
        if P0.enforce_price_bounds, p=max(min(p,P0.p_max),P0.p_min); end
        if P0.enforce_entry_share_floor, entry_shares=max(entry_shares,P0.entry_share_floor); end
        entry_shares=max(entry_shares,0);
        em=sum(entry_shares);
        if ~(em>0), entry_shares=ones(P0.I,1)/P0.I;
        else, entry_shares=entry_shares/em; end
        P0.entry_shares=entry_shares; P0.entry_by_loc=E_total*entry_shares;
    end
    p_eq=best_p;
    P=P0; P.entry_shares=best_entry; P.entry_by_loc=E_total*best_entry;
    fprintf('\n  [TIME] Full Bellman: %d calls, %.3fs (%.0f ms/call)\n', n_full, t_full, 1000*t_full/max(n_full,1));
    fprintf('  [TIME] Eval Bellman: %d calls, %.3fs (%.0f ms/call)\n', n_eval, t_eval, 1000*t_eval/max(n_eval,1));
    fprintf('  [TIME] Distribution: %d calls, %.3fs (%.0f ms/call)\n', n_dist, t_dist, 1000*t_dist/max(n_dist,1));
    % Compute full statistics once from best snapshot
    S = best_snap;
    [~,full_stats] = forward_distribution(S.bp_pol,S.hR_pol,S.tc,S.lp_j,S.fp,S.r,S.p,P,b_grid,SD,false);
    sol = pack_solution(S.V,S.c_pol,S.hR_pol,S.bp_pol,S.tc,S.lp_j,S.fp,S.fv,S.g,full_stats,P.w_hat,S.p);
end


%% ========================================================================
%%   PRECOMPUTE SHARED DATA (grid-dependent, price-independent)
%%
%%   Key optimization: identify the N_TYPES distinct (c_bar, h_bar, psi)
%%   tuples across the 21 (nn,cs) states. Utility only depends on the type,
%%   so we compute the expensive power operation once per type (not per nc).
%% ========================================================================
function SD = precompute_shared(P, b_grid)
    Nb = length(b_grid);
    nc = P.n_parity * P.n_child_states;
    K = P.n_child_stages; csm1 = K+2;

    c_bar = zeros(P.n_parity, P.n_child_states);
    h_bar = zeros(P.n_parity, P.n_child_states);
    psi_v = zeros(P.n_parity, P.n_child_states);
    for nn=1:P.n_parity
        nk=nn-1;
        for cs=1:P.n_child_states
            kp=(cs>1)&&(cs<csm1);
            if kp
                c_bar(nn,cs)=P.c_bar_0+P.c_bar_n*nk;
                switch lower(char(string(P.child_housing_spec)))
                    case "linear_only"
                        h_bar(nn,cs)=P.h_bar_0+P.h_bar_n*nk;
                    otherwise
                        h_bar(nn,cs)=P.h_bar_0+P.h_bar_jump+P.h_bar_n*nk;
                end
                psi_v(nn,cs)=P.psi_child*nk;
            else
                c_bar(nn,cs)=P.c_bar_0; h_bar(nn,cs)=P.h_bar_0;
            end
        end
    end

    SD.c_bar = c_bar;
    SD.h_bar = h_bar;
    SD.psi_v = psi_v;
    SD.cb_flat = reshape(c_bar, 1, nc);
    SD.hb_flat = reshape(h_bar, 1, nc);
    SD.psi_flat = reshape(psi_v, 1, nc);
    SD.nc = nc;
    SD.b = b_grid(:);
    SD.bp = b_grid(:)';
    SD.phi_state = get_phi_state_matrix(P);
    SD.phi_choice = get_phi_choice_tensor(P);

    % ---- TYPE FACTORING ----
    % Each nc column maps to a "type" defined by (c_bar + r*h_bar, psi).
    % For renters, the subsistence cost is d(nc) = c_bar(nc) + r * h_bar(nc),
    % and utility is K_r * (surplus - d(nc))^oms / oms + psi(nc).
    % With the same d and psi, the utility is identical → same grid search.
    %
    % Identify unique types by (c_bar, h_bar, psi) triple.
    triples = [c_bar(:), h_bar(:), psi_v(:)];
    [unique_triples, ~, type_map] = unique(triples, 'rows');
    SD.n_types = size(unique_triples, 1);
    SD.type_map = type_map;  % nc×1: which type each (nn,cs) uses
    SD.type_cb = unique_triples(:, 1);   % n_types×1
    SD.type_hb = unique_triples(:, 2);   % n_types×1
    SD.type_psi = unique_triples(:, 3);  % n_types×1

    % For each type, which nc columns belong to it
    SD.type_members = cell(SD.n_types, 1);
    for t = 1:SD.n_types
        SD.type_members{t} = find(type_map == t);
    end

    % Precompute birth-related ownership entry policies.
    nt = 1 + P.n_house;
    SD.birth_dp = false(P.n_parity, P.n_child_states, nt, nt);
    for nn=1:P.n_parity, for cs=1:P.n_child_states, for to=1:nt, for tn=1:nt
        SD.birth_dp(nn,cs,to,tn) = has_birth_dp_grant(P,nn,cs,to,tn);
    end, end, end, end
    SD.birth_entry_grant = get_birth_entry_grant_tensor(P);
end


%% ========================================================================
%%   FULL BELLMAN SOLVE (grid search — expensive, ~2 sec)
%% ========================================================================
function [V,c_pol,hR_pol,bp_pol,tenure_choice,loc_probs,fert_probs,fert_value] = ...
    solve_bellman_full(r_hat,p_hat,P,b_grid,SD)

    [V,c_pol,hR_pol,bp_pol,tenure_choice,loc_probs,fert_probs,fert_value] = ...
        solve_bellman_core(r_hat,p_hat,P,b_grid,SD, [], false);
end


%% ========================================================================
%%   POLICY EVALUATION (Howard — cheap, ~0.5 sec)
%%   Reuses stored bp_pol for savings; recomputes V, tenure, loc, fert.
%% ========================================================================
function [V,c_pol,hR_pol,bp_pol,tenure_choice,loc_probs,fert_probs,fert_value] = ...
    solve_bellman_eval(stored_bp,r_hat,p_hat,P,b_grid,SD)

    [V,c_pol,hR_pol,bp_pol,tenure_choice,loc_probs,fert_probs,fert_value] = ...
        solve_bellman_core(r_hat,p_hat,P,b_grid,SD, stored_bp, true);
end


%% ========================================================================
%%   CORE BELLMAN (shared by full and eval modes)
%% ========================================================================
function [V,c_pol,hR_pol,bp_pol,tenure_choice,loc_probs,fert_probs,fert_value] = ...
    solve_bellman_core(r_hat,p_hat,P,b_grid,SD, stored_bp, eval_mode)

    J=P.J; I=P.I; Nb=length(b_grid);
    nh=P.n_house; nt=1+nh; np=P.n_parity; ncs=P.n_child_states;
    nc=SD.nc; beta=P.beta; Rg=P.R_gross;
    sigma=P.sigma; alpha=P.alpha_cons; oms=1-sigma;
    b=SD.b; bp=SD.bp;

    V=zeros(Nb,nt,I,J,np,ncs);
    c_pol=zeros(Nb,nt,I,J,np,ncs);
    hR_pol=zeros(Nb,nt,I,J,np,ncs);
    bp_pol=ones(Nb,nt,I,J,np,ncs);
    tenure_choice=zeros(Nb,nt,I,J,np,ncs);
    loc_probs=zeros(Nb,nt,I,I,J,np,ncs);
    fert_probs=zeros(Nb,nt,I,J,np);
    fert_value=zeros(Nb,nt,I,J);

    % Price-dependent precomputation
    phi_choice=SD.phi_choice;
    birth_entry_grant=SD.birth_entry_grant;
    hcost=zeros(I,nt); heq=zeros(I,nt); dp_arr=zeros(I,nt,np,ncs); bmo=zeros(I,nt,np,ncs);
    hsrv=zeros(I,nt); ocst=zeros(I,nt);
    for i=1:I, for ten=2:nt
        hs=P.H_own(ten-1);
        hcost(i,ten)=p_hat(i)*hs; heq(i,ten)=(1-P.psi)*p_hat(i)*hs;
        hsrv(i,ten)=P.chi*hs; ocst(i,ten)=(P.delta+P.tau_H)*p_hat(i)*hs;
        for nn=1:np, for cs=1:ncs
            phi_ncs=phi_choice(i,ten,nn,cs);
            dp_arr(i,ten,nn,cs)=(1-phi_ncs)*hcost(i,ten);
            bmo(i,ten,nn,cs)=-phi_ncs*hcost(i,ten);
        end, end
    end, end

    % Bequest
    Vbq=zeros(Nb,nt,I,np,ncs);
    for i=1:I, for ten=1:nt
        hv=0; if ten>1, hv=p_hat(i)*P.H_own(ten-1); end
        for nn=1:np, for cs=1:ncs
            nk=get_completed_fertility(nn,cs,P);
            Vbq(:,ten,i,nn,cs)=bequest_utility_vec(b_grid+hv,nk,P);
        end, end
    end, end

    % DUE-style additive residential amenity and moving cost shifters
    loc_shift=zeros(I,I);
    for io=1:I, for id=1:I
        if id==io, move_cost=P.mu_stay;
        else, move_cost=P.mu_move; end
        loc_shift(io,id)=P.E_loc(id)-move_cost;
    end, end

    % Interpolation indices for movers (fast_lookup or discretize)
    iidx=zeros(Nb,I,nt); iwt=zeros(Nb,I,nt);
    for io=1:I, for to=1:nt
        ba=max(min(b_grid+heq(io,to),b_grid(end)),b_grid(1));
        [idx,wt] = fast_interp_idx(b_grid, ba, P.has_lookup_mex);
        iidx(:,io,to)=idx; iwt(:,io,to)=wt;
    end, end

    % Pre-extract type factoring data
    n_types = SD.n_types;
    type_map = SD.type_map;   % nc×1
    type_cb = SD.type_cb;     % n_types×1
    type_hb = SD.type_hb;
    type_psi = SD.type_psi;
    type_members = SD.type_members;

    % Pre-allocate work arrays outside the age loop
    Vd=zeros(Nb,nt,I,np,ncs);
    cd=zeros(Nb,nt,I,np,ncs);
    hd=zeros(Nb,nt,I,np,ncs);
    bd=zeros(Nb,nt,I,np,ncs);  % continuous b' values (not indices)
    Vo_nc = zeros(Nb, nc);
    bp_nc = zeros(Nb, nc);
    co_nc = zeros(Nb, nc);
    ho_nc = zeros(Nb, nc);

    % Golden section constants
    gs_tol = 1e-3;
    gs_alpha1 = (3 - sqrt(5))/2;
    gs_alpha2 = (sqrt(5) - 1)/2;

    % Column offset for vectorized linear indexing into Nb×nc arrays
    col_offset = Nb*(0:nc-1);  % 1×nc
    b_lo = b_grid(1); b_hi = b_grid(end);

    % ==================== BACKWARD INDUCTION ====================
    for j=J:-1:1
        in_fert=(j>=P.A_f_start)&&(j<=P.A_f_end);

        % Continuation (after child aging)
        if j==J, Vnr=Vbq;
        else, Vnr=reshape(V(:,:,:,j+1,:,:),Nb,nt,I,np,ncs); end
        Vc=apply_child_aging(Vnr,P,Nb,nt,I,np,ncs);

        % ============ STEP 1: SAVINGS (golden section on interpolated V) ============
        for i=1:I
            yj=P.income(i,j); ri=r_hat(i);
            Rv=Rg*b+yj;  % Nb×1

            if ~eval_mode
                % ---- FULL SOLVE: column-by-column GS with griddedInterpolant ----
                % (Midrigan architecture: Vbar constructed once, GS on Nb vectors)
                hRmax = P.hR_max;
                imeth = P.interp_method;

                Kr = (alpha^alpha * ((1-alpha)/ri)^(1-alpha))^oms;
                aoms = alpha*oms;

                Vcr = reshape(Vc(:,1,i,:,:), Nb, nc);
                d_nc = SD.cb_flat + ri * SD.hb_flat;
                cap_nc = ri * (hRmax - SD.hb_flat) / (1-alpha);

                % Warm-start bounds
                if j < J
                    bp_prev_r = reshape(bd(:,1,i,:,:), Nb, nc);
                end

                for c = 1:nc
                    Vbar = griddedInterpolant(b_grid, Vcr(:,c), imeth, 'linear');
                    dc = d_nc(c); pc = SD.psi_flat(c); cc = cap_nc(c);
                    cb_c = SD.cb_flat(c); hb_c = SD.hb_flat(c);
                    ht_cap_c = max(hRmax - hb_c, 1e-10);

                    lo = max(zeros(Nb,1), b_grid(1));
                    hi = max(Rv - dc - 1e-6, 0);
                    if j < J
                        lo = max(lo, bp_prev_r(:,c) - 2.0);
                        hi = min(hi, bp_prev_r(:,c) + 2.0);
                        lo = max(lo, 0); hi = max(hi, lo);
                    end

                    d = hi - lo;
                    x1 = lo + gs_alpha1*d; x2 = lo + gs_alpha2*d;
                    % Inline renter eval (Nb-vector, ~640 bytes — fits in L1)
                    s1 = Rv - dc - x1; f1 = Kr*max(s1,1e-10).^oms/oms + pc + beta*Vbar(max(min(x1,b_grid(end)),b_grid(1)));
                    cm1 = s1 > cc; if any(cm1); ct1=max(Rv(cm1)-cb_c-ri*hRmax-x1(cm1),1e-10); f1(cm1)=(ct1.^alpha.*ht_cap_c^(1-alpha)).^oms/oms+pc+beta*Vbar(max(min(x1(cm1),b_grid(end)),b_grid(1))); end
                    f1(s1<=1e-10) = -1e10;
                    s2 = Rv - dc - x2; f2 = Kr*max(s2,1e-10).^oms/oms + pc + beta*Vbar(max(min(x2,b_grid(end)),b_grid(1)));
                    cm2 = s2 > cc; if any(cm2); ct2=max(Rv(cm2)-cb_c-ri*hRmax-x2(cm2),1e-10); f2(cm2)=(ct2.^alpha.*ht_cap_c^(1-alpha)).^oms/oms+pc+beta*Vbar(max(min(x2(cm2),b_grid(end)),b_grid(1))); end
                    f2(s2<=1e-10) = -1e10;
                    d = gs_alpha1*gs_alpha2*d;

                    while any(d > gs_tol)
                        bt = f2>=f1;
                        xe = bt.*(x2+d) + (~bt).*(x1-d);
                        se = Rv - dc - xe; fe = Kr*max(se,1e-10).^oms/oms + pc + beta*Vbar(max(min(xe,b_grid(end)),b_grid(1)));
                        cme = se > cc; if any(cme); cte=max(Rv(cme)-cb_c-ri*hRmax-xe(cme),1e-10); fe(cme)=(cte.^alpha.*ht_cap_c^(1-alpha)).^oms/oms+pc+beta*Vbar(max(min(xe(cme),b_grid(end)),b_grid(1))); end
                        fe(se<=1e-10) = -1e10;
                        x1n=bt.*x2+(~bt).*xe; f1n=bt.*f2+(~bt).*fe;
                        x2n=bt.*xe+(~bt).*x1; f2n=bt.*fe+(~bt).*f1;
                        d=d*gs_alpha2; x1=x1n;x2=x2n;f1=f1n;f2=f2n;
                    end
                    bt=f2>=f1; bp_nc(:,c)=bt.*x2+(~bt).*x1; Vo_nc(:,c)=max(f1,f2);
                end

                % Policies from optimal bp
                surplus_nc = Rv - d_nc - bp_nc;
                ct_nc = alpha * max(surplus_nc, 1e-10);
                ht_nc = (1-alpha)/ri * max(surplus_nc, 1e-10);
                cm = (SD.hb_flat + ht_nc) > hRmax;
                if any(cm(:))
                    ct_cap = max(Rv - SD.cb_flat - ri*hRmax - bp_nc, 1e-10);
                    ct_nc(cm) = ct_cap(cm);
                    hRmax_rep = repmat(max(hRmax - SD.hb_flat, 1e-10), Nb, 1);
                    ht_nc(cm) = hRmax_rep(cm);
                end
                co_nc = SD.cb_flat + max(ct_nc, P.c_min);
                ho_nc = SD.hb_flat + max(ht_nc, 0.01);
                bad = surplus_nc <= 1e-10;
                co_nc(bad) = P.c_bar_0+P.c_min; ho_nc(bad) = P.h_bar_0+0.01;

                Vd(:,1,i,:,:) = reshape(Vo_nc, Nb,1,1,np,ncs);
                bd(:,1,i,:,:) = reshape(bp_nc, Nb,1,1,np,ncs);
                cd(:,1,i,:,:) = reshape(co_nc, Nb,1,1,np,ncs);
                hd(:,1,i,:,:) = reshape(ho_nc, Nb,1,1,np,ncs);

                % OWNERS: column-by-column GS with griddedInterpolant
                for ten=2:nt
                    oc=ocst(i,ten); hsv=hsrv(i,ten);
                    Vco = reshape(Vc(:,ten,i,:,:), Nb, nc);

                    if j < J
                        bp_prev_o = reshape(bd(:,ten,i,:,:), Nb, nc);
                    end

                    for c = 1:nc
                        Vbar = griddedInterpolant(b_grid, Vco(:,c), imeth, 'linear');
                        cb_c = SD.cb_flat(c); pc = SD.psi_flat(c);
                        ht_c = max(hsv - SD.hb_flat(c), 1e-10);
                        Ko_c = ht_c^((1-alpha)*oms);
                        nn_c = ceil(c/ncs); cs_c = c-(nn_c-1)*ncs;
                        bf_c = bmo(i,ten,nn_c,cs_c);

                        lo = max(bf_c, b_grid(1)) * ones(Nb,1);
                        hi = max(Rv - oc - cb_c - 1e-6, lo);
                        if j < J
                            lo = max(lo, bp_prev_o(:,c)-2.0);
                            hi = min(hi, bp_prev_o(:,c)+2.0);
                            lo = max(lo, bf_c); hi = max(hi, lo);
                        end

                        d = hi - lo;
                        x1 = lo + gs_alpha1*d; x2 = lo + gs_alpha2*d;
                        ct1=max(Rv-oc-cb_c-x1,1e-10); f1=Ko_c*ct1.^aoms/oms+pc+beta*Vbar(max(min(x1,b_grid(end)),b_grid(1))); f1(ct1<=1e-10)=-1e10;
                        ct2=max(Rv-oc-cb_c-x2,1e-10); f2=Ko_c*ct2.^aoms/oms+pc+beta*Vbar(max(min(x2,b_grid(end)),b_grid(1))); f2(ct2<=1e-10)=-1e10;
                        d = gs_alpha1*gs_alpha2*d;

                        while any(d > gs_tol)
                            bt=f2>=f1;
                            xe=bt.*(x2+d)+(~bt).*(x1-d);
                            cte=max(Rv-oc-cb_c-xe,1e-10); fe=Ko_c*cte.^aoms/oms+pc+beta*Vbar(max(min(xe,b_grid(end)),b_grid(1))); fe(cte<=1e-10)=-1e10;
                            x1n=bt.*x2+(~bt).*xe; f1n=bt.*f2+(~bt).*fe;
                            x2n=bt.*xe+(~bt).*x1; f2n=bt.*fe+(~bt).*f1;
                            d=d*gs_alpha2; x1=x1n;x2=x2n;f1=f1n;f2=f2n;
                        end
                        bt=f2>=f1; bp_nc(:,c)=bt.*x2+(~bt).*x1; Vo_nc(:,c)=max(f1,f2);
                    end
                    co_nc = SD.cb_flat + max(Rv - oc - SD.cb_flat - bp_nc, P.c_min);

                    Vd(:,ten,i,:,:) = reshape(Vo_nc, Nb,1,1,np,ncs);
                    bd(:,ten,i,:,:) = reshape(bp_nc, Nb,1,1,np,ncs);
                    cd(:,ten,i,:,:) = reshape(co_nc, Nb,1,1,np,ncs);
                end
            else
                % ---- HOWARD EVALUATION ----
                % Renters: vectorized batch interpolation (no col loop)
                bpv=reshape(stored_bp(:,1,i,j,:,:),Nb,nc);
                Vcr_nc=reshape(Vc(:,1,i,:,:),Nb,nc);
                bpv_c = max(min(bpv, b_hi), b_lo);
                if P.has_lookup_mex
                    ev_idx = lookup(b_grid, bpv_c, 3);
                else
                    ev_idx = max(min(discretize(bpv_c, [b_grid; inf]), Nb-1), 1);
                end
                ev_w = (bpv_c - b_grid(ev_idx)) ./ (b_grid(ev_idx+1) - b_grid(ev_idx));
                ev_w = max(min(ev_w, 1), 0);
                Vcbp = (1-ev_w).*Vcr_nc(ev_idx + col_offset) + ev_w.*Vcr_nc(ev_idx+1 + col_offset);

                % Vectorized surplus and utility (Kr shortcut: 1 power op instead of 2)
                Kr_ev = (alpha^alpha * ((1-alpha)/ri)^(1-alpha))^oms;
                d_nc = SD.cb_flat + ri * SD.hb_flat;
                cap_nc_ev = ri * (P.hR_max - SD.hb_flat) / (1-alpha);
                surplus_nc = Rv - d_nc - bpv;
                ss = max(surplus_nc, 1e-10);
                Vo_nc = Kr_ev * ss.^oms / oms + SD.psi_flat + beta*Vcbp;
                % Cap correction
                cm = surplus_nc > cap_nc_ev;
                if any(cm(:))
                    ct_cap = max(Rv - SD.cb_flat - ri*P.hR_max - bpv, 1e-10);
                    ht_cap_v = repmat(max(P.hR_max - SD.hb_flat, 1e-10), Nb, 1);
                    comp_cap = ct_cap.^alpha .* ht_cap_v.^(1-alpha);
                    u_cap = comp_cap.^oms / oms + SD.psi_flat;
                    Vo_nc(cm) = u_cap(cm) + beta*Vcbp(cm);
                end
                Vo_nc(surplus_nc <= 1e-10) = -1e10;
                ct_nc = alpha * max(surplus_nc, 1e-10);
                ht_nc = (1-alpha)/ri * max(surplus_nc, 1e-10);
                if any(cm(:))
                    ct_nc(cm) = ct_cap(cm);
                    ht_nc(cm) = ht_cap_v(cm);
                end
                co_nc = SD.cb_flat + max(ct_nc, P.c_min);
                ho_nc = SD.hb_flat + max(ht_nc, 0.01);
                bad = surplus_nc <= 1e-10;
                co_nc(bad) = P.c_bar_0+P.c_min; ho_nc(bad) = P.h_bar_0+0.01;

                Vd(:,1,i,:,:)=reshape(Vo_nc,Nb,1,1,np,ncs);
                bd(:,1,i,:,:)=reshape(bpv,Nb,1,1,np,ncs);
                cd(:,1,i,:,:)=reshape(co_nc,Nb,1,1,np,ncs);
                hd(:,1,i,:,:)=reshape(ho_nc,Nb,1,1,np,ncs);

                % Owners: vectorized batch interpolation (no col loop)
                for ten=2:nt
                    oc=ocst(i,ten); hsv=hsrv(i,ten);
                    bpv_o=reshape(stored_bp(:,ten,i,j,:,:),Nb,nc);
                    Vco_nc=reshape(Vc(:,ten,i,:,:),Nb,nc);
                    bpv_c = max(min(bpv_o, b_hi), b_lo);
                    if P.has_lookup_mex
                        ev_idx = lookup(b_grid, bpv_c, 3);
                    else
                        ev_idx = max(min(discretize(bpv_c, [b_grid; inf]), Nb-1), 1);
                    end
                    ev_w = (bpv_c - b_grid(ev_idx)) ./ (b_grid(ev_idx+1) - b_grid(ev_idx));
                    ev_w = max(min(ev_w, 1), 0);
                    Vcbpo = (1-ev_w).*Vco_nc(ev_idx + col_offset) + ev_w.*Vco_nc(ev_idx+1 + col_offset);
                    ct_o = max(Rv - oc - SD.cb_flat - bpv_o, 1e-10);
                    ht_o = max(hsv - SD.hb_flat, 1e-10);
                    Ko_ev = ht_o.^((1-alpha)*oms);
                    Vo_nc = Ko_ev .* ct_o.^(alpha*oms) / oms + SD.psi_flat + beta*Vcbpo;
                    Vo_nc(ct_o <= 1e-10) = -1e10;
                    co_nc = SD.cb_flat + max(ct_o, P.c_min);
                    Vd(:,ten,i,:,:)=reshape(Vo_nc,Nb,1,1,np,ncs);
                    bd(:,ten,i,:,:)=reshape(bpv_o,Nb,1,1,np,ncs);
                    cd(:,ten,i,:,:)=reshape(co_nc,Nb,1,1,np,ncs);
                end
            end
        end
        c_pol(:,:,:,j,:,:)=cd;
        hR_pol(:,:,:,j,:,:)=hd;
        bp_pol(:,:,:,j,:,:)=bd;

        % ============ STEP 2: TENURE (same for full and eval) ============
        VH=zeros(Nb,nt,I,np,ncs); tcj=zeros(Nb,nt,I,np,ncs);
        for id=1:I, for to=1:nt
            sp=0; if to>1, sp=heq(id,to); end
            Vopt=zeros(Nb,np,ncs,nt);
            if to==1, Vopt(:,:,:,1)=reshape(Vd(:,1,id,:,:),Nb,np,ncs);
            else
                ba=max(b_grid+sp,0);
                Vopt(:,:,:,1)=interp1_fast(b_grid,reshape(Vd(:,1,id,:,:),Nb,np,ncs),ba,P.has_lookup_mex);
            end
            for tn=2:nt
                hc=hcost(id,tn);
                Vow=reshape(Vd(:,tn,id,:,:),Nb,np,ncs);
                if to==tn, Vopt(:,:,:,tn)=Vow;
                elseif to==1
                    bab=b_grid-hc;
                    Vb=interp1_fast(b_grid,Vow,bab,P.has_lookup_mex);
                    for nn=1:np, for cs=1:ncs
                        dpn=dp_arr(id,tn,nn,cs); bmn=bmo(id,tn,nn,cs);
                        if SD.birth_dp(nn,cs,to,tn)
                            bag=max(bab,bmn);
                            Vb(:,nn,cs)=interp1_fast(b_grid,Vow(:,nn,cs),bag,P.has_lookup_mex);
                        elseif birth_entry_grant(id,tn,nn,cs)>0
                            gfix=birth_entry_grant(id,tn,nn,cs);
                            babg=bab+gfix;
                            Vg=interp1_fast(b_grid,Vow(:,nn,cs),babg,P.has_lookup_mex);
                            inf_m=((b_grid+gfix)<dpn)|(babg<bmn);
                            Vg(inf_m)=-1e10;
                            Vb(:,nn,cs)=Vg;
                        else
                            inf_m=(b_grid<dpn)|(bab<bmn); Vb(inf_m,nn,cs)=-1e10;
                        end
                    end, end
                    Vopt(:,:,:,tn)=Vb;
                else
                    bar=b_grid+sp-hc;
                    Vrs=interp1_fast(b_grid,Vow,bar,P.has_lookup_mex);
                    for nn=1:np, for cs=1:ncs
                        dpn=dp_arr(id,tn,nn,cs); bmn=bmo(id,tn,nn,cs);
                        dpc=dpn-sp; inf_m=(b_grid<dpc)|(bar<bmn);
                        Vrs(inf_m,nn,cs)=-1e10;
                    end, end
                    Vopt(:,:,:,tn)=Vrs;
                end
            end
            [Vm,to2]=max(Vopt,[],4); VH(:,to,id,:,:)=Vm; tcj(:,to,id,:,:)=to2;
        end, end
        tenure_choice(:,:,:,j,:,:)=tcj;

        % ============ STEP 3: LOCATION (logit) ============
        VI=zeros(Nb,nt,I,np,ncs); lpj=zeros(Nb,nt,I,I,np,ncs);
        kl=P.kappa_loc;
        for io=1:I, for to=1:nt
            Va=zeros(Nb,I,np,ncs);
            Va(:,io,:,:)=reshape(VH(:,to,io,:,:),Nb,np,ncs);
            idx=iidx(:,io,to); wt=iwt(:,io,to);
            for id=1:I
                if id==io, continue; end
                Vdst=reshape(VH(:,1,id,:,:),Nb,np,ncs);
                Va(:,id,:,:)=(1-wt).*Vdst(idx,:,:)+wt.*Vdst(idx+1,:,:);
            end
            la=Va;
            for id=1:I, la(:,id,:,:)=la(:,id,:,:)+loc_shift(io,id); end
            la=la/kl;
            lm=max(la,[],2); se=sum(exp(la-lm),2); ls=lm+log(se);
            pr=exp(la-ls);
            VI(:,to,io,:,:)=reshape(kl*ls,Nb,np,ncs);
            lpj(:,to,io,:,:,:)=pr;
        end, end
        loc_probs(:,:,:,:,j,:,:)=lpj;

        % ============ STEP 4: FERTILITY (logit) ============
        if in_fert
            Vfa=zeros(Nb,nt,I,np);
            Vfa(:,:,:,1)=VI(:,:,:,1,1);
            for nn=2:np, Vfa(:,:,:,nn)=VI(:,:,:,nn,2); end
            kf=P.kappa_fert;
            lf=Vfa/kf; lm=max(lf,[],4); se=sum(exp(lf-lm),4); ls=lm+log(se);
            fert_probs(:,:,:,j,:)=exp(lf-ls);
            fert_value(:,:,:,j)=kf*ls;
            V(:,:,:,j,1,1)=fert_value(:,:,:,j);
            V(:,:,:,j,2:end,:)=VI(:,:,:,2:end,:);
            V(:,:,:,j,1,2:end)=VI(:,:,:,1,2:end);
        else
            V(:,:,:,j,:,:)=VI;
        end
    end
end


%% ========================================================================
%%   FORWARD DISTRIBUTION
%% ========================================================================
function [g,stats] = forward_distribution(bp_pol,hR_pol,tenure_choice,...
    loc_probs,fert_probs,r_hat,p_hat,P,b_grid,SD,fast_stats)

    J=P.J; I=P.I; Nb=length(b_grid); nt=1+P.n_house;
    np=P.n_parity; ncs=P.n_child_states; nc=SD.nc;
    bmin=b_grid(1); bmax=b_grid(end);
    g=zeros(Nb,nt,I,J,np,ncs);
    be=P.b_entry_fixed; ie=find(b_grid>=be,1,'first'); if isempty(ie),ie=1;end
    for i=1:I, g(ie,1,i,1,1,1)=P.entry_by_loc(i); end
    tb=0; blc=zeros(I,1); eml=zeros(I,1); emt=0;

    hc=zeros(I,nt); he=zeros(I,nt);
    for i=1:I, for ten=2:nt
        hs=P.H_own(ten-1); hc(i,ten)=p_hat(i)*hs; he(i,ten)=(1-P.psi)*p_hat(i)*hs;
    end, end
    phi_choice=SD.phi_choice;
    use_lin=strcmpi(P.kfe_wealth_interp,'linear');

    % Precompute matrices
    lmm=cell(I,nt);
    for io=1:I, for to=1:nt
        ba=max(min(b_grid+he(io,to),bmax),bmin);
        if use_lin, lmm{io,to}=make_lin_redist(b_grid,ba);
        else, [~,di]=min(abs(b_grid-ba'),[],1);
            lmm{io,to}=sparse(di(:),(1:Nb)',ones(Nb,1),Nb,Nb); end
    end, end

    tmx=cell(I,nt,nt,np,ncs);
    for nn=1:np, for cs=1:ncs
        for id=1:I, for to=1:nt, sp=he(id,to);
            for tn=1:nt
                pn=1.0;
                if tn > 1
                    pn=phi_choice(id,tn,nn,cs);
                end
                if tn==to, tmx{id,to,tn,nn,cs}=speye(Nb);
                elseif to==1&&tn>1
                    bf=max(min(max(b_grid-hc(id,tn),-pn*hc(id,tn)),bmax),bmin);
                    tmx{id,to,tn,nn,cs}=make_redist(b_grid,bf,use_lin);
                elseif to>1&&tn==1
                    bf=max(min(max(b_grid+sp,0),bmax),bmin);
                    tmx{id,to,tn,nn,cs}=make_redist(b_grid,bf,use_lin);
                else
                    bf=max(min(max(b_grid+sp-hc(id,tn),-pn*hc(id,tn)),bmax),bmin);
                    tmx{id,to,tn,nn,cs}=make_redist(b_grid,bf,use_lin);
                end
            end
        end, end
    end, end

    K=P.n_child_stages; csm1=K+2; csm2=K+3;
    ust=P.use_stochastic_aging&&isfield(P,'Pi_child');
    if ust, Pia=P.Pi_child; end

    % Event-study style housing moments around first birth.
    % The current one-shot fertility structure can support a clean
    % first-birth +3 housing response, but not a literal second-birth event.
    birth_es3_pre_sum = 0;
    birth_es3_post_sum = 0;
    birth_es3_mass = 0;
    addchild_es3_one_sum = 0;
    addchild_es3_one_mass = 0;
    addchild_es3_two_plus_sum = 0;
    addchild_es3_two_plus_mass = 0;
    onechild_es3_pre_sum = 0;
    onechild_es3_post_sum = 0;
    onechild_es3_mass = 0;
    twoplus_es3_pre_sum = 0;
    twoplus_es3_post_sum = 0;
    twoplus_es3_mass = 0;
    event_horizon = 3;

    for j=1:J-1
        % Fertility
        if (j>=P.A_f_start)&&(j<=P.A_f_end)
            gc=g(:,:,:,j,1,1);
            pa=reshape(fert_probs(:,:,:,j,:),Nb,nt,I,np);
            pv=reshape(0:np-1,1,1,1,np);
            ba_j=gc.*sum(pa.*pv,4);
            tb=tb+sum(ba_j,'all');
            for i=1:I, blc(i)=blc(i)+sum(ba_j(:,:,i),'all'); end
            mbp=gc.*pa;
            gpf=zeros(Nb,nt,I,np,ncs);
            gpf(:,:,:,1,1)=mbp(:,:,:,1); gpf(:,:,:,2:end,2)=mbp(:,:,:,2:end);
            g(:,:,:,j,1,1)=0;
            for nn=1:np, for cs=1:ncs
                g(:,:,:,j,nn,cs)=g(:,:,:,j,nn,cs)+gpf(:,:,:,nn,cs);
            end, end

            % First-birth event-study analog:
            % compare pre-birth housing of treated childless states to mean
            % housing three periods after the birth cohort is created.
            birth_mass = sum(mbp(:,:,:,2:end), 4);
            birth_mass_total = sum(birth_mass, 'all');
            if birth_mass_total > 1e-12 && (j + event_horizon) <= J
                pre_h = mean_housing_childless_weighted(birth_mass, j, hR_pol, P);
                birth_cohort = zeros(Nb,nt,I,np,ncs);
                birth_cohort(:,:,:,2:end,2) = mbp(:,:,:,2:end);
                birth_cohort = advance_cohort_horizon(birth_cohort, j, event_horizon, ...
                    loc_probs, tenure_choice, bp_pol, P, b_grid, SD, lmm, tmx, ust, Pia);
                post_h = mean_housing_distribution(birth_cohort, j + event_horizon, hR_pol, P);

                birth_es3_pre_sum = birth_es3_pre_sum + birth_mass_total * pre_h;
                birth_es3_post_sum = birth_es3_post_sum + birth_mass_total * post_h;
                birth_es3_mass = birth_es3_mass + birth_mass_total;

                one_child = zeros(Nb,nt,I,np,ncs);
                one_child(:,:,:,2,:) = birth_cohort(:,:,:,2,:);
                one_mass = sum(one_child, 'all');
                if one_mass > 1e-12
                    one_child_birth_mass = mbp(:,:,:,2);
                    one_child_birth_mass_total = sum(one_child_birth_mass, 'all');
                    addchild_es3_one_sum = addchild_es3_one_sum + ...
                        one_mass * mean_housing_distribution(one_child, j + event_horizon, hR_pol, P);
                    addchild_es3_one_mass = addchild_es3_one_mass + one_mass;
                    if one_child_birth_mass_total > 1e-12
                        one_child_pre_h = mean_housing_childless_weighted(one_child_birth_mass, j, hR_pol, P);
                        onechild_es3_pre_sum = onechild_es3_pre_sum + ...
                            one_child_birth_mass_total * one_child_pre_h;
                        onechild_es3_post_sum = onechild_es3_post_sum + ...
                            one_mass * mean_housing_distribution(one_child, j + event_horizon, hR_pol, P);
                        onechild_es3_mass = onechild_es3_mass + one_child_birth_mass_total;
                    end
                end

                if np >= 3
                    two_plus_birth_mass = sum(mbp(:,:,:,3:end), 4);
                    two_plus_birth_mass_total = sum(two_plus_birth_mass, 'all');
                    two_plus = zeros(Nb,nt,I,np,ncs);
                    two_plus(:,:,:,3:end,:) = birth_cohort(:,:,:,3:end,:);
                    two_plus_mass = sum(two_plus, 'all');
                    if two_plus_mass > 1e-12
                        addchild_es3_two_plus_sum = addchild_es3_two_plus_sum + ...
                            two_plus_mass * mean_housing_distribution(two_plus, j + event_horizon, hR_pol, P);
                        addchild_es3_two_plus_mass = addchild_es3_two_plus_mass + two_plus_mass;
                        if two_plus_birth_mass_total > 1e-12
                            two_plus_pre_h = mean_housing_childless_weighted(two_plus_birth_mass, j, hR_pol, P);
                            twoplus_es3_pre_sum = twoplus_es3_pre_sum + ...
                                two_plus_birth_mass_total * two_plus_pre_h;
                            twoplus_es3_post_sum = twoplus_es3_post_sum + ...
                                two_plus_mass * mean_housing_distribution(two_plus, j + event_horizon, hR_pol, P);
                            twoplus_es3_mass = twoplus_es3_mass + two_plus_birth_mass_total;
                        end
                    end
                end
            end
        end

        gj=g(:,:,:,j,:,:);
        % Location
        gpl=zeros(Nb,nt,I,np,ncs);
        for io=1:I, for to=1:nt
            go=reshape(gj(:,to,io,:,:),Nb,np*ncs);
            if sum(go,'all')<1e-15, continue; end
            po=reshape(loc_probs(:,to,io,:,j,:,:),Nb,I,np*ncs);
            sp=reshape(po(:,io,:),Nb,np*ncs);
            gpl(:,to,io,:,:)=gpl(:,to,io,:,:)+reshape(go.*sp,Nb,1,1,np,ncs);
            R=lmm{io,to};
            for id=1:I
                if id==io, continue; end
                mp=reshape(po(:,id,:),Nb,np*ncs);
                gpl(:,1,id,:,:)=gpl(:,1,id,:,:)+reshape(R*(go.*mp),Nb,1,1,np,ncs);
            end
        end, end

        % Tenure
        gpt=zeros(Nb,nt,I,np,ncs);
        for nn=1:np, for id=1:I, for to=1:nt
            gs=reshape(gpl(:,to,id,nn,:),Nb,ncs);
            if sum(gs,'all')<1e-15, continue; end
            tcs=reshape(tenure_choice(:,to,id,j,nn,:),Nb,ncs);
            for tn=1:nt
                mk=(tcs==tn); if ~any(mk,'all'), continue; end
                mt=gs.*mk; if sum(mt,'all')<1e-15, continue; end
                rd=zeros(Nb,ncs);
                for cs=1:ncs, rd(:,cs)=tmx{id,to,tn,nn,cs}*mt(:,cs); end
                gpt(:,tn,id,nn,:)=gpt(:,tn,id,nn,:)+reshape(rd,Nb,1,1,1,ncs);
            end
        end, end, end

        % Savings advance (linear redistribution for continuous b')
        gps=zeros(Nb,nt,I,np,ncs);
        for i=1:I, for ten=1:nt
            gf=reshape(gpt(:,ten,i,:,:),Nb,nc);
            bpv=reshape(bp_pol(:,ten,i,j,:,:),Nb,nc);
            % Clamp to grid bounds
            bpv = max(min(bpv, b_grid(end)), b_grid(1));
            % Find bracketing grid points and weights
            if P.has_lookup_mex
                idx = lookup(b_grid, bpv, 3);
            else
                idx = max(min(discretize(bpv, [b_grid; inf]), Nb-1), 1);
            end
            w = (bpv - b_grid(idx)) ./ (b_grid(idx+1) - b_grid(idx));
            w = max(min(w, 1), 0);
            % Distribute mass
            g_new = zeros(Nb, nc);
            for col=1:nc
                g_new(:,col) = accumarray(idx(:,col), (1-w(:,col)).*gf(:,col), [Nb,1]) ...
                             + accumarray(idx(:,col)+1, w(:,col).*gf(:,col), [Nb,1]);
            end
            gps(:,ten,i,:,:)=reshape(g_new,Nb,1,1,np,ncs);
        end, end

        % Child aging
        for nn=1:np, for cs=1:ncs
            gp=gps(:,:,:,nn,cs);
            if ust
                Pn=Pia(:,:,nn);
                if (cs==K+1)&&(nn>=2)
                    if nn==2,pm=Pn(cs,csm1);else,pm=Pn(cs,csm2);end
                    if pm>0, nk=nn-1;
                        for im=1:I, fi=nk*pm*sum(gp(:,:,im),'all');
                            eml(im)=eml(im)+fi; emt=emt+fi; end
                    end
                end
                for csn=1:ncs, wt=Pn(cs,csn);
                    if wt>0, g(:,:,:,j+1,nn,csn)=g(:,:,:,j+1,nn,csn)+wt*gp; end
                end
            else
                if cs==1,csn=1; elseif cs>=csm1,csn=cs; elseif cs<K+1,csn=cs+1;
                else, if nn==1,csn=1;elseif nn==2,csn=csm1;else,csn=csm2;end, end
                if (cs==K+1)&&(csn>=csm1)&&(nn>=2), nk=nn-1;
                    for im=1:I, fi=nk*sum(gp(:,:,im),'all');
                        eml(im)=eml(im)+fi; emt=emt+fi; end
                end
                g(:,:,:,j+1,nn,csn)=g(:,:,:,j+1,nn,csn)+gp;
            end
        end, end
    end

    tm=sum(g,'all');
    if tm>1e-12
        sc=P.N_target/tm; g=g*sc; tb=tb*sc; blc=blc*sc; eml=eml*sc; emt=emt*sc;
    end
    if nargin >= 11 && fast_stats
        stats=compute_eq_stats(g,P,b_grid,p_hat,hR_pol);
    else
        stats=compute_statistics(g,fert_probs,loc_probs,P,b_grid,p_hat,hR_pol);
    end
    stats.total_births_kfe=tb; stats.births_by_loc=blc;
    stats.entry_rate=sum(g(:,:,:,1,:,:),'all'); stats.total_mass=sum(g,'all');
    stats.entrants_mature_by_loc=eml; stats.entrants_mature_total=emt;
    stats.mature_entry_shares=eml/max(emt,1e-12);
    if birth_es3_mass > 1e-12
        stats.housing_increment_0to1_eventstudy_t3 = ...
            birth_es3_post_sum / birth_es3_mass - birth_es3_pre_sum / birth_es3_mass;
    else
        stats.housing_increment_0to1_eventstudy_t3 = 0;
    end
    if addchild_es3_one_mass > 1e-12 && addchild_es3_two_plus_mass > 1e-12
        stats.housing_increment_1to2_proxy_t3 = ...
            addchild_es3_two_plus_sum / addchild_es3_two_plus_mass - ...
            addchild_es3_one_sum / addchild_es3_one_mass;
    else
        stats.housing_increment_1to2_proxy_t3 = 0;
    end
    if onechild_es3_mass > 1e-12
        stats.housing_increment_0to1_onechild_eventstudy_t3 = ...
            onechild_es3_post_sum / onechild_es3_mass - onechild_es3_pre_sum / onechild_es3_mass;
    else
        stats.housing_increment_0to1_onechild_eventstudy_t3 = 0;
    end
    if twoplus_es3_mass > 1e-12
        stats.housing_increment_0to2plus_eventstudy_t3 = ...
            twoplus_es3_post_sum / twoplus_es3_mass - twoplus_es3_pre_sum / twoplus_es3_mass;
    else
        stats.housing_increment_0to2plus_eventstudy_t3 = 0;
    end
    stats.housing_event_horizon = event_horizon;
end


%% ========================================================================
%%   FAST INTERPOLATION (uses CompEcon lookup MEX when available)
%% ========================================================================
function [idx, wt] = fast_interp_idx(b_grid, bq, has_mex)
    Nb = length(b_grid);
    bq = max(min(bq, b_grid(end)), b_grid(1));
    if has_mex
        idx = lookup(b_grid, bq, 3);
    else
        idx = max(min(discretize(bq, [b_grid; inf]), Nb-1), 1);
    end
    wt = max(min((bq - b_grid(idx)) ./ (b_grid(idx+1) - b_grid(idx)), 1), 0);
end

function Vi = interp1_fast(bg, V, bq, has_mex)
    Nb=length(bg); sz=size(V);
    [idx, wt] = fast_interp_idx(bg, bq, has_mex);
    nt=prod(sz(2:end)); Vf=reshape(V,Nb,nt);
    Vi=reshape((1-wt).*Vf(idx,:)+wt.*Vf(idx+1,:), sz);
end


%% ========================================================================
%%   FAST STATISTICS (equilibrium loop only — 5 fields)
%% ========================================================================
function stats = compute_eq_stats(g,P,bg,ph,hR)
    J=P.J;I=P.I;Nb=length(bg);nt=1+P.n_house;np=P.n_parity;ncs=P.n_child_states;
    tm=sum(g,'all');
    stats.own_rate=sum(g(:,2:end,:,:,:,:),'all')/max(tm,1e-12);
    stats.pop_share=zeros(I,1); stats.housing_demand=zeros(I,1);
    for i=1:I
        pi=sum(g(:,:,i,:,:,:),'all'); stats.pop_share(i)=pi/max(tm,1e-12);
        Hd=0;
        for j=1:J,for nn=1:np,for cs=1:ncs
            Hd=Hd+sum(g(:,1,i,j,nn,cs).*hR(:,1,i,j,nn,cs));
            for ten=2:nt,Hd=Hd+sum(g(:,ten,i,j,nn,cs))*P.H_own(ten-1);end
        end,end,end
        stats.housing_demand(i)=Hd/max(P.N_target,1e-12);
    end
    mp=sum(g(:,:,:,P.A_f_end+1:end,:,:),'all');
    stats.parity_dist=zeros(np,1);
    for nn=1:np,stats.parity_dist(nn)=sum(g(:,:,:,P.A_f_end+1:end,nn,:),'all')/max(mp,1e-12);end
    stats.mean_parity=sum((0:np-1)'.*stats.parity_dist);
end


%% ========================================================================
%%   STATISTICS (full — called once at convergence)
%% ========================================================================
function stats = compute_statistics(g,fp,lp,P,bg,ph,hR)
    J=P.J;I=P.I;Nb=length(bg);nt=1+P.n_house;np=P.n_parity;ncs=P.n_child_states;
    tm=sum(g,'all');
    stats.own_rate=sum(g(:,2:end,:,:,:,:),'all')/max(tm,1e-12);
    stats.pop_share=zeros(I,1); stats.own_by_loc=zeros(I,1); stats.housing_demand=zeros(I,1);
    for i=1:I
        pi=sum(g(:,:,i,:,:,:),'all'); stats.pop_share(i)=pi/max(tm,1e-12);
        stats.own_by_loc(i)=sum(g(:,2:end,i,:,:,:),'all')/max(pi,1e-12);
        Hd=0;
        for j=1:J,for nn=1:np,for cs=1:ncs
            Hd=Hd+sum(g(:,1,i,j,nn,cs).*hR(:,1,i,j,nn,cs));
            for ten=2:nt,Hd=Hd+sum(g(:,ten,i,j,nn,cs))*P.H_own(ten-1);end
        end,end,end
        stats.housing_demand(i)=Hd/max(P.N_target,1e-12);
    end
    stats.worker_mass_by_loc=zeros(I,1);stats.retiree_mass_by_loc=zeros(I,1);
    for i=1:I
        stats.worker_mass_by_loc(i)=sum(g(:,:,i,1:P.J_R,:,:),'all');
        stats.retiree_mass_by_loc(i)=sum(g(:,:,i,P.J_R+1:J,:,:),'all');
    end
    stats.worker_mass_total=sum(stats.worker_mass_by_loc);
    stats.retiree_mass_total=sum(stats.retiree_mass_by_loc);
    stats.parity_dist=zeros(np,1);
    mp=sum(g(:,:,:,P.A_f_end+1:end,:,:),'all');
    for nn=1:np,stats.parity_dist(nn)=sum(g(:,:,:,P.A_f_end+1:end,nn,:),'all')/max(mp,1e-12);end
    stats.mean_parity=sum((0:np-1)'.*stats.parity_dist);
    stats.own_by_parity=zeros(np,1);
    for nn=1:np,mn=sum(g(:,:,:,:,nn,:),'all');
        stats.own_by_parity(nn)=sum(g(:,2:end,:,:,nn,:),'all')/max(mn,1e-12);end
    stats.mean_parity_by_loc=zeros(I,1);stats.frac_childless_by_loc=zeros(I,1);
    for i=1:I,mip=sum(g(:,:,i,P.A_f_end+1:end,:,:),'all');
        if mip>1e-12,mn=0;
            for nn=1:np,mn=mn+(nn-1)*sum(g(:,:,i,P.A_f_end+1:end,nn,:),'all')/mip;end
            stats.mean_parity_by_loc(i)=mn;
            stats.frac_childless_by_loc(i)=sum(g(:,:,i,P.A_f_end+1:end,1,:),'all')/mip;
        end,end
    stats.own_by_age=zeros(J,1);
    for jj=1:J,gj=g(:,:,:,jj,:,:);mj=sum(gj,'all');
        if mj>1e-12,stats.own_by_age(jj)=sum(gj(:,2:end,:,:,:,:),'all')/mj;end,end
    stats.child_state_dist=zeros(J,ncs);
    for jj=1:J,mj=sum(g(:,:,:,jj,:,:),'all');if mj>1e-12
        for cs=1:ncs,stats.child_state_dist(jj,cs)=sum(g(:,:,:,jj,:,cs),'all')/mj;end,end,end
    stats.fert_by_age=zeros(J,1);
    for j=P.A_f_start:P.A_f_end,mj=sum(g(:,:,:,j,1,1),'all');if mj>1e-12
        En=0;for i=1:I,for ten=1:nt,gs=g(:,ten,i,j,1,1);nz=find(gs>1e-15);
            if isempty(nz),continue;end;pr=reshape(fp(:,ten,i,j,:),Nb,np);
            En=En+sum(gs(nz).*(pr(nz,:)*(0:np-1)'));end,end
        stats.fert_by_age(j)=En/mj;end,end
    a22s=max(1,round(22-P.age_start+1));a25s=max(1,round(25-P.age_start+1));a45e=min(J,round(45-P.age_start+1));
    a30s=max(1,round(30-P.age_start+1));a55e=min(J,round(55-P.age_start+1));
    a65s=max(1,round(65-P.age_start+1));a75e=min(J,round(75-P.age_start+1));
    asw=max(1,round(45-P.age_start+1));aew=min(J,round(55-P.age_start+1));
    newparent_cs=2:min(P.n_child_stages+1,3);
    ti=0;tmw=0;
    for i=1:I,ti=ti+P.income(i,1)*stats.worker_mass_by_loc(i);tmw=tmw+stats.worker_mass_by_loc(i);end
    mi=ti/max(tmw,1e-12);
    tw=0;tm4=0;
    for jj=asw:min(aew,J),for i=1:I,for ten=1:nt
        heq=0;if ten>1,heq=(1-P.psi)*ph(i)*P.H_own(ten-1);end
        for nn=1:np,for cs=1:ncs,gs=g(:,ten,i,jj,nn,cs);
            tw=tw+sum(gs.*(bg+heq));tm4=tm4+sum(gs);end,end
    end,end,end
    stats.mean_wealth_4555=tw/max(tm4,1e-12);stats.mean_income=mi;
    stats.wealth_to_income=stats.mean_wealth_4555/max(mi,1e-12);
    stats.entry_mass_by_loc=zeros(I,1);
    for i=1:I,stats.entry_mass_by_loc(i)=sum(g(:,:,i,1,:,:),'all');end
    stats=append_pension_budget_stats(stats,g,P);
    tmv=0;tas=0;
    for j=1:J-1,for io=1:I,for ten=1:nt,for nn=1:np,for cs=1:ncs
        gs=g(:,ten,io,j,nn,cs);mh=sum(gs);if mh<1e-15,continue;end
        sp=lp(:,ten,io,io,j,nn,cs);tmv=tmv+sum(gs.*(1-sp));tas=tas+mh;
    end,end,end,end,end
    stats.migration_rate=tmv/max(tas,1e-12);
    tmv2=0;tas2=0;
    for j=a22s:min(a45e,J),for io=1:I,for ten=1:nt,for nn=1:np,for cs=1:ncs
        gs=g(:,ten,io,j,nn,cs);mh=sum(gs);if mh<1e-15,continue;end
        sp=lp(:,ten,io,io,j,nn,cs);tmv2=tmv2+sum(gs.*(1-sp));tas2=tas2+mh;
    end,end,end,end,end
    stats.migration_rate_2245=tmv2/max(tas2,1e-12);
    tl=0;tml=0;
    for jj=asw:min(aew,J),for i=1:I,for ten=1:nt,for nn=1:np,for cs=1:ncs
        gs=g(:,ten,i,jj,nn,cs);tl=tl+sum(gs.*bg);tml=tml+sum(gs);
    end,end,end,end,end
    stats.liquid_wealth_4555=tl/max(tml,1e-12);
    stats.liquid_wealth_to_income=stats.liquid_wealth_4555/max(mi,1e-12);
    ays=max(1,round(25-P.age_start+1));aye=min(J,round(35-P.age_start+1));
    yo=0;yt=0;
    for jj=ays:aye,gj=g(:,:,:,jj,:,:);yt=yt+sum(gj,'all');
        yo=yo+sum(gj(:,2:end,:,:,:,:),'all');end
    stats.young_own_rate=yo/max(yt,1e-12);
    ylw=0;yinc=0;ycm=0;
    for jj=ays:aye,for i=1:I,gs=g(:,1,i,jj,1,1);mh=sum(gs);if mh<1e-15,continue;end
        ylw=ylw+sum(gs.*bg);yinc=yinc+P.income(i,jj)*mh;ycm=ycm+mh;end,end
    stats.young_liquid_wealth=ylw/max(ycm,1e-12);
    stats.young_childless_renter_income=yinc/max(ycm,1e-12);
    stats.young_liquid_wealth_to_income=ylw/max(yinc,1e-12);
    tba=0;tfb=0;
    for j=P.A_f_start:P.A_f_end,ra=P.age_start+(j-1)*P.da;
        for i=1:I,for ten=1:nt,gs=g(:,ten,i,j,1,1);nz=find(gs>1e-15);
            if isempty(nz),continue;end;pb=1-fp(nz,ten,i,j,1);
            wb=sum(gs(nz).*pb);tba=tba+ra*wb;tfb=tfb+wb;end,end,end
    stats.mean_age_first_birth=tba/max(tfb,1e-12);
    mg1=sum(g(:,:,:,P.A_f_end+1:end,2:end,:),'all');
    stats.parity_progression_1to2=sum(g(:,:,:,P.A_f_end+1:end,3:end,:),'all')/max(mg1,1e-12);
    ic=2;mcy=0;mct=0;
    for jj=ays:aye,mcy=mcy+sum(g(:,:,ic,jj,1,:),'all');mct=mct+sum(g(:,:,:,jj,1,:),'all');end
    stats.center_share_childless_young=mcy/max(mct,1e-12);
    stats.center_share_parents=sum(g(:,:,ic,:,2:end,:),'all')/max(sum(g(:,:,:,:,2:end,:),'all'),1e-12);
    mpo=sum(g(:,2:end,:,:,2:end,:),'all');mpa=sum(g(:,:,:,:,2:end,:),'all');
    stats.own_rate_parents=mpo/max(mpa,1e-12);
    mco=sum(g(:,2:end,:,:,1,:),'all');mca=sum(g(:,:,:,:,1,:),'all');
    stats.own_rate_childless=mco/max(mca,1e-12);
    stats.own_family_gap=stats.own_rate_parents-stats.own_rate_childless;
    a25s=max(1,round(25-P.age_start+1));
    a34e=min(J,round(34-P.age_start+1));
    a35s=max(1,round(35-P.age_start+1));
    a44e=min(J,round(44-P.age_start+1));
    prime_mass=sum(g(:,:,:,a30s:min(a55e,J),:,:),'all');
    prime_owner=sum(g(:,2:end,:,a30s:min(a55e,J),:,:),'all');
    stats.own_rate_3055=prime_owner/max(prime_mass,1e-12);
    early_mass=sum(g(:,:,:,a25s:min(a34e,J),:,:),'all');
    early_owner=sum(g(:,2:end,:,a25s:min(a34e,J),:,:),'all');
    stats.own_rate_2534=early_owner/max(early_mass,1e-12);
    mid_mass=sum(g(:,:,:,a35s:min(a44e,J),:,:),'all');
    mid_owner=sum(g(:,2:end,:,a35s:min(a44e,J),:,:),'all');
    stats.own_rate_3544=mid_owner/max(mid_mass,1e-12);
    prime_mass_p=sum(g(:,:,1,a30s:min(a55e,J),:,:),'all');
    prime_mass_c=sum(g(:,:,2,a30s:min(a55e,J),:,:),'all');
    prime_owner_p=sum(g(:,2:end,1,a30s:min(a55e,J),:,:),'all');
    prime_owner_c=sum(g(:,2:end,2,a30s:min(a55e,J),:,:),'all');
    stats.own_gradient_3055=prime_owner_p/max(prime_mass_p,1e-12)-prime_owner_c/max(prime_mass_c,1e-12);
    nonparent_mass_2245=sum(g(:,:,:,a22s:min(a45e,J),1,1),'all');
    stats.center_share_nonparents_2245=sum(g(:,:,2,a22s:min(a45e,J),1,1),'all')/max(nonparent_mass_2245,1e-12);
    nonparent_mass_3055=sum(g(:,:,:,a30s:min(a55e,J),1,1),'all');
    nonparent_owner_3055=sum(g(:,2:end,:,a30s:min(a55e,J),1,1),'all');
    stats.own_rate_nonparents_3055=nonparent_owner_3055/max(nonparent_mass_3055,1e-12);
    if isempty(newparent_cs)
        stats.center_share_newparents_2245=0;
        stats.own_rate_newparents_3055=0;
    else
        newparent_mass_2245=sum(g(:,:,:,a22s:min(a45e,J),2:end,newparent_cs),'all');
        stats.center_share_newparents_2245=sum(g(:,:,2,a22s:min(a45e,J),2:end,newparent_cs),'all')/max(newparent_mass_2245,1e-12);
        newparent_mass_3055=sum(g(:,:,:,a30s:min(a55e,J),2:end,newparent_cs),'all');
        newparent_owner_3055=sum(g(:,2:end,:,a30s:min(a55e,J),2:end,newparent_cs),'all');
        stats.own_rate_newparents_3055=newparent_owner_3055/max(newparent_mass_3055,1e-12);
    end
    stats.own_gap_newparent_nonparent_3055=stats.own_rate_newparents_3055-stats.own_rate_nonparents_3055;
    old_mass=sum(g(:,:,:,a65s:min(a75e,J),:,:),'all');
    old_owner=sum(g(:,2:end,:,a65s:min(a75e,J),:,:),'all');
    stats.old_age_own_rate_6575=old_owner/max(old_mass,1e-12);
    old_parent_mass=sum(g(:,:,:,a65s:min(a75e,J),2:end,:),'all');
    old_parent_owner=sum(g(:,2:end,:,a65s:min(a75e,J),2:end,:),'all');
    old_childless_mass=sum(g(:,:,:,a65s:min(a75e,J),1,1),'all');
    old_childless_owner=sum(g(:,2:end,:,a65s:min(a75e,J),1,1),'all');
    stats.old_age_own_rate_parents_6575=old_parent_owner/max(old_parent_mass,1e-12);
    stats.old_age_own_rate_childless_6575=old_childless_owner/max(old_childless_mass,1e-12);
    stats.old_age_parent_childless_gap_6575=stats.old_age_own_rate_parents_6575-stats.old_age_own_rate_childless_6575;
    neg_owner_mass_2545=sum(g(bg<0,2:end,:,a25s:min(a45e,J),:,:),'all');
    owner_mass_2545=sum(g(:,2:end,:,a25s:min(a45e,J),:,:),'all');
    stats.owner_neg_liquid_share_2545=neg_owner_mass_2545/max(owner_mass_2545,1e-12);
    neg_owner_mass_2534=sum(g(bg<0,2:end,:,a25s:min(a34e,J),:,:),'all');
    owner_mass_2534=sum(g(:,2:end,:,a25s:min(a34e,J),:,:),'all');
    stats.owner_neg_liquid_share_2534=neg_owner_mass_2534/max(owner_mass_2534,1e-12);
    % Prime-age childless room medians are the baseline housing-level moments
    % used in the Stage A repair of the live calibration surface.
    dep_last=P.n_child_stages+1;
    renter_room_vals={};renter_room_wts={};owner_room_vals={};owner_room_wts={};
    for j=a25s:min(a45e,J),for i=1:I,for nn=1:np,for cs=1:ncs
        if current_child_bin_dt(nn,cs,dep_last)~=2,continue;end
        gr=g(:,1,i,j,nn,cs);hr=hR(:,1,i,j,nn,cs);kr=(gr>0)&isfinite(hr)&(hr>0);
        if any(kr)
            renter_room_vals{end+1,1}=hr(kr); %#ok<AGROW>
            renter_room_wts{end+1,1}=gr(kr); %#ok<AGROW>
        end
        for ten=2:nt
            go=g(:,ten,i,j,nn,cs);ko=go>0;
            if ~any(ko),continue;end
            owner_room_vals{end+1,1}=P.H_own(ten-1)*ones(nnz(ko),1); %#ok<AGROW>
            owner_room_wts{end+1,1}=go(ko); %#ok<AGROW>
        end
    end,end,end,end
    stats.prime_childless_renter_median_rooms=weighted_median_from_cells_dt(renter_room_vals,renter_room_wts);
    stats.prime_childless_owner_median_rooms=weighted_median_from_cells_dt(owner_room_vals,owner_room_wts);
    stats.mean_housing_by_parity=zeros(np,1);
    for nn=1:np,th=0;mn=0;for i=1:I,for ten=1:nt,for j=1:J,for cs=1:ncs
        gs=g(:,ten,i,j,nn,cs);mh=sum(gs);if mh<1e-15,continue;end
        if ten==1,th=th+sum(gs.*hR(:,ten,i,j,nn,cs));else,th=th+mh*P.H_own(ten-1);end
        mn=mn+mh;end,end,end,end
        stats.mean_housing_by_parity(nn)=th/max(mn,1e-12);end
    if np>=3,stats.housing_increment_1to2=stats.mean_housing_by_parity(3)-stats.mean_housing_by_parity(2);
    else,stats.housing_increment_1to2=0;end
end

function stats=append_pension_budget_stats(stats,g,P)
    if isfield(P,'income_age_profile') && ~isempty(P.income_age_profile)
        income_profile=P.income_age_profile(:);
    else
        income_profile=get_income_age_profile(P);
    end
    payroll_tax_revenue=0;
    for i=1:P.I
        for j=1:P.J_R
            mass_ij=sum(g(:,:,i,j,:,:),'all');
            payroll_tax_revenue=payroll_tax_revenue+P.tau_pay*P.w_hat(i)*income_profile(j)*mass_ij;
        end
    end
    pension_outlays=P.pension*stats.retiree_mass_total;
    stats.payroll_tax_revenue=payroll_tax_revenue;
    stats.pension_outlays=pension_outlays;
    stats.pension_budget_residual=payroll_tax_revenue-pension_outlays;
    stats.implied_balanced_pension=payroll_tax_revenue/max(stats.retiree_mass_total,1e-12);
    stats.pension=P.pension;
end

function g_out = advance_cohort_horizon(g_in, start_age, horizon, loc_probs, tenure_choice, bp_pol, P, b_grid, SD, lmm, tmx, ust, Pia)
    g_out = g_in;
    for step = 1:horizon
        age_idx = start_age + step - 1;
        if age_idx >= P.J
            break;
        end
        g_out = advance_cohort_one_period(g_out, age_idx, loc_probs, tenure_choice, bp_pol, P, b_grid, SD, lmm, tmx, ust, Pia);
    end
end

function child_code = current_child_bin_dt(nn,cs,dep_last)
    if cs==1 || cs>dep_last
        child_code=2;
        return;
    end
    current_n=max(nn-1,0);
    if current_n<=0
        child_code=2;
    elseif current_n==1
        child_code=3;
    else
        child_code=4;
    end
end

function q = weighted_median_from_cells_dt(value_cells,weight_cells)
    if isempty(value_cells)
        q=NaN;
        return;
    end
    values=vertcat(value_cells{:});
    weights=vertcat(weight_cells{:});
    q=weighted_quantile_dt(values,weights,0.5);
end

function q = weighted_quantile_dt(values,weights,probs)
    values=values(:);weights=weights(:);
    keep=isfinite(values)&isfinite(weights)&(weights>0);
    values=values(keep);weights=weights(keep);
    if isempty(values)
        q=NaN(size(probs));
        return;
    end
    [values,ord]=sort(values);
    weights=weights(ord);
    cw=cumsum(weights)/sum(weights);
    q=zeros(size(probs));
    for k=1:numel(probs)
        idx=find(cw>=probs(k),1,'first');
        if isempty(idx),idx=numel(values);end
        q(k)=values(idx);
    end
end

function g_next = advance_cohort_one_period(gj, j, loc_probs, tenure_choice, bp_pol, P, b_grid, SD, lmm, tmx, ust, Pia)
    Nb = length(b_grid); nt = 1 + P.n_house; I = P.I; np = P.n_parity; ncs = P.n_child_states; nc = SD.nc;
    bmin = b_grid(1); bmax = b_grid(end);
    K = P.n_child_stages; csm1 = K+2; csm2 = K+3;

    gpl=zeros(Nb,nt,I,np,ncs);
    for io=1:I, for to=1:nt
        go=reshape(gj(:,to,io,:,:),Nb,np*ncs);
        if sum(go,'all')<1e-15, continue; end
        po=reshape(loc_probs(:,to,io,:,j,:,:),Nb,I,np*ncs);
        sp=reshape(po(:,io,:),Nb,np*ncs);
        gpl(:,to,io,:,:)=gpl(:,to,io,:,:)+reshape(go.*sp,Nb,1,1,np,ncs);
        R=lmm{io,to};
        for id=1:I
            if id==io, continue; end
            mp=reshape(po(:,id,:),Nb,np*ncs);
            gpl(:,1,id,:,:)=gpl(:,1,id,:,:)+reshape(R*(go.*mp),Nb,1,1,np,ncs);
        end
    end, end

    gpt=zeros(Nb,nt,I,np,ncs);
    for nn=1:np, for id=1:I, for to=1:nt
        gs=reshape(gpl(:,to,id,nn,:),Nb,ncs);
        if sum(gs,'all')<1e-15, continue; end
        tcs=reshape(tenure_choice(:,to,id,j,nn,:),Nb,ncs);
        for tn=1:nt
            mk=(tcs==tn); if ~any(mk,'all'), continue; end
            mt=gs.*mk; if sum(mt,'all')<1e-15, continue; end
            rd=zeros(Nb,ncs);
            for cs=1:ncs, rd(:,cs)=tmx{id,to,tn,nn,cs}*mt(:,cs); end
            gpt(:,tn,id,nn,:)=gpt(:,tn,id,nn,:)+reshape(rd,Nb,1,1,1,ncs);
        end
    end, end, end

    gps=zeros(Nb,nt,I,np,ncs);
    for i=1:I, for ten=1:nt
        gf=reshape(gpt(:,ten,i,:,:),Nb,nc);
        bpv=reshape(bp_pol(:,ten,i,j,:,:),Nb,nc);
        bpc=max(min(bpv,bmax),bmin);
        if P.has_lookup_mex
            idx = lookup(b_grid, bpc, 3);
        else
            idx = max(min(discretize(bpc, [b_grid; inf]), Nb-1), 1);
        end
        w = (bpc - b_grid(idx)) ./ (b_grid(idx+1) - b_grid(idx));
        w = max(min(w, 1), 0);
        g_new = zeros(Nb, nc);
        for col=1:nc
            if sum(gf(:,col)) < 1e-15, continue; end
            g_new(:,col) = accumarray(idx(:,col), (1-w(:,col)).*gf(:,col), [Nb,1]) ...
                         + accumarray(idx(:,col)+1, w(:,col).*gf(:,col), [Nb,1]);
        end
        gps(:,ten,i,:,:)=reshape(g_new,Nb,1,1,np,ncs);
    end, end

    g_next=zeros(Nb,nt,I,np,ncs);
    if ust
        for nn=1:np
            Pi=Pia(:,:,nn);
            for i=1:I, for ten=1:nt
                gin=reshape(gps(:,ten,i,nn,:),Nb,ncs);
                gout=gin*Pi;
                g_next(:,ten,i,nn,:)=reshape(gout,Nb,1,1,1,ncs);
            end, end
        end
    else
        for nn=1:np, for cs=1:ncs
            if cs==1,csn=1;elseif cs>=csm1,csn=cs;elseif cs<K+1,csn=cs+1;
            else,if nn==1,csn=1;elseif nn==2,csn=csm1;else,csn=csm2;end,end
            g_next(:,:,:,nn,csn)=g_next(:,:,:,nn,csn)+gps(:,:,:,nn,cs);
        end, end
    end
end

function mh = mean_housing_childless_weighted(weight_dist, j, hR_pol, P)
    [Nb,nt,I] = size(weight_dist);
    th = 0; mn = 0;
    for i=1:I, for ten=1:nt
        gs = weight_dist(:,ten,i);
        mh_state = sum(gs);
        if mh_state < 1e-15, continue; end
        if ten==1
            th = th + sum(gs .* hR_pol(:,ten,i,j,1,1));
        else
            th = th + mh_state * P.H_own(ten-1);
        end
        mn = mn + mh_state;
    end, end
    mh = th / max(mn, 1e-12);
end

function mh = mean_housing_distribution(g_dist, j, hR_pol, P)
    [Nb,nt,I,np,ncs] = size(g_dist);
    th = 0; mn = 0;
    for nn=1:np, for cs=1:ncs, for i=1:I, for ten=1:nt
        gs = g_dist(:,ten,i,nn,cs);
        mh_state = sum(gs);
        if mh_state < 1e-15, continue; end
        if ten==1
            th = th + sum(gs .* hR_pol(:,ten,i,j,nn,cs));
        else
            th = th + mh_state * P.H_own(ten-1);
        end
        mn = mn + mh_state;
    end, end, end, end
    mh = th / max(mn, 1e-12);
end


%% ========================================================================
%%   HELPERS
%% ========================================================================
function Vc=apply_child_aging(Vn,P,Nb,nt,I,np,ncs)
    Vc=zeros(Nb,nt,I,np,ncs); K=P.n_child_stages;
    if P.use_stochastic_aging&&isfield(P,'Pi_child')
        Pa=P.Pi_child;
        for nn=1:np, Pi=Pa(:,:,nn);
            Vnn=reshape(Vn(:,:,:,nn,:),[],ncs);
            Vc(:,:,:,nn,:)=reshape(Vnn*Pi',Nb,nt,I,1,ncs);
        end
    else
        csm1=K+2;csm2=K+3;
        for nn=1:np,for cs=1:ncs
            if cs==1,csn=1;elseif cs>=csm1,csn=cs;elseif cs<K+1,csn=cs+1;
            else,if nn==1,csn=1;elseif nn==2,csn=csm1;else,csn=csm2;end,end
            Vc(:,:,:,nn,cs)=Vn(:,:,:,nn,csn);
        end,end
    end
end

function b_grid = make_grid(P)
    Nb=P.Nb; N1=round(Nb*0.15);N2=round(Nb*0.45);N3=round(Nb*0.15);N4=Nb-N1-N2-N3;
    s1=linspace(P.b_min,-3,N1+1)';s1=s1(1:end-1);
    s2=linspace(-3,6,N2+1)';s2=s2(1:end-1);
    s3=linspace(6,20,N3+1)';s3=s3(1:end-1);
    u4=linspace(0,1,N4+1)';u4=u4(2:end);s4=20+(P.b_max-20)*u4.^P.b_grid_power;
    b_grid=[s1;s2;s3;s4];
    [~,iz]=min(abs(b_grid));b_grid(iz)=0;
    [~,ie]=min(abs(b_grid-P.b_entry_fixed));b_grid(ie)=P.b_entry_fixed;
end

function sol=pack_solution(V,c,h,bp,tc,lp,fp,fv,g,st,w,p)
    sol=struct('V',V,'c_pol',c,'hR_pol',h,'bp_pol',bp,'tenure_choice',tc,...
        'loc_probs',lp,'fert_probs',fp,'fert_value',fv,'g',g,'w_hat',w,'p_eq',p);
    fn=fieldnames(st); for k=1:length(fn),sol.(fn{k})=st.(fn{k});end
end

function P=apply_overrides(P,Po)
    fn=fieldnames(Po); for k=1:length(fn),P.(fn{k})=Po.(fn{k});end
    if isfield(Po,'eps_loc'),P.kappa_loc=Po.eps_loc;P.eps_loc=P.kappa_loc;end
    if isfield(Po,'eps_fert'),P.kappa_fert=Po.eps_fert;P.eps_fert=P.kappa_fert;end
    if isfield(Po,'kappa_loc'),P.eps_loc=P.kappa_loc;end
    if isfield(Po,'kappa_fert'),P.eps_fert=P.kappa_fert;end
    if isfield(Po,'H_bar'),P.H0=Po.H_bar;end
    if isfield(Po,'n_parity')
        P.phi=0.80*ones(P.n_parity,1);
        % Rebuild child transition matrix for new n_parity
        P.Pi_child=make_child_transition_matrix_with_matured(P.stage_durations,P.n_parity);
    end
    pension_refresh_fields = {'J','J_R','tau_pay','w_hat','entry_shares', ...
        'income_age_breaks','income_age_values','pension_mode','pension'};
    if any(isfield(Po, pension_refresh_fields)),P=set_income_given_w_and_pension(P);end
    if isfield(Po,'J'),P.entry_by_loc=P.N_0/sum(P.N_0)/P.J;end
    % Rebuild housing ladder if overridden
    if isfield(Po,'H_own')
        P.n_house=length(P.H_own);
        P.h_own_min=P.H_own(1); P.h_own_max=P.H_own(end);
    end
end

function P=finalize_location_choice_spec(P)
    if ~isfield(P,'location_choice_form') || isempty(P.location_choice_form)
        P.location_choice_form='legacy_multiplicative';
    else
        P.location_choice_form=lower(char(string(P.location_choice_form)));
    end

    switch P.location_choice_form
        case {'additive_due','additive','due'}
            P.location_choice_form='additive_due';
            P.location_choice_legacy_converted=false;
        case {'legacy_multiplicative','multiplicative'}
            if any(P.E_loc <= 0)
                error('Legacy multiplicative location amenities must be strictly positive.');
            end
            if ~isfield(P,'mu_stay') || isempty(P.mu_stay), P.mu_stay = 1.0; end
            if P.mu_stay <= 0 || P.mu_move <= 0
                error('Legacy multiplicative moving wedges must be strictly positive.');
            end
            P.E_loc = P.kappa_loc * log(P.E_loc(:));
            P.mu_stay = -P.kappa_loc * log(P.mu_stay);
            P.mu_move = -P.kappa_loc * log(P.mu_move);
            if isfield(P,'mu_move_parent') && ~isempty(P.mu_move_parent) && P.mu_move_parent > 0
                P.mu_move_parent = -P.kappa_loc * log(P.mu_move_parent);
            end
            P.location_choice_form='additive_due';
            P.location_choice_legacy_converted=true;
        otherwise
            error('Unknown location_choice_form: %s', P.location_choice_form);
    end

    if ~isfield(P,'mu_stay') || isempty(P.mu_stay), P.mu_stay = 0.0; end
end

function n=get_completed_fertility(nn,cs,P)
    K=P.n_child_stages;
    if cs==1,n=0;elseif cs==K+2,n=1;elseif cs==K+3,n=2;else,n=nn-1;end
end
function v=bequest_utility_vec(b,nk,P)
    b=max(b,0);
    scale = P.theta0 * max(1 + P.theta_n*nk, 0);
    if abs(P.sigma-1)<1e-6,v=scale.*log(P.theta1+b);
    else,v=scale.*(P.theta1+b).^(1-P.sigma)/(1-P.sigma);end
end
function R=make_lin_redist(bg,bq)
    Nb=length(bg);bq=max(min(bq(:),bg(end)),bg(1));
    idx=max(min(discretize(bq,[bg;inf]),Nb-1),1);
    w=max(min((bq-bg(idx))./(bg(idx+1)-bg(idx)),1),0);
    s=(1:Nb)';r=[idx;idx+1];c=[s;s];v=[1-w;w];
    k=v>0;R=sparse(r(k),c(k),v(k),Nb,Nb);
end
function R=make_redist(bg,bf,use_lin)
    Nb=length(bg);
    if use_lin, R=make_lin_redist(bg,bf);
    else, [~,di]=min(abs(bg-bf'),[],1);
        R=sparse(di(:),(1:Nb)',ones(Nb,1),Nb,Nb);end
end
function tf=has_birth_dp_grant(P,nn,cs,to,tn)
    tf=false;if ~(isfield(P,'birth_dp_grant')&&P.birth_dp_grant),return;end
    if to~=1||tn<=1,return;end;tf=(nn>=2)&&(cs==2);
end
%% ========================================================================
%%   GOLDEN SECTION OBJECTIVE FUNCTIONS (no anonymous handles, no closures)
%%   Pre-clamped grid bounds and col_offset passed directly.
%% ========================================================================
function f = gs_eval_rent(bp, Rv, ri, alpha, oms, cb, hb, psi, hRmax, cap, hRmax_rep, beta, Vcont, bgrid, Kr, Nb, bg_end, bg_1, col_offset, has_mex)
    surplus = Rv - cb - ri*hb - bp;
    ss = max(surplus, 1e-10);
    u = Kr * ss.^oms / oms + psi;
    cm = surplus > cap;
    if any(cm(:))
        ct_cap = max(Rv - cb - ri*hRmax - bp, 1e-10);
        comp_cap = ct_cap.^alpha .* hRmax_rep.^(1-alpha);
        u_cap = comp_cap.^oms / oms + psi;
        u(cm) = u_cap(cm);
    end
    bpc = max(min(bp, bg_end), bg_1);
    if has_mex
        idx = lookup(bgrid, bpc, 3);
    else
        idx = max(min(discretize(bpc, [bgrid; inf]), Nb-1), 1);
    end
    w = (bpc - bgrid(idx)) ./ (bgrid(idx+1) - bgrid(idx));
    lin = idx + col_offset;
    Vbp = (1-w).*Vcont(lin) + w.*Vcont(lin+1);
    f = u + beta * Vbp;
    f(surplus <= 1e-10) = -1e10;
end

function f = gs_eval_own(bp, Rv, oc, alpha, oms, cb, ht_own, psi, beta, Vcont, bgrid, Ko, Nb, bg_end, bg_1, col_offset, has_mex)
    ct = max(Rv - oc - cb - bp, 1e-10);
    u = Ko .* ct.^(alpha*oms) / oms + psi;
    bpc = max(min(bp, bg_end), bg_1);
    if has_mex
        idx = lookup(bgrid, bpc, 3);
    else
        idx = max(min(discretize(bpc, [bgrid; inf]), Nb-1), 1);
    end
    w = (bpc - bgrid(idx)) ./ (bgrid(idx+1) - bgrid(idx));
    lin = idx + col_offset;
    Vbp = (1-w).*Vcont(lin) + w.*Vcont(lin+1);
    f = u + beta * Vbp;
    f(ct <= 1e-10) = -1e10;
end

function Pa=make_child_transition_matrix_with_matured(sd,np)
    K=numel(sd);nc=K+3;Pa=zeros(nc,nc,np);
    for nn=1:np,Pi=zeros(nc);Pi(1,1)=1;
        for k=1:K,cs=k+1;pa=min(max(1/sd(k),0),1);Pi(cs,cs)=1-pa;
            if k<K,Pi(cs,cs+1)=pa;else
                if nn==1,Pi(cs,1)=pa;elseif nn==2,Pi(cs,K+2)=pa;
                else,Pi(cs,K+3)=pa;end,end,end
        Pi(K+2,K+2)=1;Pi(K+3,K+3)=1;Pa(:,:,nn)=Pi;end
end
function P=setup_parameters()
    P=struct();P.J=60;P.da=1.0;P.age_start=18;P.J_R=47;
    P.A_f_start=1;P.A_f_end=25;P.A_m=18;P.n_parity=3;
    P.use_stochastic_aging=true;P.stage_durations=[2,4,6,6];
    P.n_child_stages=length(P.stage_durations);P.n_child_states=P.n_child_stages+3;
    P.Pi_child=make_child_transition_matrix_with_matured(P.stage_durations,P.n_parity);
    P.kappa_fert=4.5;P.kappa_loc=2.0;P.eps_fert=P.kappa_fert;P.eps_loc=P.kappa_loc;
    P.alpha_cons=0.70;P.c_bar_0=0.10;P.c_bar_n=0.12;
    P.h_bar_0=4.0;P.h_bar_jump=0.75;P.h_bar_n=0.60;
    P.child_housing_spec='jump_plus_linear';
    P.kappa_h_base=0.40;P.kappa_h_slope=0;P.psi_child=0.07;
    P.theta0=0.53;P.theta_n=0.25;P.theta1=0.01;P.u_bar=0.0;P.b_entry_fixed=0.0;
    P.beta=0.96;P.rho=1/P.beta-1;P.gamma=0;P.sigma=2.0;
    P.q=0.04;P.delta=0.02;P.tau_H=0.01;P.R_gross=1+P.q;
    P.chi=1.10;P.psi=0.06;P.phi=0.80*ones(P.n_parity,1);
    P.parent_dp_waiver=false;P.parent_dp_waiver_phi=1.0;
    P.parent_dp_waiver_locations=[];P.parent_dp_waiver_owner_rungs=[];
    P.parent_dp_waiver_birth_state_only=false;
    P.birth_dp_grant=false;P.c_min=0.01;
    P.birth_entry_grant=false;P.birth_entry_grant_amount=0.0;
    P.birth_entry_grant_locations=[];P.birth_entry_grant_owner_rungs=[];
    P.n_house=6;P.h_own_min=4.0;P.h_own_max=11.0;
    P.H_own=linspace(P.h_own_min,P.h_own_max,P.n_house)';P.hR_max=8.0;
    P.I=2;P.w_hat=[1;1];P.E_loc=[0.0;0.30];
    P.income_age_breaks=[18;25;35;45;55];
    P.income_age_values=[0.565;0.838;1.000;0.985;0.935];
    P.mu_stay=0.0;P.mu_move=5.0;P.mu_move_parent=5.0;P.mu_age_decay=0;
    P.location_choice_form='additive_due';
    P.N_0=[0.5;0.5];P.entry_shares=P.N_0/sum(P.N_0);
    P.N_target=1.0;P.E_total=1/P.J;P.entry_by_loc=P.E_total*P.entry_shares;
    P.r_bar=[0.04;0.08];P.H0=[6.2;5.3];P.eta_supply=[2;0.8];P.xi_supply=P.eta_supply;
    P.p_min=0.01;P.p_max=30;P.alpha_price=0.35;
    P.rho_hat=P.rho;P.user_cost_rate=P.q+P.delta+P.tau_H;
    P.tau_pay=0.179;P.pension_mode='balanced_stationary';
    P.pension=NaN;P.pension_by_loc=NaN(P.I,1);
    P=set_income_given_w_and_pension(P);
    P.Nb=80;P.b_min=-35;P.b_max=100;P.b_grid_power=1.5;P.n_sub=1;
    P.max_iter_eq=200;P.tol_eq=1e-4;P.lambda_eq=0.30;
    P.adaptive_price_damping=true;
    P.lambda_price_min=0.005;P.lambda_price_max=0.35;
    P.lambda_entry_min=0.005;P.lambda_entry_max=0.35;
    P.adaptive_damping_decay=0.70;P.adaptive_damping_grow=1.03;
    P.cycle_guard_factor=0.85;P.target_filter_weight=0.35;
    P.enforce_price_bounds=true;P.enforce_entry_share_floor=true;
    P.entry_share_floor=1e-4;P.housing_demand_floor_for_supply=1e-6;
    P.report_clamp_hits=true;P.kfe_wealth_interp='linear';
    P.interp_method='linear';  % 'linear', 'spline', or 'pchip'
end
function P=set_income_given_w_and_pension(P)
    income_profile=get_income_age_profile(P);
    P.income_age_profile=income_profile;
    P.pension=resolve_pension_value(P,income_profile);
    P.pension_by_loc=P.pension*ones(P.I,1);
    P.income=zeros(P.I,P.J);
    for i=1:P.I,for j=1:P.J
        if j<=P.J_R,P.income(i,j)=(1-P.tau_pay)*P.w_hat(i)*income_profile(j);
        else,P.income(i,j)=P.pension;end
    end,end
end

function pension=resolve_pension_value(P,income_profile)
    if ~isfield(P,'pension_mode') || isempty(P.pension_mode)
        mode='balanced_stationary';
    else
        mode=lower(char(string(P.pension_mode)));
    end

    switch mode
        case {'balanced_stationary','balanced','paygo','fixed'}
            avg_worker_income=compute_stationary_worker_income(P,income_profile);
            retiree_ratio=P.J_R/max(P.J-P.J_R,1);
            pension=P.tau_pay*avg_worker_income*retiree_ratio;
        case {'legacy_fixed','legacy'}
            pension=P.tau_pay*1.0*(P.J_R/max(P.J-P.J_R,1));
        case {'manual','custom'}
            if ~isfield(P,'pension') || isempty(P.pension) || ~isfinite(P.pension)
                error('Manual pension_mode requires a finite P.pension.');
            end
            pension=P.pension;
        otherwise
            error('Unknown pension_mode: %s', P.pension_mode);
    end
end

function avg_worker_income=compute_stationary_worker_income(P,income_profile)
    worker_profile=income_profile(1:P.J_R);
    if isempty(worker_profile)
        avg_worker_income=0;
        return;
    end
    if isfield(P,'entry_shares') && numel(P.entry_shares)==numel(P.w_hat)
        loc_weights=max(P.entry_shares(:),0);
        if sum(loc_weights)<=0
            loc_weights=ones(numel(P.w_hat),1);
        end
    else
        loc_weights=ones(numel(P.w_hat),1);
    end
    loc_weights=loc_weights/sum(loc_weights);
    avg_wage=sum(loc_weights.*P.w_hat(:));
    avg_worker_income=avg_wage*mean(worker_profile);
end

function income_profile=get_income_age_profile(P)
    income_profile=ones(P.J,1);
    if isfield(P,'income_age_breaks') && isfield(P,'income_age_values') && ...
            numel(P.income_age_breaks)==numel(P.income_age_values) && ...
            ~isempty(P.income_age_breaks)
        age_breaks=P.income_age_breaks(:);
        age_values=P.income_age_values(:);
    else
        age_breaks=[18;25;35;45;55];
        age_values=[0.565;0.838;1.000;0.985;0.935];
    end
    age_vec=P.age_start+(0:P.J-1)'*P.da;
    for k=1:numel(age_values)
        if k < numel(age_values)
            age_hi=age_breaks(k+1);
        else
            age_hi=P.age_start+P.J_R*P.da;
        end
        mask=(age_vec>=age_breaks(k)) & (age_vec<age_hi);
        income_profile(mask)=age_values(k);
    end
end

function ps=get_phi_state_matrix(P)
    ps=repmat(P.phi(:),1,P.n_child_states);
    if ~(isfield(P,'parent_dp_waiver')&&P.parent_dp_waiver),return;end
    kp=get_parent_target_child_states(P);
    po=1.0;if isfield(P,'parent_dp_waiver_phi'),po=P.parent_dp_waiver_phi;end
    if P.n_parity>=2,ps(2:end,kp)=max(ps(2:end,kp),po);end
end

function pc=get_phi_choice_tensor(P)
    np=P.n_parity;ncs=P.n_child_states;nt=1+P.n_house;I=P.I;
    base_phi=repmat(reshape(P.phi(:),1,1,np,1),I,nt,1,ncs);
    pc=base_phi;
    if ~(isfield(P,'parent_dp_waiver')&&P.parent_dp_waiver),return;end
    if np < 2,return;end
    po=1.0;if isfield(P,'parent_dp_waiver_phi'),po=P.parent_dp_waiver_phi;end
    loc_idx=get_parent_target_locations(P);
    ten_idx=get_parent_target_owner_tenures(P);
    cs_idx=get_parent_target_child_states(P);
    if isempty(loc_idx)||isempty(ten_idx)||~any(cs_idx),return;end
    pc(loc_idx,ten_idx,2:end,cs_idx)=max(pc(loc_idx,ten_idx,2:end,cs_idx),po);
end

function loc_idx=get_parent_target_locations(P)
    loc_idx=1:P.I;
    if isfield(P,'parent_dp_waiver_locations')&&~isempty(P.parent_dp_waiver_locations)
        loc_idx=P.parent_dp_waiver_locations(:)';
    end
    loc_idx=loc_idx(loc_idx>=1 & loc_idx<=P.I);
end

function ten_idx=get_parent_target_owner_tenures(P)
    ten_idx=2:(1+P.n_house);
    if isfield(P,'parent_dp_waiver_owner_rungs')&&~isempty(P.parent_dp_waiver_owner_rungs)
        ten_idx=1+P.parent_dp_waiver_owner_rungs(:)';
    end
    ten_idx=ten_idx(ten_idx>=2 & ten_idx<=1+P.n_house);
end

function kp=get_parent_target_child_states(P)
    kp=false(1,P.n_child_states);
    if isfield(P,'parent_dp_waiver_birth_state_only')&&P.parent_dp_waiver_birth_state_only
        if P.n_child_states >= 2
            kp(2)=true;
        end
    else
        K=P.n_child_stages;
        kp(2:K+1)=true;
    end
end

function bg=get_birth_entry_grant_tensor(P)
    np=P.n_parity;ncs=P.n_child_states;nt=1+P.n_house;I=P.I;
    bg=zeros(I,nt,np,ncs);
    if ~(isfield(P,'birth_entry_grant')&&P.birth_entry_grant),return;end
    g=0.0;
    if isfield(P,'birth_entry_grant_amount'),g=P.birth_entry_grant_amount;end
    if ~(isscalar(g)&&isfinite(g)&&g>0),return;end
    if np<2||ncs<2,return;end

    loc_idx=1:I;
    if isfield(P,'birth_entry_grant_locations')&&~isempty(P.birth_entry_grant_locations)
        loc_idx=P.birth_entry_grant_locations(:)';
    end
    loc_idx=loc_idx(loc_idx>=1 & loc_idx<=I);

    ten_idx=2:nt;
    if isfield(P,'birth_entry_grant_owner_rungs')&&~isempty(P.birth_entry_grant_owner_rungs)
        ten_idx=1+P.birth_entry_grant_owner_rungs(:)';
    end
    ten_idx=ten_idx(ten_idx>=2 & ten_idx<=nt);

    if isempty(loc_idx)||isempty(ten_idx),return;end
    bg(loc_idx,ten_idx,2:end,2)=g;
end
