%% run_model.m
%
%  Spatial Lifecycle Fertility Model - STONE-GEARY with LOGIT SHOCKS
%
%  ========================================================================
%  CHANGES FROM v9:
%
%  NEW SHOCK STRUCTURE: Type-I Extreme Value (Gumbel) / Standard Logit
%
%  Location choice:
%    π^{i'} = (E_{i'} * μ_{ii'}) * exp(ν_ℓ * V^{i'}) / Σ_j (E_j * μ_{ij}) * exp(ν_ℓ * V^j)
%    V^I = (1/ν_ℓ) * log(Σ_j (E_j * μ_{ij}) * exp(ν_ℓ * V^j))
%
%  Fertility choice:
%    π^{n'} = exp(ν^n * V^I(n')) / Σ_m exp(ν^n * V^I(m))
%    V^fert = (1/ν^n) * log(Σ_m exp(ν^n * V^I(m)))
%
%  Previous v9: used Fréchet with V^ν (power formulation)
%  Current: uses Gumbel with exp(ν*V) (standard logit formulation)
%
%  ========================================================================
%  TIMING (per Model_notes_spring26.tex):
%
%  At shock ages, the sequence is: fertility → location → tenure/size
%
%  For BACKWARD INDUCTION (HJB):
%    1. Drift: solve within-period HJB using V(j+1) as continuation
%    2. Tenure choice: V_H = max over h' (includes resizing while staying)
%    3. Location logit: V_I = (1/ν_ℓ) log Σ exp(ν_ℓ V^{i'})
%    4. Fertility logit: V^n = (1/ν^n) log Σ exp(ν^n V^I(n'))
%
%  For FORWARD ITERATION (KFE):
%    1. Fertility jump (if in fertile window)
%    2. Location jump
%    3. Tenure jump (includes resizing while staying)
%    4. Drift (KFE step)
%    5. Apply demographic factor exp(-n_pop)
%    6. Advance to j+1
%
%  ========================================================================

function [sol, P, p_eq] = run_model(P_override)

    %% Start timer
    t_start = tic;

    %% Setup logging - always save a timestamped log
    log_dir = fullfile(fileparts(mfilename('fullpath')), 'logs');
    if ~exist(log_dir, 'dir')
        mkdir(log_dir);
    end
    timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
    log_file = fullfile(log_dir, sprintf('run_model_%s.log', timestamp));
    diary(log_file);
    fprintf('Log file: %s\n\n', log_file);

    %% Setup parameters
    P = setup_parameters();

    if nargin >= 1 && ~isempty(P_override)
        fields = fieldnames(P_override);
        for k = 1:length(fields)
            P.(fields{k}) = P_override.(fields{k});
        end
        if isfield(P_override, 'b_entry_base')
            P.b_entry_loc = P.b_entry_base * ones(3, 1);
        end
    end

    fprintf('============================================================\n');
    fprintf('  SPATIAL LIFECYCLE FERTILITY MODEL - LOGIT SHOCKS\n');
    fprintf('  J=%d periods, Nb=%d wealth grid points\n', P.J, P.Nb);
    fprintf('  \n');
    fprintf('  STONE-GEARY PREFERENCES:\n');
    fprintf('    c_bar(n) = %.2f + %.2f * n  (subsistence consumption)\n', P.c_bar_0, P.c_bar_n);
    fprintf('    h_bar(n) = %.2f + %.2f*(n>=1) + %.2f*n  (subsistence housing)\n', P.h_bar_0, P.h_bar_jump, P.h_bar_n);
    fprintf('    kappa_h = %.2f  (constant, slope=0)\n', P.kappa_h_base);
    fprintf('    psi(n) = %.2f * n\n', P.psi_child);
    fprintf('    theta_n = %.2f (terminal warm-glow)\n', P.theta_n);
    fprintf('  \n');
    fprintf('  SHOCK STRUCTURE (Type-I EV / Logit):\n');
    fprintf('    nu_loc = %.2f (location choice dispersion)\n', P.nu_loc);
    fprintf('    nu_fert = %.2f (fertility choice dispersion)\n', P.nu_fert);
    fprintf('  \n');
    fprintf('  NUMERICAL:\n');
    fprintf('    u_bar = %.1f (keeps V > 0)\n', P.u_bar);
    fprintf('============================================================\n\n');

    %% Setup grids
    Nb = P.Nb;
    u = linspace(0, 1, Nb)';
    b_grid = P.b_min + (P.b_max - P.b_min) * u.^P.b_grid_power;

    [~, idx_zero] = min(abs(b_grid));
    b_grid(idx_zero) = 0;

    db = diff(b_grid);
    db_f = [db; db(end)];
    db_b = [db(1); db];

    %% Pre-compute finite difference denominators
    inv_db_f = 1 ./ db_f;
    inv_db_b = 1 ./ db_b;

    %% Initial prices
    p_init = P.r_bar / P.user_cost_rate;

    %% Solve BGP equilibrium with BISECTION
    fprintf('=== Solving BGP Equilibrium (Bisection) ===\n\n');

    [n_pop_eq, p_eq, sol, P] = find_bgp_bisection(p_init, P, b_grid, db_f, db_b, inv_db_f, inv_db_b, idx_zero);

    gamma_eq = P.alpha * n_pop_eq;

    %% Display results
    fprintf('\n==================== BGP RESULTS ====================\n');
    fprintf('Population growth rate n: %.4f per period (%.4f%% annual)\n', ...
        n_pop_eq, 100*n_pop_eq/P.da);
    fprintf('Aggregate growth rate gamma: %.4f per period (%.4f%% annual)\n', ...
        gamma_eq, 100*gamma_eq/P.da);
    fprintf('Mean parity: %.2f\n', sol.mean_parity);
    fprintf('Ownership rate: %.1f%%\n', 100*sol.own_rate);

    fprintf('\nEquilibrium prices: [%.3f, %.3f, %.3f]\n', p_eq(1), p_eq(2), p_eq(3));
    fprintf('Equilibrium rents: [%.3f, %.3f, %.3f]\n', ...
        P.user_cost_rate*p_eq(1), P.user_cost_rate*p_eq(2), P.user_cost_rate*p_eq(3));

    fprintf('\n*** SPATIAL FERTILITY GRADIENT ***\n');
    loc_names = {'Peripheral', 'Secondary', 'Superstar'};
    for i = 1:P.I
        fprintf('  %s: mean_n=%.3f, childless=%.1f%%\n', loc_names{i}, ...
            sol.mean_parity_by_loc(i), 100*sol.frac_childless_by_loc(i));
    end
    gradient = sol.mean_parity_by_loc(1) - sol.mean_parity_by_loc(3);
    fprintf('  Gradient (Peripheral - Superstar): %.3f\n', gradient);

    % Report subsistence costs by location
    fprintf('\n*** SUBSISTENCE COSTS BY LOCATION (for n=2) ***\n');
    for i = 1:P.I
        r_i = P.user_cost_rate * p_eq(i);
        c_bar_2 = P.c_bar_0 + P.c_bar_n * 2;
        h_bar_2 = P.h_bar_0 + P.h_bar_jump + P.h_bar_n * 2;
        min_cost = c_bar_2 + r_i * h_bar_2;
        fprintf('  %s: r=%.3f, min_expenditure=%.3f (c_bar=%.2f + r*h_bar=%.3f)\n', ...
            loc_names{i}, r_i, min_cost, c_bar_2, r_i * h_bar_2);
    end
    fprintf('=====================================================\n');

    %% Verify BGP consistency
    verify_bgp_consistency(sol, P, n_pop_eq);

    %% Report timing
    elapsed = toc(t_start);
    fprintf('\n========================================\n');
    fprintf('Total elapsed time: %.1f seconds (%.1f minutes)\n', elapsed, elapsed/60);
    fprintf('========================================\n');

    diary off;
end


%% ========================================================================
%%                         SETUP PARAMETERS
%% ========================================================================
function P = setup_parameters()

    P = struct();

    %% Lifecycle
    P.J = 10;
    P.da = 6.0;
    P.age_start = 25;
    P.J_R = 6;

    %% Fertility
    P.A_f_start = 1;
    P.A_f_end = 3;
    P.A_m = 3;
    P.n_parity = 4;
    P.n_child_states = P.A_m + 2;

    %% ===== SHOCK PARAMETERS (Type-I EV / Logit) =====
    % These are SCALE parameters (inverse dispersion)
    % Higher nu = less noise = more deterministic
    P.nu_fert = 8.0;      % Fertility choice scale
    P.nu_loc = 2.0;       % Location choice scale

    %% ===== STONE-GEARY PARAMETERS =====

    % Subsistence consumption: c_bar(n) = c_bar_0 + c_bar_n * n
    P.c_bar_0 = 0.05;    % Baseline subsistence consumption
    P.c_bar_n = 0.05;    % Additional consumption per child

    % Subsistence housing: h_bar(n) = h_bar_0 + h_bar_jump*(n>=1) + h_bar_n * n
    % THIS IS THE KEY PARAMETER FOR SPATIAL CHILDLESSNESS
    P.h_bar_0 = 0.03;      % Baseline subsistence housing
    P.h_bar_jump = 0.50;   % Jump for first child (need a bedroom!)
    P.h_bar_n = 0.20;      % Additional housing per child

    %% ===== PREFERENCES =====

    % Housing preference: kappa_h(n) = kappa_h_base + kappa_h_slope * n
    P.kappa_h_base = 0.40;
    P.kappa_h_slope = 0;  % Stone-Geary h_bar_n does this job

    % Direct child utility: psi(n) = psi_child * n
    P.psi_child = 0.0;

    % Terminal warm-glow: theta_n * n
    P.theta_n = 1.5;

    % Bequest utility parameters
    P.theta0 = 1.5;
    P.theta1 = 0.01;

    %% ===== NUMERICAL =====
    P.u_bar = 8.0;

    %% Estate recycling
    P.phi_estate = 0.80;
    P.b_entry_base = 0.5;

    %% Discount and growth
    P.rho = 0.02 * P.da;
    P.gamma_init = 0.00;
    P.sigma = 2.0;

    %% Interest and housing costs
    P.q = 0.04 * P.da;
    P.delta = 0.02 * P.da;
    P.tau_H = 0.01 * P.da;

    %% Other preferences
    P.chi = 1.05;           % Owner-occupied housing service flow
    P.psi = 0.08;           % Transaction cost
    P.phi = 0.65 * ones(P.n_parity, 1);  % LTV limit
    P.dV_floor = 1e-4;
    P.c_min = 0.05;         % Numerical floor (below subsistence)

    %% Housing
    P.n_house = 3;
    P.h_own_min = 1.5;
    P.h_own_max = 4.0;
    P.H_own = linspace(P.h_own_min, P.h_own_max, P.n_house)';
    P.hR_max = 2.0;

    %% Locations
    P.I = 3;
    P.w_bar = [0.85; 1.0; 1.20];
    P.w_hat = P.w_bar;
    P.eta_agglom = 0.10;

    % Location amenities E_i (in notes: script E_i')
    P.E_loc = [0.95; 1.0; 1.10];

    %% Mobility (moving wedges mu_{ii'})
    P.mu_stay = 1.0;        % mu_{ii} = 1 (no penalty for staying)
    P.mu_move = 0.70;       % mu_{ij} for j != i (moving cost as wedge)
    P.mu_move_parent = 0.70;
    P.mu_age_decay = 0;

    %% Supply
    P.N_0 = [0.25; 0.35; 0.40];
    P.entry_by_loc = P.N_0 / sum(P.N_0) / P.J;
    P.b_entry_loc = P.b_entry_base * ones(3, 1);
    P.r_bar = [0.025; 0.050; 0.180] * P.da;
    P.alpha = 0.3;
    P.xi_supply = (1/P.alpha) * ones(P.I, 1);
    P.p_min = 0.1;
    P.p_max = 10.0;
    P.alpha_price = 0.15;

    %% Derived
    P.gamma = P.gamma_init;
    P.rho_hat = P.rho - (1 - P.sigma) * P.gamma;
    P.user_cost_rate = P.q + P.delta + P.tau_H - P.gamma;

    P.pension_replacement = 0.40;
    P.pension = P.pension_replacement * mean(P.w_bar);

    P.income = zeros(P.I, P.J);
    for i = 1:P.I
        for j = 1:P.J
            if j <= P.J_R
                P.income(i, j) = P.w_hat(i);
            else
                P.income(i, j) = P.pension;
            end
        end
    end

    %% Numerical
    P.Nb = 30;
    P.b_min = -6.0;
    P.b_max = 100.0;
    P.b_grid_power = 1.5;
    P.n_sub = 6;

    % ===== OUTER ROOT SOLVER (gamma bisection) =====
    P.max_iter_bisect = 80;     % Max bisection iterations
    P.tol_F = 5e-5;             % Convergence tolerance on |F(gamma)|
                                % F(gamma) = gamma - alpha * n_impl(gamma)
                                % This is the PRIMARY stopping criterion

    % ===== INNER EQUILIBRIUM TOLERANCES (tightened for stable F) =====
    % These must be tight enough that B(gamma) and E(gamma) are stable
    P.max_iter_entry = 60;      % Entry/bequest fixed point (was 40)
    P.tol_entry = 0.005;        % Tighter entry convergence (was 0.01)
    P.alpha_entry = 0.4;        % Faster mixing (was 0.3)
    P.alpha_beq = 0.15;         % Faster bequest update (was 0.1)

    P.max_iter_price = 150;     % Price fixed point (was 100)
    P.tol_price = 0.005;        % Tighter price convergence (was 0.01)

    % ===== GAMMA SEARCH RANGE =====
    P.gamma_lo = -0.05;
    P.gamma_hi = 0.15;
end


%% ========================================================================
%%                   BGP SOLVER: OUTER ROOT PROBLEM F(gamma) = 0
%%
%%   Define: F(gamma) = gamma - alpha * n_impl(gamma)
%%   where:  n_impl(gamma) = (1/A_m) * log(B(gamma)/E(gamma))
%%
%%   B(gamma) and E(gamma) come from fully solved inner equilibrium at gamma.
%%   This is a THEORETICALLY CORRECT formulation where we solve F(gamma) = 0.
%%
%%   KEY REQUIREMENTS:
%%   1. evaluate_equilibrium_given_gamma is SIDE-EFFECT FREE (no mutation of P0)
%%   2. Stopping criterion is |F(gamma)| < tol_F (NOT bracket width)
%%   3. Inner equilibrium is solved tightly enough that B, E are stable
%% ========================================================================
function [n_pop_eq, p_eq, sol, P] = find_bgp_bisection(p_init, P, b_grid, db_f, db_b, inv_db_f, inv_db_b, idx_zero)

    % =====================================================================
    % SETUP: Store immutable baseline state P0
    % =====================================================================
    P0 = P;  % This is NEVER modified - each evaluation starts from P0

    % Initial warmstart state (prices, entry shares, bequest locations)
    warmstart = struct();
    warmstart.p = p_init;
    warmstart.entry_shares = P0.N_0 / sum(P0.N_0);  % Initial guess
    warmstart.b_entry_loc = P0.b_entry_base * ones(P0.I, 1);

    gamma_lo = P0.gamma_lo;
    gamma_hi = P0.gamma_hi;

    fprintf('  ========================================================\n');
    fprintf('  OUTER ROOT PROBLEM: F(gamma) = gamma - alpha * n_impl(gamma) = 0\n');
    fprintf('  where n_impl = (1/A_m) * log(B/E)\n');
    fprintf('  Tolerance: |F| < %.2e\n', P0.tol_F);
    fprintf('  ========================================================\n\n');

    fprintf('  Evaluating bracket endpoints...\n');
    fprintf('  (gamma = alpha * n_pop, alpha = %.2f)\n\n', P0.alpha);

    % =====================================================================
    % BRACKET FINDING
    % =====================================================================
    [f_lo, ~, warmstart] = evaluate_equilibrium_given_gamma(gamma_lo, P0, warmstart, ...
        b_grid, db_f, db_b, inv_db_f, inv_db_b, idx_zero, true);
    fprintf('    F(gamma=%.5f) = %.6f\n', gamma_lo, f_lo);

    [f_hi, ~, warmstart] = evaluate_equilibrium_given_gamma(gamma_hi, P0, warmstart, ...
        b_grid, db_f, db_b, inv_db_f, inv_db_b, idx_zero, true);
    fprintf('    F(gamma=%.5f) = %.6f\n', gamma_hi, f_hi);

    % Expand brackets if needed
    max_expand = 8;
    expand_count = 0;
    while f_lo * f_hi > 0 && expand_count < max_expand
        expand_count = expand_count + 1;
        if f_lo > 0 && f_hi > 0
            gamma_lo = gamma_lo - 0.03;
            [f_lo, ~, warmstart] = evaluate_equilibrium_given_gamma(gamma_lo, P0, warmstart, ...
                b_grid, db_f, db_b, inv_db_f, inv_db_b, idx_zero, true);
            fprintf('    Expanded lo: F(gamma=%.5f) = %.6f\n', gamma_lo, f_lo);
        elseif f_lo < 0 && f_hi < 0
            gamma_hi = gamma_hi + 0.02;
            [f_hi, ~, warmstart] = evaluate_equilibrium_given_gamma(gamma_hi, P0, warmstart, ...
                b_grid, db_f, db_b, inv_db_f, inv_db_b, idx_zero, true);
            fprintf('    Expanded hi: F(gamma=%.5f) = %.6f\n', gamma_hi, f_hi);
        end
    end

    if f_lo * f_hi > 0
        warning('Could not bracket root after %d expansions! F_lo=%.4f, F_hi=%.4f', ...
            max_expand, f_lo, f_hi);
    end

    fprintf('\n  Starting bisection on gamma in [%.5f, %.5f]\n', gamma_lo, gamma_hi);
    fprintf('  (corresponds to n_pop in [%.5f, %.5f])\n', gamma_lo/P0.alpha, gamma_hi/P0.alpha);
    fprintf('  Convergence requires: |F(gamma)| < %.2e\n\n', P0.tol_F);

    % =====================================================================
    % BISECTION WITH RESIDUAL-BASED STOPPING
    % =====================================================================
    gamma_star = (gamma_lo + gamma_hi) / 2;
    sol_star = [];
    result_star = struct();
    converged = false;
    best_F = inf;
    best_gamma = gamma_star;
    best_result = struct();

    for iter = 1:P0.max_iter_bisect
        gamma_mid = (gamma_lo + gamma_hi) / 2;

        % PURE evaluation: P0 is never modified, warmstart is passed explicitly
        [f_mid, result_mid, warmstart] = evaluate_equilibrium_given_gamma(gamma_mid, P0, warmstart, ...
            b_grid, db_f, db_b, inv_db_f, inv_db_b, idx_zero, false);

        n_mid = gamma_mid / P0.alpha;
        fprintf('  Bisect %2d: gamma=%.6f (n=%.5f), F=%.2e, n_impl=%.5f, own=%.1f%%, parity=%.2f\n', ...
            iter, gamma_mid, n_mid, f_mid, result_mid.n_impl, ...
            100*result_mid.sol.own_rate, result_mid.sol.mean_parity);

        % Track best solution found
        if abs(f_mid) < best_F
            best_F = abs(f_mid);
            best_gamma = gamma_mid;
            best_result = result_mid;
        end

        % CONVERGENCE CHECK: RESIDUAL ONLY (not bracket width!)
        if abs(f_mid) < P0.tol_F
            gamma_star = gamma_mid;
            result_star = result_mid;
            converged = true;
            fprintf('\n  *** CONVERGED: |F| = %.2e < tol = %.2e ***\n', abs(f_mid), P0.tol_F);
            break;
        end

        % Update brackets (standard bisection)
        if f_lo * f_mid < 0
            gamma_hi = gamma_mid;
            f_hi = f_mid;
        else
            gamma_lo = gamma_mid;
            f_lo = f_mid;
        end

        gamma_star = gamma_mid;
        result_star = result_mid;
    end

    if ~converged
        fprintf('\n  WARNING: Bisection did not converge to |F| < %.2e\n', P0.tol_F);
        fprintf('           Best |F| found: %.2e at gamma = %.6f\n', best_F, best_gamma);
        fprintf('           Bracket width: %.2e\n', gamma_hi - gamma_lo);

        % Use best solution found
        gamma_star = best_gamma;
        result_star = best_result;
    end

    % =====================================================================
    % FINAL CONSISTENCY CHECK AND OUTPUT
    % =====================================================================
    sol = result_star.sol;
    p_eq = result_star.p;

    % Reconstruct P with the converged gamma (but don't modify P0!)
    P = P0;
    P.gamma = gamma_star;
    P.rho_hat = P0.rho - (1 - P0.sigma) * gamma_star;
    P.user_cost_rate = P0.q + P0.delta + P0.tau_H - gamma_star;
    P.entry_by_loc = result_star.entry_by_loc;
    P.b_entry_loc = result_star.b_entry_loc;
    P.w_hat = sol.w_hat;

    n_pop_eq = gamma_star / P.alpha;
    sol.n_pop = n_pop_eq;
    sol.gamma = gamma_star;
    sol.b_entry_loc = P.b_entry_loc;
    sol.entry_by_loc = P.entry_by_loc;

    % Final diagnostics
    fprintf('\n  ========================================================\n');
    fprintf('  FINAL CONSISTENCY CHECK:\n');
    fprintf('    gamma_star  = %.8f\n', gamma_star);
    fprintf('    n_used      = %.8f (= gamma_star / alpha)\n', n_pop_eq);
    fprintf('    n_impl      = %.8f (= log(B/E) / A_m)\n', result_star.n_impl);
    fprintf('    B           = %.8f\n', result_star.B);
    fprintf('    E           = %.8f\n', result_star.E);
    fprintf('    F(gamma)    = %.2e (should be < %.2e)\n', ...
        gamma_star - P.alpha * result_star.n_impl, P0.tol_F);
    fprintf('    n_used - n_impl = %.2e\n', n_pop_eq - result_star.n_impl);
    fprintf('  --------------------------------------------------------\n');
    fprintf('  MATURITY FLOW VERIFICATION (kids follow parents):\n');
    fprintf('    M_total (maturing children)   = %.8f\n', result_star.M_total);
    fprintf('    E (entry rate)                = %.8f\n', result_star.E);
    fprintf('    B*exp(-n*A_m) (theory)        = %.8f\n', result_star.B * exp(-n_pop_eq * P0.A_m));
    fprintf('    M_total / E                   = %.6f (should be ~1.0)\n', result_star.M_total / max(result_star.E, 1e-12));
    fprintf('    E / [B*exp(-n*A_m)]           = %.6f (should be ~1.0)\n', result_star.E / max(result_star.B * exp(-n_pop_eq * P0.A_m), 1e-12));
    fprintf('  --------------------------------------------------------\n');
    fprintf('  ENTRY SHARES (maturity-based):\n');
    for i = 1:P.I
        fprintf('    Location %d: mature_share=%.4f, entry_share=%.4f\n', i, result_star.mature_shares(i), result_star.entry_shares(i));
    end
    fprintf('  ========================================================\n');
end


%% ========================================================================
%%   EVALUATE EQUILIBRIUM GIVEN GAMMA - PURE, SIDE-EFFECT FREE
%%
%%   This function NEVER modifies P0. It creates a local copy P_local
%%   and returns all state changes via the output warmstart struct.
%%
%%   PERFORMANCE OPTIMIZATION (Option C):
%%   HJB depends on prices/wages, NOT on entry distribution or bequests.
%%   So we solve: Price loop { HJB once, then Entry/KFE loop to convergence }
%%   This reduces HJB solves from ~1200 to ~400 (3x speedup on HJB).
%%
%%   Returns:
%%     F_gamma   - residual F(gamma) = gamma - alpha * n_impl
%%     result    - struct with sol, p, B, E, n_impl, entry_by_loc, b_entry_loc
%%     warmstart - updated warmstart for next evaluation
%% ========================================================================
function [F_gamma, result, warmstart_out] = evaluate_equilibrium_given_gamma(gamma, P0, warmstart, ...
    b_grid, db_f, db_b, inv_db_f, inv_db_b, idx_zero, verbose)

    % =====================================================================
    % CREATE LOCAL COPY OF PARAMETERS (P0 is NEVER modified)
    % =====================================================================
    P_local = P0;
    P_local.gamma = gamma;
    P_local.rho_hat = P0.rho - (1 - P0.sigma) * gamma;
    P_local.user_cost_rate = P0.q + P0.delta + P0.tau_H - gamma;

    % Safety bounds
    if P_local.user_cost_rate <= 0.01
        P_local.user_cost_rate = 0.01;
    end
    if P_local.rho_hat <= 0.01
        P_local.rho_hat = 0.01;
    end

    n_pop = gamma / P0.alpha;  % The n we're ASSUMING for this gamma

    % =====================================================================
    % INITIALIZE FROM WARMSTART (explicit, not hidden state)
    % =====================================================================
    p = warmstart.p;
    entry_shares = warmstart.entry_shares;
    b_entry_loc = warmstart.b_entry_loc;

    total_entry = 1 / P0.J;  % Normalized entry rate
    P_local.entry_by_loc = total_entry * entry_shares;
    P_local.b_entry_loc = b_entry_loc;

    % =====================================================================
    % RESTRUCTURED FIXED POINT: Price { HJB once, Entry/KFE loop }
    %
    % Key insight: HJB policies depend on (p, r, w) but NOT on entry shares
    % or bequest distribution. So within a price iteration:
    %   1. Solve HJB once -> get policies
    %   2. Iterate KFE + entry/bequest to convergence using fixed policies
    %   3. Check price residual -> update prices
    % =====================================================================

    best_err = inf;
    best_p = p;
    best_sol = [];

    for iter_price = 1:P0.max_iter_price
        r = P_local.user_cost_rate * p;

        % =================================================================
        % STEP 1: SOLVE HJB ONCE FOR THESE PRICES
        % =================================================================
        [V, c_pol, hR_pol, tenure_choice, loc_probs, fert_probs, fert_value] = ...
            solve_hjb_logit(r, p, P_local, b_grid, db_f, db_b, inv_db_f, inv_db_b);

        % =================================================================
        % STEP 2: ITERATE KFE + ENTRY/BEQUEST TO CONVERGENCE (fixed HJB)
        % =================================================================
        for iter_entry = 1:P0.max_iter_entry
            % Solve KFE with current entry distribution
            [g, stats] = solve_kfe(c_pol, hR_pol, tenure_choice, loc_probs, fert_probs, ...
                                   r, p, P_local, b_grid, db_f, db_b, idx_zero, n_pop);

            % Compute maturity shares from KFE
            mature_shares = stats.mature_entry_shares;
            b_entry_implied = stats.bequest_per_entrant;

            % Convergence check for entry/bequest
            err_shares = max(abs(mature_shares - entry_shares));
            err_beq = max(abs(b_entry_implied - b_entry_loc) ./ max(abs(b_entry_loc), 0.1));
            err_entry = max(err_shares, err_beq);

            if err_entry < P0.tol_entry
                break;
            end

            % Update entry distribution (Picard iteration)
            entry_shares = (1 - P0.alpha_entry) * entry_shares + P0.alpha_entry * mature_shares;
            P_local.entry_by_loc = total_entry * entry_shares;

            b_entry_loc = (1 - P0.alpha_beq) * b_entry_loc + P0.alpha_beq * b_entry_implied;
            b_entry_loc = max(b_entry_loc, 0.1);
            P_local.b_entry_loc = b_entry_loc;
        end

        % =================================================================
        % STEP 3: UPDATE WAGES AND CHECK PRICE CONVERGENCE
        % =================================================================
        w_new = zeros(P0.I, 1);
        for i = 1:P0.I
            pop_ratio = max(stats.pop_share(i), 0.05) / P0.N_0(i);
            w_new(i) = P0.w_bar(i) * pop_ratio^P0.eta_agglom;
        end
        P_local.w_hat = w_new;
        for i = 1:P0.I
            for j = 1:P0.J
                if j <= P0.J_R
                    P_local.income(i, j) = P_local.w_hat(i);
                else
                    P_local.income(i, j) = P0.pension;
                end
            end
        end

        % Compute price residual
        r_supply = zeros(P0.I, 1);
        for i = 1:P0.I
            pop_ratio = max(stats.pop_share(i), 0.05) / P0.N_0(i);
            r_supply(i) = P0.r_bar(i) * pop_ratio^(1/P0.xi_supply(i));
        end
        p_supply = r_supply / P_local.user_cost_rate;

        excess = (p - p_supply) ./ p_supply;
        max_price_err = max(abs(excess));

        % Track best solution
        if max_price_err < best_err
            best_err = max_price_err;
            best_p = p;
            best_sol = struct('V', V, 'c_pol', c_pol, 'hR_pol', hR_pol, ...
                'tenure_choice', tenure_choice, 'loc_probs', loc_probs, ...
                'fert_probs', fert_probs, 'fert_value', fert_value, 'g', g, ...
                'own_rate', stats.own_rate, 'pop_share', stats.pop_share, ...
                'own_by_loc', stats.own_by_loc, 'own_by_age', stats.own_by_age, ...
                'parity_dist', stats.parity_dist, 'mean_parity', stats.mean_parity, ...
                'mean_parity_by_loc', stats.mean_parity_by_loc, ...
                'frac_childless_by_loc', stats.frac_childless_by_loc, ...
                'death_wealth_by_loc', stats.death_wealth_by_loc, ...
                'entry_mass_by_loc', stats.entry_mass_by_loc, ...
                'bequest_per_entrant', stats.bequest_per_entrant, ...
                'own_by_parity', stats.own_by_parity, ...
                'child_state_dist', stats.child_state_dist, ...
                'fert_by_age', stats.fert_by_age, ...
                'total_births_kfe', stats.total_births_kfe, ...
                'births_by_loc', stats.births_by_loc, ...
                'entry_rate', stats.entry_rate, ...
                'entrants_mature_by_loc', stats.entrants_mature_by_loc, ...
                'entrants_mature_total', stats.entrants_mature_total, ...
                'mature_entry_shares', stats.mature_entry_shares, ...
                'w_hat', P_local.w_hat, ...
                'p_eq', p);
        end

        if max_price_err < P0.tol_price
            break;
        end

        % Update prices
        p = (1 - P0.alpha_price) * p + P0.alpha_price * p_supply;
        p = max(p, P0.p_min);
        p = min(p, P0.p_max);
    end

    % Use best solution found
    p = best_p;
    sol = best_sol;

    % =====================================================================
    % COMPUTE n_impl FROM FLOW IDENTITY
    % =====================================================================
    B = sol.total_births_kfe;
    E = sol.entry_rate;
    M_total = sol.entrants_mature_total;  % Maturing children flow

    if E > 1e-12 && B > 1e-12
        n_impl = log(B / E) / P0.A_m;
    else
        n_impl = 0;
    end

    gamma_impl = P0.alpha * n_impl;

    % THE RESIDUAL: F(gamma) = gamma - alpha * n_impl(gamma)
    F_gamma = gamma - gamma_impl;

    % Theoretical check: M_total should equal E which should equal B*exp(-n*A_m)
    % This verifies we're counting CHILDREN not PARENTS
    n_used = gamma / P0.alpha;
    E_theory = B * exp(-n_used * P0.A_m);

    if verbose
        fprintf('      [Flow] B=%.6f, E=%.6f, B/E=%.6f\n', B, E, B/max(E,1e-12));
        fprintf('      [BGP] n_impl=%.6f, gamma_impl=%.6f, F(gamma)=%.6f\n', n_impl, gamma_impl, F_gamma);
        fprintf('      [Maturity check] M_total=%.6f, E=%.6f, B*exp(-n*A_m)=%.6f\n', M_total, E, E_theory);
        fprintf('      [Maturity gaps] M_total/E=%.4f, E/E_theory=%.4f\n', M_total/max(E,1e-12), E/max(E_theory,1e-12));
    end

    % =====================================================================
    % PACK OUTPUT
    % =====================================================================
    result = struct();
    result.sol = sol;
    result.p = p;
    result.B = B;
    result.E = E;
    result.M_total = M_total;
    result.n_impl = n_impl;
    result.gamma_impl = gamma_impl;
    result.entry_by_loc = P_local.entry_by_loc;
    result.b_entry_loc = b_entry_loc;
    result.entry_shares = entry_shares;
    result.mature_shares = sol.mature_entry_shares;

    % Update warmstart for next call (pass forward the converged state)
    warmstart_out = struct();
    warmstart_out.p = p;
    warmstart_out.entry_shares = entry_shares;
    warmstart_out.b_entry_loc = b_entry_loc;
end


%% ========================================================================
%%                    HJB SOLVER WITH LOGIT SHOCKS
%%                    CORRECT TIMING: DRIFT -> TENURE -> LOCATION -> FERTILITY
%% ========================================================================
function [V, c_pol, hR_pol, tenure_choice, loc_probs, fert_probs, fert_value] = ...
    solve_hjb_logit(r_hat, p_hat, P, b_grid, db_f, db_b, inv_db_f, inv_db_b)

    % TIMING (for backward induction at age j):
    %   1. Start with V(j+1) - continuation values at next shock age
    %   2. DRIFT: Solve within-period HJB to get V_drift at age j
    %   3. TENURE: V_H = max over h' of V_drift (after accounting for buy/sell)
    %   4. LOCATION: V_I = logit aggregation over destinations using V_H
    %   5. FERTILITY: V^n = logit aggregation over parities using V_I (if fertile)
    %   6. Store V^n (or V_I if not fertile) as V(j)

    J = P.J;
    I = P.I;
    Nb = length(b_grid);
    n_h = P.n_house;
    n_tenure = 1 + n_h;
    n_parity = P.n_parity;
    n_child = P.n_child_states;
    Delta = P.da;

    % Pre-allocate outputs
    % V is the "pre-shock" value (what enters the shock sequence)
    V = zeros(Nb, n_tenure, I, J, n_parity, n_child);
    c_pol = zeros(Nb, n_tenure, I, J, n_parity, n_child);
    hR_pol = zeros(Nb, n_tenure, I, J, n_parity, n_child);

    % tenure_choice now has ten_orig dimension to allow resizing while staying
    % Dimensions: (Nb, n_tenure, I, J, n_parity, n_child) -> optimal ten' given (b, ten_orig, i', ...)
    tenure_choice = zeros(Nb, n_tenure, I, J, n_parity, n_child);

    loc_probs = zeros(Nb, n_tenure, I, I, J, n_parity, n_child);
    fert_probs = zeros(Nb, n_tenure, I, J, n_parity);
    fert_value = zeros(Nb, n_tenure, I, J);

    % Pre-compute house costs for each location and tenure
    house_costs = zeros(I, n_tenure);
    down_payments = zeros(I, n_tenure, n_parity);
    b_min_owner = zeros(I, n_tenure, n_parity);
    house_equity = zeros(I, n_tenure);

    for i = 1:I
        for ten = 2:n_tenure
            h_size = P.H_own(ten - 1);
            house_costs(i, ten) = p_hat(i) * h_size;
            house_equity(i, ten) = (1 - P.psi) * p_hat(i) * h_size;
            for nn = 1:n_parity
                phi_n = P.phi(nn);
                down_payments(i, ten, nn) = (1 - phi_n) * house_costs(i, ten);
                b_min_owner(i, ten, nn) = -phi_n * house_costs(i, ten);
            end
        end
    end

    % Pre-compute segment solver parameters for each (loc, tenure)
    seg_params = cell(n_tenure, I);
    for i = 1:I
        for ten = 1:n_tenure
            sp = struct();
            sp.r_hat = r_hat(i);
            sp.p_hat = p_hat(i);
            sp.y = zeros(J, 1);
            for j = 1:J
                sp.y(j) = P.income(i, j);
            end
            if ten == 1
                sp.b_constraint = 0;
                sp.h_size = 0;
                sp.is_owner = false;
                sp.owner_cost = 0;
                sp.h_services = 0;
            else
                h_size = P.H_own(ten - 1);
                sp.h_size = h_size;
                sp.is_owner = true;
                sp.owner_cost = (P.delta + P.tau_H) * p_hat(i) * h_size;
                sp.h_services = P.chi * h_size;
            end
            seg_params{ten, i} = sp;
        end
    end

    b_min_grid = b_grid(1);
    b_max_grid = b_grid(end);

    % Terminal condition at age J (bequest utility)
    for i = 1:I
        for ten = 1:n_tenure
            if ten == 1
                house_value = 0;
            else
                house_value = house_equity(i, ten);
            end
            for nn = 1:n_parity
                n_children = nn - 1;
                for cs = 1:n_child
                    V(:, ten, i, J, nn, cs) = bequest_utility_vec(b_grid + house_value, n_children, P);
                end
            end
        end
    end

    % Backward induction
    for j = J-1:-1:1
        in_fert_window = (j >= P.A_f_start) && (j <= P.A_f_end);

        % ================================================================
        % STEP 1: DRIFT - Solve within-period HJB for ALL states
        % This uses V(j+1) as continuation and produces V_drift at age j
        % ================================================================
        V_drift = zeros(Nb, n_tenure, I, n_parity, n_child);
        c_drift = zeros(Nb, n_tenure, I, n_parity, n_child);
        hR_drift = zeros(Nb, n_tenure, I, n_parity, n_child);

        for nn = 1:n_parity
            n_children = nn - 1;
            phi_n = P.phi(nn);

            for cs = 1:n_child
                kids_present = (cs > 1);

                % kappa_h(n) from draft
                kappa_h_eff = P.kappa_h_base + P.kappa_h_slope * n_children;

                % STONE-GEARY: Subsistence levels depend on children present
                if kids_present
                    c_bar = P.c_bar_0 + P.c_bar_n * n_children;
                    h_bar = P.h_bar_0 + P.h_bar_jump + P.h_bar_n * n_children;
                    psi_bonus = P.psi_child * n_children;
                else
                    c_bar = P.c_bar_0;
                    h_bar = P.h_bar_0;
                    psi_bonus = 0;
                end

                % Child state transition
                if cs == 1
                    cs_next = 1;
                elseif cs < n_child
                    cs_next = cs + 1;
                else
                    cs_next = 1;
                end

                % Continuation value at j+1
                V_cont = V(:, :, :, j+1, nn, cs_next);

                % Solve HJB segment (drift) with V(j+1) as boundary
                [V_seg, c_seg, hR_seg] = solve_hjb_segment_stone_geary(...
                    V_cont, seg_params, j, kappa_h_eff, phi_n, P, b_grid, ...
                    db_f, db_b, inv_db_f, inv_db_b, Delta, psi_bonus, c_bar, h_bar);

                V_drift(:, :, :, nn, cs) = V_seg;
                c_drift(:, :, :, nn, cs) = c_seg;
                hR_drift(:, :, :, nn, cs) = hR_seg;

                % Store policies
                c_pol(:, :, :, j, nn, cs) = c_seg;
                hR_pol(:, :, :, j, nn, cs) = hR_seg;
            end
        end

        % ================================================================
        % STEP 2: TENURE CHOICE - V_H = max over h' of V_drift
        % Now includes resizing while staying (ten_orig matters)
        % Per notes eq (127): V_H(b,h,i',a,n) = max_{h'} V(b_after, h', i', a, n)
        % ================================================================
        V_H = zeros(Nb, n_tenure, I, n_parity, n_child);  % V_H(b, ten_orig, i_dest, nn, cs)
        tenure_choice_j = zeros(Nb, n_tenure, I, n_parity, n_child);

        for nn = 1:n_parity
            phi_n = P.phi(nn);

            for cs = 1:n_child
                % Create interpolants for V_drift at this (nn, cs)
                F_drift = cell(n_tenure, I);
                for i = 1:I
                    for ten = 1:n_tenure
                        F_drift{ten, i} = griddedInterpolant(b_grid, V_drift(:, ten, i, nn, cs), 'linear', 'linear');
                    end
                end

                for i_dest = 1:I
                    for ten_orig = 1:n_tenure
                        % Arriving wealth depends on whether we're staying or moving
                        % (handled in location stage, here we assume we're AT i_dest)
                        % For stayers (i_dest = i_orig): b_arrival = b
                        % For movers: b_arrival = b + house_equity(i_orig, ten_orig)
                        % But tenure choice happens AFTER location is decided
                        % So here b is already the post-location-move wealth

                        % Current house equity if selling
                        if ten_orig > 1
                            sell_proceeds = house_equity(i_dest, ten_orig);  % (1-psi)*p*h
                        else
                            sell_proceeds = 0;
                        end

                        % Value of each tenure option
                        V_tenure_options = zeros(Nb, n_tenure);

                        % Option 1: Rent (ten_new = 1)
                        if ten_orig == 1
                            % Already renting, no transaction
                            V_tenure_options(:, 1) = V_drift(:, 1, i_dest, nn, cs);
                        else
                            % Sell house, become renter
                            % IMPORTANT: Renters must have b >= 0 (borrowing constraint)
                            b_after_sell = max(b_grid + sell_proceeds, 0);
                            b_query = max(min(b_after_sell, b_max_grid), b_min_grid);
                            V_tenure_options(:, 1) = F_drift{1, i_dest}(b_query);
                        end

                        % Options 2+: Own (ten_new = 2, 3, ...)
                        for ten_new = 2:n_tenure
                            hc_new = house_costs(i_dest, ten_new);
                            dp_new = down_payments(i_dest, ten_new, nn);
                            bmo_new = b_min_owner(i_dest, ten_new, nn);

                            if ten_orig == ten_new
                                % Keep same house - no transaction
                                V_tenure_options(:, ten_new) = V_drift(:, ten_new, i_dest, nn, cs);
                            elseif ten_orig == 1
                                % Buy from renting
                                b_after_buy = b_grid - hc_new;
                                feasible = (b_grid >= dp_new) & (b_after_buy >= bmo_new);
                                b_query = max(min(b_after_buy, b_max_grid), b_min_grid);
                                V_own = F_drift{ten_new, i_dest}(b_query);
                                V_tenure_options(:, ten_new) = -1e10;
                                V_tenure_options(feasible, ten_new) = V_own(feasible);
                            else
                                % Resize: sell old, buy new
                                b_after_sell = b_grid + sell_proceeds;
                                b_after_resize = b_after_sell - hc_new;
                                dp_check = dp_new - sell_proceeds;  % net down payment needed
                                feasible = (b_grid >= dp_check) & (b_after_resize >= bmo_new);
                                b_query = max(min(b_after_resize, b_max_grid), b_min_grid);
                                V_own = F_drift{ten_new, i_dest}(b_query);
                                V_tenure_options(:, ten_new) = -1e10;
                                V_tenure_options(feasible, ten_new) = V_own(feasible);
                            end
                        end

                        % Optimal tenure choice
                        [V_H_opt, ten_opt] = max(V_tenure_options, [], 2);
                        V_H(:, ten_orig, i_dest, nn, cs) = V_H_opt;
                        tenure_choice_j(:, ten_orig, i_dest, nn, cs) = ten_opt;
                    end
                end
            end
        end

        % Store tenure choice (now with ten_orig dimension)
        tenure_choice(:, :, :, j, :, :) = tenure_choice_j;

        % ================================================================
        % STEP 3: LOCATION CHOICE with LOGIT - V_I
        % π^{i'} = (E_{i'} * μ_{ii'}) * exp(ν_ℓ * V^{i'}) / Σ
        % V^I = (1/ν_ℓ) * log(Σ exp(ν_ℓ * V^{i'}))
        % ================================================================
        V_I = zeros(Nb, n_tenure, I, n_parity, n_child);
        loc_probs_j = zeros(Nb, n_tenure, I, I, n_parity, n_child);

        mu_move_j = P.mu_move;

        for nn = 1:n_parity
            for cs = 1:n_child
                % Create interpolants for V_H at this (nn, cs)
                % V_H is indexed by (b, ten_orig, i_dest, nn, cs)
                F_H = cell(n_tenure, I);  % F_H{ten_orig, i_dest}
                for ten_orig = 1:n_tenure
                    for i_dest = 1:I
                        F_H{ten_orig, i_dest} = griddedInterpolant(b_grid, V_H(:, ten_orig, i_dest, nn, cs), 'linear', 'linear');
                    end
                end

                for i_orig = 1:I
                    for ten_orig = 1:n_tenure
                        house_value_orig = house_equity(i_orig, ten_orig);

                        % Compute log(E*mu) + nu_loc * V^H for each destination
                        log_Q_all = zeros(Nb, I);
                        for i_dest = 1:I
                            if i_dest == i_orig
                                % Stay: no moving cost, keep house
                                mu = P.mu_stay;
                                b_arrival = b_grid;
                                % V_H with same ten_orig at same location
                                V_H_dest = V_H(:, ten_orig, i_dest, nn, cs);
                            else
                                % Move: sell house, arrive as renter
                                mu = mu_move_j;
                                b_arrival = b_grid + house_value_orig;
                                % Movers arrive as renters (ten_orig=1 for V_H lookup)
                                b_query = max(min(b_arrival, b_max_grid), b_min_grid);
                                V_H_dest = F_H{1, i_dest}(b_query);  % ten_orig=1 (renter)
                            end

                            % LOGIT kernel
                            log_Q_all(:, i_dest) = log(P.E_loc(i_dest) * mu) + P.nu_loc * V_H_dest;
                        end

                        % Log-sum-exp for numerical stability
                        log_Q_max = max(log_Q_all, [], 2);
                        sum_exp = sum(exp(log_Q_all - log_Q_max), 2);
                        log_sum_Q = log_Q_max + log(sum_exp);

                        % Choice probabilities (softmax)
                        probs = exp(log_Q_all - log_sum_Q);

                        % Inclusive value
                        V_inclusive = log_sum_Q / P.nu_loc;

                        V_I(:, ten_orig, i_orig, nn, cs) = V_inclusive;
                        loc_probs_j(:, ten_orig, i_orig, :, nn, cs) = probs;
                    end
                end
            end
        end

        loc_probs(:, :, :, :, j, :, :) = loc_probs_j;

        % ================================================================
        % STEP 4: FERTILITY CHOICE with LOGIT (if in fertile window)
        % π^{n'} = exp(ν^n * V^I(n')) / Σ_m exp(ν^n * V^I(m))
        % V^fert = (1/ν^n) * log(Σ_m exp(ν^n * V^I(m)))
        % ================================================================
        if in_fert_window
            for i = 1:I
                for ten = 1:n_tenure
                    % Collect V_I for each parity option
                    V_options = zeros(Nb, n_parity);
                    for nn_opt = 1:n_parity
                        n_kids = nn_opt - 1;
                        if n_kids == 0
                            % Childless: cs=1
                            V_options(:, nn_opt) = V_I(:, ten, i, 1, 1);
                        else
                            % With children: cs=2 (kids present)
                            V_options(:, nn_opt) = V_I(:, ten, i, nn_opt, 2);
                        end
                    end

                    % LOGIT: exp(nu_fert * V) formulation
                    log_V_scaled = P.nu_fert * V_options;

                    % Log-sum-exp for numerical stability
                    log_V_max = max(log_V_scaled, [], 2);
                    sum_exp = sum(exp(log_V_scaled - log_V_max), 2);
                    log_sum_V = log_V_max + log(sum_exp);

                    % Choice probabilities (softmax)
                    fert_probs(:, ten, i, j, :) = exp(log_V_scaled - log_sum_V);

                    % Inclusive value (post-fertility value)
                    fert_value(:, ten, i, j) = log_sum_V / P.nu_fert;
                end
            end

            % For childless agents at fertile ages, V = V^fert
            % This is the value BEFORE the fertility shock is realized
            for ten = 1:n_tenure
                for i = 1:I
                    V(:, ten, i, j, 1, 1) = fert_value(:, ten, i, j);
                end
            end

            % For agents who already have children, V = V_I (no fertility choice)
            for nn = 2:n_parity
                for cs = 1:n_child
                    V(:, :, :, j, nn, cs) = V_I(:, :, :, nn, cs);
                end
            end
            % Also for childless agents in non-kid state
            for cs = 2:n_child
                V(:, :, :, j, 1, cs) = V_I(:, :, :, 1, cs);
            end
        else
            % Not in fertile window: V = V_I (location is the last shock)
            for nn = 1:n_parity
                for cs = 1:n_child
                    V(:, :, :, j, nn, cs) = V_I(:, :, :, nn, cs);
                end
            end
        end
    end
end


%% ========================================================================
%%                    STONE-GEARY HJB SEGMENT SOLVER
%% ========================================================================
function [V_all, c_all, hR_all] = solve_hjb_segment_stone_geary(...
    V_cont, seg_params, j, kappa_h_eff, phi_n, P, b_grid, ...
    db_f, db_b, inv_db_f, inv_db_b, Delta, psi_bonus, c_bar, h_bar)
    % STONE-GEARY VERSION:
    % c* = c_bar + V_b^(-1/sigma)
    % h* = h_bar + (kappa_h / (r * V_b))^(1/sigma)

    Nb = length(b_grid);
    I = size(V_cont, 3);
    n_tenure = size(V_cont, 2);
    n_sub = P.n_sub;
    dt = Delta / n_sub;

    V_all = V_cont;
    c_all = zeros(Nb, n_tenure, I);
    hR_all = zeros(Nb, n_tenure, I);

    % Pre-extract parameters
    r_hat_vec = zeros(I, 1);
    p_hat_vec = zeros(I, 1);
    y_vec = zeros(I, 1);
    for i = 1:I
        r_hat_vec(i) = seg_params{1, i}.r_hat;
        p_hat_vec(i) = seg_params{1, i}.p_hat;
        y_vec(i) = seg_params{1, i}.y(j);
    end

    b_constraint = zeros(n_tenure, I);
    owner_cost = zeros(n_tenure, I);
    h_services_owner = zeros(n_tenure, I);
    is_owner = false(n_tenure, I);

    for i = 1:I
        for ten = 1:n_tenure
            sp = seg_params{ten, i};
            if sp.is_owner
                is_owner(ten, i) = true;
                b_constraint(ten, i) = -phi_n * sp.p_hat * sp.h_size;
                owner_cost(ten, i) = sp.owner_cost;
                h_services_owner(ten, i) = sp.h_services;
            else
                is_owner(ten, i) = false;
                b_constraint(ten, i) = 0;
                owner_cost(ten, i) = 0;
                h_services_owner(ten, i) = 0;
            end
        end
    end

    r_hat_3d = reshape(r_hat_vec, 1, 1, I);
    y_3d = reshape(y_vec, 1, 1, I);
    b_constraint_3d = reshape(b_constraint, 1, n_tenure, I);
    owner_cost_3d = reshape(owner_cost, 1, n_tenure, I);
    h_services_owner_3d = reshape(h_services_owner, 1, n_tenure, I);

    b_grid_3d = b_grid;

    for sub = 1:n_sub
        dV_f = zeros(Nb, n_tenure, I);
        dV_b = zeros(Nb, n_tenure, I);

        dV_f(1:Nb-1, :, :) = (V_all(2:Nb, :, :) - V_all(1:Nb-1, :, :)) .* inv_db_f(1:Nb-1);
        dV_f(Nb, :, :) = P.dV_floor;

        dV_b(2:Nb, :, :) = (V_all(2:Nb, :, :) - V_all(1:Nb-1, :, :)) .* inv_db_b(2:Nb);
        dV_b(1, :, :) = P.dV_floor;

        dV_f = max(dV_f, P.dV_floor);
        dV_b = max(dV_b, P.dV_floor);

        dV = max(dV_f, dV_b);

        % ================================================================
        % STONE-GEARY OPTIMAL POLICIES
        % c* = c_bar + V_b^(-1/sigma)
        % h* = h_bar + (kappa_h / (r * V_b))^(1/sigma)
        % ================================================================

        % Supernumerary consumption
        c_tilde_opt = dV.^(-1/P.sigma);
        c_opt = c_bar + c_tilde_opt;
        c_opt = max(c_opt, c_bar + P.c_min);

        % Supernumerary rental housing (only for renters)
        h_tilde_opt = (kappa_h_eff ./ (r_hat_3d .* dV)).^(1/P.sigma);
        hR_opt = h_bar + h_tilde_opt;
        hR_opt = max(hR_opt, h_bar + 0.01);
        hR_opt = min(hR_opt, P.hR_max);

        h_services = zeros(Nb, n_tenure, I);
        flow_cost = zeros(Nb, n_tenure, I);

        % Renters: housing services = hR_opt
        h_services(:, 1, :) = hR_opt(:, 1, :);
        for i = 1:I
            flow_cost(:, 1, i) = r_hat_vec(i) * hR_opt(:, 1, i);
        end

        % Owners: housing services = chi * H_own (fixed)
        for ten = 2:n_tenure
            for i = 1:I
                h_services(:, ten, i) = h_services_owner(ten, i);
                flow_cost(:, ten, i) = owner_cost(ten, i);
            end
            hR_opt(:, ten, :) = 0;
        end

        % Resource constraint: ensure c + flow_cost <= resources
        resources = y_3d + (P.q - P.gamma) * b_grid_3d;
        % Maximum consumption given resources and required housing
        c_max = resources - flow_cost;
        c_max = max(c_max, c_bar + P.c_min);
        c_opt = min(c_opt, c_max);

        % Handle infeasibility for owners
        for ten = 2:n_tenure
            for i = 1:I
                infeasible = b_grid < b_constraint(ten, i) - 1e-6;
                c_opt(infeasible, ten, i) = c_bar + P.c_min;
            end
        end

        % Drift (same as before, but with new c_opt)
        mu = y_3d + (P.q - P.gamma) * b_grid_3d - c_opt - flow_cost;

        for ten = 1:n_tenure
            for i = 1:I
                at_constraint = (b_grid <= b_constraint(ten, i) + 1e-6) & (mu(:, ten, i) < 0);
                mu(at_constraint, ten, i) = 0;
            end
        end

        x = max(mu, 0) .* inv_db_f;
        y_drift = max(-mu, 0) .* inv_db_b;

        % ================================================================
        % STONE-GEARY UTILITY
        % u = (c - c_bar)^(1-σ)/(1-σ) + κ_h (h - h_bar)^(1-σ)/(1-σ) + u_bar + ψ
        % ================================================================
        u_flow = utility_flow_stone_geary(c_opt, h_services, kappa_h_eff, P, psi_bonus, c_bar, h_bar);

        % Thomas algorithm solve
        n_systems = n_tenure * I;
        x_flat = reshape(x, Nb, n_systems);
        y_flat = reshape(y_drift, Nb, n_systems);
        u_flat = reshape(u_flow, Nb, n_systems);
        V_flat = reshape(V_all, Nb, n_systems);

        b_diag = (1/dt + P.rho_hat) + x_flat + y_flat;
        c_diag = -x_flat(1:Nb-1, :);
        a_diag = -y_flat(2:Nb, :);
        rhs = u_flat + V_flat / dt;

        V_flat = thomas_solve_vec(a_diag, b_diag, c_diag, rhs);

        V_all = reshape(V_flat, Nb, n_tenure, I);

        c_all = c_opt;
        hR_all = hR_opt;
    end
end


%% ========================================================================
%%                    STONE-GEARY UTILITY FUNCTION
%% ========================================================================
function u = utility_flow_stone_geary(c, h, kappa_h, P, psi_bonus, c_bar, h_bar)
    % STONE-GEARY UTILITY:
    % u = (c - c_bar)^(1-σ)/(1-σ) + κ_h (h - h_bar)^(1-σ)/(1-σ) + u_bar + ψ

    % Supernumerary consumption and housing
    c_tilde = max(c - c_bar, 1e-10);
    h_tilde = max(h - h_bar, 1e-10);

    if abs(P.sigma - 1) < 1e-6
        u = log(c_tilde) + kappa_h * log(h_tilde) + P.u_bar + psi_bonus;
    else
        u = c_tilde.^(1-P.sigma)/(1-P.sigma) + kappa_h * h_tilde.^(1-P.sigma)/(1-P.sigma) + P.u_bar + psi_bonus;
    end
end


%% ========================================================================
%%                    VECTORIZED THOMAS ALGORITHM
%% ========================================================================
function x = thomas_solve_vec(a, b, c, d)
    [n, m] = size(b);

    c_star = zeros(n-1, m);
    d_star = zeros(n, m);

    c_star(1, :) = c(1, :) ./ b(1, :);
    d_star(1, :) = d(1, :) ./ b(1, :);

    for i = 2:n-1
        denom = b(i, :) - a(i-1, :) .* c_star(i-1, :);
        c_star(i, :) = c(i, :) ./ denom;
        d_star(i, :) = (d(i, :) - a(i-1, :) .* d_star(i-1, :)) ./ denom;
    end

    denom = b(n, :) - a(n-1, :) .* c_star(n-1, :);
    d_star(n, :) = (d(n, :) - a(n-1, :) .* d_star(n-1, :)) ./ denom;

    x = zeros(n, m);
    x(n, :) = d_star(n, :);
    for i = n-1:-1:1
        x(i, :) = d_star(i, :) - c_star(i, :) .* x(i+1, :);
    end
end


%% ========================================================================
%%                    KFE SOLVER
%%                    CORRECT TIMING: FERTILITY -> LOCATION -> TENURE -> DRIFT -> ADVANCE
%% ========================================================================
function [g, stats] = solve_kfe(c_pol, hR_pol, tenure_choice, loc_probs, fert_probs, ...
                                r_hat, p_hat, P, b_grid, db_f, db_b, idx_zero, n_pop)

    % TIMING (for forward iteration at age j):
    %   1. FERTILITY: Apply fertility shock (if in fertile window)
    %   2. LOCATION: Apply location transition using loc_probs
    %   3. TENURE: Apply tenure choice (with resizing while staying)
    %   4. DRIFT: Evolve wealth distribution using KFE
    %   5. Apply demographic factor exp(-n_pop)
    %   6. ADVANCE: Move to age j+1

    J = P.J;
    I = P.I;
    Nb = length(b_grid);
    n_h = P.n_house;
    n_tenure = 1 + n_h;
    n_parity = P.n_parity;
    n_child = P.n_child_states;
    Delta = P.da;

    b_min_grid = b_grid(1);
    b_max_grid = b_grid(end);

    g = zeros(Nb, n_tenure, I, J, n_parity, n_child);

    % Entry: age 1, renters, n=0, no kids
    for i = 1:I
        idx_entry_i = find(b_grid >= P.b_entry_loc(i), 1, 'first');
        if isempty(idx_entry_i), idx_entry_i = Nb; end
        g(idx_entry_i, 1, i, 1, 1, 1) = P.entry_by_loc(i);
    end

    demographic_factor = exp(-n_pop);

    total_births_kfe = 0;
    births_by_loc = zeros(I, 1);

    % MATURITY FLOW ACCUMULATORS
    % Track where children mature (reach adulthood) after parents have moved
    % This is the correct object for entry shares: "kids follow parents"
    entrants_mature_by_loc = zeros(I, 1);
    entrants_mature_total = 0;

    % Pre-compute house costs and equity
    house_costs = zeros(I, n_tenure);
    house_equity = zeros(I, n_tenure);
    for i = 1:I
        for ten = 2:n_tenure
            h_size = P.H_own(ten - 1);
            house_costs(i, ten) = p_hat(i) * h_size;
            house_equity(i, ten) = (1 - P.psi) * p_hat(i) * h_size;
        end
    end

    diag_idx = (1:Nb)';
    upper_rows = (2:Nb)';
    upper_cols = (1:Nb-1)';
    lower_rows = (1:Nb-1)';
    lower_cols = (2:Nb)';

    inv_db_f = 1 ./ db_f;
    inv_db_b = 1 ./ db_b;

    for j = 1:J-1
        in_fert_window = (j >= P.A_f_start) && (j <= P.A_f_end);

        % ================================================================
        % STEP 1: FERTILITY - Apply fertility shock (if in fertile window)
        % ================================================================
        if in_fert_window
            g_childless = g(:, :, :, j, 1, 1);

            % Count births for BGP calculation
            for i = 1:I
                for ten = 1:n_tenure
                    g_slice = g_childless(:, ten, i);
                    mass_total = sum(g_slice);
                    if mass_total < 1e-15, continue; end

                    probs_mat = squeeze(fert_probs(:, ten, i, j, :));
                    parity_vec = (0:n_parity-1)';
                    expected_parity = probs_mat * parity_vec;
                    births_this = sum(g_slice .* expected_parity);
                    total_births_kfe = total_births_kfe + births_this;
                    births_by_loc(i) = births_by_loc(i) + births_this;
                end
            end

            % Redistribute mass according to fertility probabilities
            g_post_fert = zeros(Nb, n_tenure, I, n_parity, n_child);

            for i = 1:I
                for ten = 1:n_tenure
                    g_slice = g_childless(:, ten, i);
                    if sum(g_slice) < 1e-15, continue; end

                    probs_fert = squeeze(fert_probs(:, ten, i, j, :));

                    for nn = 1:n_parity
                        n_kids = nn - 1;
                        cs_new = (n_kids == 0) * 1 + (n_kids > 0) * 2;

                        mass_nn = g_slice .* probs_fert(:, nn);
                        g_post_fert(:, ten, i, nn, cs_new) = ...
                            g_post_fert(:, ten, i, nn, cs_new) + mass_nn;
                    end
                end
            end

            % Update g with post-fertility distribution
            g(:, :, :, j, 1, 1) = 0;
            for nn = 1:n_parity
                for cs = 1:n_child
                    g(:, :, :, j, nn, cs) = g(:, :, :, j, nn, cs) + g_post_fert(:, :, :, nn, cs);
                end
            end
        end

        % ================================================================
        % Process each (parity, child_state) combination
        % ================================================================
        for nn = 1:n_parity
            phi_n = P.phi(nn);

            for cs = 1:n_child
                % Child state transition
                if cs == 1
                    cs_next = 1;
                elseif cs < n_child
                    cs_next = cs + 1;
                else
                    cs_next = 1;
                end

                % ============================================================
                % STEP 2: LOCATION - Apply location transition
                % ============================================================
                g_post_loc = zeros(Nb, n_tenure, I);

                for i_orig = 1:I
                    for ten_orig = 1:n_tenure
                        g_slice = g(:, ten_orig, i_orig, j, nn, cs);
                        nz_idx = find(g_slice > 1e-15);
                        if isempty(nz_idx), continue; end

                        house_value_orig = house_equity(i_orig, ten_orig);
                        probs = squeeze(loc_probs(:, ten_orig, i_orig, :, j, nn, cs));

                        for k = 1:length(nz_idx)
                            ib = nz_idx(k);
                            mass = g_slice(ib);
                            b = b_grid(ib);

                            for i_dest = 1:I
                                mass_dest = mass * probs(ib, i_dest);
                                if mass_dest < 1e-15, continue; end

                                if i_dest ~= i_orig
                                    % Movers: sell house, arrive as renter
                                    % IMPORTANT: Renters must have b >= 0 (borrowing constraint)
                                    b_after_loc = max(b + house_value_orig, 0);
                                    arriving_tenure = 1;  % Movers arrive as renters
                                else
                                    % Stayers: keep current tenure
                                    b_after_loc = b;
                                    arriving_tenure = ten_orig;
                                end

                                % Place mass at arriving wealth and tenure
                                [~, ib_after] = min(abs(b_grid - b_after_loc));
                                g_post_loc(ib_after, arriving_tenure, i_dest) = ...
                                    g_post_loc(ib_after, arriving_tenure, i_dest) + mass_dest;
                            end
                        end
                    end
                end

                % ============================================================
                % STEP 3: TENURE - Apply tenure choice (with resizing)
                % tenure_choice has dimensions: (Nb, n_tenure, I, J, n_parity, n_child)
                % tenure_choice(ib, ten_orig, i_dest, j, nn, cs) = optimal ten'
                % ============================================================
                g_post_tenure = zeros(Nb, n_tenure, I);

                for i_dest = 1:I
                    for ten_orig = 1:n_tenure
                        g_slice = g_post_loc(:, ten_orig, i_dest);
                        nz_idx = find(g_slice > 1e-15);
                        if isempty(nz_idx), continue; end

                        % Current house equity if selling
                        if ten_orig > 1
                            sell_proceeds = house_equity(i_dest, ten_orig);
                        else
                            sell_proceeds = 0;
                        end

                        for k = 1:length(nz_idx)
                            ib = nz_idx(k);
                            mass = g_slice(ib);
                            b = b_grid(ib);

                            % Get optimal tenure choice
                            ten_new = tenure_choice(ib, ten_orig, i_dest, j, nn, cs);

                            % Compute wealth after tenure transition
                            if ten_new == ten_orig
                                % No change
                                b_final = b;
                            elseif ten_orig == 1 && ten_new > 1
                                % Buy from renting
                                hc_new = house_costs(i_dest, ten_new);
                                dp_new = (1 - phi_n) * hc_new;
                                bmo_new = -phi_n * hc_new;

                                if b >= dp_new
                                    b_final = max(b - hc_new, bmo_new);
                                else
                                    % Can't afford, stay renter
                                    ten_new = 1;
                                    b_final = max(b, 0);
                                end
                            elseif ten_orig > 1 && ten_new == 1
                                % Sell to rent
                                % IMPORTANT: Renters must have b >= 0 (borrowing constraint)
                                b_final = max(b + sell_proceeds, 0);
                            else
                                % Resize: sell old, buy new
                                hc_new = house_costs(i_dest, ten_new);
                                dp_new = (1 - phi_n) * hc_new;
                                bmo_new = -phi_n * hc_new;

                                b_after_sell = b + sell_proceeds;
                                b_after_resize = b_after_sell - hc_new;
                                dp_check = dp_new - sell_proceeds;

                                if b >= dp_check && b_after_resize >= bmo_new
                                    b_final = b_after_resize;
                                else
                                    % Can't afford, keep original or become renter
                                    if ten_orig > 1
                                        ten_new = ten_orig;
                                        b_final = b;
                                    else
                                        ten_new = 1;
                                        b_final = max(b, 0);
                                    end
                                end
                            end

                            % Place mass
                            b_final = max(min(b_final, b_max_grid), b_min_grid);
                            [~, ib_final] = min(abs(b_grid - b_final));
                            g_post_tenure(ib_final, ten_new, i_dest) = ...
                                g_post_tenure(ib_final, ten_new, i_dest) + mass;
                        end
                    end
                end

                % ============================================================
                % STEP 4: DRIFT - Evolve wealth distribution using KFE
                % ============================================================
                g_post_drift = zeros(Nb, n_tenure, I);

                for i = 1:I
                    y_j = P.income(i, j);

                    for ten = 1:n_tenure
                        g_curr = g_post_tenure(:, ten, i);

                        if sum(g_curr) < 1e-15
                            g_post_drift(:, ten, i) = g_curr;
                            continue;
                        end

                        if ten == 1
                            c_j = c_pol(:, 1, i, j, nn, cs);
                            hR_j = hR_pol(:, 1, i, j, nn, cs);
                            mu = y_j + (P.q - P.gamma) * b_grid - c_j - r_hat(i) * hR_j;
                            b_constraint = 0;
                        else
                            h_size = P.H_own(ten - 1);
                            c_j = c_pol(:, ten, i, j, nn, cs);
                            owner_cost = (P.delta + P.tau_H) * p_hat(i) * h_size;
                            mu = y_j + (P.q - P.gamma) * b_grid - c_j - owner_cost;
                            b_constraint = -phi_n * p_hat(i) * h_size;
                        end

                        at_constraint = (b_grid <= b_constraint + 1e-6) & (mu < 0);
                        mu(at_constraint) = 0;

                        x = max(mu, 0) .* inv_db_f;
                        y_kfe = max(-mu, 0) .* inv_db_b;
                        x(Nb) = 0;
                        y_kfe(1) = 0;

                        diag_vals = -(x + y_kfe);

                        A_kfe = sparse([diag_idx; upper_rows; lower_rows], ...
                                      [diag_idx; upper_cols; lower_cols], ...
                                      [diag_vals; x(1:end-1); y_kfe(2:end)], Nb, Nb);

                        g_post_drift(:, ten, i) = max((speye(Nb) - Delta * A_kfe) \ g_curr, 0);
                    end
                end

                % ============================================================
                % COUNT MATURITY EVENTS
                % Maturity occurs when cs_next == 1 and cs > 1 (kids finish growing)
                % Only for parents with children (nn >= 2)
                % g_post_drift contains parents' distribution AFTER they've moved
                % This is the correct location for "kids follow parents"
                % ============================================================
                if (cs_next == 1) && (cs > 1) && (nn >= 2)
                    n_kids = nn - 1;  % Number of children maturing

                    for i_m = 1:I
                        % Flow of maturing children = n_kids * parent mass in location i
                        flow_i = n_kids * sum(g_post_drift(:, :, i_m), 'all');
                        entrants_mature_by_loc(i_m) = entrants_mature_by_loc(i_m) + flow_i;
                        entrants_mature_total = entrants_mature_total + flow_i;
                    end
                end

                % ============================================================
                % STEP 5 & 6: Apply demographic factor and advance to j+1
                % ============================================================
                g(:, :, :, j+1, nn, cs_next) = g(:, :, :, j+1, nn, cs_next) + demographic_factor * g_post_drift;
            end
        end
    end

    % Renormalize
    total_mass = sum(g, 'all');
    if total_mass > 1e-12
        g = g / total_mass;
        total_births_kfe = total_births_kfe / total_mass;
        births_by_loc = births_by_loc / total_mass;
    end

    entry_rate = sum(g(:, :, :, 1, :, :), 'all');

    stats = compute_statistics(g, fert_probs, P, b_grid, p_hat);
    stats.total_births_kfe = total_births_kfe;
    stats.births_by_loc = births_by_loc;
    stats.entry_rate = entry_rate;

    % Maturity flow statistics (for "kids follow parents" entry allocation)
    % Normalize by total mass (same as births)
    if total_mass > 1e-12
        entrants_mature_by_loc = entrants_mature_by_loc / total_mass;
        entrants_mature_total = entrants_mature_total / total_mass;
    end

    stats.entrants_mature_by_loc = entrants_mature_by_loc;
    stats.entrants_mature_total = entrants_mature_total;
    stats.mature_entry_shares = entrants_mature_by_loc / max(entrants_mature_total, 1e-12);
end


%% ========================================================================
%%                  COMPUTE STATISTICS
%% ========================================================================
function stats = compute_statistics(g, fert_probs, P, b_grid, p_hat)

    J = P.J;
    I = P.I;
    Nb = length(b_grid);
    n_tenure = 1 + P.n_house;
    n_parity = P.n_parity;
    n_child = P.n_child_states;

    total_mass = sum(g, 'all');

    own_mass = sum(g(:, 2:end, :, :, :, :), 'all');
    stats.own_rate = own_mass / max(total_mass, 1e-12);

    stats.pop_share = zeros(I, 1);
    stats.own_by_loc = zeros(I, 1);
    for i = 1:I
        pop_i = sum(g(:, :, i, :, :, :), 'all');
        stats.pop_share(i) = pop_i / max(total_mass, 1e-12);
        own_i = sum(g(:, 2:end, i, :, :, :), 'all');
        stats.own_by_loc(i) = own_i / max(pop_i, 1e-12);
    end

    stats.parity_dist = zeros(n_parity, 1);
    mass_post_fert = sum(g(:, :, :, P.A_f_end+1:end, :, :), 'all');
    for nn = 1:n_parity
        mass_nn = sum(g(:, :, :, P.A_f_end+1:end, nn, :), 'all');
        stats.parity_dist(nn) = mass_nn / max(mass_post_fert, 1e-12);
    end
    stats.mean_parity = sum((0:n_parity-1)' .* stats.parity_dist);

    stats.own_by_parity = zeros(n_parity, 1);
    for nn = 1:n_parity
        mass_nn = sum(g(:, :, :, :, nn, :), 'all');
        own_nn = sum(g(:, 2:end, :, :, nn, :), 'all');
        stats.own_by_parity(nn) = own_nn / max(mass_nn, 1e-12);
    end

    stats.mean_parity_by_loc = zeros(I, 1);
    stats.frac_childless_by_loc = zeros(I, 1);
    for i = 1:I
        mass_i_post = sum(g(:, :, i, P.A_f_end+1:end, :, :), 'all');
        if mass_i_post > 1e-12
            mean_n = 0;
            for nn = 1:n_parity
                mass_inn = sum(g(:, :, i, P.A_f_end+1:end, nn, :), 'all');
                mean_n = mean_n + (nn-1) * mass_inn / mass_i_post;
            end
            stats.mean_parity_by_loc(i) = mean_n;
            mass_childless = sum(g(:, :, i, P.A_f_end+1:end, 1, :), 'all');
            stats.frac_childless_by_loc(i) = mass_childless / mass_i_post;
        end
    end

    stats.own_by_age = zeros(J, 1);
    for jj = 1:J
        g_j = g(:, :, :, jj, :, :);
        mass_j = sum(g_j, 'all');
        if mass_j > 1e-12
            stats.own_by_age(jj) = sum(g_j(:, 2:end, :, :, :, :), 'all') / mass_j;
        end
    end

    stats.child_state_dist = zeros(J, n_child);
    for jj = 1:J
        mass_j = sum(g(:, :, :, jj, :, :), 'all');
        if mass_j > 1e-12
            for cs = 1:n_child
                stats.child_state_dist(jj, cs) = sum(g(:, :, :, jj, :, cs), 'all') / mass_j;
            end
        end
    end

    stats.fert_by_age = zeros(J, 1);
    for j = P.A_f_start:P.A_f_end
        mass_j = sum(g(:, :, :, j, 1, 1), 'all');
        if mass_j > 1e-12
            E_n = 0;
            for i = 1:I
                for ten = 1:n_tenure
                    g_slice = g(:, ten, i, j, 1, 1);
                    nz_idx = find(g_slice > 1e-15);
                    if isempty(nz_idx), continue; end

                    probs = squeeze(fert_probs(:, ten, i, j, :));
                    parity_vec = (0:n_parity-1)';
                    E_n_vec = probs * parity_vec;
                    E_n = E_n + sum(g_slice(nz_idx) .* E_n_vec(nz_idx));
                end
            end
            stats.fert_by_age(j) = E_n / mass_j;
        end
    end

    stats.death_wealth_by_loc = zeros(I, 1);
    stats.entry_mass_by_loc = zeros(I, 1);
    stats.bequest_per_entrant = zeros(I, 1);

    for i = 1:I
        death_wealth_i = 0;
        for ten = 1:n_tenure
            if ten > 1
                H_size = P.H_own(ten - 1);
                house_eq = (1 - P.psi) * p_hat(i) * H_size;
            else
                house_eq = 0;
            end
            for nn = 1:n_parity
                for cs = 1:n_child
                    g_slice = g(:, ten, i, J, nn, cs);
                    death_wealth_i = death_wealth_i + sum(g_slice .* (b_grid + house_eq));
                end
            end
        end
        stats.death_wealth_by_loc(i) = death_wealth_i;
        stats.entry_mass_by_loc(i) = sum(g(:, :, i, 1, :, :), 'all');

        T_i = stats.entry_mass_by_loc(i);
        if T_i > 1e-12
            B_i = P.phi_estate * death_wealth_i;
            stats.bequest_per_entrant(i) = P.b_entry_base + B_i / T_i;
        else
            stats.bequest_per_entrant(i) = P.b_entry_base;
        end
    end
end


%% ========================================================================
%%                    VECTORIZED UTILITY FUNCTIONS
%% ========================================================================
function v = bequest_utility_vec(b, n_children, P)
    b = max(b, 0);
    if abs(P.sigma - 1) < 1e-6
        v = P.theta0 * log(P.theta1 + b);
    else
        v = P.theta0 * (P.theta1 + b).^(1-P.sigma) / (1-P.sigma);
    end
    v = v + P.theta_n * n_children;
end


%% ========================================================================
%%                    BGP CONSISTENCY VERIFICATION
%% ========================================================================
function verify_bgp_consistency(sol, P, n_pop)

    fprintf('\n=== BGP CONSISTENCY VERIFICATION ===\n');

    g = sol.g;
    J = P.J;
    A_m = P.A_m;

    fprintf('\n1. Age distribution (should decay at rate exp(-n) = %.4f):\n', exp(-n_pop));
    mass_by_age = zeros(J, 1);
    for jj = 1:J
        mass_by_age(jj) = sum(g(:, :, :, jj, :, :), 'all');
    end

    fprintf('   Age   Mass       Ratio(j/j-1)  Expected\n');
    for jj = 1:min(J, 6)
        if jj == 1
            ratio = 1.0;
        else
            ratio = mass_by_age(jj) / max(mass_by_age(jj-1), 1e-12);
        end
        fprintf('   %2d    %.4f     %.4f        %.4f\n', jj, mass_by_age(jj), ratio, exp(-n_pop));
    end

    birth_rate = sol.total_births_kfe;
    entry_rate = sol.entry_rate;
    death_rate = sum(g(:, :, :, J, :, :), 'all');

    fprintf('\n2. KEY BGP CHECK: E = B × exp(-n × A_m):\n');
    implied_entry = birth_rate * exp(-n_pop * A_m);
    fprintf('   Birth rate:     %.6f\n', birth_rate);
    fprintf('   Entry rate:     %.6f\n', entry_rate);
    fprintf('   E_implied:      %.6f  (= B × exp(-n × A_m))\n', implied_entry);
    fprintf('   Gap (E - E_impl): %.6f (%.2f%%)\n', ...
        entry_rate - implied_entry, 100*abs(entry_rate - implied_entry)/max(entry_rate, 1e-12));
    if abs(entry_rate - implied_entry) / max(entry_rate, 1e-12) < 0.05
        fprintf('   STATUS: *** BGP CONSISTENT ***\n');
    else
        fprintf('   STATUS: *** NOT AT BGP ***\n');
    end

    fprintf('\n3. Implied n from B/E:\n');
    n_from_BE = log(birth_rate / entry_rate) / A_m;
    fprintf('   n from log(B/E)/A_m: %.6f\n', n_from_BE);
    fprintf('   n used:              %.6f\n', n_pop);
    fprintf('   Error: %.6f\n', abs(n_from_BE - n_pop));

    fprintf('\n4. Fertility accounting:\n');
    fprintf('   Mean parity (completed): %.3f\n', sol.mean_parity);
    fprintf('   Parity distribution: [%.1f%%, %.1f%%, %.1f%%, %.1f%%]\n', ...
        100*sol.parity_dist(1), 100*sol.parity_dist(2), ...
        100*sol.parity_dist(3), 100*sol.parity_dist(4));
    fprintf('   Maturation lag A_m: %d periods\n', A_m);

    fprintf('\n5. Scale invariance check:\n');
    total_mass = sum(g, 'all');
    fprintf('   Total mass = %.6f (should be 1.0)\n', total_mass);

    fprintf('\n========================================\n');
end
