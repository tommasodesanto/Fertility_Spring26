%% test_annual_v2.m
%
% Test annual periods (da=1) with properly adjusted parameters.
%
% Key insight: With more periods, there are more fertility decision points.
% To get similar aggregate fertility, we need to adjust nu_fert.
%
% Timeline:
%   Ages 18-78 (60 years of adult life)
%   J = 60 periods
%   A_m = 18 (children mature at 18)
%   Fertile ages: 18-35 (17 years = 17 periods)
%
% ========================================================================

clear; clc;
addpath(fileparts(mfilename('fullpath')));

fprintf('=============================================================\n');
fprintf('  ANNUAL PERIODS TEST (da=1.0, J=60)\n');
fprintf('=============================================================\n\n');

%% Parameters for annual model
P = struct();

% Timeline
P.J = 60;           % 60 years of adult life (ages 18-78)
P.da = 1.0;         % Annual periods
P.J_R = 47;         % Retire at age 65 (65-18 = 47 periods)
P.age_start = 18;

% Fertility
P.A_m = 18;         % Children mature at 18
P.A_f_start = 1;    % Fertile from age 18
P.A_f_end = 17;     % Fertile until age 35 (17 periods)

% Wealth grid
P.Nb = 30;          % Keep same grid size for speed

% ========================================================================
% SCALING FROM BASE (da=6 equivalent to 5 periods per 30 years)
% Base model uses da=6, so da_ratio = 1/6 ≈ 0.167
% ========================================================================
da_base = 6.0;
da_ratio = P.da / da_base;

fprintf('Scaling from base (da=%.1f) to annual (da=%.1f)\n', da_base, P.da);
fprintf('da_ratio = %.4f\n\n', da_ratio);

% SUBSISTENCE PARAMETERS (scale with da_ratio)
% These are per-period requirements
P.c_bar_0 = 0.05 * da_ratio;     % 0.05 / 6 ≈ 0.0083
P.c_bar_n = 0.05 * da_ratio;
P.h_bar_0 = 0.03 * da_ratio;
P.h_bar_jump = 0.50 * da_ratio;
P.h_bar_n = 0.20 * da_ratio;
P.c_min = 0.05 * da_ratio;

% SHOCK PARAMETERS
% With more fertility decision points, need higher nu_fert to get similar
% aggregate fertility. Rough scaling: inversely with number of fertile periods.
% Base: 3 fertile periods with nu_fert=8
% Annual: 17 fertile periods
% Scale factor for nu: (17/3) * da_ratio ≈ 5.67 * 0.167 ≈ 0.94
% But this is approximate - may need tuning

% Alternative approach: keep nu at reasonable levels that work
P.nu_fert = 3.0;    % Moderate fertility dispersion (tuned for annual)
P.nu_loc = 0.5;     % Location dispersion (scaled down)

% UTILITY OFFSET
% With logit shocks, only V differences matter - u_bar not needed!
P.u_bar = 0.0;      % Can be zero with Gumbel/logit

% TERMINAL VALUES (don't scale - these are lump-sum at death)
P.theta_n = 1.5;    % Keep same
P.theta0 = 1.5;
P.theta1 = 0.01;

% BEQUEST AND ENTRY
P.b_entry_base = 0.5;  % Initial wealth (stock, don't scale)
P.phi_estate = 0.80;

fprintf('Scaled subsistence parameters:\n');
fprintf('  c_bar_0 = %.4f, c_bar_n = %.4f\n', P.c_bar_0, P.c_bar_n);
fprintf('  h_bar_0 = %.4f, h_bar_jump = %.4f, h_bar_n = %.4f\n', ...
    P.h_bar_0, P.h_bar_jump, P.h_bar_n);
fprintf('\n');
fprintf('Shock parameters:\n');
fprintf('  nu_fert = %.2f, nu_loc = %.2f\n', P.nu_fert, P.nu_loc);
fprintf('  u_bar = %.1f\n', P.u_bar);
fprintf('\n');

%% Run model
fprintf('Running annual model...\n');
fprintf('State space: %d locations × %d wealth × %d tenure × %d periods × %d parity × %d child_states\n', ...
    3, P.Nb, 2, P.J, 4, P.A_m + 2);
total_states = 3 * P.Nb * 2 * P.J * 4 * (P.A_m + 2);
fprintf('Total states: %d (%.1f million)\n\n', total_states, total_states/1e6);

tic;
[sol, P_out, p_eq] = run_model_fast(P);
elapsed = toc;

%% Results
fprintf('\n========================================\n');
fprintf('ANNUAL MODEL RESULTS\n');
fprintf('========================================\n');
fprintf('Total elapsed time: %.1f seconds (%.1f minutes)\n', elapsed, elapsed/60);
fprintf('\n');
fprintf('Population growth (annual): %.4f%%\n', 100*sol.n_pop);
fprintf('Mean parity: %.3f\n', sol.mean_parity);
fprintf('Ownership rate: %.1f%%\n', 100*sol.own_rate);
fprintf('Childless rate: %.1f%%\n', 100*mean(sol.frac_childless_by_loc));
fprintf('\n');
fprintf('Prices: [%.3f, %.3f, %.3f]\n', p_eq(1), p_eq(2), p_eq(3));
fprintf('Rents:  [%.4f, %.4f, %.4f]\n', ...
    P_out.user_cost_rate*p_eq(1), P_out.user_cost_rate*p_eq(2), P_out.user_cost_rate*p_eq(3));
fprintf('\n');
fprintf('Fertility by location:\n');
loc_names = {'Peripheral', 'Secondary', 'Superstar'};
for i = 1:3
    fprintf('  %s: mean_n=%.3f, childless=%.1f%%\n', ...
        loc_names{i}, sol.mean_parity_by_loc(i), 100*sol.frac_childless_by_loc(i));
end
