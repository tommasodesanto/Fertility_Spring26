%% test_parfor_speedup.m
%
% Test BATCHED HJB solver optimization.
% The batched solver processes all (nn, cs) combinations in one call
% instead of 40 separate calls per age.
%
% ========================================================================

clear; clc;
addpath(fileparts(mfilename('fullpath')));

fprintf('=============================================================\n');
fprintf('  BATCHED HJB SOLVER TEST\n');
fprintf('=============================================================\n\n');

%% Test with quick annual model (J=20)
fprintf('Running quick test (J=20, annual periods)...\n');
fprintf('This tests the parfor optimization is working.\n\n');

P = struct();
P.J = 20;
P.da = 1.0;
P.J_R = 16;
P.A_m = 8;
P.A_f_start = 1;
P.A_f_end = 8;
P.Nb = 30;

% Use zero u_bar (logit doesn't need it)
P.u_bar = 0.0;

% Moderate shock parameters
P.nu_fert = 3.0;
P.nu_loc = 0.5;

fprintf('State space: 3 loc × %d wealth × 2 tenure × %d periods × 4 parity × %d child_states\n', ...
    P.Nb, P.J, P.A_m + 2);
total_states = 3 * P.Nb * 2 * P.J * 4 * (P.A_m + 2);
fprintf('Total states: %d (%.2f million)\n\n', total_states, total_states/1e6);

% Run model
tic;
[sol, P_out, p_eq] = run_model_fast(P);
elapsed = toc;

fprintf('\n========================================\n');
fprintf('RESULTS\n');
fprintf('========================================\n');
fprintf('Total elapsed time: %.1f seconds (%.1f minutes)\n', elapsed, elapsed/60);
fprintf('\n');
fprintf('Population growth (annual): %.4f%%\n', 100*sol.n_pop);
fprintf('Mean parity: %.3f\n', sol.mean_parity);
fprintf('Ownership rate: %.1f%%\n', 100*sol.own_rate);
fprintf('Childless rate: %.1f%%\n', 100*mean(sol.frac_childless_by_loc));
fprintf('\n');
fprintf('Prices: [%.3f, %.3f, %.3f]\n', p_eq(1), p_eq(2), p_eq(3));

% Basic sanity checks
fprintf('\n========================================\n');
fprintf('SANITY CHECKS\n');
fprintf('========================================\n');

passed = true;

if sol.mean_parity < 0 || sol.mean_parity > 10
    fprintf('[FAIL] Mean parity out of range: %.2f\n', sol.mean_parity);
    passed = false;
else
    fprintf('[PASS] Mean parity in reasonable range: %.2f\n', sol.mean_parity);
end

if sol.own_rate < 0 || sol.own_rate > 1
    fprintf('[FAIL] Ownership rate out of range: %.2f%%\n', 100*sol.own_rate);
    passed = false;
else
    fprintf('[PASS] Ownership rate in [0,1]: %.1f%%\n', 100*sol.own_rate);
end

if any(p_eq <= 0)
    fprintf('[FAIL] Non-positive prices detected\n');
    passed = false;
else
    fprintf('[PASS] All prices positive\n');
end

if abs(sol.n_pop) > 0.5
    fprintf('[WARN] Population growth seems extreme: %.4f\n', sol.n_pop);
else
    fprintf('[PASS] Population growth reasonable: %.4f\n', sol.n_pop);
end

fprintf('\n');
if passed
    fprintf('*** ALL CHECKS PASSED ***\n');
else
    fprintf('*** SOME CHECKS FAILED ***\n');
end

fprintf('\n========================================\n');
fprintf('TIMING SUMMARY\n');
fprintf('========================================\n');
fprintf('Previous baseline (40 calls/age): ~5.5 min for J=20\n');
fprintf('Current time (batched): %.1f min\n', elapsed/60);
fprintf('Speedup: %.1fx\n', 5.5 / (elapsed/60));
