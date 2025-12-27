%% test_drift_optimization.m
%
% Quick test of the in-place drift optimization in run_model_fast.m
%

clear; clc;

fprintf('==========================================================\n');
fprintf('TESTING IN-PLACE DRIFT OPTIMIZATION\n');
fprintf('==========================================================\n\n');

P_test = struct();
P_test.J = 10;
P_test.da = 6.0;
P_test.J_R = 6;
P_test.A_f_start = 1;
P_test.A_f_end = 3;
P_test.A_m = 3;

fprintf('Test configuration: J=%d, Nb=30\n', P_test.J);
fprintf('Running model with drift precomputation optimization...\n\n');

t_start = tic;
[sol, P, p_eq] = run_model_fast(P_test);
elapsed = toc(t_start);

fprintf('\n==========================================================\n');
fprintf('RESULTS\n');
fprintf('==========================================================\n');
fprintf('Elapsed time: %.1f seconds (%.1f minutes)\n', elapsed, elapsed/60);
fprintf('gamma: %.6f\n', sol.gamma);
fprintf('n_pop: %.6f\n', sol.n_pop);
fprintf('mean_parity: %.4f\n', sol.mean_parity);
fprintf('own_rate: %.4f\n', sol.own_rate);
fprintf('prices: [%.3f, %.3f, %.3f]\n', p_eq(1), p_eq(2), p_eq(3));
fprintf('==========================================================\n');
