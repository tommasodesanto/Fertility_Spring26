%% test_production_size.m
%
% Test the optimized run_model_fast.m with production-size parameters
%

clear; clc;

fprintf('==========================================================\n');
fprintf('PRODUCTION SIZE TEST - run_model_fast.m with drift optimization\n');
fprintf('==========================================================\n\n');

% Larger model - closer to production
P_test = struct();
P_test.J = 20;       % More age periods
P_test.da = 3.0;     % Smaller period length
P_test.J_R = 14;     % Retirement
P_test.A_f_start = 1;
P_test.A_f_end = 5;
P_test.A_m = 6;      % Kids take 6 periods to mature

fprintf('Configuration: J=%d, da=%.1f, Nb=30\n', P_test.J, P_test.da);
fprintf('Running optimized model...\n\n');

t_start = tic;
[sol, P, p_eq] = run_model_fast(P_test);
elapsed = toc(t_start);

fprintf('\n==========================================================\n');
fprintf('RESULTS\n');
fprintf('==========================================================\n');
fprintf('Total time: %.1f seconds (%.1f minutes)\n', elapsed, elapsed/60);
fprintf('gamma: %.6f (%.4f%% annual)\n', sol.gamma, 100*sol.gamma/P.da);
fprintf('n_pop: %.6f (%.4f%% annual)\n', sol.n_pop, 100*sol.n_pop/P.da);
fprintf('mean_parity: %.4f\n', sol.mean_parity);
fprintf('own_rate: %.4f\n', sol.own_rate);
fprintf('prices: [%.3f, %.3f, %.3f]\n', p_eq(1), p_eq(2), p_eq(3));
fprintf('==========================================================\n');
