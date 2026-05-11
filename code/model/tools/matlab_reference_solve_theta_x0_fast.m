function matlab_reference_solve_theta_x0_fast()
% MATLAB_REFERENCE_SOLVE_THETA_X0_FAST
% Direct theta -> GE solve benchmark matching:
%   python -m dt_cp_model.cli solve-theta --setup fast --theta x0 --max-iter-eq 120
%
% No SMM loss and no geography inversion.

    this_dir = fileparts(mfilename('fullpath'));
    port_dir = fileparts(this_dir);
    dt_dir = fileparts(port_dir);
    cd(dt_dir);
    addpath(dt_dir);

    [~, ~, ~, P_base, names, ~, ~, x0] = build_calibration_setup('fast');
    P_base.max_iter_eq = 120;

    for k = 1:numel(names)
        P_base.(names{k}) = x0(k);
    end
    P_base.eps_fert = P_base.kappa_fert;
    P_base.eps_loc = P_base.kappa_loc;

    tic;
    [sol, ~, p_eq] = run_model_cp_dt(P_base);
    elapsed = toc;

    report = struct();
    report.elapsed_sec = elapsed;
    report.p_eq = p_eq(:)';
    report.tfr = 2 * sol.mean_parity;
    report.own_rate = sol.own_rate;
    report.pop_share = sol.pop_share(:)';
    report.mean_age_first_birth = sol.mean_age_first_birth;
    report.migration_rate_2245 = sol.migration_rate_2245;
    report.prime_childless_renter_median_rooms = sol.prime_childless_renter_median_rooms;
    report.prime_childless_owner_median_rooms = sol.prime_childless_owner_median_rooms;
    report.housing_increment_0to1 = sol.housing_increment_0to1_eventstudy_t3;
    report.housing_increment_1to2 = sol.housing_increment_1to2_proxy_t3;
    report.total_mass = sol.total_mass;
    report.theta = x0(:)';

    out_dir = fullfile(port_dir, 'benchmarks');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    fid = fopen(fullfile(out_dir, 'matlab_solve_theta_x0_fast.json'), 'w');
    fwrite(fid, jsonencode(report, 'PrettyPrint', true));
    fclose(fid);

    fprintf('%s\n', jsonencode(report, 'PrettyPrint', true));
end

