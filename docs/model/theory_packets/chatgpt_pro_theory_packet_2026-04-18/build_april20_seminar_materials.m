function out = build_april20_seminar_materials(anchor_name)
% BUILD_APRIL20_SEMINAR_MATERIALS
% Generate seminar-grade benchmark figures and a counterfactual block from
% the live April 20 presentation anchor.

    if nargin < 1 || isempty(anchor_name)
        anchor_name = 'stageA_bench_main';
    end

    this_dir = fileparts(mfilename('fullpath'));
    addpath(this_dir);
    out_dir = fullfile(this_dir, 'output');
    fig_dir = fullfile(out_dir, 'figures');
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    old_vis = get(0, 'DefaultFigureVisible');
    cleanup_vis = onCleanup(@() set(0, 'DefaultFigureVisible', old_vis));
    set(0, 'DefaultFigureVisible', 'off');

    % Refresh the standard battery and scorecard from the live anchor.
    generate_candidate_figure_battery(anchor_name, 'current_candidate');

    anchor_file = fullfile(out_dir, sprintf('%s.mat', anchor_name));
    S = load(anchor_file, 'sol', 'P_out', 'p_eq');
    sol = S.sol;
    P = S.P_out;
    p_eq = S.p_eq;
    save(fullfile(out_dir, 'last_run_dt.mat'), 'sol', 'P', 'p_eq', '-v7');

    style = figure_style();
    b_grid = rebuild_grid(P);
    [targets, ~, inversion_targets] = build_calibration_setup('benchmark'); %#ok<ASGLU>

    export_benchmark_figures(sol, P, p_eq, b_grid, fig_dir, style, targets, inversion_targets);
    build_april20_temp_equilibrium_figures(anchor_name);

    counterfactuals = run_counterfactual_block(anchor_name, P);
    export_counterfactual_figures(counterfactuals, fig_dir, style);
    write_counterfactual_table(counterfactuals, fullfile(fig_dir, 'seminar_counterfactual_table.tex'));
    write_summary_file(sol, P, p_eq, counterfactuals, fullfile(fig_dir, 'seminar_materials_summary.txt'), anchor_name);
    write_live_manifest(anchor_name, fig_dir);
    organize_slide_figures(fig_dir);

    out = struct();
    out.anchor_name = anchor_name;
    out.anchor_file = anchor_file;
    out.figure_dir = fig_dir;
    out.counterfactuals = counterfactuals;
end

function write_live_manifest(anchor_name, fig_dir)
    fid = fopen(fullfile(fig_dir, 'manifest.txt'), 'w');
    if fid < 0
        return;
    end
    fprintf(fid, 'live seminar figure bundle\n');
    fprintf(fid, 'anchor=%s\n', anchor_name);
    fprintf(fid, 'built=%s\n', datestr(now, 31));
    fclose(fid);
end


function export_benchmark_figures(sol, P, p_eq, b_grid, fig_dir, style, targets, inversion_targets)
    age_focus = 15;  % age 32
    age_alt = 10;    % age 27
    age_vec = P.age_start + (0:P.J-1);

    export_spatial_equilibrium(sol, P, p_eq, fig_dir, style, targets, inversion_targets);
    export_lifecycle_sorting(sol, P, age_vec, fig_dir, style);
    export_housing_policy(sol, P, p_eq, b_grid, age_focus, fig_dir, style);
    export_savings_tenure(sol, P, p_eq, b_grid, age_focus, fig_dir, style);
    export_fertility_state(sol, P, b_grid, age_alt, age_focus, fig_dir, style);
    export_h01_h12_mechanism(sol, P, b_grid, age_vec, fig_dir, style, targets);
    export_owner_ladder_thresholds(P, fig_dir, style);
end


function export_spatial_equilibrium(sol, P, p_eq, fig_dir, style, targets, inversion_targets)
    fig = make_figure([90 90 1500 900], style);
    tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    colors = style.loc_colors;
    rent_eq = P.user_cost_rate * p_eq;

    nexttile(tl, 1);
    cats = categorical({'Price', 'Annual user cost', 'Population share'});
    vals = [p_eq(:), rent_eq(:), sol.pop_share(:)]';
    bh = bar(cats, vals, 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = colors(1,:);
    bh(2).FaceColor = colors(2,:);
    ylabel('Level');
    title('Equilibrium Prices and Population');
    legend({'Periphery','Center'}, 'Location', 'northwest');
    grid on;
    ylim([0, max(vals(:)) * 1.22]);

    nexttile(tl, 2);
    metric_labels = categorical({'Own 30-55', 'Own family gap', 'Center share new parents', 'Old-age parent gap'});
    model_vals = [sol.own_rate_3055; sol.own_gap_newparent_nonparent_3055; ...
        sol.center_share_newparents_2245; sol.old_age_parent_childless_gap_6575];
    target_vals = [targets.own_rate; targets.own_family_gap; ...
        targets.center_share_newparents; targets.old_age_parent_childless_gap];
    bh = bar(metric_labels, [model_vals, target_vals], 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.main_blue;
    bh(2).FaceColor = style.light_gray;
    ylabel('Share');
    title('Benchmark vs Targeted Housing and Tenure Moments');
    legend({'Model','Target'}, 'Location', 'northwest');
    grid on;

    nexttile(tl, 3);
    gap_labels = categorical({'H01: +3 first birth', 'H12: +3 second-birth proxy', ...
        'Young liquid wealth / income', 'Center pop share', 'Rent ratio C/P'});
    model_vals = [sol.housing_increment_0to1_eventstudy_t3; sol.housing_increment_1to2_proxy_t3; ...
        sol.young_liquid_wealth_to_income; sol.pop_share(2); ...
        (P.user_cost_rate * p_eq(2)) / max(P.user_cost_rate * p_eq(1), 1e-12)];
    target_vals = [targets.housing_increment_0to1; targets.housing_increment_1to2; ...
        targets.young_liquid_wealth_to_income; inversion_targets.pop_share_C; inversion_targets.rent_ratio];
    bh = bar(gap_labels, [model_vals, target_vals], 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.gold;
    bh(2).FaceColor = style.light_gray;
    ylabel('Level');
    title('Housing and Geography Anchors');
    legend({'Model','Target'}, 'Location', 'northwest');
    grid on;

    nexttile(tl, 4);
    own_by_loc = [sum(sol.g(:,2:end,1,1:P.J_R,:,:), 'all') / max(sum(sol.g(:,:,1,1:P.J_R,:,:), 'all'), 1e-12); ...
        sum(sol.g(:,2:end,2,1:P.J_R,:,:), 'all') / max(sum(sol.g(:,:,2,1:P.J_R,:,:), 'all'), 1e-12)];
    bar(categorical({'Periphery','Center'}), [own_by_loc, sol.mean_parity_by_loc], 'grouped', 'LineWidth', 0.8);
    ax = gca;
    ax.Children(1).FaceColor = style.soft_green;
    ax.Children(2).FaceColor = style.deep_red;
    ylabel('Rate');
    title('Prime-Age Ownership and Completed Fertility by Location');
    legend({'Own rate (workers)','Mean parity'}, 'Location', 'northwest');
    grid on;

    title(tl, 'Benchmark Anchor: Spatial Equilibrium and Target Fit', ...
        'FontWeight', 'bold', 'FontSize', style.title_size);
    exportgraphics(fig, fullfile(fig_dir, 'seminar01_spatial_equilibrium.png'), 'Resolution', 200);
    close(fig);
end


function export_lifecycle_sorting(sol, P, age_vec, fig_dir, style)
    fig = make_figure([110 110 1500 900], style);
    tl = tiledlayout(fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    g = sol.g;
    np = P.n_parity;
    b_grid = rebuild_grid(P);

    nexttile(tl, 1);
    own_loc_age = zeros(P.J, 2);
    for j = 1:P.J
        for i = 1:2
            gi = g(:,:,i,j,:,:);
            own_loc_age(j,i) = sum(gi(:,2:end,:,:,:,:), 'all') / max(sum(gi, 'all'), 1e-12);
        end
    end
    plot(age_vec, own_loc_age(:,1), '-', 'Color', style.loc_colors(1,:), 'LineWidth', 2.5); hold on;
    plot(age_vec, own_loc_age(:,2), '-', 'Color', style.loc_colors(2,:), 'LineWidth', 2.5);
    hold off; grid on;
    ylim([0 1]);
    xlabel('Age'); ylabel('Ownership rate');
    title('Ownership over the Lifecycle');
    legend({'Periphery','Center'}, 'Location', 'northwest');

    nexttile(tl, 2);
    par_age = zeros(P.J, np);
    for j = 1:P.J
        gj = g(:,:,:,j,:,:);
        mj = sum(gj, 'all');
        for nn = 1:np
            par_age(j, nn) = sum(gj(:,:,:,:,nn,:), 'all') / max(mj, 1e-12);
        end
    end
    ha = area(age_vec, par_age);
    parity_colors = style.parity_colors(1:np,:);
    for k = 1:np
        ha(k).FaceColor = parity_colors(k,:);
        ha(k).EdgeColor = 'none';
    end
    grid on; ylim([0 1]);
    xlabel('Age'); ylabel('Share');
    title('Family Size Composition');
    legend(arrayfun(@(x) sprintf('n = %d', x), 0:np-1, 'UniformOutput', false), ...
        'Location', 'eastoutside');

    nexttile(tl, 3);
    fert_rate = sol.fert_by_age(:);
    plot(age_vec, 2 * cumsum(fert_rate), '-', 'Color', style.main_blue, 'LineWidth', 2.5); hold on;
    yyaxis right;
    bar(age_vec(P.A_f_start:P.A_f_end), fert_rate(P.A_f_start:P.A_f_end), 0.8, ...
        'FaceColor', style.gold, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    yyaxis left;
    ylabel('Cumulative births');
    yyaxis right;
    ylabel('Birth flow');
    xlabel('Age');
    title('Fertility Timing');
    grid on;

    nexttile(tl, 4);
    center_share_status = NaN(P.J, 3);
    dep_states = 2:(P.n_child_stages + 1);
    for j = 1:P.J
        mass0 = sum(g(:,:,:,j,1,1), 'all');
        if mass0 > 1e-12
            center_share_status(j,1) = sum(g(:,:,2,j,1,1), 'all') / mass0;
        end
        if np >= 2
            mass1 = sum(g(:,:,:,j,2,dep_states), 'all');
            if mass1 > 1e-12
                center_share_status(j,2) = sum(g(:,:,2,j,2,dep_states), 'all') / mass1;
            end
        end
        if np >= 3
            mass2 = sum(g(:,:,:,j,3:np,dep_states), 'all');
            if mass2 > 1e-12
                center_share_status(j,3) = sum(g(:,:,2,j,3:np,dep_states), 'all') / mass2;
            end
        end
    end
    plot(age_vec, center_share_status(:,1), '-', 'Color', style.parity_colors(1,:), 'LineWidth', 2.5); hold on;
    plot(age_vec, center_share_status(:,2), '-', 'Color', style.parity_colors(2,:), 'LineWidth', 2.5);
    plot(age_vec, center_share_status(:,3), '-', 'Color', style.parity_colors(3,:), 'LineWidth', 2.5);
    hold off; grid on; ylim([0 1]);
    xlabel('Age'); ylabel('Center share');
    title('Where Families Live');
    legend({'Childless','One child','Two-plus'}, 'Location', 'southwest');

    nexttile(tl, 5);
    mean_b_age = zeros(P.J, 1);
    for j = 1:P.J
        gj = g(:,:,:,j,:,:);
        mass_by_b = squeeze(sum(gj, [2, 3, 5, 6]));
        mean_b_age(j) = sum(mass_by_b(:) .* b_grid(:)) / max(sum(gj, 'all'), 1e-12);
    end
    plot(age_vec, mean_b_age, '-', 'Color', style.deep_red, 'LineWidth', 2.5);
    xlabel('Age'); ylabel('Liquid wealth');
    title('Financial Wealth');
    grid on;

    nexttile(tl, 6);
    own_by_parity = zeros(np, 1);
    for nn = 1:np
        own_by_parity(nn) = sum(g(:,2:end,:,:,nn,:), 'all') / max(sum(g(:,:,:,:,nn,:), 'all'), 1e-12);
    end
    bh = bar(0:np-1, own_by_parity, 'FaceColor', 'flat', 'LineWidth', 0.8);
    bh.CData = style.parity_colors(1:np,:);
    xlabel('Completed fertility');
    ylabel('Ownership rate');
    title('Ownership by Family Size');
    ylim([0 1]); grid on;

    title(tl, 'Benchmark Anchor: Lifecycle and Spatial Sorting', ...
        'FontWeight', 'bold', 'FontSize', style.title_size);
    exportgraphics(fig, fullfile(fig_dir, 'seminar02_lifecycle_sorting.png'), 'Resolution', 200);
    close(fig);
end


function export_housing_policy(sol, P, p_eq, b_grid, j, fig_dir, style)
    fig = make_figure([120 120 1500 900], style);
    tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    ib = find(b_grid >= -2 & b_grid <= 18);
    labels = {'Current renter state', sprintf('Current owner state: H_1 = %.1f rooms', P.H_own(1))};
    from_tenure = [1, 2];

    for idx = 1:2
        to = from_tenure(idx);
        for iloc = 1:2
            nexttile(tl, (idx-1) * 2 + iloc);
            hold on;
            for nn = 1:min(3, P.n_parity)
                cs = child_state_for_plot(P, nn);
                h_act = realized_housing_from_state(sol, P, ib, to, iloc, j, nn, cs);
                plot(b_grid(ib), h_act, '-', 'Color', style.parity_colors(nn,:), 'LineWidth', 2.4);
            end

            if to == 1
                yline(P.hR_max, '--', 'Color', style.dark_gray, 'LineWidth', 1.2, 'HandleVisibility', 'off');
                dp = (1 - P.phi(1)) * p_eq(iloc) * P.H_own(:);
                for k = 1:numel(dp)
                    xline(dp(k), ':', 'Color', [0.82 0.82 0.82], 'LineWidth', 1.0, 'HandleVisibility', 'off');
                end
            end

            hold off;
            grid on;
            xlabel('Liquid wealth');
            ylabel('Housing services');
            title(sprintf('%s, %s', labels{idx}, loc_name(iloc)));
            if idx == 1 && iloc == 1
                legend({'Childless', 'One child', 'Two-plus'}, 'Location', 'northwest');
            end
            xlim([b_grid(ib(1)), b_grid(ib(end))]);
        end
    end

    title(tl, sprintf('Housing Policy Functions at Age %d', P.age_start + j - 1), ...
        'FontWeight', 'bold', 'FontSize', style.title_size);
    exportgraphics(fig, fullfile(fig_dir, 'seminar03_housing_policy.png'), 'Resolution', 200);
    close(fig);
end


function export_savings_tenure(sol, P, p_eq, b_grid, j, fig_dir, style)
    fig = make_figure([130 130 1500 900], style);
    tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    ib = find(b_grid >= -2 & b_grid <= 16);
    parities = 1:min(3, P.n_parity);

    for iloc = 1:2
        nexttile(tl, iloc);
        hold on;
        plot(b_grid(ib), b_grid(ib), '--', 'Color', style.dark_gray, 'LineWidth', 1.2, 'HandleVisibility', 'off');
        for nn = parities
            cs = child_state_for_plot(P, nn);
            bp = squeeze(sol.bp_pol(ib,1,iloc,j,nn,cs));
            plot(b_grid(ib), bp, '-', 'Color', style.parity_colors(nn,:), 'LineWidth', 2.4);
        end
        hold off;
        grid on;
        xlabel('Liquid wealth');
        ylabel('Next-period liquid wealth');
        title(sprintf('Renter savings, %s', loc_name(iloc)));
        if iloc == 1
            legend({'Childless','One child','Two-plus'}, 'Location', 'northwest');
        end
    end

    for iloc = 1:2
        nexttile(tl, 2 + iloc);
        hold on;
        for nn = parities
            cs = child_state_for_plot(P, nn);
            tc = squeeze(sol.tenure_choice(ib,1,iloc,j,nn,cs));
            stairs(b_grid(ib), tc - 1, '-', 'Color', style.parity_colors(nn,:), 'LineWidth', 2.4);
        end
        dp1 = (1 - P.phi(1)) * p_eq(iloc) * P.H_own(1);
        xline(dp1, '--', 'Color', style.dark_gray, 'LineWidth', 1.2, 'HandleVisibility', 'off');
        hold off;
        grid on;
        xlabel('Liquid wealth');
        ylabel('Chosen owner rung (0 = rent)');
        title(sprintf('Tenure thresholds, %s', loc_name(iloc)));
        yticks(0:P.n_house);
        ylim([-0.2, P.n_house + 0.2]);
    end

    title(tl, sprintf('Savings and Tenure Thresholds at Age %d', P.age_start + j - 1), ...
        'FontWeight', 'bold', 'FontSize', style.title_size);
    exportgraphics(fig, fullfile(fig_dir, 'seminar04_savings_tenure.png'), 'Resolution', 200);
    close(fig);
end


function export_fertility_state(sol, P, b_grid, j1, j2, fig_dir, style)
    fig = make_figure([140 140 1500 900], style);
    tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    ib = find(b_grid >= 0 & b_grid <= 16);
    ages = [j1, j2];
    tenures = [1, 2];
    title_bits = {'Current renter', 'Current owner H_1'};

    for a = 1:2
        for t = 1:2
            nexttile(tl, (a-1) * 2 + t);
            hold on;
            for iloc = 1:2
                pr = squeeze(sol.fert_probs(ib, tenures(t), iloc, ages(a), 2));
                plot(b_grid(ib), pr, '-', 'Color', style.loc_colors(iloc,:), 'LineWidth', 2.6);
            end
            hold off;
            grid on;
            xlabel('Liquid wealth');
            ylabel('Pr(first birth)');
            title(sprintf('Age %d, %s', P.age_start + ages(a) - 1, title_bits{t}));
            ylim([0, 0.22]);
            if a == 1 && t == 1
                legend({'Periphery','Center'}, 'Location', 'northwest');
            end
        end
    end

    title(tl, 'Fertility Policy by Housing State', ...
        'FontWeight', 'bold', 'FontSize', style.title_size);
    exportgraphics(fig, fullfile(fig_dir, 'seminar05_fertility_state.png'), 'Resolution', 200);
    close(fig);
end


function export_h01_h12_mechanism(sol, P, b_grid, age_vec, fig_dir, style, targets)
    fig = make_figure([150 150 1500 950], style);
    tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    g = sol.g;
    dep_states = 2:(P.n_child_stages + 1);

    nexttile(tl, 1);
    cats = categorical({'H01', 'H12'});
    model_vals = [sol.housing_increment_0to1_eventstudy_t3; sol.housing_increment_1to2_proxy_t3];
    target_vals = [targets.housing_increment_0to1; targets.housing_increment_1to2];
    bh = bar(cats, [model_vals, target_vals], 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.main_blue;
    bh(2).FaceColor = style.light_gray;
    ylabel('Rooms, +3 coefficient');
    title('Housing Responses at the Target Horizon');
    legend({'Model','Target'}, 'Location', 'northwest');
    grid on;

    nexttile(tl, 2);
    own_status = NaN(P.J, 3);
    for j = 1:P.J
        mass0 = sum(g(:,:,:,j,1,1), 'all');
        if mass0 > 1e-12
            own_status(j,1) = sum(g(:,2:end,:,j,1,1), 'all') / mass0;
        end
        if P.n_parity >= 2
            mass1 = sum(g(:,:,:,j,2,dep_states), 'all');
            if mass1 > 1e-12
                own_status(j,2) = sum(g(:,2:end,:,j,2,dep_states), 'all') / mass1;
            end
        end
        if P.n_parity >= 3
            mass2 = sum(g(:,:,:,j,3:P.n_parity,dep_states), 'all');
            if mass2 > 1e-12
                own_status(j,3) = sum(g(:,2:end,:,j,3:P.n_parity,dep_states), 'all') / mass2;
            end
        end
    end
    plot(age_vec, own_status(:,1), '-', 'Color', style.parity_colors(1,:), 'LineWidth', 2.5); hold on;
    plot(age_vec, own_status(:,2), '-', 'Color', style.parity_colors(2,:), 'LineWidth', 2.5);
    plot(age_vec, own_status(:,3), '-', 'Color', style.parity_colors(3,:), 'LineWidth', 2.5);
    hold off; grid on; ylim([0 1]);
    xlabel('Age'); ylabel('Ownership rate');
    title('Most of the Ownership Adjustment Happens at First Birth');
    legend({'Childless','One child','Two-plus'}, 'Location', 'northwest');

    nexttile(tl, 3);
    center_status = NaN(P.J, 3);
    for j = 1:P.J
        mass0 = sum(g(:,:,:,j,1,1), 'all');
        if mass0 > 1e-12
            center_status(j,1) = sum(g(:,:,2,j,1,1), 'all') / mass0;
        end
        if P.n_parity >= 2
            mass1 = sum(g(:,:,:,j,2,dep_states), 'all');
            if mass1 > 1e-12
                center_status(j,2) = sum(g(:,:,2,j,2,dep_states), 'all') / mass1;
            end
        end
        if P.n_parity >= 3
            mass2 = sum(g(:,:,:,j,3:P.n_parity,dep_states), 'all');
            if mass2 > 1e-12
                center_status(j,3) = sum(g(:,:,2,j,3:P.n_parity,dep_states), 'all') / mass2;
            end
        end
    end
    plot(age_vec, center_status(:,1), '-', 'Color', style.parity_colors(1,:), 'LineWidth', 2.5); hold on;
    plot(age_vec, center_status(:,2), '-', 'Color', style.parity_colors(2,:), 'LineWidth', 2.5);
    plot(age_vec, center_status(:,3), '-', 'Color', style.parity_colors(3,:), 'LineWidth', 2.5);
    hold off; grid on; ylim([0 1]);
    xlabel('Age'); ylabel('Center share');
    title('Families Are Already More Peripheral by the Second Child');

    nexttile(tl, 4);
    ib = find(b_grid >= 0 & b_grid <= 18);
    cap_share = zeros(2, 3);
    for iloc = 1:2
        g0 = squeeze(sol.g(ib,1,iloc,15,1,1));
        h0 = squeeze(sol.hR_pol(ib,1,iloc,15,1,1));
        cap_share(iloc,1) = sum(g0(h0 >= P.hR_max - 1e-8)) / max(sum(g0), 1e-12);
        if P.n_parity >= 2
            cs1 = child_state_for_plot(P, 2);
            g1 = squeeze(sol.g(ib,1,iloc,15,2,cs1));
            h1 = squeeze(sol.hR_pol(ib,1,iloc,15,2,cs1));
            cap_share(iloc,2) = sum(g1(h1 >= P.hR_max - 1e-8)) / max(sum(g1), 1e-12);
        end
        if P.n_parity >= 3
            cs2 = child_state_for_plot(P, 3);
            g2 = squeeze(sol.g(ib,1,iloc,15,3,cs2));
            h2 = squeeze(sol.hR_pol(ib,1,iloc,15,3,cs2));
            cap_share(iloc,3) = sum(g2(h2 >= P.hR_max - 1e-8)) / max(sum(g2), 1e-12);
        end
    end
    bh = bar(categorical({'Periphery','Center'}), cap_share, 'grouped', 'LineWidth', 0.8);
    for k = 1:size(cap_share,2)
        bh(k).FaceColor = style.parity_colors(k,:);
    end
    ylabel('Share of renter branch at cap');
    title('Extra Housing for Larger Families Comes Through the Owner Ladder');
    legend({'Childless','One child','Two-plus'}, 'Location', 'northwest');
    grid on; ylim([0 1]);

    title(tl, 'Why H01 Is Close While H12 Remains Low', ...
        'FontWeight', 'bold', 'FontSize', style.title_size);
    exportgraphics(fig, fullfile(fig_dir, 'seminar06_h01_h12_mechanism.png'), 'Resolution', 200);
    close(fig);
end


function export_owner_ladder_thresholds(P, fig_dir, style)
    fig = make_figure([155 155 1400 760], style);
    tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    owner_rungs = 1:P.n_house;
    owner_rooms = P.H_own(:);
    owner_services = P.chi * owner_rooms;
    need0 = P.h_bar_0;
    need1 = P.h_bar_0 + P.h_bar_jump + P.h_bar_n;
    need2 = P.h_bar_0 + P.h_bar_jump + 2 * P.h_bar_n;

    nexttile(tl, 1);
    bh = bar(owner_rungs, owner_rooms, 'FaceColor', style.main_blue, 'LineWidth', 0.8);
    %#ok<NASGU>
    hold on;
    yline(P.hR_max, '--', 'Color', style.deep_red, 'LineWidth', 1.8);
    yline(need0, ':', 'Color', style.dark_gray, 'LineWidth', 1.6);
    yline(need1, ':', 'Color', style.parity_colors(2,:), 'LineWidth', 1.8);
    yline(need2, ':', 'Color', style.parity_colors(3,:), 'LineWidth', 1.8);
    hold off;
    xlabel('Owner rung');
    ylabel('Rooms');
    title('Owner Ladder in Raw Rooms');
    xticks(owner_rungs);
    xticklabels(compose('H_%d', owner_rungs));
    legend({'Owner rung','Rental cap', 'Need: childless', 'Need: first child', 'Need: second child'}, ...
        'Location', 'northwest');
    grid on;

    nexttile(tl, 2);
    bh = bar(owner_rungs, owner_services, 'FaceColor', style.soft_green, 'LineWidth', 0.8);
    %#ok<NASGU>
    hold on;
    yline(need0, ':', 'Color', style.dark_gray, 'LineWidth', 1.6);
    yline(need1, ':', 'Color', style.parity_colors(2,:), 'LineWidth', 1.8);
    yline(need2, ':', 'Color', style.parity_colors(3,:), 'LineWidth', 1.8);
    hold off;
    xlabel('Owner rung');
    ylabel('Effective housing services');
    title('Effective Owner Services: \chi H_k');
    xticks(owner_rungs);
    xticklabels(compose('H_%d', owner_rungs));
    grid on;

    title(tl, 'Benchmark Calibration: Owner Ladder versus Family-Space Needs', ...
        'FontWeight', 'bold', 'FontSize', style.title_size);
    exportgraphics(fig, fullfile(fig_dir, 'seminar11_owner_ladder_thresholds.png'), 'Resolution', 200);
    close(fig);
end


function counterfactuals = run_counterfactual_block(anchor_name, P_anchor)
    this_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(this_dir, 'output');
    cases = cell(12,1);

    cases{1} = build_cf_case('Benchmark', 'baseline', struct(), [], true);
    cases{2} = build_cf_case('10% all', 'own_access', ...
        struct('phi', 0.90 * ones(P_anchor.n_parity, 1)), ...
        'Easier ownership access for all households', false);
    cases{3} = build_cf_case('0% all', 'own_access_zero', ...
        struct('phi', 1.00 * ones(P_anchor.n_parity, 1)), ...
        'Zero down payment for all households', false);
    cases{4} = build_cf_case('10% parents', 'parent_access_10', ...
        struct('parent_dp_waiver', true, 'parent_dp_waiver_phi', 0.90), ...
        'Parent-only 10 percent down payment', false);
    cases{5} = build_cf_case('0% parents', 'parent_access_zero', ...
        struct('parent_dp_waiver', true, 'parent_dp_waiver_phi', 1.00), ...
        'Parent-only zero down payment', false);
    cases{6} = build_cf_case('Birth grant', 'birth_grant', ...
        struct('birth_dp_grant', true), ...
        'Birth-triggered down-payment grant', false);
    cases{7} = build_cf_case('10% starter', 'starter_access_10', ...
        struct('parent_dp_waiver', true, 'parent_dp_waiver_phi', 0.90, ...
            'parent_dp_waiver_birth_state_only', true, ...
            'parent_dp_waiver_owner_rungs', 1:2), ...
        'Starter-home access at birth', false);
    cases{8} = build_cf_case('0% starter', 'starter_access_zero', ...
        struct('parent_dp_waiver', true, 'parent_dp_waiver_phi', 1.00, ...
            'parent_dp_waiver_birth_state_only', true, ...
            'parent_dp_waiver_owner_rungs', 1:2), ...
        'Zero-down starter-home access at birth', false);
    cases{9} = build_cf_case('10% center starter', 'center_starter_10', ...
        struct('parent_dp_waiver', true, 'parent_dp_waiver_phi', 0.90, ...
            'parent_dp_waiver_birth_state_only', true, ...
            'parent_dp_waiver_owner_rungs', 1:2, ...
            'parent_dp_waiver_locations', 2), ...
        'Center-only starter-home access at birth', false);
    cases{10} = build_cf_case('0% center starter', 'center_starter_zero', ...
        struct('parent_dp_waiver', true, 'parent_dp_waiver_phi', 1.00, ...
            'parent_dp_waiver_birth_state_only', true, ...
            'parent_dp_waiver_owner_rungs', 1:2, ...
            'parent_dp_waiver_locations', 2), ...
        'Zero-down center starter-home access at birth', false);
    cases{11} = build_cf_case('15% supply', 'family_supply', ...
        struct('H0', 1.15 * P_anchor.H0), ...
        'More family-sized housing supply', false);
    cases{12} = build_cf_case('30% center supply', 'core_family', ...
        struct('H0', [P_anchor.H0(1); 1.30 * P_anchor.H0(2)]), ...
        'Family space in the core', false);

    counterfactuals = cell(numel(cases), 1);
    for i = 1:numel(cases)
        case_i = cases{i};
        if case_i.is_baseline
            out_file = fullfile(out_dir, sprintf('%s.mat', anchor_name));
            S = load(out_file, 'sol', 'P_out', 'p_eq');
            counterfactuals{i} = extract_case_summary(case_i, S.sol, S.P_out, S.p_eq, out_file, '');
        else
            save_file = fullfile(out_dir, sprintf('%s_%s_%s.mat', 'benchmark', anchor_name, case_i.save_label));
            summary_file = fullfile(out_dir, sprintf('%s_%s_%s.txt', 'benchmark', anchor_name, case_i.save_label));
            if exist(save_file, 'file')
                S = load(save_file, 'sol', 'P_out', 'p_eq');
                counterfactuals{i} = extract_case_summary(case_i, S.sol, S.P_out, S.p_eq, save_file, summary_file);
            else
                out = confirm_override_from_saved_anchor(anchor_name, case_i.override, case_i.save_label, 'benchmark');
                S = load(out.save_file, 'sol', 'P_out', 'p_eq');
                counterfactuals{i} = extract_case_summary(case_i, S.sol, S.P_out, S.p_eq, out.save_file, out.summary_file);
            end
        end
    end
    counterfactuals = vertcat(counterfactuals{:});
end


function export_counterfactual_figures(counterfactuals, fig_dir, style)
    export_access_figure(counterfactuals(1:6), fig_dir, style);
    export_targeted_access_figure(counterfactuals([1, 7, 8, 9, 10]), fig_dir, style);
    export_supply_figure(counterfactuals([1, 11, 12]), fig_dir, style);
    export_summary_figure(counterfactuals([1:6, 11, 12]), fig_dir, style);
end


function export_access_figure(cases, fig_dir, style)
    fig = make_figure([160 160 1750 820], style);
    tl = tiledlayout(fig, 1, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
    labels = {cases.display_label};
    xcat = categorical(labels, labels);

    nexttile(tl, 1);
    vals = [[cases.dp_small_P]', [cases.dp_small_C]'];
    bh = bar(xcat, vals, 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.loc_colors(1,:);
    bh(2).FaceColor = style.loc_colors(2,:);
    ylabel('Required liquid wealth for H_1');
    title('Down-Payment Threshold');
    legend({'Periphery','Center'}, 'Location', 'northwest');
    grid on;
    ax = gca;
    ax.XTickLabelRotation = 25;

    nexttile(tl, 2);
    vals = [[cases.own_nonparent_3055]', [cases.own_newparent_3055]'];
    bh = bar(xcat, vals, 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.parity_colors(1,:);
    bh(2).FaceColor = style.parity_colors(2,:);
    ylabel('Ownership rate');
    title('Prime-Age Ownership');
    legend({'Nonparents','New parents'}, 'Location', 'northwest');
    grid on; ylim([0 1]);
    ax = gca;
    ax.XTickLabelRotation = 25;

    nexttile(tl, 3);
    vals = [[cases.h01]', [cases.h12]'];
    bh = bar(xcat, vals, 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.main_blue;
    bh(2).FaceColor = style.gold;
    ylabel('Rooms, +3 coefficient');
    title('Family Housing Response');
    legend({'H01','H12'}, 'Location', 'northwest');
    grid on;
    ax = gca;
    ax.XTickLabelRotation = 25;

    nexttile(tl, 4);
    vals = [[cases.center_np]', [cases.center_newparent]'];
    bh = bar(xcat, vals, 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.light_gray;
    bh(2).FaceColor = style.deep_red;
    ylabel('Center share');
    title('Within-Metro Sorting');
    legend({'Nonparents','New parents'}, 'Location', 'northwest');
    grid on; ylim([0 1]);
    ax = gca;
    ax.XTickLabelRotation = 25;

    title(tl, 'Counterfactual 1: Easier Ownership Access', ...
        'FontWeight', 'bold', 'FontSize', style.title_size);
    exportgraphics(fig, fullfile(fig_dir, 'seminar07_cf_access.png'), 'Resolution', 200);
    close(fig);
end


function export_supply_figure(cases, fig_dir, style)
    fig = make_figure([170 170 1500 820], style);
    tl = tiledlayout(fig, 1, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
    labels = {cases.display_label};
    xcat = categorical(labels, labels);

    nexttile(tl, 1);
    vals = [[cases.p_P]', [cases.p_C]'];
    bh = bar(xcat, vals, 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.loc_colors(1,:);
    bh(2).FaceColor = style.loc_colors(2,:);
    ylabel('Equilibrium price');
    title('Housing Prices');
    legend({'Periphery','Center'}, 'Location', 'northwest');
    grid on;

    nexttile(tl, 2);
    vals = [[cases.h01]', [cases.h12]'];
    bh = bar(xcat, vals, 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.main_blue;
    bh(2).FaceColor = style.gold;
    ylabel('Rooms, +3 coefficient');
    title('Housing Responses');
    legend({'H01','H12'}, 'Location', 'northwest');
    grid on;

    nexttile(tl, 3);
    vals = [[cases.center_newparent]', [cases.own_3055]'];
    bh = bar(xcat, vals, 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.deep_red;
    bh(2).FaceColor = style.soft_green;
    ylabel('Share');
    title('Parent Centrality and Ownership');
    legend({'Center share new parents','Own 30-55'}, 'Location', 'northwest');
    grid on; ylim([0 1]);

    nexttile(tl, 4);
    vals = [cases.tfr]';
    bh = bar(xcat, vals, 'FaceColor', style.main_blue, 'LineWidth', 0.8);
    %#ok<NASGU>
    ylabel('Total fertility rate');
    title('Aggregate Fertility');
    grid on;

    title(tl, 'Counterfactual 2-3: More Family Space and Family Space in the Core', ...
        'FontWeight', 'bold', 'FontSize', style.title_size);
    exportgraphics(fig, fullfile(fig_dir, 'seminar08_cf_supply.png'), 'Resolution', 200);
    close(fig);
end


function export_targeted_access_figure(cases, fig_dir, style)
    fig = make_figure([175 175 1650 820], style);
    tl = tiledlayout(fig, 1, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
    labels = {cases.display_label};
    xcat = categorical(labels, labels);

    nexttile(tl, 1);
    vals = [[cases.dp_small_P]', [cases.dp_small_C]'];
    bh = bar(xcat, vals, 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.loc_colors(1,:);
    bh(2).FaceColor = style.loc_colors(2,:);
    ylabel('Required liquid wealth for H_1');
    title('Starter-Home Threshold');
    legend({'Periphery','Center'}, 'Location', 'northwest');
    grid on;
    ax = gca;
    ax.XTickLabelRotation = 25;

    nexttile(tl, 2);
    vals = [[cases.h01]', [cases.h12]'];
    bh = bar(xcat, vals, 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.main_blue;
    bh(2).FaceColor = style.gold;
    ylabel('Rooms, +3 coefficient');
    title('Birth-Related Housing Response');
    legend({'H01','H12'}, 'Location', 'northwest');
    grid on;
    ax = gca;
    ax.XTickLabelRotation = 25;

    nexttile(tl, 3);
    vals = [cases.tfr]';
    bh = bar(xcat, vals, 'FaceColor', style.main_blue, 'LineWidth', 0.8);
    %#ok<NASGU>
    ylabel('TFR');
    title('Aggregate Fertility');
    grid on;
    ax = gca;
    ax.XTickLabelRotation = 25;

    nexttile(tl, 4);
    vals = [[cases.center_newparent]', [cases.own_3055]'];
    bh = bar(xcat, vals, 'grouped', 'LineWidth', 0.8);
    bh(1).FaceColor = style.deep_red;
    bh(2).FaceColor = style.soft_green;
    ylabel('Share');
    title('Parent Centrality and Ownership');
    legend({'Center share new parents','Own 30-55'}, 'Location', 'northwest');
    grid on; ylim([0 1]);
    ax = gca;
    ax.XTickLabelRotation = 25;

    title(tl, 'Counterfactual 1B: Targeted Starter-Home and Center Access', ...
        'FontWeight', 'bold', 'FontSize', style.title_size);
    exportgraphics(fig, fullfile(fig_dir, 'seminar10_cf_targeted_access.png'), 'Resolution', 200);
    close(fig);
end


function export_summary_figure(cases, fig_dir, style)
    fig = make_figure([180 180 1550 920], style);
    tl = tiledlayout(fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    labels = {cases.display_label};
    x = 1:numel(cases);

    plot_metric_panel(nexttile(tl, 1), x, [cases.tfr], labels, 'TFR', style);
    plot_metric_panel(nexttile(tl, 2), x, [cases.h01], labels, 'H01: +3 first-birth rooms', style);
    plot_metric_panel(nexttile(tl, 3), x, [cases.h12], labels, 'H12: +3 additional-child proxy', style);
    plot_metric_panel(nexttile(tl, 4), x, [cases.own_3055], labels, 'Ownership, ages 30-55', style);
    plot_metric_panel(nexttile(tl, 5), x, [cases.center_newparent], labels, 'Center share, new parents', style);
    plot_metric_panel(nexttile(tl, 6), x, [cases.age_first_birth], labels, 'Mean age at first birth', style);

    title(tl, 'Counterfactual Summary: Benchmark-to-Policy Comparisons', ...
        'FontWeight', 'bold', 'FontSize', style.title_size);
    exportgraphics(fig, fullfile(fig_dir, 'seminar09_cf_summary.png'), 'Resolution', 200);
    close(fig);
end


function write_counterfactual_table(cases, path)
    fid = fopen(path, 'w');
    if fid < 0
        error('Could not open %s for writing.', path);
    end
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '\\begin{tabular}{lcccccc}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Scenario & TFR & H01 & H12 & Own 30--55 & Center share new parents & Age first birth \\\\\n');
    fprintf(fid, '\\midrule\n');
    for i = 1:numel(cases)
        fprintf(fid, '%s & %.3f & %.3f & %.3f & %.3f & %.3f & %.2f \\\\\n', ...
            latex_escape(cases(i).display_label), cases(i).tfr, cases(i).h01, cases(i).h12, ...
            cases(i).own_3055, cases(i).center_newparent, cases(i).age_first_birth);
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
end


function write_summary_file(sol, P, p_eq, cases, path, anchor_name)
    fid = fopen(path, 'w');
    if fid < 0
        error('Could not open %s for writing.', path);
    end
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, 'April 20 seminar materials\n\n');
    fprintf(fid, 'Benchmark anchor\n');
    fprintf(fid, 'label = %s\n', char(string(anchor_name)));
    fprintf(fid, 'TFR = %.6f\n', 2 * sol.mean_parity);
    fprintf(fid, 'H01 = %.6f\n', sol.housing_increment_0to1_eventstudy_t3);
    fprintf(fid, 'H12 = %.6f\n', sol.housing_increment_1to2_proxy_t3);
    fprintf(fid, 'Own 30-55 = %.6f\n', sol.own_rate_3055);
    fprintf(fid, 'Own family gap = %.6f\n', sol.own_gap_newparent_nonparent_3055);
    fprintf(fid, 'yW/I = %.6f\n', sol.young_liquid_wealth_to_income);
    fprintf(fid, 'Center share new parents = %.6f\n', sol.center_share_newparents_2245);
    fprintf(fid, 'Old-age parent-childless gap = %.6f\n', sol.old_age_parent_childless_gap_6575);
    fprintf(fid, 'p_P = %.6f\n', p_eq(1));
    fprintf(fid, 'p_C = %.6f\n\n', p_eq(2));

    fprintf(fid, 'Counterfactuals\n');
    for i = 1:numel(cases)
        fprintf(fid, '[%s]\n', cases(i).display_label);
        fprintf(fid, 'summary_file = %s\n', cases(i).summary_file);
        fprintf(fid, 'TFR = %.6f\n', cases(i).tfr);
        fprintf(fid, 'H01 = %.6f\n', cases(i).h01);
        fprintf(fid, 'H12 = %.6f\n', cases(i).h12);
        fprintf(fid, 'Own 30-55 = %.6f\n', cases(i).own_3055);
        fprintf(fid, 'Center share new parents = %.6f\n', cases(i).center_newparent);
        fprintf(fid, 'Age first birth = %.6f\n\n', cases(i).age_first_birth);
    end
end


function s = build_cf_case(display_label, save_label, override, mechanism_label, is_baseline)
    if nargin < 4 || isempty(mechanism_label)
        mechanism_label = display_label;
    end
    if nargin < 5
        is_baseline = false;
    end
    s = struct();
    s.display_label = display_label;
    s.save_label = save_label;
    s.override = override;
    s.mechanism_label = mechanism_label;
    s.is_baseline = is_baseline;
end


function out = extract_case_summary(case_i, sol, P, p_eq, save_file, summary_file)
    out = case_i;
    out.save_file = save_file;
    out.summary_file = summary_file;
    out.p_P = p_eq(1);
    out.p_C = p_eq(2);
    out.rent_ratio = (P.user_cost_rate * p_eq(2)) / max(P.user_cost_rate * p_eq(1), 1e-12);
    out.dp_small_P = seminar_required_wealth_for_small_owner(P, p_eq, 1);
    out.dp_small_C = seminar_required_wealth_for_small_owner(P, p_eq, 2);
    out.tfr = 2 * sol.mean_parity;
    out.childless = sol.parity_dist(1);
    out.age_first_birth = sol.mean_age_first_birth;
    out.h01 = sol.housing_increment_0to1_eventstudy_t3;
    out.h12 = sol.housing_increment_1to2_proxy_t3;
    out.own_3055 = sol.own_rate_3055;
    out.own_gap = sol.own_gap_newparent_nonparent_3055;
    out.own_nonparent_3055 = sol.own_rate_nonparents_3055;
    out.own_newparent_3055 = sol.own_rate_newparents_3055;
    out.own_nonparent = sol.own_rate_nonparents_3055;
    out.own_newparent = sol.own_rate_newparents_3055;
    out.center_np = sol.center_share_nonparents_2245;
    out.center_newparent = sol.center_share_newparents_2245;
    out.ywi = sol.young_liquid_wealth_to_income;
    out.old_gap = sol.old_age_parent_childless_gap_6575;
end


function dp = seminar_required_wealth_for_small_owner(P, p_eq, loc_idx)
    if P.n_parity >= 2
        nn = 2;
        cs = 2;
    else
        nn = 1;
        cs = 1;
    end
    phi_eff = seminar_effective_phi(P, loc_idx, 1, nn, cs);
    dp = (1 - phi_eff) * p_eq(loc_idx) * P.H_own(1);
end


function phi_eff = seminar_effective_phi(P, loc_idx, owner_rung, nn, cs)
    phi_eff = P.phi(min(nn, numel(P.phi)));
    if ~(isfield(P, 'parent_dp_waiver') && P.parent_dp_waiver)
        return;
    end
    if nn < 2
        return;
    end
    if isfield(P, 'parent_dp_waiver_birth_state_only') && P.parent_dp_waiver_birth_state_only
        if cs ~= 2
            return;
        end
    else
        K = P.n_child_stages;
        if cs < 2 || cs > K + 1
            return;
        end
    end
    if isfield(P, 'parent_dp_waiver_locations') && ~isempty(P.parent_dp_waiver_locations)
        if ~ismember(loc_idx, P.parent_dp_waiver_locations)
            return;
        end
    end
    if isfield(P, 'parent_dp_waiver_owner_rungs') && ~isempty(P.parent_dp_waiver_owner_rungs)
        if ~ismember(owner_rung, P.parent_dp_waiver_owner_rungs)
            return;
        end
    end
    phi_eff = max(phi_eff, P.parent_dp_waiver_phi);
end


function plot_metric_panel(ax, x, vals, labels, ttl, style)
    axes(ax); %#ok<LAXES>
    plot(x, vals, '-o', 'Color', style.main_blue, 'LineWidth', 2.5, ...
        'MarkerFaceColor', style.main_blue, 'MarkerSize', 7);
    grid on;
    title(ttl);
    xticks(x);
    xticklabels(labels);
    xtickangle(30);
    set(gca, 'FontSize', max(style.axis_size - 1, 10));
end


function fig = make_figure(pos, style)
    fig = figure('Color', 'w', 'Position', pos);
    set(fig, 'DefaultAxesFontName', style.font_name);
    set(fig, 'DefaultTextFontName', style.font_name);
    set(fig, 'DefaultAxesFontSize', style.axis_size);
    set(fig, 'DefaultTextFontSize', style.axis_size);
end


function style = figure_style()
    style = struct();
    style.font_name = 'Helvetica';
    style.axis_size = 12;
    style.title_size = 18;
    style.main_blue = [0.11 0.34 0.67];
    style.gold = [0.82 0.55 0.08];
    style.deep_red = [0.70 0.18 0.12];
    style.soft_green = [0.28 0.56 0.28];
    style.dark_gray = [0.30 0.30 0.30];
    style.light_gray = [0.80 0.82 0.84];
    style.loc_colors = [0.18 0.42 0.74; 0.80 0.30 0.14];
    style.parity_colors = [0.35 0.39 0.44; 0.08 0.47 0.56; 0.86 0.48 0.08; 0.62 0.28 0.55];
end


function b_grid = rebuild_grid(P)
    Nb = P.Nb;
    Nb1 = round(Nb * 0.15);
    Nb2 = round(Nb * 0.45);
    Nb3 = round(Nb * 0.15);
    Nb4 = Nb - Nb1 - Nb2 - Nb3;
    s1 = linspace(P.b_min, -3, Nb1 + 1)'; s1 = s1(1:end-1);
    s2 = linspace(-3, 6, Nb2 + 1)'; s2 = s2(1:end-1);
    s3 = linspace(6, 20, Nb3 + 1)'; s3 = s3(1:end-1);
    u4 = linspace(0, 1, Nb4 + 1)'; u4 = u4(2:end);
    s4 = 20 + (P.b_max - 20) * u4.^P.b_grid_power;
    b_grid = [s1; s2; s3; s4];
    [~, iz] = min(abs(b_grid)); b_grid(iz) = 0;
    [~, ie] = min(abs(b_grid - P.b_entry_fixed)); b_grid(ie) = P.b_entry_fixed;
end


function h = realized_housing_from_state(sol, P, ib, to, iloc, j, nn, cs)
    tc = squeeze(sol.tenure_choice(ib, to, iloc, j, nn, cs));
    h = NaN(numel(ib), 1);
    rent_mask = (tc == 1);
    if any(rent_mask)
        h(rent_mask) = squeeze(sol.hR_pol(ib(rent_mask), 1, iloc, j, nn, cs));
    end
    own_mask = (tc > 1);
    if any(own_mask)
        h(own_mask) = P.chi * P.H_own(tc(own_mask) - 1);
    end
end


function cs = child_state_for_plot(P, nn)
    if nn <= 1
        cs = 1;
    else
        cs = 2;
    end
end


function name = loc_name(i)
    if i == 1
        name = 'Periphery';
    else
        name = 'Center';
    end
end


function out = latex_escape(str)
    out = strrep(str, '%', '\%');
end
