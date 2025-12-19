% run_profiler.m - Run model with profiling
% Execute from MATLAB command line or batch mode

profile on;
run_model;
profile off;
p = profile('info');

% Print summary sorted by total time
fprintf('\n\n========================================\n');
fprintf('PROFILER RESULTS (Top 30 by Total Time)\n');
fprintf('========================================\n');

% Check available fields
ft1 = p.FunctionTable(1);
field_names = fieldnames(ft1);
fprintf('Available fields: %s\n', strjoin(field_names, ', '));
fprintf('%s\n', repmat('-', 1, 85));

% Determine correct field name for self time
if isfield(ft1, 'SelfTime')
    self_time_field = 'SelfTime';
elseif isfield(ft1, 'SelfTimeSec')
    self_time_field = 'SelfTimeSec';
else
    self_time_field = '';
end

fprintf('%-50s %10s %10s %10s\n', 'Function', 'Calls', 'TotalTime', 'SelfTime');
fprintf('%s\n', repmat('-', 1, 85));

% Sort by total time
total_times = arrayfun(@(x) x.TotalTime, p.FunctionTable);
[~, idx] = sort(total_times, 'descend');

for i = 1:min(30, length(idx))
    ft = p.FunctionTable(idx(i));
    name = ft.FunctionName;
    if length(name) > 50
        name = ['...' name(end-46:end)];
    end

    if ~isempty(self_time_field)
        self_time = ft.(self_time_field);
    else
        self_time = 0;
    end

    fprintf('%-50s %10d %10.2f %10.2f\n', name, ft.NumCalls, ft.TotalTime, self_time);
end

% Also show by self time (where actual computation happens)
if ~isempty(self_time_field)
    fprintf('\n\n========================================\n');
    fprintf('PROFILER RESULTS (Top 30 by Self Time)\n');
    fprintf('========================================\n');
    fprintf('%-50s %10s %10s %10s\n', 'Function', 'Calls', 'TotalTime', 'SelfTime');
    fprintf('%s\n', repmat('-', 1, 85));

    self_times = arrayfun(@(x) x.(self_time_field), p.FunctionTable);
    [~, idx] = sort(self_times, 'descend');

    for i = 1:min(30, length(idx))
        ft = p.FunctionTable(idx(i));
        name = ft.FunctionName;
        if length(name) > 50
            name = ['...' name(end-46:end)];
        end
        fprintf('%-50s %10d %10.2f %10.2f\n', name, ft.NumCalls, ft.TotalTime, ft.(self_time_field));
    end
end

% Save profile data for detailed inspection in GUI if needed
save('profile_data.mat', 'p');
fprintf('\nProfile data saved to profile_data.mat\n');
fprintf('To view in GUI: load profile_data.mat; profview(0, p);\n');
