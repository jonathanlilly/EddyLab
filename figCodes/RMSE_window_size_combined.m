%% Comprehensive RMSE Analysis for Multiple Models and Datasets
window_size_array = [366, 100, 50, 25, 10, 5, 1];

% Define all combinations to test
test_cases = {
    {'composite', 'alongtrack', 'Alongtrack_{composite}'};
    {'Gaussian', 'alongtrack', 'Alongtrack_{Gaussian}'};
    {'composite', 'multiMission', 'Multi_{composite}'};
    {'Gaussian', 'multiMission', 'Multi_{Gaussian}'};
    {'composite', 'fullfield', 'Truth'}
    {'composite', 'j3nMission', 'j3n_{composite}'};
    {'Gaussian', 'j3nMission', 'j3n_{Gaussian}'};
    %{'elliptical', 'alongtrack', 'Alongtrack_elliptical'};
};

% Initialize storage
n_cases = length(test_cases);
n_windows = length(window_size_array);
rmse_results = zeros(n_cases, n_windows);
rmse_daily_results = cell(n_cases, 1);
getModelForDay_results = cell(n_cases, n_windows);
window_center_results = cell(n_cases, n_windows);
case_names = cell(n_cases, 1);

% Compute RMSE for all cases
fprintf('Computing RMSE for all model-data combinations...\n');
%%
for case_idx = 1:n_cases
    model_type = test_cases{case_idx}{1};
    data_type = test_cases{case_idx}{2};
    case_name = test_cases{case_idx}{3};
    case_names{case_idx} = case_name;
    
    fprintf('Processing %s (model: %s, data: %s)\n', case_name, model_type, data_type);
    
    % Select appropriate dataset
    switch data_type
        case 'alongtrack'
            current_data = alongtrack;
        case 'multiMission'
            current_data = multiMission;
        case 'fullfield'
            current_data = fullfield;
        case 'j3nMission'
            current_data = j3n_mission;
        otherwise
            error('Unknown data type: %s', data_type);
    end
    
    % Initialize arrays for this case
    rmse_array_case = zeros(1, n_windows);
    rmse_daily_case = zeros(totalDays, n_windows);
    
    % Loop through window sizes
    for w = 1:n_windows
        window_size = window_size_array(w);
        
        try
            [rmse_daily, rmse, getModelForDay, window_center] = compute_model_error(current_data, eddy_field, ...
                eddyPath_fun_t, window_size, model_type);
            
            rmse_array_case(w) = rmse;
            rmse_daily_case(:, w) = rmse_daily;
            
            % Store the model functions and window centers
            getModelForDay_results{case_idx, w} = getModelForDay;
            window_center_results{case_idx, w} = window_center;
            
            fprintf('  Window size %d days: RMSE = %.3f cm\n', window_size, rmse*1e2);
            
        catch ME
            warning('Failed for %s with window size %d: %s', case_name, window_size, ME.message);
            rmse_array_case(w) = NaN;
            rmse_daily_case(:, w) = NaN;
            getModelForDay_results{case_idx, w} = [];
            window_center_results{case_idx, w} = [];
        end
    end
    
    % Store results
    rmse_results(case_idx, :) = rmse_array_case;
    rmse_daily_results{case_idx} = rmse_daily_case;
end

%% figure function that plots in test_case order
function createfigure_RMSE_windowSize(window_size_array, rmse_results, case_names)
% Create figure
figure('Position', [100, 100, 900, 600]);

% Create axes
axes1 = axes;
hold(axes1,'on');

% Plot each case in order with specific styling rules
plot_handles = [];
for case_idx = 1:length(case_names)
    case_name = case_names{case_idx};
    
    % Determine color based on model type
    if contains(case_name, 'Gauss')
        color = [0.96078431372549 0.466666666666667 0.16078431372549]; % Orange
    elseif contains(case_name, 'composite')
        color = [0.0666666666666667 0.443137254901961 0.745098039215686]; % Blue
    elseif contains(case_name, 'Truth')
        color = [0.501960784313725 0.501960784313725 0.501960784313725]; % Gray
    else
        color = [0, 0, 0]; % Black for any other case
    end
    
    % Determine line style based on data type
    if contains(case_name, 'Alongtrack')
        line_style = '--'; % Dashed for alongtrack
    elseif contains(case_name, 'Multi') 
        line_style = '-';  % Solid for multi-mission and j3n
    elseif contains(case_name, 'Truth')
        line_style = ':';  % Dotted for Truth
    else
        line_style = '-';  % Default solid
    end
    
    % Determine marker
    if contains(case_name, 'Truth')
        marker = 'none';
        marker_face_color = 'none';
    else
        marker = 'o';
        marker_face_color = color;
    end
    
    h = loglog(window_size_array, rmse_results(case_idx, :) * 1e2, ...
        'DisplayName', case_name, ...
        'Color', color, ...
        'LineStyle', line_style, ...
        'Marker', marker, ...
        'MarkerFaceColor', marker_face_color, ...
        'MarkerSize', 8, ...
        'LineWidth', 2);
    
    plot_handles = [plot_handles, h];
end

% Create ylabel and xlabel
ylabel('RMSE (cm)', 'FontName', 'times', 'FontSize', 16);
xlabel('Time window size (days)', 'FontName', 'times', 'FontSize', 16);

xlim([1,366])

box(axes1, 'on');
hold(axes1, 'off');

% Set axes properties
set(axes1, 'FontName', 'times', 'FontSize', 16, 'XMinorTick', 'on', 'XScale', 'log', ...
    'YMinorTick', 'on', 'YScale', 'log');

% Create legend in the same order as test cases
legend1 = legend(plot_handles, case_names, 'Location', 'northeast', ...
    'FontName', 'times', 'FontSize', 12, 'NumColumns', 2,'orientation','horizontal');

grid(axes1, 'off');
end

%% Call the updated figure function
createfigure_RMSE_windowSize(window_size_array, rmse_results, case_names);

%% Save results
save('RMSE_comprehensive_analysis.mat', 'rmse_results', 'rmse_daily_results', ...
     'getModelForDay_results', 'window_center_results', ...
     'window_size_array', 'case_names', 'test_cases');

%% Display summary table
fprintf('\n=== RMSE Summary Table (cm) ===\n');
fprintf('Window Size: ');
fprintf('%8d ', window_size_array);
fprintf('\n');
for i = 1:n_cases
    fprintf('%-20s: ', case_names{i});
    fprintf('%8.2f ', rmse_results(i, :) * 1e2);
    fprintf('\n');
end

%% Display information about saved model functions
fprintf('\n=== Saved Model Functions Summary ===\n');
fprintf('Saved getModelForDay functions and window_center data for:\n');
for i = 1:n_cases
    fprintf('  %s: %d window sizes\n', case_names{i}, n_windows);
end
fprintf('\nAccess examples:\n');
fprintf('  Model function for case 1, window size 1: getModelForDay_results{1,1}\n');
fprintf('  Window centers for case 2, window size 3: window_center_results{2,3}\n');

%% Optional: Create individual time series plots for each case
figure('Position', [100, 100, 1200, 800]);
colors = lines(n_cases);

for case_idx = 1:n_cases
    subplot(2, 3, case_idx);
    hold on;
    
    rmse_daily_case = rmse_daily_results{case_idx};
    
    for w = 1:min(6, n_windows)  % Plot up to 6 window sizes
        plot([1:totalDays]-1, rmse_daily_case(:, w)*1e2, ...
             'LineWidth', 2, 'Color', colors(w, :));
    end
    
    xlabel('Time (day)', 'FontName', 'times');
    ylabel('RMSE (cm)', 'FontName', 'times');
    title(case_names{case_idx}, 'FontName', 'times', 'FontSize', 14);
    set(gca, 'FontName', 'times', 'FontSize', 12);
    xlim([0, totalDays-1]);
    
    if case_idx == 1
        legend_labels = arrayfun(@(x) sprintf('$W_t=%d$ days', x), ...
                               window_size_array(1:min(6, n_windows)), 'UniformOutput', false);
        legend(legend_labels, 'Location', 'best', 'FontSize', 16,'interpreter','latex');
    end
end

sgtitle('RMSE Time Series for All Model-Data Combinations', ...
        'FontName', 'times', 'FontSize', 16);