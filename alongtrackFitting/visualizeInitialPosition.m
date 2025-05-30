function visualizeInitialPosition(eddyPath_fun_t, eddyParams, paramsCell, trueParamsCell)
    % Create figure
    figure('Position', [100, 100, 1200, 600]);
    
    totalTimeWindows = length(paramsCell);
    colors = jet(totalTimeWindows);
    
    % Panel 1: Full trajectory with all initial positions markefunction plotEnhancedEddyFit(alongtrack, eddyPath_fun_t, paramsCell, trueParamsCell)
    % Create a multi-panel figure to better visualize eddy fit results
    totalTimeWindows = length(paramsCell);
    
    % Create figure with appropriate size
    figure('Position', [100, 100, 1200, 900]);
    
    % Panel 1: Spatial view with all windows
    subplot(2, 2, 1);
    hold on;
    
    % Define color progression for windows
    colors = jet(totalTimeWindows);
    
    % Initialize arrays to store errors for later plots
    position_errors = zeros(1, totalTimeWindows);
    amplitude_errors = zeros(1, totalTimeWindows);
    scale_errors = zeros(1, totalTimeWindows);
    velocity_errors = zeros(1, totalTimeWindows);
    window_centers = zeros(1, totalTimeWindows);
    
    % Plot spatial representation with eddy trajectories
    for j = 1:totalTimeWindows
        params = paramsCell{j};
        trueParams = trueParamsCell{j};
        
        % Get unique days in window
        window_days = unique(params.elapsed_time + params.t0);
        
        % Get window center time for x-axis in time series plots
        window_centers(j) = params.t0 + mean(unique(params.elapsed_time));
        
        % True positions through window - directly from eddyPath_fun_t
        xo_true_full = eddyPath_fun_t.xe(window_days);
        yo_true_full = eddyPath_fun_t.ye(window_days);
        
        % True amplitudes and length scales - assuming these match window days
        A_true_full = eddyParams.A(1:length(window_days));
        L_true_full = eddyParams.L(1:length(window_days));
        
        % Generate time points for model trajectory
        t_unique = unique(params.elapsed_time);
        
        % Fit positions through window - using fixed initial position and velocity
        % This respects your model where x0 is the position at t=0 in the window
        xo_fit_full = params.x0 + params.cx * t_unique;
        yo_fit_full = params.y0 + params.cy * t_unique;
        
        % Plot trajectories
        plot(xo_true_full/1e3, yo_true_full/1e3, 'Color', colors(j,:), 'LineWidth', 1.5);
        plot(xo_fit_full/1e3, yo_fit_full/1e3, '--', 'Color', colors(j,:), 'LineWidth', 1.5);
        
        % Plot start and end points of trajectories
        plot(xo_true_full(1)/1e3, yo_true_full(1)/1e3, 'o', 'Color', colors(j,:), 'MarkerFaceColor', colors(j,:));
        plot(xo_true_full(end)/1e3, yo_true_full(end)/1e3, '*', 'Color', colors(j,:), 'MarkerSize', 10);
        
        % Calculate errors for time series
        position_errors(j) = sqrt(mean((xo_fit_full - xo_true_full).^2 + (yo_fit_full - yo_true_full).^2))/1e3;
        amplitude_errors(j) = abs(params.A - mean(A_true_full));
        scale_errors(j) = abs(params.L - mean(L_true_full))/1e3;
        velocity_errors(j) = sqrt((params.cx - trueParams.cx)^2 + (params.cy - trueParams.cy)^2)/1e3;
    end
    
    grid on;
    xlabel('x (km)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('y (km)', 'FontName', 'Times', 'FontSize', 12);
    title('Eddy Trajectories', 'FontName', 'Times', 'FontSize', 14);
    legend('True', 'Fitted', 'Location', 'best');
    axis equal;
    
    % Panel 2: Error metrics time series
    subplot(2, 2, 2);
    hold on;
    
    % Add text explaining initial position definition
    text(0.05, 0.95, 'Note: x0/y0 represent position at t=0 of each window', 'Units', 'normalized', ...
         'FontName', 'Times', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
    
    % Plot position errors over time
    yyaxis left;
    plot(window_centers, position_errors, 'o-', 'LineWidth', 1.5);
    ylabel('Position Error (km)', 'FontName', 'Times', 'FontSize', 12);
    
    % Plot amplitude errors on secondary y-axis
    yyaxis right;
    plot(window_centers, amplitude_errors*100, 's-', 'LineWidth', 1.5);
    ylabel('Amplitude Error (cm)', 'FontName', 'Times', 'FontSize', 12);
    
    grid on;
    title('Fit Errors Over Time', 'FontName', 'Times', 'FontSize', 14);
    xlabel('Time (days)', 'FontName', 'Times', 'FontSize', 12);
    legend('Position Error', 'Amplitude Error', 'Location', 'best');
    
    % Panel 3: Amplitude and Length Scale Comparison
    subplot(2, 2, 3);
    hold on;
    
    % Initialize arrays for plotting parameters vs time
    all_times = [];
    all_A_true = [];
    all_L_true = [];
    all_window_indices = [];
    
    % Collect all parameter values across windows
    for j = 1:totalTimeWindows
        params = paramsCell{j};
        trueParams = trueParamsCell{j};
        
        % Get days in this window
        window_days = unique(params.elapsed_time + params.t0);
        
        % Extract parameters
        A_true_values = eddyParams.A(1:length(window_days));
        L_true_values = eddyParams.L(1:length(window_days));
        
        % Store window index for each time point for coloring
        window_idx = j * ones(size(window_days));
        
        % Accumulate data
        all_times = [all_times, window_days];
        all_A_true = [all_A_true, A_true_values];
        all_L_true = [all_L_true, L_true_values];
        all_window_indices = [all_window_indices, window_idx];
        
        % Plot fitted values as horizontal lines spanning window
        plot([min(window_days), max(window_days)], [params.A*100, params.A*100], '-', 'Color', colors(j,:), 'LineWidth', 2);
        plot([min(window_days), max(window_days)], [params.L/1e3, params.L/1e3], '--', 'Color', colors(j,:), 'LineWidth', 2);
    end
    
    % Plot true parameter values
    scatter(all_times, all_A_true*100, 30, all_window_indices, 'o', 'filled');
    scatter(all_times, all_L_true/1e3, 30, all_window_indices, 's');
    
    grid on;
    xlabel('Time (days)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('Parameter Value', 'FontName', 'Times', 'FontSize', 12);
    title('Amplitude (cm) and Length Scale (km) Comparison', 'FontName', 'Times', 'FontSize', 14);
    legend('Fitted Amplitude', 'Fitted Length Scale', 'True Amplitude', 'True Length Scale', 'Location', 'best');
    
    % Panel 4: Error distribution with colormap
    subplot(2, 2, 4);
    
    % Create a 2D view of parameter space with errors
    scatter(scale_errors, amplitude_errors*100, 100, velocity_errors, 'filled');
    colorbar;
    
    xlabel('Length Scale Error (km)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('Amplitude Error (cm)', 'FontName', 'Times', 'FontSize', 12);
    title('Parameter Error Space (color = velocity error)', 'FontName', 'Times', 'FontSize', 14);
    grid on;
    
    % Add overall title
    sgtitle('Eddy Gaussian Fit Analysis', 'FontName', 'Times', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Adjust spacing
    set(gcf, 'Color', 'w');
    colormap(jet);
    
    % Create second figure for detailed SSH visualization
    figure('Position', [100, 100, 1200, 500]);
    
    % Choose a representative window to visualize in detail
    j = min(3, totalTimeWindows); % Third window or last if fewer
    params = paramsCell{j};
    trueParams = trueParamsCell{j};
    
    % Extract window data
    window_days = params.elapsed_time + params.t0;
    first_day = window_days(1);
    last_day = window_days(end);
    [alongtrack_window, ~] = extractAlongtrackWindow(alongtrack, first_day, last_day);
    
    % Panel 1: SSH data with fitted contours (first day)
    subplot(1, 2, 1);
    
    % Get or calculate SSH grid for first day
    [alongtrack_day, ~] = extractAlongtrackWindow(alongtrack, first_day, first_day);
    
    % Plot SSH data
    scatter(alongtrack_day.x/1e3, alongtrack_day.y/1e3, 30, alongtrack_day.ssh*100, 'filled');
    colorbar;
    hold on;
    
    % Add eddy contours from model (true and fitted)
    th = 0:pi/50:2*pi;
    
    % True contours (first day)
    xo_true = eddyPath_fun_t.xe(first_day);
    yo_true = eddyPath_fun_t.ye(first_day);
    A_true = trueParams.A(1);
    L_true = trueParams.L(1);
    plot(xo_true/1e3 + L_true/1e3*cos(th), yo_true/1e3 + L_true/1e3*sin(th), 'r-', 'LineWidth', 2);
    
    % Fitted contours (first day)
    xo_fit = params.x0;
    yo_fit = params.y0;
    L_fit = params.L;
    plot(xo_fit/1e3 + L_fit/1e3*cos(th), yo_fit/1e3 + L_fit/1e3*sin(th), 'b--', 'LineWidth', 2);
    
    % Plot centers
    plot(xo_true/1e3, yo_true/1e3, 'r*', 'MarkerSize', 10);
    plot(xo_fit/1e3, yo_fit/1e3, 'b*', 'MarkerSize', 10);
    
    grid on;
    xlabel('x (km)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('y (km)', 'FontName', 'Times', 'FontSize', 12);
    title(['SSH and Eddy Fit (Day ', num2str(first_day), ')'], 'FontName', 'Times', 'FontSize', 14);
    axis equal;
    colormap(brewermap([], '-Spectral'));
    
    % Panel 2: SSH data with fitted contours (last day)
    subplot(1, 2, 2);
    
    % Get or calculate SSH grid for last day
    [alongtrack_day, ~] = extractAlongtrackWindow(alongtrack, last_day, last_day);
    
    % Plot SSH data
    scatter(alongtrack_day.x/1e3, alongtrack_day.y/1e3, 30, alongtrack_day.ssh*100, 'filled');
    colorbar;
    hold on;
    
    % True contours (last day)
    xo_true = eddyPath_fun_t.xe(last_day);
    yo_true = eddyPath_fun_t.ye(last_day);
    A_true = trueParams.A(end);
    L_true = trueParams.L(end);
    plot(xo_true/1e3 + L_true/1e3*cos(th), yo_true/1e3 + L_true/1e3*sin(th), 'r-', 'LineWidth', 2);
    
    % Fitted contours (last day)
    t_elapsed = params.elapsed_time(end);
    xo_fit = params.x0 + params.cx * t_elapsed;
    yo_fit = params.y0 + params.cy * t_elapsed;
    L_fit = params.L;
    plot(xo_fit/1e3 + L_fit/1e3*cos(th), yo_fit/1e3 + L_fit/1e3*sin(th), 'b--', 'LineWidth', 2);
    
    % Plot centers
    plot(xo_true/1e3, yo_true/1e3, 'r*', 'MarkerSize', 10);
    plot(xo_fit/1e3, yo_fit/1e3, 'b*', 'MarkerSize', 10);
    
    grid on;
    xlabel('x (km)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('y (km)', 'FontName', 'Times', 'FontSize', 12);
    title(['SSH and Eddy Fit (Day ', num2str(last_day), ')'], 'FontName', 'Times', 'FontSize', 14);
    axis equal;
    colormap(brewermap([], '-Spectral'));
    
    % Add overall title
    sgtitle('SSH Data with Fitted Eddy Model', 'FontName', 'Times', 'FontSize', 16, 'FontWeight', 'bold');
end

% Function to visualize evolution of model parameters over time
function plotParameterEvolution(paramsCell, trueParamsCell)
    % Create figure
    figure('Position', [100, 100, 1200, 800]);
    
    totalTimeWindows = length(paramsCell);
    
    % Extract parameter values and times
    window_centers = zeros(1, totalTimeWindows);
    A_fit = zeros(1, totalTimeWindows);
    L_fit = zeros(1, totalTimeWindows);
    cx_fit = zeros(1, totalTimeWindows);
    cy_fit = zeros(1, totalTimeWindows);
    x0_fit = zeros(1, totalTimeWindows);  % Initial positions
    y0_fit = zeros(1, totalTimeWindows);
    
    A_true_mean = zeros(1, totalTimeWindows);
    L_true_mean = zeros(1, totalTimeWindows);
    cx_true = zeros(1, totalTimeWindows);
    cy_true = zeros(1, totalTimeWindows);
    x0_true = zeros(1, totalTimeWindows);  % Initial positions (from trueParams)
    y0_true = zeros(1, totalTimeWindows);
    
    for j = 1:totalTimeWindows
        params = paramsCell{j};
        trueParams = trueParamsCell{j};
        
        % Get window center time
        window_centers(j) = params.t0 + mean(params.elapsed_time);
        
        % Fitted parameters
        A_fit(j) = params.A;
        L_fit(j) = params.L;
        cx_fit(j) = params.cx;
        cy_fit(j) = params.cy;
        x0_fit(j) = params.x0;  % Initial position (t=0 of window)
        y0_fit(j) = params.y0;
        
        % True parameters (means over window)
        A_true_mean(j) = mean(trueParams.A);
        L_true_mean(j) = mean(trueParams.L);
        cx_true(j) = trueParams.cx;
        cy_true(j) = trueParams.cy;
        x0_true(j) = trueParams.x0;  % Initial position from trueParams
        y0_true(j) = trueParams.y0;
    end
    
    % Plot amplitude evolution
    subplot(2, 2, 1);
    plot(window_centers, A_fit*100, 'bo-', 'LineWidth', 1.5);
    hold on;
    plot(window_centers, A_true_mean*100, 'ro-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (days)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('Amplitude (cm)', 'FontName', 'Times', 'FontSize', 12);
    title('Amplitude Evolution', 'FontName', 'Times', 'FontSize', 14);
    legend('Fitted', 'True (mean)', 'Location', 'best');
    
    % Plot length scale evolution
    subplot(2, 2, 2);
    plot(window_centers, L_fit/1e3, 'bo-', 'LineWidth', 1.5);
    hold on;
    plot(window_centers, L_true_mean/1e3, 'ro-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (days)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('Length Scale (km)', 'FontName', 'Times', 'FontSize', 12);
    title('Length Scale Evolution', 'FontName', 'Times', 'FontSize', 14);
    legend('Fitted', 'True (mean)', 'Location', 'best');
    
    % Plot x-velocity evolution
    subplot(2, 2, 3);
    plot(window_centers, cx_fit*86400, 'bo-', 'LineWidth', 1.5);
    hold on;
    plot(window_centers, cx_true*86400, 'ro-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (days)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('x-velocity (m/day)', 'FontName', 'Times', 'FontSize', 12);
    title('x-velocity Evolution', 'FontName', 'Times', 'FontSize', 14);
    legend('Fitted', 'True', 'Location', 'best');
    
    % Plot y-velocity evolution
    subplot(2, 2, 4);
    plot(window_centers, cy_fit*86400, 'bo-', 'LineWidth', 1.5);
    hold on;
    plot(window_centers, cy_true*86400, 'ro-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (days)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('y-velocity (m/day)', 'FontName', 'Times', 'FontSize', 12);
    title('y-velocity Evolution', 'FontName', 'Times', 'FontSize', 14);
    legend('Fitted', 'True', 'Location', 'best');
    
    % Create a third figure for initial positions comparison
    figure('Position', [100, 100, 800, 400]);
    
    % Plot x0 evolution
    subplot(1, 2, 1);
    plot(window_centers, x0_fit/1e3, 'bo-', 'LineWidth', 1.5);
    hold on;
    plot(window_centers, x0_true/1e3, 'ro-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (days)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('Initial x-position (km)', 'FontName', 'Times', 'FontSize', 12);
    title('x₀ Evolution (First position in window)', 'FontName', 'Times', 'FontSize', 14);
    legend('Fitted', 'True', 'Location', 'best');
    
    % Plot y0 evolution
    subplot(1, 2, 2);
    plot(window_centers, y0_fit/1e3, 'bo-', 'LineWidth', 1.5);
    hold on;
    plot(window_centers, y0_true/1e3, 'ro-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (days)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('Initial y-position (km)', 'FontName', 'Times', 'FontSize', 12);
    title('y₀ Evolution (First position in window)', 'FontName', 'Times', 'FontSize', 14);
    legend('Fitted', 'True', 'Location', 'best');
    
    % Add overall title
    sgtitle('Initial Position Parameters (x₀, y₀)', 'FontName', 'Times', 'FontSize', 16, 'FontWeight', 'bold');
end

% Function to create model-data misfit visualization
function visualizeModelDataMisfit(alongtrack, eddyPath_fun_t, paramsCell, trueParamsCell, window_index)
    % Select window to visualize
    params = paramsCell{window_index};
    trueParams = trueParamsCell{window_index};
    
    % Extract window data
    window_days = unique(params.elapsed_time + params.t0);
    [alongtrack_window, ~] = extractAlongtrackWindow(alongtrack, min(window_days), max(window_days));
    
    % Create figure
    figure('Position', [100, 100, 1200, 800]);
    
    % Compute model SSH using fitted parameters
    ssh_model = zeros(size(alongtrack_window.ssh));
    for i = 1:length(alongtrack_window.ssh)
        % Find elapsed time for this data point
        current_day = alongtrack_window.t(i);
        if current_day >= min(window_days) && current_day <= max(window_days)
            t_elapsed = current_day - params.t0;  % Calculate elapsed time
            ssh_model(i) = params.A * exp(-((alongtrack_window.x(i) - params.x0 - params.cx*t_elapsed)^2 + ...
                              (alongtrack_window.y(i) - params.y0 - params.cy*t_elapsed)^2) / params.L^2);
        end
    end
    
    % Panel 1: Model SSH
    subplot(2, 2, 1);
    scatter(alongtrack_window.x/1e3, alongtrack_window.y/1e3, 30, ssh_model*100, 'filled');
    colorbar;
    hold on;
    
    % Add eddy trajectory
    xo_fit = params.x0 + params.cx * params.elapsed_time;
    yo_fit = params.y0 + params.cy * params.elapsed_time;
    plot(xo_fit/1e3, yo_fit/1e3, 'k-', 'LineWidth', 2);
    plot(xo_fit(1)/1e3, yo_fit(1)/1e3, 'ko', 'MarkerFaceColor', 'k');
    plot(xo_fit(end)/1e3, yo_fit(end)/1e3, 'k*', 'MarkerSize', 10);
    
    grid on;
    xlabel('x (km)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('y (km)', 'FontName', 'Times', 'FontSize', 12);
    title('Model SSH (cm)', 'FontName', 'Times', 'FontSize', 14);
    axis equal;
    colormap(brewermap([], '-Spectral'));
    
    % Panel 2: Data SSH
    subplot(2, 2, 2);
    scatter(alongtrack_window.x/1e3, alongtrack_window.y/1e3, 30, alongtrack_window.ssh*100, 'filled');
    colorbar;
    hold on;
    
    % Add true eddy trajectory - directly from eddyPath_fun_t
    xo_true = eddyPath_fun_t.xe(window_days);
    yo_true = eddyPath_fun_t.ye(window_days);
    plot(xo_true/1e3, yo_true/1e3, 'k-', 'LineWidth', 2);
    plot(xo_true(1)/1e3, yo_true(1)/1e3, 'ko', 'MarkerFaceColor', 'k');
    plot(xo_true(end)/1e3, yo_true(end)/1e3, 'k*', 'MarkerSize', 10);
    
    % Generate unique time points for fitted trajectory
    t_unique = 0:1:(max(window_days)-min(window_days));
    
    % Calculate fitted trajectory
    xo_fit = params.x0 + params.cx * t_unique;
    yo_fit = params.y0 + params.cy * t_unique;
    
    % Add text showing the first position (x0, y0)
    text(0.05, 0.95, ['x₀ true: ', num2str(xo_true(1)/1e3, '%.1f'), ' km, y₀ true: ', num2str(yo_true(1)/1e3, '%.1f'), ' km'], ...
         'Units', 'normalized', 'FontName', 'Times', 'FontSize', 10, 'HorizontalAlignment', 'left');
    text(0.05, 0.90, ['x₀ fit: ', num2str(params.x0/1e3, '%.1f'), ' km, y₀ fit: ', num2str(params.y0/1e3, '%.1f'), ' km'], ...
         'Units', 'normalized', 'FontName', 'Times', 'FontSize', 10, 'HorizontalAlignment', 'left');
    
    grid on;
    xlabel('x (km)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('y (km)', 'FontName', 'Times', 'FontSize', 12);
    title('Data SSH (cm)', 'FontName', 'Times', 'FontSize', 14);
    axis equal;
    colormap(brewermap([], '-Spectral'));
    
    % Panel 3: Misfit (Model - Data)
    subplot(2, 2, 3);
    misfit = ssh_model - alongtrack_window.ssh;
    scatter(alongtrack_window.x/1e3, alongtrack_window.y/1e3, 30, misfit*100, 'filled');
    colorbar;
    hold on;
    
    % Add both trajectories
    plot(xo_true/1e3, yo_true/1e3, 'r-', 'LineWidth', 2);
    plot(xo_fit/1e3, yo_fit/1e3, 'b--', 'LineWidth', 2);
    
    grid on;
    xlabel('x (km)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('y (km)', 'FontName', 'Times', 'FontSize', 12);
    title('Misfit: Model - Data (cm)', 'FontName', 'Times', 'FontSize', 14);
    axis equal;
    colormap(brewermap([], 'RdBu'));
    caxis([-max(abs(misfit*100)), max(abs(misfit*100))]);
    
    % Panel 4: Statistical analysis of misfit
    subplot(2, 2, 4);
    
    % Calculate misfit statistics
    mean_misfit = mean(misfit*100);
    std_misfit = std(misfit*100);
    rms_misfit = sqrt(mean((misfit*100).^2));
    
    % Create histogram of misfit
    histogram(misfit*100, 20, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'k');
    hold on;
    
    % Add vertical lines for mean and std
    xline(mean_misfit, 'r-', 'LineWidth', 2);
    xline(mean_misfit + std_misfit, 'r--', 'LineWidth', 1.5);
    xline(mean_misfit - std_misfit, 'r--', 'LineWidth', 1.5);
    
    grid on;
    xlabel('Misfit (cm)', 'FontName', 'Times', 'FontSize', 12);
    ylabel('Frequency', 'FontName', 'Times', 'FontSize', 12);
    title(['Misfit Statistics: Mean = ', num2str(mean_misfit, '%.2f'), ' cm, RMS = ', num2str(rms_misfit, '%.2f'), ' cm'], ...
          'FontName', 'Times', 'FontSize', 14);
    
    % Add text with misfit statistics
    text(0.05, 0.9, ['Mean: ', num2str(mean_misfit, '%.2f'), ' cm'], 'Units', 'normalized', ...
         'FontName', 'Times', 'FontSize', 12);
    text(0.05, 0.8, ['Std: ', num2str(std_misfit, '%.2f'), ' cm'], 'Units', 'normalized', ...
         'FontName', 'Times', 'FontSize', 12);
    text(0.05, 0.7, ['RMS: ', num2str(rms_misfit, '%.2f'), ' cm'], 'Units', 'normalized', ...
         'FontName', 'Times', 'FontSize', 12);
    text(0.05, 0.6, ['Max Abs: ', num2str(max(abs(misfit*100)), '%.2f'), ' cm'], 'Units', 'normalized', ...
         'FontName', 'Times', 'FontSize', 12);
    
    % Add overall title
    sgtitle(['Model-Data Misfit Analysis (Window ', num2str(window_index), ')'], ...
            'FontName', 'Times', 'FontSize', 16, 'FontWeight', 'bold');
end