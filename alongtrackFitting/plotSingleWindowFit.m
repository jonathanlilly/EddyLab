% Function to visualize a specific window fit
function plotSingleWindowFit(alongtrack, eddyPath_fun_t, paramsCell, trueParamsCell, window_index)
j=window_index;
    if iscell(paramsCell)
        trueParams=trueParamsCell{j};
        params=paramsCell{j};
        t0=paramsCell{1}.t0;
        current_window_t0=paramsCell{j}.t0;
        window_start_day=current_window_t0-t0;
    else
        trueParams=trueParamsCell;
        params=paramsCell;
        window_start_day=0;
    end

    % Extract window days
    window_days = unique(params.elapsed_time)+window_start_day;
    % window_center=round(mean(window_days));
    % Extract alongtrack data for this window
    [alongtrack_window, ~] = extractAlongtrackWindow(alongtrack, max(window_days)+t0, max(window_days)+t0);
    
    % Create figure
    figure('Position', [100, 100, 900, 400]);
    
    % Plot SSH data
    % subplot(1, 2, 1);
    scatter(alongtrack_window.x/1e3, alongtrack_window.y/1e3, 30, alongtrack_window.ssh*100, 'filled');
    cmap = brewermap(256, '-Spectral');
    colormap(cmap)
    colorbar;
    hold on;
    
    % Add true and fitted trajectories
    xo_true = eddyPath_fun_t.xe(window_days);
    yo_true = eddyPath_fun_t.ye(window_days);
    
    % Generate time points for fitted trajectory
    t_unique = 0:1:(max(window_days)-min(window_days));
    xo_fit = params.x0 + params.cx * t_unique;
    yo_fit = params.y0 + params.cy * t_unique;
    
    % Plot trajectories
    plot(xo_true/1e3, yo_true/1e3, 'r-', 'LineWidth', 2);
    plot(xo_fit/1e3, yo_fit/1e3, 'b--', 'LineWidth', 2);
    
    % Add circles at characteristic radius
    th = 0:pi/50:2*pi;
    
    % Circle at initial position (true)
    plot(xo_true(end)/1e3 + mean(trueParams.L)/1e3*cos(th), yo_true(end)/1e3 + mean(trueParams.L)/1e3*sin(th), 'r-', 'LineWidth', 1.5);
    
    % Circle at initial position (fitted)
    plot(xo_fit(end)/1e3 + params.L/1e3*cos(th), yo_fit(end)/1e3 + params.L/1e3*sin(th), 'b--', 'LineWidth', 1.5);
    
    % Mark centers
    plot(xo_true(1)/1e3, yo_true(1)/1e3, 'r*', 'MarkerSize', 10);
    plot(params.x0/1e3, params.y0/1e3, 'b*', 'MarkerSize', 10);
    
    grid on;
    xlabel('x (km)', 'FontName', 'times', 'FontSize', 12);
    ylabel('y (km)', 'FontName', 'times', 'FontSize', 12);
    % title('SSH Data with Fitted Model', 'FontName', 'times', 'FontSize', 14);
    lg=legend('SSH Data', 'True Path', 'Fitted Path', 'True Radius', 'Fitted Radius', 'True Center', 'Fitted Center','location','eastoutside');
    axis equal;
    
    
    % % Plot parameter comparison
    % subplot(1, 2, 2);
    % 
    % % Define parameters to compare
    % param_names = {'A (m)', 'L (km)', 'x₀ (km)', 'y₀ (km)', 'c_x (m/s)', 'c_y (m/s)'};
    % fitted_vals = [params.A, params.L/1e3, params.x0/1e3, params.y0/1e3, params.cx, params.cy];
    % true_vals = [mean(trueParams.A), mean(trueParams.L)/1e3, xo_true(1)/1e3, yo_true(1)/1e3, trueParams.cx, trueParams.cy];
    % 
    % % Calculate percent differences
    % percent_diff = 100 * (fitted_vals - true_vals) ./ true_vals;
    % 
    % % Create bar chart
    % bar([fitted_vals; true_vals]', 'grouped');
    % grid on;
    % set(gca, 'XTickLabel', param_names, 'XTickLabelRotation', 45);
    % legend('Fitted', 'True', 'Location', 'best');
    % ylabel('Parameter Value', 'FontName', 'times', 'FontSize', 12);
    % title('Parameter Comparison', 'FontName', 'times', 'FontSize', 14);
    % 
    % % Add text with percent differences
    % for i = 1:length(param_names)
    %     text(i, max(fitted_vals(i), true_vals(i))*1.1, sprintf('%.1f%%', percent_diff(i)), ...
    %          'HorizontalAlignment', 'center', 'FontName', 'times', 'FontSize', 10);
    % end
    % 
    % % Adjust y-limits to accommodate text
    % ylim([min(min(fitted_vals), min(true_vals))*0.9, max(max(fitted_vals), max(true_vals))*1.3]);
    
    % Add overall title
    sgtitle(sprintf('Window %d Analysis (Days %d-%d)', window_index, min(window_days), max(window_days)), ...
            'FontName', 'times', 'FontSize', 16);
end