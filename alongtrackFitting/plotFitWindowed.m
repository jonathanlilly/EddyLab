function plotFitWindowed(paramsCell, trueParamsCell, eddyPath_fun_t)
    totalTimeWindows = length(paramsCell);
    
    % Get window size from first window
    if ~isempty(paramsCell) && isfield(paramsCell{1}, 'elapsed_time')
        window_size = max(unique(paramsCell{1}.elapsed_time)) + 1;
    else
        window_size = 10; % Default window size if not determinable
    end
    
    % Define parameter labels and variable names
    param_label = {'A (m)','L (km)','x_o (km)','y_o (km)','v_x (m/s)','v_y (m/s)'};
    param_var = {'A','L','x0','y0','cx','cy'};
    scale_factors = [1e2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3]; % For display units: A in m, L/x0/y0 in km
    
    % Extract values from each window
    paramValues = zeros(totalTimeWindows, 6);
    trueParamValues = zeros(totalTimeWindows, 6);
    
    for j = 1:totalTimeWindows
        params = paramsCell{j};
        trueParams = trueParamsCell{j};
        t0=paramsCell{1}.t0;
        current_window_t0=paramsCell{j}.t0;
        window_start_day=current_window_t0-t0;
        
        % Get unique days in this window
        % window_days = unique(params.elapsed_time + params.t0);
        window_days = unique(params.elapsed_time)+window_start_day;
        % Get true x0, y0 from the eddy path at the start of this window
        x0_true = eddyPath_fun_t.xe(min(window_days));
        y0_true = eddyPath_fun_t.ye(min(window_days));
        
        % Store parameter values
        paramValues(j, 1) = params.A;
        paramValues(j, 2) = params.L;
        paramValues(j, 3) = params.x0;
        paramValues(j, 4) = params.y0;
        paramValues(j, 5) = params.cx;
        paramValues(j, 6) = params.cy;
        
        % Store true parameter values
        trueParamValues(j, 1) = mean(trueParams.A);
        trueParamValues(j, 2) = mean(trueParams.L);
        trueParamValues(j, 3) = x0_true;
        trueParamValues(j, 4) = y0_true;
        trueParamValues(j, 5) = trueParams.cx;
        trueParamValues(j, 6) = trueParams.cy;
    end
    
    % Create figure
    figure('Position', [100, 100, 800, 800]);
    
    % Plot each parameter
    for i = 1:6
        subplot(6, 1, i)
        hold on
        
        % Plot fitted values
        plot(1:totalTimeWindows, paramValues(:, i).*scale_factors(i), 'k.-', 'MarkerSize', 15, 'LineWidth', 1.5);
        
        % Plot true values
        plot(1:totalTimeWindows, trueParamValues(:, i).*scale_factors(i), 'r.:', 'LineWidth', 1.5);
        
        % Labels and formatting
        ylabel(param_label{i}, 'FontName', 'times', 'FontSize', 12);
        grid on
        
        % Only add x-label to bottom plot
        if i == 6
            xlabel(sprintf('Window Number (window size = %d days)', window_size), 'FontName', 'times', 'FontSize', 12);
            legend('Fitted', 'True', 'Location', 'best');
        else
            set(gca, 'XTickLabel', []);
        end
        
        xlim([0.5, totalTimeWindows+0.5]);
        
        % Adjust spacing between subplots
        pos = get(gca, 'Position');
        set(gca, 'Position', [pos(1), pos(2)+0.01*(6-i), pos(3), pos(4)]);
    end
    
    % Add overall title
    sgtitle(sprintf('Parameter Evolution Across %d Windows', totalTimeWindows), 'FontName', 'times', 'FontSize', 14);
end
