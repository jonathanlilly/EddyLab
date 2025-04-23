function [mse,min_MSE] = compute_min_MSE_center(alongtrack, eddyPath, fullfield, window_size_array, model, options)
arguments
    alongtrack struct
    eddyPath struct
    fullfield struct
    window_size_array (:,1) {mustBeNumeric}
    model string
    options.bin_size = 12.5* 1e3; % in meters
    options.max_r = 250 *1e3 % in meters
    options.showplot = false; %control whether to display plots
end

% ensure that the arrays are in ascending order in time before windowing
[alongtrack.t,sort_idx]=sort(alongtrack.t,'ascend');
alongtrack.x=alongtrack.x(sort_idx);
alongtrack.y=alongtrack.y(sort_idx);
alongtrack.ssh=alongtrack.ssh(sort_idx);

totalDays = max(alongtrack.t)-min(alongtrack.t);

[fullfield.t,sort_idx]=sort(fullfield.t,'ascend');
fullfield.x=fullfield.x(sort_idx);
fullfield.y=fullfield.y(sort_idx);
fullfield.ssh=fullfield.ssh(sort_idx);

T = length(window_size_array); %number of time variation
it_options = optimset('TolX',1e-3,'TolFun',1e-3);

% Calculate the center of the time series
center_time = min(alongtrack.t) + totalDays / 2;

% Process each window size
for j = 1:T
    window_days = window_size_array(j);
    
    % Calculate the start and end times for this window centered on the eddy
    window_start_day = center_time - window_days / 2;
    window_end_day = center_time + window_days / 2;

    % Ensure the window does not exceed the bounds of the time series
    if window_start_day < min(alongtrack.t)
        window_start_day = min(alongtrack.t);
    end
    if window_end_day > max(alongtrack.t)
        window_end_day = max(alongtrack.t);
    end

    % Find indices that correspond to times within this window
    fullfield_window_indices = find(fullfield.t >= window_start_day & fullfield.t <= window_end_day);

    % Extract time window data
    fullfield_window.x = fullfield.x(fullfield_window_indices);
    fullfield_window.y = fullfield.y(fullfield_window_indices);
    fullfield_window.t = fullfield.t(fullfield_window_indices);
    fullfield_window.ssh = fullfield.ssh(fullfield_window_indices);

    % Compute the composite for the current window size
    eddyPath_window.xe = eddyPath.xe(fullfield_window.t - min(fullfield.t) + 1);
    eddyPath_window.ye = eddyPath.ye(fullfield_window.t - min(fullfield.t) + 1);
    
    % Time-averaged eddy composite from full field
    [mz_true, xmid_true, ymid_true, ~, ~] = composite2D(fullfield_window, eddyPath_window, bin_size=options.bin_size,showplot=0);
    
    % Define spatial window
    spatial_window = [fliplr(-options.bin_size/2:-options.bin_size:-options.max_r), options.bin_size/2:options.bin_size:options.max_r]';
    
    % Extract values within spatial window
    mask_x = xmid_true >= min(spatial_window) & xmid_true <= max(spatial_window);
    mask_y = ymid_true >= min(spatial_window) & ymid_true <= max(spatial_window);
    ssh_true = mz_true(mask_x, mask_y);
    

    % Compute the SSH model based on the selected model type
    switch model
        case 'composite'
        
        alongtrack_window_indices = find(alongtrack.t >= window_start_day & alongtrack.t <= window_end_day);
        alongtrack_window.x = alongtrack.x(alongtrack_window_indices);
        alongtrack_window.y = alongtrack.y(alongtrack_window_indices);
        alongtrack_window.t = alongtrack.t(alongtrack_window_indices);
        alongtrack_window.ssh = alongtrack.ssh(alongtrack_window_indices);  
        % Compute eddy path for alongtrack data
        eddyPath_window_alongtrack.xe = eddyPath.xe(alongtrack_window.t - min(alongtrack.t) + 1);
        eddyPath_window_alongtrack.ye = eddyPath.ye(alongtrack_window.t - min(alongtrack.t) + 1);
        
        % Time-averaged eddy composite from alongtrack
        [mz_xy, xmid_xy, ymid_xy, ~, ~] = composite2D(alongtrack_window, eddyPath_window_alongtrack, bin_size=options.bin_size,showplot=0);
        
        % Extract model values within same spatial window
        mask_x_model = xmid_xy >= min(spatial_window) & xmid_xy <= max(spatial_window);
        mask_y_model = ymid_xy >= min(spatial_window) & ymid_xy <= max(spatial_window);
        ssh_model = mz_xy(mask_x_model, mask_y_model);
        
        % Calculate MSE for this window size
        mse(j) = mean((ssh_true(:) - ssh_model(:)).^2, 'omitnan');

        % Calculate percentage of filled cells within the window
        total_cells = length(spatial_window)^2;
        filled_cells = sum(~isnan(ssh_model(:)));
        coverage_percentage(j) = filled_cells / total_cells;
            
        case 'Gaussian'
            % % Fit a Gaussian model
            % [paramsCell, initParamsCell] = FitAlongTrackXYToEddyModelWindowed(alongtrack, eddyFit_fun, initParams, eddyPath, it_options,window=window_days);
            % for i=1:totalTimeWindows
            % % Calculate the start and end times for this window in days
            % window_start_day = min(alongtrack.t) + (i-1)*time_step;
            % window_end_day = window_start_day + window_days;
            % 
            % % Find indices that correspond to times within this window
            % window_indices = find(alongtrack.t >= window_start_day & alongtrack.t <= window_end_day);
            % 
            % % Extract time window
            % alongtrack_window.x = alongtrack.x(window_indices);
            % alongtrack_window.y = alongtrack.y(window_indices);
            % alongtrack_window.t = alongtrack.t(window_indices);
            % alongtrack_window.ssh = alongtrack.ssh(window_indices);
            % 
            % eddy_model = analyticalEddyModel(eddyPath,paramsCell{i});
            % ssh_model = eddy_model(alongtrack_window.x,alongtrack_window.y,alongtrack_window.t);        
            % % MSE per time window
            % mse_t(i) = mean((ssh_true(:,i) - ssh_model).^2);
            % end
        case 'Elliptical'
            % Fit an elliptical model
            % [paramsCell, initParamsCell] = FitAlongTrackXYToEddyModelWindowed(alongtrack_window, eddyFit_fun, initParams, eddyPath_fun_t, it_options);
        otherwise
            error('Invalid model type');
    end
    % ssh_true_interp = interp2(ssh_true, xmid, ymid); % not necessary
    % because everything will be in the bin size resolution
    % Compute the MSE between the true SSH and the SSH model
end
    figure
    plot(window_size_array,mse*1e2,'LineWidth',2);
    xlabel('Window size', 'FontName', 'times');ylabel('MSE (cm^2)', 'FontName', 'times')
    set(gca, 'fontname', 'times','fontsize', 16)
    % Store the absolute minimum MSE for the current model and time window
    min_MSE = min(mse);
    figure
    plot(window_size_array,coverage_percentage,'LineWidth',2);
end