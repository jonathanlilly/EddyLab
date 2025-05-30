function [mse,min_MSE] = compute_min_MSE(alongtrack, eddyPath, fullfield, window_size_array, model, options)
arguments
    alongtrack struct
    eddyPath struct
    fullfield struct
    window_size_array (:,1) {mustBeNumeric}
    model string
    options.bin_size = 12.5 * 1e3; % in meters
    options.max_r = 250 *1e3 % in meters
    options.showplot = false; %control whether to display plots
end

overlap = 0; %50% overlap between time_window

% ensure that the arrays are in ascending order in time before windowing
[alongtrack.t,sort_idx]=sort(alongtrack.t,'ascend');
alongtrack.x=alongtrack.x(sort_idx);
alongtrack.y=alongtrack.y(sort_idx);
alongtrack.ssh=alongtrack.ssh(sort_idx);

[fullfield.t,sort_idx]=sort(fullfield.t,'ascend');
fullfield.x=fullfield.x(sort_idx);
fullfield.y=fullfield.y(sort_idx);
fullfield.ssh=fullfield.ssh(sort_idx);

totalDays = max(fullfield.t)-min(fullfield.t);

T = length(window_size_array); %number of time variation
it_options = optimset('TolX',1e-3,'TolFun',1e-3);

for j = 1:T
    window_days = window_size_array(j);
    
    time_step=floor(window_days*(1-overlap));
    totalTimeWindows=floor((totalDays - window_days) / time_step) + 1;

    % True SSH from full field composite at each time window
    % rmid_max = [0:options.bin_size:options.max_r]';
    % ssh_true = nan(length(rmid_max),totalTimeWindows);
    spatial_window = [fliplr(-options.bin_size/2:-options.bin_size:-options.max_r),options.bin_size/2:options.bin_size:options.max_r]';
    ssh_true = nan(length(spatial_window),length(spatial_window),totalTimeWindows);

    for i=totalTimeWindows:-1:1
    % Calculate the start and end times for this window in days
    % window_start_day = min(fullfield.t) + (i-1)*time_step;
    % window_end_day = window_start_day + window_days-1;
    window_end_day = max(fullfield.t)-(i-1)*time_step;
    window_start_day = window_end_day - window_days;


    % Find indices that correspond to times within this window
    window_indices = find(fullfield.t >= window_start_day & fullfield.t <= window_end_day);

    % Extract time window
    fullfield_window.x = fullfield.x(window_indices);
    fullfield_window.y = fullfield.y(window_indices);
    fullfield_window.t = fullfield.t(window_indices);
    fullfield_window.ssh = fullfield.ssh(window_indices);

    % Compute the 2D composite for the current time window
    eddyPath_window.xe=eddyPath.xe(fullfield_window.t-min(fullfield.t));
    eddyPath_window.ye=eddyPath.ye(fullfield_window.t-min(fullfield.t));
    % [mz_true, rmid_true, numz_true, stdz_true] = radialProfile(fullfield_window,eddyPath_window,showplot=0);
    % ssh_true(:,i) = mz_true(1:length(rmid_max));
    % time-averaged eddy composite from full field
    [mz_true, xmid_true, ymid_true, numz_true, stdz_true] = composite2D(fullfield_window,eddyPath_window,showplot=0);
    ssh_true(:,:,i) = mz_true(xmid_true>=min(spatial_window)&xmid_true<=max(spatial_window),ymid_true>=min(spatial_window)&ymid_true<=max(spatial_window));
    end

    % Compute the SSH model based on the selected model type
    switch model
        case 'composite'
        mse_t=zeros(totalTimeWindows,1);
        coverage_t=zeros(totalTimeWindows,1);

        for i=totalTimeWindows:-1:1
        % Calculate the start and end times for this window in days
        % window_start_day = min(alongtrack.t) + (i-1)*time_step;
        % window_end_day = window_start_day + window_days-1;
        window_end_day = max(alongtrack.t)-(i-1)*time_step;
        window_start_day = window_end_day - window_days;
    
        % Find indices that correspond to times within this window
        window_indices = find(alongtrack.t >= window_start_day & alongtrack.t <= window_end_day);
    
        % Extract time window
        alongtrack_window.x = alongtrack.x(window_indices);
        alongtrack_window.y = alongtrack.y(window_indices);
        alongtrack_window.t = alongtrack.t(window_indices);
        alongtrack_window.ssh = alongtrack.ssh(window_indices);
    
        % Compute the 2D composite for the current time window
        eddyPath_window.xe=eddyPath.xe(alongtrack_window.t-min(alongtrack.t));
        eddyPath_window.ye=eddyPath.ye(alongtrack_window.t-min(alongtrack.t));
        % [mz_r, rmid_r, numz_r, stdz_r] = radialProfile(alongtrack_window,eddyPath_window,showplot=0);
        % ssh_model = mz_r(1:length(rmid_max));
        % time-averaged eddy composite from full field
        [mz_xy, xmid_xy, ymid_xy, numz_xy, stdz_xy] = composite2D(alongtrack_window,eddyPath_window,showplot=0,interp=true);% options: bin_size=12.5*1e3
        ssh_model = mz_xy(xmid_xy>=min(spatial_window)&xmid_xy<=max(spatial_window),ymid_xy>=min(spatial_window)&ymid_xy<=max(spatial_window));
        % MSE per time window
        mse_t(i) = mean((ssh_true(:,:,i) - ssh_model).^2,'all','omitmissing');
        
        % Calculate percentage of filled cells within the window
        total_cells = length(spatial_window)^2;
        filled_cells = sum(~isnan(ssh_model(:)));
        coverage_t(i) = filled_cells / total_cells;    
        end

        case 'Gaussian'
            % Fit a Gaussian model
            [paramsCell, initParamsCell] = FitAlongTrackXYToEddyModelWindowed(alongtrack, eddyFit_fun, initParams, eddyPath, it_options,window=window_days);
            for i=1:totalTimeWindows
            % Calculate the start and end times for this window in days
            window_start_day = min(alongtrack.t) + (i-1)*time_step;
            window_end_day = window_start_day + window_days;
        
            % Find indices that correspond to times within this window
            window_indices = find(alongtrack.t >= window_start_day & alongtrack.t <= window_end_day);
        
            % Extract time window
            alongtrack_window.x = alongtrack.x(window_indices);
            alongtrack_window.y = alongtrack.y(window_indices);
            alongtrack_window.t = alongtrack.t(window_indices);
            alongtrack_window.ssh = alongtrack.ssh(window_indices);
    
            eddy_model = analyticalEddyModel(eddyPath,paramsCell{i});
            ssh_model = eddy_model(alongtrack_window.x,alongtrack_window.y,alongtrack_window.t);        
            % MSE per time window
            mse_t(i) = mean((ssh_true(:,i) - ssh_model).^2);
            end
        case 'Elliptical'
            % Fit an elliptical model
            % [paramsCell, initParamsCell] = FitAlongTrackXYToEddyModelWindowed(alongtrack_window, eddyFit_fun, initParams, eddyPath_fun_t, it_options);
        otherwise
            error('Invalid model type');
    end
    % ssh_true_interp = interp2(ssh_true, xmid, ymid); % not necessary
    % because everything will be in the bin size resolution
    % Compute the MSE between the true SSH and the SSH model

    % MSE per window size
    mse(j) = mean(mse_t);
    coverage_percentage(j)=mean(coverage_t);
    
end
    figure
    plot(window_size_array,sqrt(mse*1e2),'LineWidth',2);
    xlabel('Window size', 'FontName', 'times');ylabel('MSE (cm^2)', 'FontName', 'times')
    set(gca, 'fontname', 'times','fontsize', 16)
    % Store the absolute minimum MSE for the current model and time window
    min_MSE = min(mse);


    % weight=0.5;
    % combined_score=weight*(mse - min(mse)) / (max(mse) - min(mse))+(1-weight)*(1-coverage_percentage);
    % figure,plot(window_size_array,combined_score,'LineWidth',2);
    % ylabel('Grid coverage'); xlabel('Window size')
    % set(gca, 'fontname', 'times','fontsize', 16)
end