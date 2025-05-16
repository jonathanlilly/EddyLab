function [mse_t,mse]=compute_mse_model_fun(alongtrack,eddy_field, eddyPath_fun_t, window_size, model, options)
arguments
    alongtrack struct
    eddy_field struct
    eddyPath_fun_t struct
    window_size (1,1) {mustBeNumeric}
    model string
    options.bin_size = 12.5 * 1e3; % in meters
    options.max_r = 250 *1e3 % in meters
    options.overlap = 0; %50% overlap between time_window
    options.showplot = false; %control whether to display plots
end
use options

t0=min(alongtrack.t);
% totalDays = max(fullfield.t)-min(fullfield.t)+1;
spatial_window = [fliplr(-bin_size/2:-bin_size:-max_r),bin_size/2:bin_size:max_r]';
%%

%generate model
switch model
case 'composite'
    % Calculate the start and end times for this window in days
    [window_start_day, window_end_day, totalTimeWindows] = timeWindowBounds(alongtrack.t, window_size,overlap);
    
    
    %for each time window generate a model (x,y,t) and compare to truth
    for i = 1:totalTimeWindows
    % Extract time window
    [alongtrack_window, window_indices] = extractAlongtrackWindow(alongtrack, window_start_day(i), window_end_day(i));
    
    if isempty(window_indices)
        continue
    end
    
    % Compute the 2D composite for the current time window
    eddyPath_window.xe=eddyPath_fun_t.xe(alongtrack_window.t-t0);
    eddyPath_window.ye=eddyPath_fun_t.ye(alongtrack_window.t-t0);

    % time-averaged eddy composite from full field
    [mz, xmid, ymid, numz, stdz] = composite2D(alongtrack_window,eddyPath_window,showplot=0);

    ssh_model{i} = @(x,y) findSSHmodel(x,y,xmid,ymid,mz,spatial_window);
    end
    window_center = (window_start_day + window_end_day) / 2;
case 'Gaussian'
    % Fit a Gaussian model
    it_options = optimset('TolX',1e-3,'TolFun',1e-3);
    [~,amplitude,radius] = findEddyCentroid(eddy_field.x, eddy_field.y, eddy_field.ssh,'thresholdratio',0.9,'GetBoundary', true);
    initParams.A=amplitude;
    initParams.L=radius;
    eddyFit_fun = @(x,y,t,A,L,x0,y0,cx,cy) A.*exp(-((x-x0-cx*t).^2 + (y-y0-cy*t).^2)/L^2);
    [paramsCell, initParamsCell,window_center] = FitAlongTrackXYToEddyModelWindowed(alongtrack, eddyFit_fun, initParams, eddyPath_fun_t, it_options,window=window_size);
    
    % Store the SSH model function for each window
    ssh_model = cell(totalTimeWindows,1);
    
    for i = 1:totalTimeWindows
        params = paramsCell{i};
        ssh_model{i} = @(x,y,t) eddyFit_fun(x,y,t-params.t0,params.A,params.L,params.x0,params.y0,params.cx,params.cy);
    end
    % eddy_model = analyticalEddyModel(eddyPath,paramsCell{i});
    % ssh_model{i} = eddy_model(alongtrack_window.x,alongtrack_window.y,alongtrack_window.t);
case 'Elliptical'
    % Fit an elliptical model
    % [paramsCell, initParamsCell] = FitAlongTrackXYToEddyModelWindowed(alongtrack_window, eddyFit_fun, initParams, eddyPath_fun_t, it_options);
otherwise
    error('Invalid model type');
end
end
%%
% % For final adjusted window, only evaluate the non-overlapping portion
% if isFinalWindow
%     % Calculate the start of non-overlapping days
%     eval_start_day = t0 + (totalTimeWindows-1)*time_step;
% else
%     eval_start_day = window_start_day;
% end

%% Create a function that returns an interpolated model for any day
getModelForDay = @(day) createInterpolatedModel(day, window_center, ssh_model);

% for n=1:length([eval_start_day:window_end_day])
mse_t=nan(length(eddy_field.t),1);

for n=1:totalDays
% Convert to eddy-centered coordinates
x_rel = eddy_field.x - eddyPath_fun_t.xe(n-1);
y_rel = eddy_field.y - eddyPath_fun_t.ye(n-1);

% Find points within spatial window
in_window_x = x_rel >= min(spatial_window) & x_rel <= max(spatial_window);
in_window_y = y_rel >= min(spatial_window) & y_rel <= max(spatial_window);
[x_grid,y_grid] = ndgrid(x_rel(in_window_x),y_rel(in_window_y));

ssh_true_n=eddy_field.ssh(in_window_x,in_window_y,n)';
if totalTimeWindows==1
    ssh_model_n=ssh_model{1}(x_grid',y_grid',n);
else
% Get the model function for this specific day
ssh_model_interp = getModelForDay(t0+n-1);
ssh_model_n=ssh_model_interp(x_grid',y_grid',n);
end
% MSE per day
mse_t(n) = mean((ssh_true_n - ssh_model_n).^2,'all','omitnan');
end
% MSE per window size
mse = mean(mse_t,'omitnan');

if options.showplot
    %% Plot true, model, and diff SSH
    figure;
    subplot(1, 3, 1);
    jpcolor(x_rel(in_window_x),x_rel(in_window_y),eddy_field.ssh(in_window_x,in_window_y,end)')
    shading flat;axis tight
    hold on;
    colorbar;
    title('True SSH');
    axis equal;
    
    % Plot model SSH
    subplot(1, 3, 2);
    jpcolor(x_rel(in_window_x),y_rel(in_window_y),ssh_model_n); shading flat;
    hold on;axis tight
    colorbar;
    title('Model SSH');
    axis equal;
    
    % Plot difference
    subplot(1, 3, 3);
    jpcolor(x_rel(in_window_x),y_rel(in_window_y), eddy_field.ssh(in_window_x,in_window_y,end)' - ssh_model_n); shading flat;
    hold on;axis tight
    colorbar;
    title('Difference (True - Model)');
    axis equal;

    figure;
    plot([1:totalDays]-1,mse_t(:),'linewidth',2)
    xlabel('Time (day)')
    ylabel('MSE (m^2)')
    set(gca,'fontname','times','fontsize',16)

end

end

plotFitPosition(paramsCell,initParamsCell)
xlim([min(x),max(x)]);ylim([min(y),max(y)])
plotFitWindowed(paramsCell,true_params, eddyPath_fun_t)