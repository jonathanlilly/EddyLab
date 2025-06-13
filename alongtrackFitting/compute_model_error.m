function [rmse_daily,rmse,getModelForDay,window_center]=compute_model_error(alongtrack,eddy_field, eddyPath_fun_t, window_size, model, options)
arguments
    alongtrack struct
    eddy_field struct
    eddyPath_fun_t struct
    window_size (1,1) {mustBeNumeric}
    model string
    options.bin_size = 12.5 * 1e3 % in meters
    options.max_r = 300 *1e3 % in meters
    options.overlap = 0.5 % overlap between time_window 0-1
    options.showplot = false %control whether to display plots
end
use options
%% make video
% clearvars outputVideo
% outputVideo = VideoWriter(strcat('E:\Research\myCode\alongtrack_Gaussian_window10.avi'));
% fr_rate = 365/60*5;
% outputVideo.FrameRate = fr_rate;
% open(outputVideo)
% f1=figure;
%%
t0=floor(min(alongtrack.t));
totalDays = length(eddy_field.t);
spatial_window = [fliplr(-bin_size/2:-bin_size:-max_r),bin_size/2:bin_size:max_r]';

% Calculate the start and end times for this window in days
[window_start_day, window_end_day, totalTimeWindows] = timeWindowBounds(alongtrack.t, window_size,overlap);
window_center = (window_start_day + window_end_day) / 2;

%generate model for each window
ssh_model = cell(totalTimeWindows,1);

switch model
case 'composite'
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
    [mz, xmid, ymid, numz, stdz] = composite2D(alongtrack_window,eddyPath_window,showplot=0,interp=false);
    [XGrid, YGrid] = ndgrid(xmid, ymid);
    % Create scatteredInterpolant and interpolate onto full grid
    F = scatteredInterpolant(XGrid(~isnan(mz)), YGrid(~isnan(mz)), mz(~isnan(mz)), 'natural', 'nearest');
    mz_interp = F(XGrid, YGrid);
    ssh_model{i} = @(x,y,t) findSSHmodel(x,y,xmid,ymid,mz_interp,spatial_window);
    end
    
case 'Gaussian'
    % Get eddy parameters
    [~,amplitude,radius] = findEddyCentroid(eddy_field.x, eddy_field.y, eddy_field.ssh,'thresholdratio',0.9,'GetBoundary', true);
    eddyParams.A=amplitude;
    eddyParams.L=radius;
    eddyParams.xe=eddyPath_fun_t.xe;
    eddyParams.ye=eddyPath_fun_t.ye;

    % Fit models to windows
    eddyFit_fun = @(x,y,t,A,L,x0,y0,cx,cy) A.*exp(-((x-x0-cx*t).^2 + (y-y0-cy*t).^2)/L^2);
    it_options = optimset('TolX',1e-3,'TolFun',1e-3);
    [paramsCell, trueParamsCell, ~] = FitAlongTrackXYToEddyModelWindowed(alongtrack, eddyFit_fun, eddyParams, it_options,"window", window_size,usePreviousFitConstraints=true);
    
    % Create model functions
    for i = 1:totalTimeWindows
        params = paramsCell{i};
        ssh_model{i} = @(x,y,t) eddyFit_fun(x,y,t-(params.t0-t0), params.A, params.L, params.x0, params.y0, params.cx, params.cy);
    end

case 'Elliptical'
    % Fit an elliptical model
    % [paramsCell, initParamsCell] = FitAlongTrackXYToEddyModelWindowed(alongtrack_window, eddyFit_fun, initParams, eddyPath_fun_t, it_options);
otherwise
    error('Invalid model type');
end

%% Compute MSE per day
% Create a function that returns an interpolated model for any day
getModelForDay = @(day) createInterpolatedModel(day, window_center, ssh_model);

% Initialize MSE array
mse_t = nan(totalDays, 1);

for n = 1:totalDays
    % Get model for this day
    if totalTimeWindows == 1
        ssh_model_interp = ssh_model{1};
    else
        ssh_model_interp = getModelForDay(t0 + n - 1);
    end

    % eddy-centered coordinates
    x_rel = eddy_field.x - eddyPath_fun_t.xe(n-1);
    y_rel = eddy_field.y - eddyPath_fun_t.ye(n-1);
    
    % Find points within spatial window
    in_window_x = x_rel >= min(spatial_window) & x_rel <= max(spatial_window);
    in_window_y = y_rel >= min(spatial_window) & y_rel <= max(spatial_window);
    
    ssh_true_n = eddy_field.ssh(in_window_x, in_window_y, n)';
    % Handle coordinates based on model type
    switch model
        case 'composite'
            % Composite models expect eddy-centered coordinates
            [x_grid, y_grid] = ndgrid(x_rel(in_window_x), y_rel(in_window_y));
            ssh_model_n = ssh_model_interp(x_grid', y_grid', n-1);
            
        case 'Gaussian'
            % Gaussian models expect earth fixed coordinates
            [x_grid, y_grid] = ndgrid(eddy_field.x, eddy_field.y);
            ssh_model_fulldomain = ssh_model_interp(x_grid', y_grid', n-1);
            ssh_model_n = ssh_model_fulldomain(in_window_y, in_window_x);
    end
% figure(f1);jpcolor(x_rel(in_window_x)/1e3,y_rel(in_window_y)/1e3,ssh_model_n*1e2);
% shading flat
% axis equal
% xlabel('Distance East (km)', 'FontName', 'times')
% ylabel('Distance North (km)', 'FontName', 'times')
% set(gca, 'fontname', 'times','fontsize',16)
% if exist('brewermap')
%     colormap(brewermap([], '-Spectral'))
% end
% c = colorbar('EastOutside');
% c.Label.String = 'Time-averaged SSH (cm)';
% c.Label.FontSize=16;
% c.Label.FontName='times';
% xlim([-250, 250]), ylim([-250, 250]);
% set(gcf,'position',[2872,291,482,348])
% set(gcf,"Color",'white')
% clim([-2.5,13])
% text(0,270,['Day ',num2str(n)],'FontSize',12,'fontname','times','HorizontalAlignment','center')
% F=getframe(gcf);
% writeVideo(outputVideo,F)

    % Compute MSE for this day
    mse_daily(n) = mean((ssh_true_n - ssh_model_n).^2, 'all', 'omitnan');
end

rmse_daily=sqrt(mse_daily);
% MSE per window size
rmse = sqrt(mean(mse_daily,'omitnan'));
%% 

if options.showplot
    %% Plot true, model, and diff SSH
    figure;
    subplot(1, 3, 1);
    jpcolor(x_rel(in_window_x),x_rel(in_window_y),ssh_true_n)
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
    jpcolor(x_rel(in_window_x),y_rel(in_window_y), ssh_true_n - ssh_model_n); shading flat;
    hold on;axis tight
    colorbar;
    title('True - Model');
    axis equal;

    figure;
    plot([1:totalDays]-1,rmse_daily(:),'linewidth',2)
    xlabel('Time (day)')
    ylabel('RMSE (m)')
    set(gca,'fontname','times','fontsize',16)

plotFitPosition(eddyPath_fun_t,paramsCell,trueParamsCell)
xlim([min(x/1e3),max(x/1e3)]);ylim([min(y/1e3),max(y/1e3)])
plotFitWindowed(paramsCell,trueParamsCell, eddyPath_fun_t)

% To analyze a specific window in detail (e.g., window 3):
plotSingleWindowFit(alongtrack, eddyPath_fun_t, paramsCell, trueParamsCell, 60);

end
