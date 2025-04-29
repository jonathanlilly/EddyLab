% Load QG model
current_file_path = matlab.desktop.editor.getActiveFilename;
[current_dir, ~, ~] = fileparts(current_file_path);
readdir = 'G:\My Drive\AlongTrack\';
filename = 'BetaEddyOne.nc';

x_QG = ncread([readdir, filename], 'x'); %meters
y_QG = ncread([readdir, filename], 'y'); %meters
ssh = squeeze(ncread([readdir, filename], 'ssh')); %meter;% matrix order in x,y,z
lato = 24;

totalDays = size(ssh, 3);
eddy_time = ([1:totalDays]-1)';
x = x_QG - mean(x_QG);
y = y_QG - mean(y_QG);

eddy_field.x = x;
eddy_field.y = y;
eddy_field.t = eddy_time;
eddy_field.ssh = ssh;

% Full field struct
[X,Y,T] = ndgrid(1:length(eddy_field.x), 1:length(eddy_field.y), 1:length(eddy_field.t));
fullfield.x = eddy_field.x(reshape(X,[],1));
fullfield.y = eddy_field.y(reshape(Y,[],1));
fullfield.t = eddy_field.t(reshape(T,[],1))+timeo(1); % Subtract 1 to match your indexing
fullfield.ssh = reshape(eddy_field.ssh,[],1);

overlap = 0; %50% overlap between time_window
window_size_array=[1,100,totalDays];
model='composite';
[mse_t,mse]=compute_mse_fun(fullfield,eddy_field, eddyPath_fun_t, window_size_array, model,showplot=1)

function [mse_t,mse]=compute_mse_fun(alongtrack,eddy_field, eddyPath_fun_t, window_size_array, model, options)
arguments
    alongtrack struct
    eddy_field struct
    eddyPath_fun_t struct
    window_size_array (:,1) {mustBeNumeric}
    model string
    options.bin_size = 12.5 * 1e3; % in meters
    options.max_r = 250 *1e3 % in meters
    options.showplot = false; %control whether to display plots
end

[alongtrack.t,sort_idx]=sort(alongtrack.t,'ascend');
alongtrack.x=alongtrack.x(sort_idx);
alongtrack.y=alongtrack.y(sort_idx);
alongtrack.ssh=alongtrack.ssh(sort_idx);

% totalDays = max(fullfield.t)-min(fullfield.t)+1;
%%

T = length(window_size_array); %number of time variation
it_options = optimset('TolX',1e-3,'TolFun',1e-3);
spatial_window = [fliplr(-options.bin_size/2:-options.bin_size:-options.max_r),options.bin_size/2:options.bin_size:options.max_r]';

%test various time window
for j = 1:T
    window_days = window_size_array(j);
    
    time_step=floor(window_days*(1-overlap));
    totalTimeWindows=floor((totalDays) / time_step) + 1;

    %for each time window generate a model (x,y,t) and compare to truth
    mse_t=nan(length(eddy_field.t),1);

    for i = 1:totalTimeWindows

    % Calculate the start and end times for this window in days
    window_start_day = min(alongtrack.t) + (i-1)*time_step;
    window_end_day = window_start_day + window_days-1;
    
    % if the last window is small, count backwards from the end to
    % composite the same window size
    if (i-1)*time_step + window_days > totalDays
        window_end_day = max(alongtrack.t);
        window_start_day = window_end_day - window_days+1;
        isFinalWindow = true;
    else
        isFinalWindow = false;
    end

    % Find indices that correspond to times within this window
    window_indices = find(alongtrack.t >= window_start_day & alongtrack.t <= window_end_day);

    % Extract time window
    alongtrack_window.x = alongtrack.x(window_indices);
    alongtrack_window.y = alongtrack.y(window_indices);
    alongtrack_window.t = alongtrack.t(window_indices);
    alongtrack_window.ssh = alongtrack.ssh(window_indices);

    % Compute the 2D composite for the current time window
    eddyPath_window.xe=eddyPath_fun_t.xe(alongtrack_window.t-min(alongtrack.t));
    eddyPath_window.ye=eddyPath_fun_t.ye(alongtrack_window.t-min(alongtrack.t));

    %generate model
    switch model
    case 'composite'
    % time-averaged eddy composite from full field
    [mz, xmid, ymid, numz, stdz] = composite2D(alongtrack_window,eddyPath_window,showplot=0);

    ssh_model = @(x,y) findSSHmodel(x,y,xmid,ymid,mz,spatial_window);
    % % Create a function handle for this time window using interpolation
    % ssh_model = @(x, y) interp2(xmid_true, ymid_true', mz_true, x, y, 'linear', NaN);
    case 'Gaussian'
    end
    %%
    % For final adjusted window, only evaluate the non-overlapping portion
    if isFinalWindow
        % Calculate the start of non-overlapping days
        eval_start_day = min(alongtrack.t) + (totalTimeWindows-1)*time_step;
    else
        eval_start_day = window_start_day;
    end

    for n=1:length([eval_start_day:window_end_day])
        day_idx=(eval_start_day+n)-timeo(1);
    % Convert to eddy-centered coordinates
    x_rel = eddy_field.x - eddyPath_fun_t.xe(day_idx);
    y_rel = eddy_field.y - eddyPath_fun_t.ye(day_idx);
    
    % Find points within spatial window
    in_window_x = x_rel >= min(spatial_window) & x_rel <= max(spatial_window);
    in_window_y = y_rel >= min(spatial_window) & y_rel <= max(spatial_window);
    [x_grid,y_grid] = ndgrid(x_rel(in_window_x),y_rel(in_window_y));
    
    ssh_true_n=eddy_field.ssh(in_window_x,in_window_y,day_idx)';
    ssh_model_n=ssh_model(x_grid',y_grid');
    % MSE per day
    mse_t{j}(day_idx) = mean((ssh_true_n - ssh_model_n).^2,'all','omitnan');
    end
    end
    % MSE per window size
    mse(j) = mean(mse_t{j},'omitnan');
end

    if options.showplot
        %% Plot true, model, and diff SSH
        figure;
        subplot(1, 3, 1);
        jpcolor(x_rel(in_window_x),x_rel(in_window_y),eddy_field.ssh(in_window_x,in_window_y,window_start_day+n-timeo(1))')
        shading flat;axis tight
        hold on;
        colorbar;
        title('True SSH');
        axis equal;
        
        % Plot model SSH
        subplot(1, 3, 2);
        jpcolor(x_rel(in_window_x),y_rel(in_window_y),ssh_model(x_grid',y_grid')); shading flat;
        hold on;axis tight
        colorbar;
        title('Model SSH');
        axis equal;
        
        % Plot difference
        subplot(1, 3, 3);
        jpcolor(x_rel(in_window_x),y_rel(in_window_y), eddy_field.ssh(in_window_x,in_window_y,window_start_day+n-timeo(1))' - ssh_model(x_grid',y_grid')); shading flat;
        hold on;axis tight
        colorbar;
        title('Difference (True - Model)');
        axis equal;

        figure,hold on
        for j = 1:T
            plot(eddy_time,mse_t{j},'linewidth',2)
        end
        xlabel('Time (day)')
        ylabel('MSE (m^2)')
        set(gca,'fontname','times','fontsize',16)


