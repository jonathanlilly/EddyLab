function [paramsCell, initParamsCell] = FitAlongTrackXYToEddyModelWindowed(alongtrack, eddyFit_fun, initParams, eddyPath_fun_t, it_options, options)
arguments
    alongtrack struct
    eddyFit_fun function_handle 
    initParams struct
    eddyPath_fun_t struct
    it_options struct
    options.LB (1,6) double %= [0, 50, -1000e3, -500e3, -10, -10]
    options.UB (1,6) double %= [30, 150, 1000e3, 500e3, 0, 0]
end
%alongtrackLatLon: lat,lon,time
% Deduce time window length from alongtrackLatLon
totalDays = max(alongtrack.t)-min(alongtrack.t);

min_days_per_window = 45;  % Need at least 4.5 cycles for the eddy core coverage
min_total_windows = 3;      % Need at least 3 windows to see evolution
overlap = 0.5; %50% overlap between time_window
window_size = max(min_days_per_window, floor(totalDays/(min_total_windows/overlap)));

time_step=floor(window_size*(1-overlap));
totalTimeWindows=floor((totalDays - window_size) / time_step) + 1;

% ensure that the arrays are in ascending order in time before windowing
[alongtrack.t,sort_idx]=sort(alongtrack.t,'ascend');
alongtrack.x=alongtrack.x(sort_idx);
alongtrack.y=alongtrack.y(sort_idx);
alongtrack.lon=alongtrack.lon(sort_idx);
alongtrack.lat=alongtrack.lat(sort_idx);
alongtrack.ssh=alongtrack.ssh(sort_idx);

% pre-allocate cells for number of windows
paramsCell = cell(totalTimeWindows, 1);

for i=1:totalTimeWindows
    % Calculate the start and end times for this window in days
    window_start_day = min(alongtrack.t) + (i-1)*time_step;
    window_end_day = window_start_day + window_size;

    % Find indices that correspond to times within this window
    window_indices = find(alongtrack.t >= window_start_day & alongtrack.t <= window_end_day);

% Extract time window
at_window.x = alongtrack.x(window_indices);
at_window.y = alongtrack.y(window_indices);
at_window.t = alongtrack.t(window_indices);
at_window.ssh = alongtrack.ssh(window_indices);

% Define t0 for this specific window
t0_window = min(at_window.t);
elapsed_time_window = at_window.t-t0_window;

% % latc and lonc are the center of the alongtrack domain
% lonc=(min(at_window.lon(:))+max(at_window.lon(:)))/2;
% latc=(min(at_window.lat(:))+max(at_window.lat(:)))/2;
% 
% % Project {lat,lon} -> {x,y}
% [at_window.x, at_window.y] = latlon2xy(at_window.latitude, at_window.longitude, latc, lonc);

% new inital parameters for this window
initParams_window.A = initParams.A;
initParams_window.L = initParams.L;
% initParams_window.x0 = initParams.x0+i*initParams.cx*max(elapsed_time_window);
% initParams_window.y0 = initParams.y0+i*initParams.cy*max(elapsed_time_window);
%assuming you roughly know the eddy center from eddy-tracking algorithm
%beginning of this particular window minus t0 offset of the entire eddy lifetime
initParams_window.x0 = eddyPath_fun_t.xe(t0_window-min(alongtrack.t))+rand*20e3-10e3; %random uncertainty +/- 10e3
initParams_window.y0 = eddyPath_fun_t.ye(t0_window-min(alongtrack.t))+rand*20e3-10e3;
initParams_window.cx = initParams.cx;
initParams_window.cy = initParams.cy;
initParamsCell{i,1} = initParams_window;

% bound


% Call eddy model fit in XY
params = FitAlongTrackXYToEddyModel(at_window, eddyFit_fun, initParams_window, it_options);

paramsCell{i,1} = params;
end

