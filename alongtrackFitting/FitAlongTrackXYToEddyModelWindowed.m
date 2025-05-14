function [paramsCell, initParamsCell, window_center] = FitAlongTrackXYToEddyModelWindowed(alongtrack, eddyFit_fun, initParams, eddyPath_fun_t, it_options, options)
arguments
    alongtrack struct
    eddyFit_fun function_handle 
    initParams struct
    eddyPath_fun_t struct
    it_options struct
    options.window (1,1) {mustBeNumeric}
    options.overlap = 0; %50% overlap between time_window
    options.LB (1,6) double %= [0, 50, -1000e3, -500e3, -10, -10]
    options.UB (1,6) double %= [30, 150, 1000e3, 500e3, 0, 0]
end
%alongtrackLatLon: lat,lon,time
% Deduce time window length from alongtrackLatLon
t0=min(alongtrack.t);
totalDays=max(alongtrack.t)-t0+1;

min_days_per_window = 45;  % Need at least 4.5 cycles for the eddy core coverage
min_total_windows = 3;      % Need at least 3 windows to see evolution
overlap = options.overlap; %50% overlap between time_window

if isfield(options,'window')
    window_size=options.window;
else
    window_size = max(min_days_per_window, floor(totalDays/(min_total_windows/overlap)));
end

[window_start_day, window_end_day, totalTimeWindows] = timeWindowBounds(alongtrack.t, window_size, overlap);
window_center = (window_start_day + window_end_day) / 2;

% pre-allocate cells for number of windows
paramsCell = cell(totalTimeWindows, 1);

for i=1:totalTimeWindows
    % Extract time window
    [alongtrack_window, window_indices] = extractAlongtrackWindow(alongtrack, window_start_day(i), window_end_day(i));

%calculate velocity
cx=vdiff(eddyPath_fun_t.xe(alongtrack_window.t-t0),2);
cy=vdiff(eddyPath_fun_t.ye(alongtrack_window.t-t0),2);

if length(initParams.A)==1
    A=initParams.A;
    L=initParams.L;
else
    A=initParams.A(window_start_day(i)-t0+1);
    L=initParams.L(window_start_day(i)-t0+1);
end

% Define t0 for this specific window
t0_window = min(alongtrack_window.t);
elapsed_time_window = alongtrack_window.t-t0_window;

% new inital parameters for this window
initParams_window.A = A+2*(rand-0.5)*1e-2; %random uncertainty +/- 1e-2
initParams_window.L = L+2*(rand-0.5)*1e3; %random uncertainty +/- 1e3
%assuming you roughly know the eddy center from eddy-tracking algorithm
%beginning of this particular window minus t0 offset of the entire eddy lifetime
% initParams_window.x0 = initParams.x0+i*initParams.cx*max(elapsed_time_window);
% initParams_window.y0 = initParams.y0+i*initParams.cy*max(elapsed_time_window);
initParams_window.x0 = eddyPath_fun_t.xe(t0_window-t0)+2*(rand-0.5)*10e3; %random uncertainty +/- 10e3
initParams_window.y0 = eddyPath_fun_t.ye(t0_window-t0)+2*(rand-0.5)*10e3;
initParams_window.cx = cx(1)+2*(rand-0.5)*1e2;  %random uncertainty +/- 1e2
initParams_window.cy = cy(1)+2*(rand-0.5)*1e2;
initParamsCell{i,1} = initParams_window;

% bound


% Call eddy model fit in XY
params = FitAlongTrackXYToEddyModel(alongtrack_window, eddyFit_fun, initParams_window, it_options);

paramsCell{i,1} = params;
end

