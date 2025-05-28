function [paramsCell, trueParamsCell, window_center] = FitAlongTrackXYToEddyModelWindowed(alongtrack, eddyFit_fun, eddyParams, it_options, options)
arguments
    alongtrack struct
    eddyFit_fun function_handle 
    eddyParams struct
    it_options struct
    options.window (1,1) {mustBeNumeric}
    options.overlap = 0; %50% overlap between time_window
    options.LB (1,6) double %= [0, 50, -1000e3, -500e3, -10, -10]
    options.UB (1,6) double %= [30, 150, 1000e3, 500e3, 0, 0]
    options.setOffset=false % if you want to add offset to 
end
%alongtrackLatLon: lat,lon,time
% Deduce time window length from alongtrackLatLon
t0=min(alongtrack.t);
totalDays=max(alongtrack.t)-t0+1;

% min_days_per_window = 45;  % Need at least 4.5 cycles for the eddy core coverage
% min_total_windows = 3;      % Need at least 3 windows to see evolution
overlap = options.overlap; %50% overlap between time_window

% if isfield(options,'window')
    window_size=options.window;
% else
    % window_size = max(min_days_per_window, floor(totalDays/(min_total_windows/overlap)));
% end

[window_start_day, window_end_day, totalTimeWindows] = timeWindowBounds(alongtrack.t, window_size, overlap);
window_center = (window_start_day + window_end_day) / 2;

% pre-allocate cells for number of windows
% paramsCell = cell(totalTimeWindows, 1);

for i=1:totalTimeWindows
    % Extract time window
    [alongtrack_window, window_indices] = extractAlongtrackWindow(alongtrack, window_start_day(i), window_end_day(i));

%calculate velocity
xe=eddyParams.xe((window_start_day(i):window_end_day(i))-t0);
ye=eddyParams.ye((window_start_day(i):window_end_day(i))-t0);
cx=vdiff(xe,2);
cy=vdiff(ye,2);

% Define t0 for this specific window
t0_window = min(alongtrack_window.t);
elapsed_time_window = alongtrack_window.t-t0_window;

trueParams.A=mean(eddyParams.A((window_start_day(i):window_end_day(i))-t0+1));
trueParams.L=mean(eddyParams.L((window_start_day(i):window_end_day(i))-t0+1));
trueParams.cx=mean(cx);
trueParams.cy=mean(cy);
trueParams.x0=xe(1);
trueParams.y0=ye(1);
trueParamsCell{i,1} = trueParams;

% initParams_window.A=mean(trueParams.A);
% initParams_window.L=mean(trueParams.L);
% initParams_window.cx=mean(trueParams.cx);
% initParams_window.cy=mean(trueParams.cy);

% Call eddy model fit in XY
params = FitAlongTrackXYToEddyModel(alongtrack_window, eddyFit_fun, trueParams, it_options);
params.t0 = window_start_day(i);
paramsCell{i,1} = params;
end

