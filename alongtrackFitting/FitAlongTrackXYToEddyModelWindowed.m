function [paramsCell, trueParamsCell, window_center] = FitAlongTrackXYToEddyModelWindowed(alongtrack, eddyFit_fun, eddyParams, it_options, options)
arguments
    alongtrack struct
    eddyFit_fun function_handle 
    eddyParams struct
    it_options struct
    options.window (1,1) {mustBeNumeric}
    options.overlap = 0.5; %50% overlap between time_window
    options.usePreviousFitConstraints logical = true   % Enable/disable previous fit bounds
    options.setOffset=false % if you want to add offset to 
    % Tracking-based bound options (primary)
    options.amplitudeTolerance double = 0.5      % ±50% around tracking amplitude
    options.radiusTolerance double = 0.5         % ±50% around tracking radius
    options.positionTolerance double = 1.0       % ±N radii around tracking position
    options.velocityTolerance double = 1.0       % ±N times tracking velocity
    
end
%alongtrackLatLon: lat,lon,time
% Deduce time window length from alongtrackLatLon
t0=floor(min(alongtrack.t));

% totalDays=max(alongtrack.t)-t0+1;
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

%velocity
eddyParams.cx=vdiff(eddyParams.xe((window_start_day(1):window_end_day(totalTimeWindows))-t0),2);
eddyParams.cy=vdiff(eddyParams.ye((window_start_day(1):window_end_day(totalTimeWindows))-t0),2);

for i=1:totalTimeWindows
    % Extract time window
    [alongtrack_window, window_indices] = extractAlongtrackWindow(alongtrack, window_start_day(i), window_end_day(i));

    %calculate velocity
    xe=eddyParams.xe((window_start_day(i):window_end_day(i))-t0);%function - start at 0
    ye=eddyParams.ye((window_start_day(i):window_end_day(i))-t0);%function - start at 0
    cx=eddyParams.cx((window_start_day(i):window_end_day(i))-t0+1);%index - start at 1
    cy=eddyParams.cy((window_start_day(i):window_end_day(i))-t0+1);%index - start at 1
    
    % Define t0 for this specific window
    t0_window = min(alongtrack_window.t);
    elapsed_time_window = alongtrack_window.t-t0_window;
    
    % Generate initial parameters from tracking data
    trueParams.A=mean(eddyParams.A((window_start_day(i):window_end_day(i))-t0+1));
    trueParams.L=mean(eddyParams.L((window_start_day(i):window_end_day(i))-t0+1));
    trueParams.cx=mean(cx);
    trueParams.cy=mean(cy);
    trueParams.x0=xe(1);
    trueParams.y0=ye(1);
    trueParamsCell{i,1} = trueParams;
    
    % initial guess
    if i == 1
        initParams = trueParams;  % First window: use tracking
    else
        % Later windows: project previous fit
        time_step = window_start_day(i) - window_start_day(i-1);
        initParams.x0 = paramsCell{i-1}.x0 + paramsCell{i-1}.cx * time_step;
        initParams.y0 = paramsCell{i-1}.y0 + paramsCell{i-1}.cy * time_step;
        initParams.A = paramsCell{i-1}.A;
        initParams.L = paramsCell{i-1}.L;
        initParams.cx = paramsCell{i-1}.cx;
        initParams.cy = paramsCell{i-1}.cy;
    end
    
%     radius_tolerance = 0.35;  % ±15% only
%     pos_tolerance = 0.4;      % ±0.4 radii for position
%     vel_tolerance = 0.5;     % ±25% for velocity
% 
%     LB = [
%         trueParams.A * 0.7,                                    % A: -30%
%         trueParams.L * (1 - radius_tolerance),                 % L: -15% (tight!)
%         trueParams.x0 - pos_tolerance * trueParams.L,          % x0: ±0.4 radii
%         trueParams.y0 - pos_tolerance * trueParams.L,          % y0: ±0.4 radii
%         trueParams.cx - max(abs(trueParams.cx) * vel_tolerance, 0.3),  % cx: min ±0.3 m/s
%         trueParams.cy - max(abs(trueParams.cy) * vel_tolerance, 0.3)   % cy: min ±0.3 m/s
%     ];
% 
%     UB = [
%         trueParams.A * 1.3,                                    % A: +30%
%         trueParams.L * (1 + radius_tolerance),                 % L: +15% (tight!)
%         trueParams.x0 + pos_tolerance * trueParams.L,          % x0: ±0.4 radii
%         trueParams.y0 + pos_tolerance * trueParams.L,          % y0: ±0.4 radii
%         trueParams.cx + max(abs(trueParams.cx) * vel_tolerance, 0.3),  % cx: min ±0.3 m/s
%         trueParams.cy + max(abs(trueParams.cy) * vel_tolerance, 0.3)   % cy: min ±0.3 m/s
%     ];
% 
%     % Call your existing fitting function with bounds
%     bound.lower=LB;
%     bound.upper=UB;
%     params = FitAlongTrackXYToEddyModel(alongtrack_window, eddyFit_fun, initParams, it_options, ...
%         bound=bound);
% 
%     % % velocity smoothing
%     % if i > 1
%     %     % If velocities changed too much, smooth them
%     %     if abs(params.cx - paramsCell{i-1}.cx) > 0.5  % >0.5 m/s change
%     %         params.cx = 0.7 * params.cx + 0.3 * paramsCell{i-1}.cx;
%     %     end
%     %     if abs(params.cy - paramsCell{i-1}.cy) > 0.5
%     %         params.cy = 0.7 * params.cy + 0.3 * paramsCell{i-1}.cy;
%     %     end
%     % end
% 
%     params.t0 = window_start_day(i);
%     params.elapsed_time = elapsed_time_window;
%     paramsCell{i,1} = params;
% end

% Generate bounds
tracking_bound_options = struct(...
    'amplitudeTolerance', options.amplitudeTolerance, ...
    'radiusTolerance', options.radiusTolerance, ...
    'positionTolerance', options.positionTolerance, ...
    'velocityTolerance', options.velocityTolerance);

[LB_tracking, UB_tracking] = generateTrackingBounds(initParams, window_end_day(i) - window_start_day(i), tracking_bound_options);

% STEP 3: Optionally combine with previous fit bounds using intersection
if i > 1 && options.usePreviousFitConstraints
    fprintf('\nApplying additional constraints from previous fit...\n');

    prevParams = paramsCell{i-1};
    time_diff = window_start_day(i) - window_start_day(i-1);

    prevfit_bound_options = struct(...
        'amplitudeTolerance', 0.15, ...
        'radiusTolerance', 0.2 , ...
        'velocityTolerance', 0.3, ...
        'positionTolerance', 0.5);

    [LB_previous, UB_previous] = generatePreviousFitBounds(prevParams, time_diff, prevfit_bound_options);

    % Combine bounds using intersection (&&)
    [LB_final, UB_final] = intersectBounds(LB_tracking, UB_tracking, LB_previous, UB_previous);
else
    % Use only tracking bounds
    LB_final = LB_tracking;
    UB_final = UB_tracking;
    fprintf('Using tracking bounds only (first window or previous constraints disabled)\n');
end

bound.lower=LB_final;
bound.upper=UB_final;

% Call eddy model fit in XY
params = FitAlongTrackXYToEddyModel(alongtrack_window, eddyFit_fun, trueParams, it_options,bound=bound);
params.t0 = window_start_day(i);
paramsCell{i,1} = params;
end
end

function [LB_previous, UB_previous] = generatePreviousFitBounds(prevParams, time_diff, options)
% Set default options if not provided
if nargin < 3 || isempty(options)
    options = struct();
end

% Set default values for missing fields
if ~isfield(options, 'radiusTolerance'), options.radiusTolerance = 0.2; end
if ~isfield(options, 'velocityTolerance'), options.velocityTolerance = 0.3; end
if ~isfield(options, 'amplitudeTolerance'), options.amplitudeTolerance = 0.15; end
if ~isfield(options, 'positionTolerance'), options.positionTolerance = 0.5; end

% Project previous fit to current window start time
x0_projected = prevParams.x0 + prevParams.cx * time_diff;
y0_projected = prevParams.y0 + prevParams.cy * time_diff;

% Generate conservative bounds around previous fit values
A_range = prevParams.A * [1-options.amplitudeTolerance, 1+options.amplitudeTolerance];
L_range = prevParams.L * [1-options.radiusTolerance, 1+options.radiusTolerance];

% Position bounds: around projected position
pos_tolerance = options.positionTolerance * prevParams.L;
x0_range = [x0_projected - pos_tolerance, x0_projected + pos_tolerance];
y0_range = [y0_projected - pos_tolerance, y0_projected + pos_tolerance];

% Velocity bounds: around previous velocity
cx_tolerance = abs(prevParams.cx) * options.velocityTolerance;
cy_tolerance = abs(prevParams.cy) * options.velocityTolerance;
cx_range = [prevParams.cx - cx_tolerance, prevParams.cx + cx_tolerance];
cy_range = [prevParams.cy - cy_tolerance, prevParams.cy + cy_tolerance];

% Assemble bounds
LB_previous = [A_range(1), L_range(1), x0_range(1), y0_range(1), cx_range(1), cy_range(1)];
UB_previous = [A_range(2), L_range(2), x0_range(2), y0_range(2), cx_range(2), cy_range(2)];

fprintf('Previous fit bounds (projected %.1f days forward):\n', time_diff);
fprintf('  Position projected from (%.1f, %.1f) to (%.1f, %.1f) km\n', ...
    prevParams.x0/1e3, prevParams.y0/1e3, x0_projected/1e3, y0_projected/1e3);
end

function [LB_combined, UB_combined] = intersectBounds(LB_tracking, UB_tracking, LB_previous, UB_previous)
% Combine tracking bounds with previous fit bounds using intersection (&&)
% This creates tighter bounds where both constraints must be satisfied

LB_combined = max(LB_tracking, LB_previous);  % Take the maximum of lower bounds
UB_combined = min(UB_tracking, UB_previous);  % Take the minimum of upper bounds

% Check for valid bounds (LB <= UB)
invalid_bounds = LB_combined > UB_combined;
if any(invalid_bounds)
    param_names = {'A', 'L', 'x0', 'y0', 'cx', 'cy'};
    fprintf('Warning: Incompatible bounds detected for parameters: %s\n', ...
        strjoin(param_names(invalid_bounds), ', '));
    fprintf('Using tracking bounds only for incompatible parameters.\n');

    % Revert to tracking bounds for incompatible parameters
    LB_combined(invalid_bounds) = LB_tracking(invalid_bounds);
    UB_combined(invalid_bounds) = UB_tracking(invalid_bounds);
end

fprintf('Combined bounds (tracking && previous):\n');
param_names = {'A', 'L', 'x0', 'y0', 'cx', 'cy'};
units = {'m', 'km', 'km', 'km', 'm/s', 'm/s'};
scale_factors = [1, 1e-3, 1e-3, 1e-3, 1, 1];

for i = 1:6
    fprintf('  %s: [%.3f, %.3f] %s\n', param_names{i}, ...
        LB_combined(i)*scale_factors(i), UB_combined(i)*scale_factors(i), units{i});
end
end