function [LB, UB] = generateTrackingBounds(trackingParams, window_duration, tracking_bound_options)
arguments
    trackingParams struct  % Tracking-derived parameters
    window_duration double % Length of the time window
    tracking_bound_options struct
    % tracking_bound_options.amplitudeTolerance double = 0.3    % ±30% around tracking amplitude
    % tracking_bound_options.radiusTolerance double = 0.2       % ±20% around tracking radius
    % tracking_bound_options.positionTolerance double = 0.5     % ±0.5 radii around tracking position
    % tracking_bound_options.velocityTolerance double = 0.3     % ±30% around current velocity
end

% Physical constraints (absolute min/max bounds)
tracking_bound_options.minAmplitude = 0.001;         % Minimum 1 mm
tracking_bound_options.maxAmplitude = 1.0;          % Maximum 1 m
tracking_bound_options.minRadius = 5e3;             % Minimum 5 km
tracking_bound_options.maxRadius = 200e3;           % Maximum 200 km
tracking_bound_options.maxVelocity = 8e3;           % Maximum 8 km/day

% Generate bounds based on tracking estimates with physical tolerances
A_center = trackingParams.A;
L_center = trackingParams.L;
x0_center = trackingParams.x0;
y0_center = trackingParams.y0;
cx_center = trackingParams.cx;
cy_center = trackingParams.cy;

% Amplitude bounds
A_range = A_center * [1-tracking_bound_options.amplitudeTolerance, 1+tracking_bound_options.amplitudeTolerance];
A_range(1) = max(A_range(1), tracking_bound_options.minAmplitude);
A_range(2) = min(A_range(2), tracking_bound_options.maxAmplitude);

% Radius bounds
L_range = L_center * [1-tracking_bound_options.radiusTolerance, 1+tracking_bound_options.radiusTolerance];
L_range(1) = max(L_range(1), tracking_bound_options.minRadius);
L_range(2) = min(L_range(2), tracking_bound_options.maxRadius);

% Position bounds: center ± N radii
pos_tolerance = tracking_bound_options.positionTolerance * L_center;
x0_range = [x0_center - pos_tolerance, x0_center + pos_tolerance];
y0_range = [y0_center - pos_tolerance, y0_center + pos_tolerance];

% Velocity bounds: center ± tolerance * |velocity|
cx_tolerance = tracking_bound_options.velocityTolerance * max(abs(cx_center), 0.5); % Minimum 0.5 m/s tolerance
cy_tolerance = tracking_bound_options.velocityTolerance * max(abs(cy_center), 0.5);
cx_range = [cx_center - cx_tolerance, cx_center + cx_tolerance];
cy_range = [cy_center - cy_tolerance, cy_center + cy_tolerance];

% Apply maximum velocity constraints
cx_range(1) = max(cx_range(1), -tracking_bound_options.maxVelocity);
cx_range(2) = min(cx_range(2), tracking_bound_options.maxVelocity);
cy_range(1) = max(cy_range(1), -tracking_bound_options.maxVelocity);
cy_range(2) = min(cy_range(2), tracking_bound_options.maxVelocity);

% Assemble bounds
LB = [A_range(1), L_range(1), x0_range(1), y0_range(1), cx_range(1), cy_range(1)];
UB = [A_range(2), L_range(2), x0_range(2), y0_range(2), cx_range(2), cy_range(2)];

fprintf('Tracking-based bounds (window duration: %.1f days):\n', window_duration);
fprintf('  A: [%.3f, %.3f] m (tracking: %.3f m)\n', LB(1), UB(1), A_center);
fprintf('  L: [%.1f, %.1f] km (tracking: %.1f km)\n', LB(2)/1e3, UB(2)/1e3, L_center/1e3);
fprintf('  x0: [%.1f, %.1f] km (tracking: %.1f km)\n', LB(3)/1e3, UB(3)/1e3, x0_center/1e3);
fprintf('  y0: [%.1f, %.1f] km (tracking: %.1f km)\n', LB(4)/1e3, UB(4)/1e3, y0_center/1e3);
fprintf('  cx: [%.3f, %.3f] m/s (tracking: %.3f m/s)\n', LB(5), UB(5), cx_center);
fprintf('  cy: [%.3f, %.3f] m/s (tracking: %.3f m/s)\n', LB(6), UB(6), cy_center);
end
