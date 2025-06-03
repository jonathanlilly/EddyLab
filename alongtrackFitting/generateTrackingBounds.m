function [LB, UB] = generateTrackingBounds(trackingParams, window_duration, options)
% arguments
%     trackingParams struct  % Tracking-derived parameters
%     window_duration double % Length of the time window
%     options.amplitudeTolerance double = 0.5    % ±50% around tracking amplitude
%     options.radiusTolerance double = 0.5       % ±50% around tracking radius
%     options.positionTolerance double = 1.0     % ±N radii around tracking position
%     options.velocityTolerance double = 1.0     % ±N times current velocity
% end

% Physical constraints (absolute min/max bounds)
options.minAmplitude = 0.001;         % Minimum 1 mm
options.maxAmplitude = 1.0;          % Maximum 1 m
options.minRadius = 5e3;             % Minimum 5 km
options.maxRadius = 200e3;           % Maximum 200 km
options.maxVelocity = 8e3;           % Maximum 8 km/day

% Generate bounds based on tracking estimates with physical tolerances
A_center = trackingParams.A;
L_center = trackingParams.L;
x0_center = trackingParams.x0;
y0_center = trackingParams.y0;
cx_center = trackingParams.cx;
cy_center = trackingParams.cy;

% Amplitude bounds
A_range = A_center * [1-options.amplitudeTolerance, 1+options.amplitudeTolerance];
A_range(1) = max(A_range(1), options.minAmplitude);
A_range(2) = min(A_range(2), options.maxAmplitude);

% Radius bounds
L_range = L_center * [1-options.radiusTolerance, 1+options.radiusTolerance];
L_range(1) = max(L_range(1), options.minRadius);
L_range(2) = min(L_range(2), options.maxRadius);

% Position bounds: center ± N radii
pos_tolerance = options.positionTolerance * L_center;
x0_range = [x0_center - pos_tolerance, x0_center + pos_tolerance];
y0_range = [y0_center - pos_tolerance, y0_center + pos_tolerance];

% Velocity bounds: center ± tolerance * |velocity|
cx_tolerance = options.velocityTolerance * max(abs(cx_center), 0.1); % Minimum 0.1 m/s tolerance
cy_tolerance = options.velocityTolerance * max(abs(cy_center), 0.1);
cx_range = [cx_center - cx_tolerance, cx_center + cx_tolerance];
cy_range = [cy_center - cy_tolerance, cy_center + cy_tolerance];

% Apply maximum velocity constraints
cx_range(1) = max(cx_range(1), -options.maxVelocity);
cx_range(2) = min(cx_range(2), options.maxVelocity);
cy_range(1) = max(cy_range(1), -options.maxVelocity);
cy_range(2) = min(cy_range(2), options.maxVelocity);

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
