function [center_xy, amplitude, varargout] = findEddyCentroid_fromVelocity(x, y, u3D, v3D, options)
arguments
    x (:,:) {mustBeNumeric}
    y (:,:) {mustBeNumeric}
    u3D (:,:,:) {mustBeNumeric}  % x-component of velocity
    v3D (:,:,:) {mustBeNumeric}  % y-component of velocity
    options.ThresholdRatio (1,1) {mustBeNumeric} = 0.9 %default for SSH
    options.GetBoundary (1,1) logical = false

end
% Eddy detection based on velocity field vorticity
% Inputs:
% x, y: 1D coordinate arrays [m]
% u3D: 3D array of x-velocity values (size: nx x ny x nt) [m/s]
% v3D: 3D array of y-velocity values (size: nx x ny x nt) [m/s]
% Optional inputs (Name-Value pairs):
% 'ThresholdRatio': Ratio for thresholding (default: 0.9)
% 'GetBoundary': Whether to compute boundary (default: false)
% 'epsilon': Tolerance for shield zero-crossing vorticity (default: 0.01)

% Outputs:
% center_xy: [nt x 2] array of centroid coordinates
% amplitude: Maximum vorticity amplitude over time [1/s]
% Optional outputs (if GetBoundary is true):
% radius: Equivalent radius of core [m]
% core_xy: Cell array of core boundary coordinates [x;y]

    [xMat, yMat] = ndgrid(x, y);
    nt = size(u3D, 3);
    center_xy = zeros(nt, 2);
    amplitude = zeros(nt, 1);

    % Calculate grid spacing
    dx = x(2) - x(1);
    dy = y(2) - y(1);

    if options.GetBoundary
        core_xy = cell(nt, 1);
        radius = zeros(nt, 1);
    end
    
    % Process each time step
    for t = 1:nt
        % Get current velocity fields
        u = u3D(:,:,t);
        v = v3D(:,:,t);
        
        % Calculate vorticity: zeta = dv/dx - du/dy
        [du_dy, ~] = gradient(u, dy, dx);
        [~, dv_dx] = gradient(v, dy, dx);
        zeta = dv_dx - du_dy;
        
        % Find maximum vorticity
        [maxzeta, maxind] = max(abs(zeta(:)));
        [maxx, maxy] = ind2sub(size(zeta), maxind);
        max_point = [x(maxx), y(maxy)];
        
        % Store amplitude with sign
        amplitude(t) = zeta(maxind);
        
        % Find zero vorticity contour (core boundary)
        C = contourc(x, y, zeta', [0 0]);
        
        if isempty(C)
            % No zero contour found, use centroid at maximum
            center_xy(t,:) = max_point;
            if options.GetBoundary
                core_xy{t} = [];
                radius(t) = NaN;
            end
            continue;
        end
        
        % Find contour closest to vorticity maximum
        contourArray = splitContours(C);
        
        % Get the most relevant contour (closest to vorticity maximum)
        distances = zeros(length(contourArray), 1);
        for i = 1:length(contourArray)
            distances(i) = min(sqrt((max_point(1)-contourArray{i}(1,:)').^2 + ...
                (max_point(2)-contourArray{i}(2,:)').^2));
        end
        [~, closest_idx] = min(distances);
        boundary = contourArray{closest_idx}';
        
        % Calculate centroid
        polyin = polyshape(boundary(:,1), boundary(:,2));
        [center_xy(t,1), center_xy(t,2)] = centroid(polyin);
        
        % Get boundary if requested
        if options.GetBoundary
            core_xy{t} = boundary';
            
            % Calculate area and radius
            area = polyarea(boundary(:,1), boundary(:,2));
            radius(t) = sqrt(area/pi);
        end
    end
    
    % Set output arguments
    if options.GetBoundary
        varargout{1} = radius;
        varargout{2} = core_xy;
    end
end

function contourArray = splitContours(C)
% Helper function to split contour matrix into cell array of contours
i = 1;
contourArray = {};
while i < size(C,2)
    n = C(2,i);
    contourArray{end+1} = C(:,i+1:i+n);
    i = i + n + 1;
end
end