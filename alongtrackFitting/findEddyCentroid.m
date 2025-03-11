function [center_xy, amplitude, varargout] = findEddyCentroid(x, y, ssh3D, varargin)
% following Chelton eddy tracking
% Inputs:
% x, y: 1D coordinate arrays
% ssh3D: 3D array of SSH values (size: ny x nx x nt)
% Optional inputs (Name-Value pairs):
% 'GetBoundary': Whether to compute boundary (default: false)

% Outputs:
% center_xy: [nt x 2] array of centroid coordinates
% Optional outputs (if GetBoundary is true):
% core_xy: Cell array of boundary coordinates [x;y]
% amplitude: SSH amplitude over time

    % Parse inputs
    p = inputParser;
    addParameter(p, 'ThresholdRatio', 0.9);
    addParameter(p, 'GetBoundary', false);
    parse(p, varargin{:});
    
    threshold_ratio = p.Results.ThresholdRatio;
    get_boundary = p.Results.GetBoundary;
    
    [xMat, yMat] = ndgrid(x, y);
    nt = size(ssh3D, 3);
    center_xy = zeros(nt, 2);
    amplitude = zeros(nt, 1);

    if get_boundary
        core_xy = cell(nt, 1);
        % Calculate grid spacing
        dx = x(2) - x(1);
        dy = y(2) - y(1);
    end
    
    % Process each time step
    for t = 1:nt
        % Get current SSH field
        ssh = ssh3D(:,:,t);
        
        % Find maximum
        [maxssh, maxind] = max(ssh(:));
        [maxx, maxy] = ind2sub(size(ssh), maxind);
        max_point = [x(maxx), y(maxy)];
        amplitude(t) = maxssh;
        if threshold_ratio == 1
            % If threshold is 1, just return maximum position
            center_xy(t,:) = max_point;
        else
            % Find points above threshold
            threshold = threshold_ratio * maxssh;
            above_threshold = ssh >= threshold;
            
            % Calculate weighted centroid
            x_points = xMat(above_threshold);
            y_points = yMat(above_threshold);
            ssh_weights = ssh(above_threshold);
            
            center_xy(t,1) = sum(x_points .* ssh_weights) / sum(ssh_weights);
            center_xy(t,2) = sum(y_points .* ssh_weights) / sum(ssh_weights);
        end

        % Get boundary only if requested
        if get_boundary
            % Find zero vorticity contour
            % Calculate vorticity (zeta=g/f nabla^2 eta)
            [Sx, Sy] = gradient(ssh, dx, dy);
            [Sxx, ~] = gradient(Sx, dx, dy);
            [~, Syy] = gradient(Sy, dx, dy);
            zeta = (Sxx + Syy); % without the g/f

            % Find SSH e-fold
            SSHefold=max(ssh(:)) * exp(-1);

            % C = contourc(x, y, zeta', [0 0]);
            C = contourc(x, y, ssh', [SSHefold SSHefold]);
            % Find contour closest to SSH maximum
            contourArray = splitContours(C);
            
            % Get the most relevant contour (closest to SSH maximum)
            distances = zeros(length(contourArray), 1);
            for i = 1:length(contourArray)
            distances(i) = min(sqrt((max_point(1)-contourArray{i}(1,:)').^2 + ...
            (max_point(2)-contourArray{i}(2,:)').^2)); %%min distance
            end
            [~, closest_idx] = min(distances);
            boundary = contourArray{closest_idx}';
            core_xy{t} = boundary';
            % Calculate area using polyarea
            area = polyarea(boundary(:,1), boundary(:,2));
            radius(t) = sqrt(area/pi);
        end
    end
    
    % Set output arguments
    if get_boundary
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
