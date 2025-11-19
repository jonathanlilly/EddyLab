function [center_xy, amplitude, varargout] = findEddyCentroid(x, y, ssh3D, options)
arguments
    x (:,:) {mustBeNumeric}
    y (:,:) {mustBeNumeric}
    ssh3D (:,:,:) {mustBeNumeric}
    options.ThresholdRatio (1,1) {mustBeNumeric} = 0.9 %default for SSH
    options.GetBoundary (1,1) logical = 'false'
    options.BoundaryType (1,1) string = 'SSH-e-fold'
    options.lato (1,1) {mustBeNumeric} = 24
    options.epsilon = 0.01%1e-4/(4e3).^2*g/fo/fo %Tolerance for shield zero-crossing, in case shield plateaus for zero-crossing zeta value. dSSH ~ 1e-2 m, dx ~ 1e4 m -> dSSH/dx*2=1e-11
end
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
use options

    [xMat, yMat] = ndgrid(x, y);
    nt = size(ssh3D, 3);
    center_xy = zeros(nt, 2);
    amplitude = zeros(nt, 1);

    if GetBoundary
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

        if isequal(BoundaryType,'zero-zeta')
            % Find zero vorticity contour
            [zeta] = zetaField(dx,dy,ssh,lato=lato);
            C = contourc(x, y, zeta', [0 0]);
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

            polyin=polyshape(boundary(:,1),boundary(:,2));
            [center_xy(t,1), center_xy(t,2)] = centroid(polyin);

        else %max ssh
            if ThresholdRatio == 1
                % If threshold is 1, just return maximum position
                center_xy(t,:) = max_point;
            else
                % Find points above threshold
                threshold = ThresholdRatio * maxssh;
                above_threshold = ssh >= threshold;
                
                % Calculate weighted centroid
                x_points = xMat(above_threshold);
                y_points = yMat(above_threshold);
                ssh_weights = ssh(above_threshold);
                
                center_xy(t,1) = sum(x_points .* ssh_weights) / sum(ssh_weights);
                center_xy(t,2) = sum(y_points .* ssh_weights) / sum(ssh_weights);
            end
        end

        % Get boundary only if requested
        if GetBoundary
            if isequal(BoundaryType,'zero-zeta')
                %core is already calculated above
    
                % Find the contour with vorticity close to epsilon that encloses the core
                C2 = contourc(x, y, zeta', [epsilon epsilon]);
                
                % Split the new contour matrix into separate contours
                contourArray2 = splitContours(C2);
                
                % Find the peak vorticity contour
                [max_zeta, max_idx] = max(zeta(:));
                max_contour_level = zeta(max_idx);
                C_max = contourc(x, y, zeta', [max_contour_level, max_contour_level]);
                
                % Split the peak vorticity contour matrix into separate contours
                contourArray_max = splitContours(C_max);
                
                % Find the contour with the peak vorticity closest to the core
                peak_contour = [];
                min_distance = Inf;
                for i = 1:length(contourArray_max)
                    current_contour = contourArray_max{i}';
                    distance = min(sqrt((current_contour(:,1) - center_xy(t,1)).^2 + (current_contour(:,2) - center_xy(t,2)).^2));
                    if distance < min_distance
                        peak_contour = current_contour;
                        min_distance = distance;
                    end
                end
                
                % Find the smallest contour that is bigger than the peak vorticity contour and encloses it
                second_boundary = [];
                min_area_diff = Inf;
                
                for i = 1:length(contourArray2)
                    current_boundary = contourArray2{i}';
                    
                    % Check if the current contour encloses the peak vorticity contour
                    if all(inpolygon(peak_contour(:,1), peak_contour(:,2), current_boundary(:,1), current_boundary(:,2)))
                        % Calculate the area difference between the current contour and the peak vorticity contour
                        current_area = polyarea(current_boundary(:,1), current_boundary(:,2));
                        peak_area = polyarea(peak_contour(:,1), peak_contour(:,2));
                        area_diff = current_area - peak_area;
                        
                        % Update the second boundary if the current contour is the smallest that encloses the peak vorticity contour
                        if area_diff < min_area_diff
                            second_boundary = current_boundary;
                            min_area_diff = area_diff;
                        end
                    end
                end
    
                if ~isempty(second_boundary)
                shield{t} = second_boundary'; % Store second boundary
                end
                % Calculate area using polyarea
                area = polyarea(boundary(:,1), boundary(:,2));
                radius(t) = sqrt(area/pi);

                % Set output arguments
                varargout{1} = radius;
                varargout{2} = core_xy;
                varargout{3} = shield;

            else
                % Find SSH e-fold
                SSHefold=max(ssh(:)) * exp(-1);
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

                % Set output arguments
                varargout{1} = radius;
                varargout{2} = core_xy;
            end
        end
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
