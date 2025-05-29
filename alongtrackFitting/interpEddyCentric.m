function [ssh_interp, x_out, y_out] = interpEddyCentric(x, y, t, xo, yo, ssh, bin_size, max_r)
% INTERPEDDYCENTRIC Interpolate gridded data onto eddy-centric coordinate system
%
% Inputs:
%   x, y, t    - Original coordinate vectors
%   xo, yo     - Eddy center positions at each time (arrays)
%   ssh        - 3D SSH data array (x, y, t)
%   bin_size   - Grid spacing for interpolation (in meters)
%
% Outputs:
%   ssh_interp   - Interpolated SSH data on eddy-centric grid
%   XGrid, YGrid - Meshgrid matrices of interpolated coordinates
% Determine data type if set to auto
% Automatically detect data type based on ssh structure
if isvector(ssh) && length(ssh) == length(x) && length(ssh) == length(y) && length(ssh) == length(t)
    % If ssh is a vector with same length as x, y, t, it's alongtrack data
    data_type = "alongtrack";
elseif ndims(ssh) == 3 && size(ssh,1) == length(x) && size(ssh,2) == length(y) && size(ssh,3) == length(t)
    % If ssh is a 3D array with dimensions matching x, y, t lengths, it's grid data
    data_type = "grid";
elseif issparse(ssh) || isscalar(size(ssh))
    % If ssh is sparse or has a single dimension, it's alongtrack data
    data_type = "alongtrack";
else
    warning('Could not automatically determine data type from ssh structure. Defaulting to alongtrack.');
    data_type = "alongtrack";
end

% Set up interpolation grid
max_r = (max_r/bin_size)*bin_size;
xbin = -max_r:bin_size:max_r;
xmid = (xbin(1:end-1) + xbin(2:end)) / 2;
ymid = xmid;
% Create meshgrid for interpolated coordinates
[XGrid, YGrid] = ndgrid(xmid, xmid);

if data_type == "grid"
% Loop over time to interpolate onto eddy-centric grid

% Pre-allocate for interpolated data
ssh_interp = zeros(length(xmid), length(xmid), length(t));
for i = 1:length(t)
    % Calculate eddy-relative coordinates for this time step
    xE = x - xo(i);
    yE = y - yo(i);
    
    % Create meshgrid of original coordinates
    [XE, YE] = ndgrid(xE, yE);
    
    % Interpolate SSH data onto eddy-centric grid
    ssh_interp(:,:,i) = interpn(XE, YE, ssh(:,:,i), XGrid, YGrid, 'linear', 0);
    x_out=XGrid;
    y_out=YGrid;
end

elseif data_type == "alongtrack"
    % Scattered data interpolation using scatteredInterpolant
    t_days=sort(unique(floor(t)),'ascend');
    
    % Pre-allocate for interpolated data
    ssh_interp=[];
    xE_mid=[];yE_mid=[];
    % Loop over time to interpolate scattered data
    for i = 1:length(t_days)
        % Calculate eddy-relative coordinates for this time step
        xE = x - xo(i);
        yE = y - yo(i);
        current_idx = find(floor(t) == t_days(i) & abs(xE) <= max_r & abs(yE) <= max_r);
        
        if ~isempty(current_idx)
            ssh_i = ssh(current_idx);
            xE_i = xE(current_idx);
            yE_i = yE(current_idx);
        
            % Find valid data points for this time step
            valid_idx = ~isnan(ssh_i) & ~isnan(xE_i) & ~isnan(yE_i);
            
            if sum(valid_idx) > 3 % Need at least 3 points for interpolation
                % Only process valid points
                xE_valid = xE_i(valid_idx);
                yE_valid = yE_i(valid_idx);
                ssh_valid = ssh_i(valid_idx);
                
                % Find nearest bin center for each VALID point
                xmid_i = zeros(length(xE_valid), 1);
                ymid_i = zeros(length(yE_valid), 1);
                
                for j = 1:length(xE_valid)
                    [~, xmid_idx] = min(abs(xmid - xE_valid(j)));
                    xmid_i(j) = xmid(xmid_idx);
                    [~, ymid_idx] = min(abs(xmid - yE_valid(j))); % Assuming square grid
                    ymid_i(j) = xmid(ymid_idx);
                end
                
                % Interpolate
                F = scatteredInterpolant(xE_valid, yE_valid, ssh_valid, 'natural', 'nearest');
                ssh_interp_i = F(xmid_i, ymid_i);
                
                % Store results
                xE_mid = [xE_mid; xmid_i];
                yE_mid = [yE_mid; ymid_i];
                ssh_interp = [ssh_interp; ssh_interp_i];
            end
        end
    end
    x_out=xE_mid;
    y_out=yE_mid;
end