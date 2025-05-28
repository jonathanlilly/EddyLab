function [ssh_interp, XGrid, YGrid] = interpEddyCentric(x, y, t, xo, yo, ssh, bin_size, max_r)
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
end

elseif data_type == "alongtrack"
    % Scattered data interpolation using scatteredInterpolant
    t_days=sort(unique(floor(t)),'ascend');
    
    % Pre-allocate for interpolated data
    ssh_interp = zeros(length(xmid), length(xmid), length(t_days));
    % Loop over time to interpolate scattered data
    for i = 1:length(t_days)
        % Calculate eddy-relative coordinates for this time step
        xE = x - xo(i);
        yE = y - yo(i);
        current_idx=find(floor(t)==t_days(i));
        ssh_i = ssh(current_idx);
        xE_i = xE(current_idx);
        yE_i = yE(current_idx);
        
        % Find valid data points for this time step
        valid_idx = ~isnan(ssh_i) & ~isnan(xE_i) & ~isnan(yE_i);
        
        if sum(valid_idx) > 3 % Need at least 3 points for interpolation
            F = scatteredInterpolant(xE_i(valid_idx), yE_i(valid_idx), ssh_i(valid_idx), 'linear', 'none');
            ssh_interp(:,:,i) = F(XGrid, YGrid);
        else
            % Not enough valid points, leave as zeros/NaN
            ssh_interp(:,:,i) = NaN(size(XGrid));
        end
    end
    
else
    error('data_type must be either "grid" or "scattered"');
end
end