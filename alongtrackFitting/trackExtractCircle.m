function [latt, lont, varargout] = trackExtractCircle(lat, lon, varargin)
% region = [min(long(:)), max(long(:)), min(latg(:)), max(latg(:))];

% Get region from the last argument
    region = varargin{end};
    varargin = varargin(1:end-1); % Remove region from varargin

% After calling trackextract with a rectangular region
[latt, lont, varargout{1:length(varargin)}] = trackextract(lat, lon, varargin{:}, region);

% Define circle parameters
circleCenterLat = (region(3)+region(4))/2;  % Latitude of circle center
circleCenterLon = (region(1)+region(2))/2; % Longitude of circle center
circleRadius = max(region(2)-region(1),region(4)-region(3));      % Radius in degrees

% Calculate distance from each point to circle center
% Note: For more accuracy, you should use great circle distance
distances = sqrt((latt - circleCenterLat).^2 + (lont - circleCenterLon).^2);

% Create mask for points within the circle
circularMask = distances <= circleRadius;

% Set points outside circle to NaN
latt(~circularMask) = NaN;
lont(~circularMask) = NaN;

% Apply mask to all additional outputs (time, SSH, etc.)
for i = 1:length(varargin)
    temp = varargout{i};
    for j = 1:size(temp, 3) % Handle 3D arrays like SSH
        temp(:,:,j) = temp(:,:,j);
        temp(~circularMask,:,j) = NaN;
    end
    varargout{i} = temp;
end

end