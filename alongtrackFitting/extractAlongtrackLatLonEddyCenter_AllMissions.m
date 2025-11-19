function alongtrackLatLon = extractAlongtrackLatLonEddyCenter_AllMissions(alongtrackSimulator, eddyPath, timeo, options)
% Extract alongtrack lat/lon/time from all missions for eddy tracking
%
% Inputs:
%   alongtrackSimulator - AlongTrackSimulator object
%   eddyPath - Structure with eddy path information (lat, lon fields)
%   timeo - Times corresponding to eddy positions (MATLAB datenum format)
%   options.radius - Radius in km around eddy center (default: 300)
%   options.getSSH - Whether to get SSH data (default: false, not implemented for simulator)
%
% Output:
%   alongtrackLatLon - Structure with fields:
%       .t - Time (seconds from epoch)
%       .lon - Longitude
%       .lat - Latitude
%       .mission - Mission name for each point

arguments
    alongtrackSimulator AlongTrackSimulator
    eddyPath struct
    timeo (:,:) {mustBeNumeric, mustBePositive}
    options.radius (1,1) {mustBeNumeric} = 300
    options.getSSH = false
    options.groupByBaseMission (1,1) logical = false  % Combine mission phases
    options.debugMission string = ""  % Debug a specific mission (e.g., "enn")
end

radius = options.radius;
getSSH = options.getSSH;

% Find eddy center bounds
latc = (max(eddyPath.lat) + min(eddyPath.lat)) / 2;
lonc = (max(eddyPath.lon) + min(eddyPath.lon)) / 2;
[lat_rad, lon_rad] = xy2latlon(radius, radius, latc, lonc);
lat_rad = lat_rad - latc;
lon_rad = lon_rad - lonc;

% Initialize output structure
alongtrackLatLon.t = [];
alongtrackLatLon.lon = [];
alongtrackLatLon.lat = [];
alongtrackLatLon.mission = [];
if getSSH
    alongtrackLatLon.ssh = [];
    warning('SSH retrieval not implemented for simulator-based extraction');
end

% Get all available missions
allMissions = alongtrackSimulator.missions();

% Convert timeo (datenum) to seconds
% Assuming some reference epoch - you may need to adjust this
t_min = min(timeo);
t_max = max(timeo);

% Create time vector covering the entire eddy lifetime with daily resolution
time_days = t_min:1:t_max; % One point per day

fprintf('Processing %d days across eddy lifetime\n', length(time_days));
fprintf('Eddy center: lat=%.2f, lon=%.2f\n', latc, lonc);
fprintf('Search radius: %.1f km (lat ±%.2f°, lon ±%.2f°)\n', radius, lat_rad, lon_rad);

% Loop through each day in eddy lifetime
for day_idx = 1:length(time_days)
    current_day = time_days(day_idx);
    
    % Find nearest eddy position for this day
    [~, nearest_idx] = min(abs(timeo - current_day));
    eddy_lat = eddyPath.lat(nearest_idx);
    eddy_lon = eddyPath.lon(nearest_idx);
    
    % Define search region around eddy center
    latg = eddy_lat - round(lat_rad):0.25:eddy_lat + round(lat_rad);
    long = eddy_lon - round(lon_rad):0.25:eddy_lon + round(lon_rad);
    
    % Create a region box
    region = [min(long(:)), max(long(:)), min(latg(:)), max(latg(:))];
    minLon = region(1);
    maxLon = region(2);
    minLat = region(3);
    maxLat = region(4);

    % Check if region crosses ±180° longitude boundary
    % This happens when the search box wraps around (e.g., eddy near dateline)
    % Similar to your mapped data extraction's wrap-around handling
    crosses_dateline = (maxLon - minLon) > 180;
    
    
    % Convert day to seconds for time window (24-hour window)
    time_start = (current_day - t_min) * 86400; % seconds since start
    time_end = time_start + 86400; % 24 hours later
    
    % Loop through all missions
    for mission_idx = 1:length(allMissions)
        mission_name = allMissions(mission_idx);
        
        try
            % Get ground track for this mission for the day
            % Create time vector covering several orbits to ensure coverage
            T_orbit = alongtrackSimulator.orbitalPeriodForMissionWithName(mission_name);
            N_orbits = ceil(24*60*60 / T_orbit) + 1; % Number of orbits in 24 hours + buffer
            time_vec = time_start:1:time_end;
            
            % Check if mission was active during this time
            % Handle missions with multiple phases (arrays of start/end dates)
            mission_params = alongtrackSimulator.missionParameters(mission_name);
            mission_start = datenum(mission_params.start_date);
            mission_end = datenum(mission_params.end_date);
            
            % Check if current_day falls within any phase
            is_active = false;
            if isscalar(mission_start)
                % Single phase mission
                if current_day >= mission_start && (isinf(mission_end) || current_day <= mission_end)
                    is_active = true;
                end
            else
                % Multi-phase mission - check each phase
                for phase_idx = 1:length(mission_start)
                    if current_day >= mission_start(phase_idx) && ...
                       (isinf(mission_end(phase_idx)) || current_day <= mission_end(phase_idx))
                        is_active = true;
                        break;
                    end
                end
            end
            
            if ~is_active
                continue;
            end
            
            % Get ground track
            [lat, lon, time] = alongtrackSimulator.groundTrackForMissionWithName(...
                mission_name, time=time_vec);
            
            % Ensure all arrays are column vectors
            lat = lat(:);
            lon = lon(:);
            time = time(:);
            
            % Filter points within the region (handle dateline crossing)
            if crosses_dateline
                % Eddy crosses ±180° - use OR logic for longitude
                % Points are valid if they're in eastern OR western part of wrapped region
                within_box = (lat >= minLat) & (lat <= maxLat) & ...
                            ((lon >= minLon) | (lon <= maxLon));
            else
                % Normal case - standard bounding box
                within_box = (lat >= minLat) & (lat <= maxLat) & ...
                            (lon >= minLon) & (lon <= maxLon);
            end
            
            if any(within_box)
                % Extract points within box
                lat_in = lat(within_box);
                lon_in = lon(within_box);
                time_in = time(within_box);
                
                % Further filter by distance from eddy center using haversine
                dist_km = haversine_distance(eddy_lat, eddy_lon, lat_in, lon_in);
                within_radius = dist_km <= radius;
                
                if any(within_radius)
                    % Concatenate results
                    n_points = sum(within_radius);
                    alongtrackLatLon.t = [alongtrackLatLon.t; time_in(within_radius)];
                    alongtrackLatLon.lon = [alongtrackLatLon.lon; lon_in(within_radius)];
                    alongtrackLatLon.lat = [alongtrackLatLon.lat; lat_in(within_radius)];
                    
                    % Store mission name (optionally group by base mission)
                    if options.groupByBaseMission
                        % Extract base mission name (e.g., "ERS-1" from "ERS-1 35-day cycle")
                        base_mission = extractBaseMissionName(mission_name);
                        alongtrackLatLon.mission = [alongtrackLatLon.mission; ...
                            repmat(base_mission, n_points, 1)];
                    else
                        alongtrackLatLon.mission = [alongtrackLatLon.mission; ...
                            repmat(mission_name, n_points, 1)];
                    end
                    
                    %fprintf('  Day %d/%d: %s - found %d points\n', ...
                        % day_idx, length(time_days), mission_name, n_points);
                end
            end
        catch ME
            % Skip missions that cause errors but provide more detail
            if contains(ME.message, 'Dimensions') || contains(ME.message, 'concatenat')
                fprintf('  Warning: Dimension mismatch for mission %s (day %d): %s\n', ...
                    mission_name, day_idx, ME.message);
            else
                fprintf('  Warning: Error processing mission %s (day %d): %s\n', ...
                    mission_name, day_idx, ME.message);
            end
        end
    end
end

% Sort by time
[alongtrackLatLon.t, sort_idx] = sort(alongtrackLatLon.t);
alongtrackLatLon.lon = alongtrackLatLon.lon(sort_idx);
alongtrackLatLon.lat = alongtrackLatLon.lat(sort_idx);
alongtrackLatLon.mission = alongtrackLatLon.mission(sort_idx);

fprintf('\nTotal points found: %d\n', length(alongtrackLatLon.t));
fprintf('Unique missions: %s\n', strjoin(unique(alongtrackLatLon.mission), ', '));

end


function base_name = extractBaseMissionName(mission_name)
% Extract base mission name by removing phase descriptors
% Examples:
%   "ERS-1 35-day cycle" -> "ERS-1"
%   "ERS-1 168-day cycle" -> "ERS-1"
%   "Jason-1" -> "Jason-1"

% Common patterns to remove
patterns = {
    ' 35-day cycle'
    ' 168-day cycle'
    ' 3-day cycle'
    ' multidisciplinary phase'
    ' geodetic phase'
    ' ice phase'
    ' Phase A'
    ' Phase B'
    ' Phase C'
    ' Phase D'
    ' Phase E'
    ' Phase F'
};

base_name = mission_name;
for i = 1:length(patterns)
    base_name = strrep(base_name, patterns{i}, '');
end

% Trim any remaining whitespace
base_name = strtrim(base_name);
end


function dist_km = haversine_distance(lat1, lon1, lat2, lon2)
% Calculate distance in km between points on Earth using Haversine formula
% Inputs can be scalars or arrays (lat2, lon2 can be arrays)
    R = 6371; % Earth radius in km
    lat1_rad = deg2rad(lat1);
    lat2_rad = deg2rad(lat2);
    delta_lat = deg2rad(lat2 - lat1);
    delta_lon = deg2rad(lon2 - lon1);
    
    a = sin(delta_lat/2).^2 + cos(lat1_rad) .* cos(lat2_rad) .* sin(delta_lon/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    dist_km = R * c;
end