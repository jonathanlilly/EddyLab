function plotTrackandAVISO(alongtrack, eddyPath, mapped_field, radius, options)
% plotAlongtrack Visualizes along-track data, eddy paths, and mapped AVISO field
%   This function creates plots of alongtrack data, eddy paths, and gridded AVISO.
%   Can handle both lat/lon and x/y coordinate systems.
%
% Inputs:
%   alongtrack - Struct containing along-track data with fields:
%                Either (lon, lat, t) or (x, y, t)
%                Optional: ssh field for third subplot
%   eddyPath - Struct containing eddy path:
%              Either (lon, lat) or (xe, ye) as arrays or function handles
%   mapped_field - Struct containing gridded AVISO data:
%                  (time, lat, lon) and filename
%
arguments
    alongtrack struct
    eddyPath struct
    mapped_field struct
    radius double
    options.time_window (1,:) {mustBeNumeric} = []
    options.snapshot_index (1,1) {mustBeNumeric} = 1
    options.lon_bounds (1,2) {mustBeNumeric} = [-40, 30]
    options.lat_bounds (1,2) {mustBeNumeric} = [-40, -20]
    options.eddy_id string
end

% Check if SSH data is available
has_ssh = isfield(alongtrack, 'ssh');

t0 = min(alongtrack.t);
timeo = [floor(t0):max(floor(alongtrack.t))+1];

if ~isempty(options.time_window)
    valid_time = logical(sum(floor(alongtrack.t)==timeo(options.time_window), 2));
    alongtrack.t = alongtrack.t(valid_time);
    if has_ssh
        alongtrack.ssh = alongtrack.ssh(valid_time);
    end
    if isfield(alongtrack, 'lat') && isfield(alongtrack, 'lon')
        alongtrack.lat = alongtrack.lat(valid_time);
        alongtrack.lon = alongtrack.lon(valid_time);
    end
    if isfield(alongtrack, 'x') && isfield(alongtrack, 'y')
        alongtrack.x = alongtrack.x(valid_time);
        alongtrack.y = alongtrack.y(valid_time);
    end
end

% Calculate eddy positions at each alongtrack time
elapsed_time = alongtrack.t - t0;

% Initialize eddytrack struct
eddytrack = struct();

% Determine coordinate system (lat/lon or x/y)
if isfield(alongtrack, 'lat') && isfield(alongtrack, 'lon')
    isLatLon = true;
    % Working in lat/lon coordinates
    if isfield(eddyPath, 'lat') && isfield(eddyPath, 'lon')
        % Both alongtrack and eddyPath in lat/lon
        if isa(eddyPath.lon, 'function_handle')
            eddytrack.longitude = eddyPath.lon(elapsed_time);
            eddytrack.latitude = eddyPath.lat(elapsed_time);
        else
            eddytrack.longitude = eddyPath.lon;
            eddytrack.latitude = eddyPath.lat;
        end
    else
        if isa(eddyPath.xe, 'function_handle')
            eddytrack.x = eddyPath.xe(elapsed_time);
            eddytrack.y = eddyPath.ye(elapsed_time);
        else
            eddytrack.x = eddyPath.xe;
            eddytrack.y = eddyPath.ye;
        end
        % lato and lono are the center of the alongtrack domain
        lono = (min(alongtrack.lon(:))+max(alongtrack.lon(:)))/2;
        lato = (min(alongtrack.lat(:))+max(alongtrack.lat(:)))/2;
    
        [eddytrack.latitude, eddytrack.longitude] = xy2latlon(eddytrack.x/1e3, eddytrack.y/1e3, lato, lono);
        eddytrack.longitude = deg180(eddytrack.longitude-lono) + lono; %avoids unwrapping issues
    end
    
    % Prepare for plotting in lat/lon
    plot_x = alongtrack.lon;
    plot_y = alongtrack.lat;
    xlabel_str = 'Longitude';
    ylabel_str = 'Latitude';
    eddy_x = eddytrack.longitude;
    eddy_y = eddytrack.latitude;
    
elseif isfield(alongtrack, 'x') && isfield(alongtrack, 'y')
    isLatLon = false;
    % Working in x/y coordinates
    if isfield(eddyPath, 'xe') && isfield(eddyPath, 'ye')
        % Both alongtrack and eddyPath in x/y
        if isa(eddyPath.xe, 'function_handle')
            eddytrack.x = eddyPath.xe(elapsed_time);
            eddytrack.y = eddyPath.ye(elapsed_time);
        else
            eddytrack.x = eddyPath.xe;
            eddytrack.y = eddyPath.ye;
        end
    else
        error('Coordinate system mismatch: alongtrack is in x/y but eddyPath is not');
    end
    
    % Convert meters to kilometers for plotting
    eddytrack.x_km = eddytrack.x / 1e3;
    eddytrack.y_km = eddytrack.y / 1e3;
    
    % Prepare for plotting in km
    plot_x = alongtrack.x / 1e3; % Convert to km
    plot_y = alongtrack.y / 1e3; % Convert to km
    xlabel_str = '$x$ (km)';
    ylabel_str = '$y$ (km)';
    eddy_x = eddytrack.x_km;
    eddy_y = eddytrack.y_km;
else
    error('alongtrack must contain either lat/lon or x/y coordinate fields');
end

% Create color map for missions
uniqueMissions = unique(alongtrack.mission);
n_missions = length(uniqueMissions);
colors = lines(n_missions);  % Use distinct colors for each mission

% Create figure with appropriate number of subplots
if has_ssh
    n_subplots = 3;
    figure('Position', [100, 100, 650, 1000]); 
else
    n_subplots = 2;
    figure('Position', [100, 100, 650, 550]); 
end

% First subplot: Track visualization
ax1 = subplot(n_subplots, 1, 1);
hold on
% add a title
if isfield(options, 'eddy_id')
title(['Eddy' options.eddy_id], 'interpreter', 'tex')
end

for i = 1:n_missions
    mission_mask = strcmp(alongtrack.mission, uniqueMissions(i));
    h(1+i) = scatter(alongtrack.lon(mission_mask), alongtrack.lat(mission_mask), ...
        2, colors(i,:), 'filled', 'DisplayName', char(uniqueMissions(i)));
    h(1+i).MarkerFaceAlpha = 1;
end
xlabel(xlabel_str, 'FontName', 'times', 'fontsize', 16, 'Interpreter', 'latex')
ylabel(ylabel_str, 'FontName', 'times', 'fontsize', 14, 'Interpreter', 'latex')
h(1) = plot(eddy_x, eddy_y, 'LineWidth', 3, 'Color', 0*[1 1 1]);
plot(eddy_x(1), eddy_y(1), 'ko', 'markersize', 8, 'markerfacecolor', 'k')
plot(eddy_x(1), eddy_y(1), 'wo', 'markersize', 5, 'markerfacecolor', 'w')
plot(eddy_x(end), eddy_y(end), 'wx', 'markersize', 10, 'linewidth', 2)
plot(eddy_x(end), eddy_y(end), 'kx', 'markersize', 8, 'linewidth', 1.5)
box on
axis tight

if isLatLon
    set(gca, 'dataaspectratio', [1 cosd(mean(get(gca, 'ylim'))) 1])
    xlim(options.lon_bounds), ylim(options.lat_bounds)
    topoplot
    latratio(30)
else
    set(gca, 'dataaspectratio', [1 1 1])
end

set(gca, 'fontname', 'times', 'fontsize', 16)
legend(h, ['Eddy path'; uniqueMissions], 'location', 'northoutside', 'interpreter', 'tex','NumColumns',4)

% Store position of first subplot for alignment
pos1=[0.11,0.48,0.72,0.29];
set(ax1, 'Position', pos1);
%get(ax1, 'Position');

% Second subplot: Mapped AVISO field
ax2 = subplot(n_subplots, 1, 2);
hold on

% Select time snapshot
n = options.snapshot_index;
if n > length(timeo)
    n = length(timeo);
end

% Find time index closest to selected time
[~, time_idx] = min(abs(mapped_field.time - timeo(n)));

% Find indices for lat/lon bounds
[~, start_lon] = min(abs(mapped_field.lon - options.lon_bounds(1)));
[~, last_lon] = min(abs(mapped_field.lon - options.lon_bounds(2)));
[~, start_lat] = min(abs(mapped_field.lat - options.lat_bounds(1)));
[~, last_lat] = min(abs(mapped_field.lat - options.lat_bounds(2)));

% Handle longitude wrap-around
if start_lon > last_lon
    % First part: end of longitude array
    start1 = [start_lon, start_lat, time_idx];
    count1 = [length(mapped_field.lon) - start_lon + 1, last_lat - start_lat + 1, 1];
    ssh_part1 = ncread(mapped_field.filename, 'sla', start1, count1) * 100; % Convert to cm
    lon_part1 = mapped_field.lon(start_lon:end);
    
    % Second part: beginning of longitude array
    start2 = [1, start_lat, time_idx];
    count2 = [last_lon, last_lat - start_lat + 1, 1];
    ssh_part2 = ncread(mapped_field.filename, 'sla', start2, count2) * 100; % Convert to cm
    lon_part2 = mapped_field.lon(1:last_lon);
    
    % Combine
    plot_ssh = cat(1, ssh_part1, ssh_part2);
    plot_lon = [lon_part1; lon_part2];
else
    % Normal case
    start = [start_lon, start_lat, time_idx];
    count = [last_lon - start_lon + 1, last_lat - start_lat + 1, 1];
    plot_ssh = ncread(mapped_field.filename, 'sla', start, count) * 100; % Convert to cm
    plot_lon = mapped_field.lon(start_lon:last_lon);
end
plot_lat = mapped_field.lat(start_lat:last_lat);

% Create the mapped field plot
jpcolor(plot_lon, plot_lat, plot_ssh')
c = colorbar('EastOutside');
c.Label.String = 'Mapped SSH (cm)';
c.Label.FontSize = 12;
c.Label.FontName = 'times';
colormap(gca, brewermap([], '-Spectral'))
shading flat

% Plot eddy track on top
plot(eddy_x, eddy_y, 'LineWidth', 3, 'Color', 0*[1 1 1]);
plot(eddy_x(1), eddy_y(1), 'ko', 'markersize', 8, 'markerfacecolor', 'k')
plot(eddy_x(1), eddy_y(1), 'wo', 'markersize', 5, 'markerfacecolor', 'w')
plot(eddy_x(end), eddy_y(end), 'wx', 'markersize', 10, 'linewidth', 2)
plot(eddy_x(end), eddy_y(end), 'kx', 'markersize', 8, 'linewidth', 1.5)

% Plot a circle
theta = linspace(0,2*pi,100);
circle_x = radius * cos(theta);
circle_y = radius * sin(theta);

if isLatLon
    [circle_lat, circle_lon] = xy2latlon(circle_x, circle_y, eddy_y(options.snapshot_index), eddy_x(options.snapshot_index));
    plot(circle_lon,circle_lat,'linewidth',2,'Color','k')
else
    plot(eddy_x(options.snapshot_index)+circle_x,eddy_y(options.snapshot_index)+circle_y, ...
    'linewidth',3,'Color','k')
end
set(gca, 'dataaspectratio', [1 cosd(mean(get(gca,'ylim'))) 1])
topoplot
latratio(30)
xlim(options.lon_bounds), ylim(options.lat_bounds)
clim([-20, 40])
set(gca, 'fontname', 'times', 'fontsize', 16)
xlabel(xlabel_str, 'FontName', 'times', 'Interpreter', 'latex')
ylabel(ylabel_str, 'FontName', 'times', 'Interpreter', 'latex')

% Align second subplot with first subplot6399
% pos2 = get(ax2, 'Position');
% pos2(1) = pos1(1);  % Match left edge
% pos2(3) = pos1(3);  % Match width
pos2=[0.11,0.11,0.72,0.29];
set(ax2, 'Position', pos2);

% Third subplot: Along-track SSH (only if SSH data exists)
if has_ssh
    ax3 = subplot(3, 1, 3);
    cmap = brewermap(256, '-Spectral');
    
    s = scatter3(plot_x, plot_y, alongtrack.ssh*1e2, 7, alongtrack.ssh*1e2, 'filled', ...
        'markerEdgeColor', 'none');
    set(gca, 'fontname', 'times', 'fontsize', 14)
    colormap(gca, cmap)
    c = colorbar('EastOutside');
    c.Label.String = 'SSH (cm)';
    c.Label.FontSize = 16;
    c.Label.FontName = 'times';
    
    caxis([min(alongtrack.ssh(:)*1e2), max(alongtrack.ssh(:)*1e2)])
    xlabel(xlabel_str, 'FontName', 'times', 'Interpreter', 'latex')
    ylabel(ylabel_str, 'FontName', 'times', 'Interpreter', 'latex')
    zlabel('SSH (cm)', 'FontName', 'times', 'Interpreter', 'latex')
    grid on
    xlim([min(plot_x), max(plot_x)]);
    ylim([min(plot_y), max(plot_y)]);
    view(-30, 35)
    
    % Align third subplot with first subplot
    pos3 = get(ax3, 'Position');
    pos3(1) = pos1(1);  % Match left edge
    pos3(3) = pos1(3);  % Match width
    set(ax3, 'Position', pos3);
end

end