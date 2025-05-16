function visualizeEddyFit(alongtrack, eddyPath_fun_t, paramsCell, window_center, t0)
% VISUALIZEEDDYFIT Creates comprehensive visualizations to debug eddy center detection
%
% Inputs:
%   alongtrack - Struct with alongtrack data (x, y, t, ssh)
%   eddyPath_fun_t - Struct with eddy path functions (xe, ye)
%   paramsCell - Cell array of fitted parameters for each time window
%   window_center - Array of center times for each window
%   t0 - Reference time (minimum time in the dataset)
%
% This function creates multiple visualizations to help diagnose eddy center detection problems:
%   1. Tracks of estimated eddy centers vs. ground truth
%   2. SSH data with fitted model overlays for each time window
%   3. Parameter evolution over time windows
%   4. Optimization cost function landscape around detected centers

%% 1. Create figure to compare eddy tracks
figure('Position', [100, 100, 1000, 800], 'Name', 'Eddy Center Tracking');

% Plot the alongtrack data points
subplot(2,2,1);
scatter(alongtrack.x, alongtrack.y, 10, alongtrack.ssh, 'filled');
hold on;
colorbar;
xlabel('X (m)');
ylabel('Y (m)');
title('Alongtrack Data with Eddy Centers');
colormap(gca, parula);

% Plot the ground truth eddy centers
days = 0:(max(alongtrack.t)-t0);
xe_truth = eddyPath_fun_t.xe(days);
ye_truth = eddyPath_fun_t.ye(days);
plot(xe_truth, ye_truth, 'k-', 'LineWidth', 2);

% Plot the fitted eddy centers for each window
colors = jet(length(paramsCell));
centers_x = zeros(length(paramsCell), 1);
centers_y = zeros(length(paramsCell), 1);
for i = 1:length(paramsCell)
    if ~isempty(paramsCell{i})
        params = paramsCell{i};
        % Calculate center at the window midpoint
        t_mid = window_center(i) - params.t0;
        centers_x(i) = params.x0 + params.cx * t_mid;
        centers_y(i) = params.y0 + params.cy * t_mid;
        plot(centers_x(i), centers_y(i), 'o', 'Color', colors(i,:), 'MarkerSize', 10, 'LineWidth', 2);
    end
end
legend('Alongtrack data', 'True eddy path', 'Fitted centers');

% Plot differences between true and fitted centers
subplot(2,2,2);
for i = 1:length(paramsCell)
    if ~isempty(paramsCell{i})
        % Find closest day to window center
        [~, idx] = min(abs(days - (window_center(i) - t0)));
        dx = centers_x(i) - xe_truth(idx);
        dy = centers_y(i) - ye_truth(idx);
        quiver(xe_truth(idx), ye_truth(idx), dx, dy, 0, 'Color', colors(i,:), 'LineWidth', 2);
        text(xe_truth(idx)+dx/2, ye_truth(idx)+dy/2, num2str(i), 'Color', colors(i,:));
    end
end
xlabel('X (m)');
ylabel('Y (m)');
title('Displacement Vectors: True to Fitted Centers');
grid on;

% Plot error distribution
subplot(2,2,3);
for i = 1:length(paramsCell)
    if ~isempty(paramsCell{i})
        % Find closest day to window center
        [~, idx] = min(abs(days - (window_center(i) - t0)));
        dx = centers_x(i) - xe_truth(idx);
        dy = centers_y(i) - ye_truth(idx);
        scatter(dx/1e3, dy/1e3, 100, colors(i,:), 'filled');
        hold on;
        text(dx/1e3+2, dy/1e3, num2str(i), 'Color', colors(i,:));
    end
end
xlabel('X error (km)');
ylabel('Y error (km)');
title('Center Detection Error Distribution');
grid on;
axis equal;
% Add circle for 10km error radius
th = 0:pi/50:2*pi;
r = 10; % 10km radius
plot(r*cos(th), r*sin(th), 'k--');

% Plot parameter evolution
subplot(2,2,4);
amplitudes = zeros(length(paramsCell), 1);
radii = zeros(length(paramsCell), 1);
velocities = zeros(length(paramsCell), 1);
for i = 1:length(paramsCell)
    if ~isempty(paramsCell{i})
        params = paramsCell{i};
        amplitudes(i) = params.A;
        radii(i) = params.L/1e3; % Convert to km
        velocities(i) = sqrt(params.cx^2 + params.cy^2);
    end
end

% Plot parameters
yyaxis left;
plot(window_center-t0, amplitudes, 'b-o', 'LineWidth', 2);
hold on;
plot(window_center-t0, radii, 'g-o', 'LineWidth', 2);
ylabel('Amplitude (m) / Radius (km)');

yyaxis right;
plot(window_center-t0, velocities, 'r-o', 'LineWidth', 2);
ylabel('Velocity (m/day)');

xlabel('Time (days)');
title('Parameter Evolution');
legend('Amplitude', 'Radius', 'Velocity');
grid on;

%% 2. Visualize SSH and model fit for selected windows
numWindows = min(length(paramsCell), 3); % Show at most 3 windows
figure('Position', [100, 100, 1200, 400*numWindows], 'Name', 'SSH Data vs Model Fits');

for i = 1:numWindows
    if isempty(paramsCell{i})
        continue;
    end
    
    params = paramsCell{i};
    
    % Extract data for this window
    t_window = window_center(i);
    window_indices = abs(alongtrack.t - t_window) <= 10; % Data around window center
    x_window = alongtrack.x(window_indices);
    y_window = alongtrack.y(window_indices);
    ssh_window = alongtrack.ssh(window_indices);
    
    % Get true center
    [~, idx] = min(abs(days - (t_window - t0)));
    xe_true = xe_truth(idx);
    ye_true = ye_truth(idx);
    
    % Get fitted center
    xe_fit = params.x0 + params.cx * (t_window - params.t0);
    ye_fit = params.y0 + params.cy * (t_window - params.t0);
    
    % Calculate model values on a grid for visualization
    [X, Y] = meshgrid(linspace(min(x_window)-50e3, max(x_window)+50e3, 100), ...
                     linspace(min(y_window)-50e3, max(y_window)+50e3, 100));
    t_rel = t_window - params.t0;
    SSH_model = params.A .* exp(-((X-params.x0-params.cx*t_rel).^2 + ...
                                  (Y-params.y0-params.cy*t_rel).^2)/params.L^2);
    
    % Plot data, model, and difference
    subplot(numWindows, 3, (i-1)*3+1);
    scatter(x_window-xe_true, y_window-ye_true, 30, ssh_window, 'filled');
    hold on;
    plot(0, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2); % True center
    plot(xe_fit-xe_true, ye_fit-ye_true, 'ro', 'MarkerSize', 12, 'LineWidth', 2); % Fitted center
    colorbar;
    xlabel('X - X_{true} (m)');
    ylabel('Y - Y_{true} (m)');
    title(sprintf('Window %d: Data (Day %.1f)', i, t_window-t0));
    axis equal;
    
    subplot(numWindows, 3, (i-1)*3+2);
    pcolor(X-xe_true, Y-ye_true, SSH_model);
    shading interp;
    hold on;
    plot(0, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2); % True center
    plot(xe_fit-xe_true, ye_fit-ye_true, 'ro', 'MarkerSize', 12, 'LineWidth', 2); % Fitted center
    colorbar;
    xlabel('X - X_{true} (m)');
    ylabel('Y - Y_{true} (m)');
    title(sprintf('Window %d: Model', i));
    axis equal;
    
    % Plot cost function landscape
    subplot(numWindows, 3, (i-1)*3+3);
    visualizeCostFunctionLandscape(x_window, y_window, ssh_window, t_window, params, xe_true, ye_true);
    xlabel('X - X_{true} (km)');
    ylabel('Y - Y_{true} (km)');
    title(sprintf('Window %d: Cost Function Landscape', i));
end

% Adjust colormaps for better visualization
colormaps = {flipud(brewermap(256, 'RdBu')), ...
             brewermap(256, 'YlGnBu'), ...
             brewermap(256, 'OrRd')};
for i = 1:numWindows*3
    if i <= length(colormaps)
        subplot(numWindows, 3, i);
        try
            colormap(gca, colormaps{mod(i-1,3)+1});
        catch
            colormap(gca, jet);
        end
    end
end

end

function visualizeCostFunctionLandscape(x, y, ssh, t, params, xe_true, ye_true)
% VISUALIZECOSTFUNCTIONLANDSCAPE Creates a visualization of the cost function landscape
% around the fitted eddy center to help diagnose optimization issues

% Create a grid of potential center positions
[X, Y] = meshgrid(linspace(xe_true-50e3, xe_true+50e3, 50), ...
                 linspace(ye_true-50e3, ye_true+50e3, 50));

% Calculate cost function value at each position
cost = zeros(size(X));
t_rel = t - params.t0;
for i = 1:numel(X)
    x0 = X(i);
    y0 = Y(i);
    % Calculate model values with these center coordinates
    ssh_model = params.A .* exp(-((x-x0-params.cx*t_rel).^2 + ...
                               (y-y0-params.cy*t_rel).^2)/params.L^2);
    % Calculate MSE
    cost(i) = mean((ssh - ssh_model).^2);
end

% Normalize for visualization
cost = log10(cost);

% Plot the landscape
pcolor((X-xe_true)/1e3, (Y-ye_true)/1e3, cost);
shading interp;
hold on;

% Plot the true center
plot(0, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2);

% Plot the fitted center
xe_fit = params.x0 + params.cx * t_rel;
ye_fit = params.y0 + params.cy * t_rel;
plot((xe_fit-xe_true)/1e3, (ye_fit-ye_true)/1e3, 'ro', 'MarkerSize', 12, 'LineWidth', 2);

% Add colorbar
colorbar;
title('Log10(MSE)');
end

function visualizeAlongtrackCoverage(alongtrack, eddyPath_fun_t, t0)
% VISUALIZEALONGTRACKCORVERAGE Creates visualizations to assess alongtrack data coverage
% relative to the eddy position and structure

% Get eddy positions for each day
days = 0:(max(alongtrack.t)-t0);
xe_truth = eddyPath_fun_t.xe(days);
ye_truth = eddyPath_fun_t.ye(days);

% Create figure
figure('Position', [100, 100, 1200, 800], 'Name', 'Alongtrack Data Coverage');

% Plot the density of measurements in eddy-centered coordinates
subplot(2,2,1);
x_rel = zeros(length(alongtrack.x), 1);
y_rel = zeros(length(alongtrack.y), 1);

% Calculate eddy-relative coordinates for each measurement
for i = 1:length(alongtrack.t)
    day_idx = round(alongtrack.t(i) - t0) + 1;
    if day_idx > 0 && day_idx <= length(days)
        x_rel(i) = alongtrack.x(i) - xe_truth(day_idx);
        y_rel(i) = alongtrack.y(i) - ye_truth(day_idx);
    end
end

% Create histogram of measurements in eddy-relative frame
histogram2(x_rel/1e3, y_rel/1e3, 50, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on');
xlabel('X - X_{eddy} (km)');
ylabel('Y - Y_{eddy} (km)');
title('Data Density in Eddy-Centered Frame');
colorbar;
axis equal;

% Add circle showing typical eddy radius (100km)
hold on;
th = 0:pi/50:2*pi;
r = 100; % 100km radius
plot(r*cos(th), r*sin(th), 'r--', 'LineWidth', 2);
r = 200; % 200km radius
plot(r*cos(th), r*sin(th), 'r:', 'LineWidth', 1);

% Plot coverage as a function of time
subplot(2,2,2);
days_unique = unique(floor(alongtrack.t - t0));
counts = zeros(length(days_unique), 1);
for i = 1:length(days_unique)
    counts(i) = sum(floor(alongtrack.t - t0) == days_unique(i));
end
bar(days_unique, counts);
xlabel('Day');
ylabel('Number of measurements');
title('Temporal Coverage');
grid on;

% Plot data availability in azimuth vs radius
subplot(2,2,3);
r = sqrt(x_rel.^2 + y_rel.^2)/1e3; % km
theta = atan2(y_rel, x_rel);
polarscatter(theta, r, 10, alongtrack.ssh, 'filled');
title('Data in Polar Coordinates (Eddy-Centered)');
colorbar;

% Plot average data density as a function of radius
subplot(2,2,4);
r_bins = 0:10:300; % 10km bins
counts = histcounts(r, r_bins);
area_ring = pi * diff((r_bins).^2); % area of each annular ring in km^2
density = counts ./ area_ring; % points per km^2
bar(r_bins(1:end-1)+5, density);
xlabel('Distance from Eddy Center (km)');
ylabel('Measurement Density (points/km²)');
title('Radial Data Density');
grid on;
xlim([0, 300]);

end

function visualizeModelEvolution(alongtrack, eddyPath_fun_t, paramsCell, window_center, t0)
% VISUALIZEMODELVALUES Creates an animation showing how the eddy model evolves
% compared to the data over time

% Set up the figure
figure('Position', [100, 100, 800, 800], 'Name', 'Eddy Model Evolution');

% Get eddy positions for each day
days = 0:(max(alongtrack.t)-t0);
xe_truth = eddyPath_fun_t.xe(days);
ye_truth = eddyPath_fun_t.ye(days);

% Determine the day step for the animation (show ~20 frames)
day_step = max(1, floor(length(days)/20));
days_to_show = 1:day_step:length(days);

% Create a function to get the model for any day
getModelForDay = @(day) createInterpolatedModel(day, window_center, paramsCell, t0);

% Precompute the model grid
x_range = [min(alongtrack.x), max(alongtrack.x)];
y_range = [min(alongtrack.y), max(alongtrack.y)];
buffer = 50e3; % 50km buffer
[X, Y] = meshgrid(linspace(x_range(1)-buffer, x_range(2)+buffer, 100), ...
                 linspace(y_range(1)-buffer, y_range(2)+buffer, 100));

% Animation loop
for i = 1:length(days_to_show)
    day = days_to_show(i);
    
    % Get current true eddy center
    xe = xe_truth(day);
    ye = ye_truth(day);
    
    % Get data for this day (within a 2-day window)
    day_indices = abs(alongtrack.t - (day + t0)) <= 2;
    x_day = alongtrack.x(day_indices);
    y_day = alongtrack.y(day_indices);
    ssh_day = alongtrack.ssh(day_indices);
    
    % Get the model for this day
    model_params = getModelForDay(day);
    if ~isempty(model_params)
        t_rel = day;  % Days since t0
        SSH_model = model_params.A .* exp(-((X-model_params.x0-model_params.cx*t_rel).^2 + ...
                                        (Y-model_params.y0-model_params.cy*t_rel).^2)/model_params.L^2);
        
        % Plot model
        pcolor(X-xe, Y-ye, SSH_model);
        shading interp;
        hold on;
        
        % Plot true center
        plot(0, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2);
        
        % Plot fitted center
        xe_fit = model_params.x0 + model_params.cx * t_rel;
        ye_fit = model_params.y0 + model_params.cy * t_rel;
        plot(xe_fit-xe, ye_fit-ye, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
        
        % Plot data points from this day
        if ~isempty(x_day)
            scatter(x_day-xe, y_day-ye, 30, ssh_day, 'filled', 'MarkerEdgeColor', 'k');
        end
        
        % Add informative text
        text(-0.9*max(abs(X(:)-xe)), 0.9*max(abs(Y(:)-ye)), ...
             sprintf('Day %d\nA = %.2f m\nL = %.1f km\nError = %.1f km', ...
             day, model_params.A, model_params.L/1e3, ...
             sqrt((xe_fit-xe)^2+(ye_fit-ye)^2)/1e3), ...
             'BackgroundColor', 'white', 'EdgeColor', 'black');
        
        % Format plot
        xlabel('X - X_{true} (m)');
        ylabel('Y - Y_{true} (m)');
        title(sprintf('Eddy Model vs Data (Day %d)', day));
        axis equal;
        colorbar;
        colormap(flipud(brewermap(256, 'RdBu')));
        caxis([-max(abs(SSH_model(:))), max(abs(SSH_model(:)))]);
        
        % Set consistent axis limits
        limit = 200e3;  % 200km around the eddy center
        xlim([-limit, limit]);
        ylim([-limit, limit]);
        
        hold off;
        
        % Pause briefly for animation effect
        pause(0.2);
        
        % Clear for next frame if not the last one
        if i < length(days_to_show)
            clf;
        end
    end
end
end

function model_params = createInterpolatedModel(day, window_center, paramsCell, t0)
% CREATEINTERPOLATEDMODEL Returns the model parameters for a specific day
% by interpolating between window centers

% Convert day to absolute time
t_day = day + t0;

% Find the two closest windows
[~, idx] = sort(abs(window_center - t_day));
if length(idx) < 2 || isempty(paramsCell{idx(1)}) || isempty(paramsCell{idx(2)})
    % If we don't have two valid windows, return the closest one
    if ~isempty(paramsCell{idx(1)})
        model_params = paramsCell{idx(1)};
    else
        model_params = [];
    end
    return;
end

% Get the two closest windows
idx1 = idx(1);
idx2 = idx(2);

% Calculate weights for interpolation (inverse distance weighting)
d1 = abs(window_center(idx1) - t_day);
d2 = abs(window_center(idx2) - t_day);
w1 = d2 / (d1 + d2);
w2 = d1 / (d1 + d2);

% Interpolate parameters
params1 = paramsCell{idx1};
params2 = paramsCell{idx2};

model_params = struct();
model_params.A = w1 * params1.A + w2 * params2.A;
model_params.L = w1 * params1.L + w2 * params2.L;
model_params.cx = w1 * params1.cx + w2 * params2.cx;
model_params.cy = w1 * params1.cy + w2 * params2.cy;

% The x0, y0 parameters need special handling because they're initial positions
% Calculate the actual centers at the window centers
t1_rel = window_center(idx1) - params1.t0;
x1 = params1.x0 + params1.cx * t1_rel;
y1 = params1.y0 + params1.cy * t1_rel;

t2_rel = window_center(idx2) - params2.t0;
x2 = params2.x0 + params2.cx * t2_rel;
y2 = params2.y0 + params2.cy * t2_rel;

% Interpolate the centers
x_center = w1 * x1 + w2 * x2;
y_center = w1 * y1 + w2 * y2;

% Set t0 to the day we're evaluating at
model_params.t0 = t_day;

% Calculate new x0, y0 that would give the interpolated center at t_day
model_params.x0 = x_center;
model_params.y0 = y_center;
end
