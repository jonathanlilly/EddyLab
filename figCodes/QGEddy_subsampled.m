%% subplot 1: eddy field
n=300;
cmap = brewermap(256, '-Spectral');
figure;
subplot(3,1,1)
jpcolor(eddy_field.x/1e3, eddy_field.y/1e3, eddy_field.ssh(:,:,n)'*1e2);
shading flat
axis equal tight
xlabel('$x$ (km)', 'FontName', 'times', 'Interpreter', 'latex')
ylabel('$y$ (km)', 'FontName', 'times', 'Interpreter', 'latex')
set(gca,'xticklabel',[],'xlabel',[])
set(gca, 'fontname', 'times','fontsize',12)
colormap(cmap)
c = colorbar('EastOutside');
% c.Label.String = 'Time-averaged SSH (cm)';
c.Label.FontSize=16;
c.Label.FontName='times';

clim([-2.5, 12])

% Add eddy path
hold on;
eddy_x_km = eddyPath.xe(1:n)/1000;
eddy_y_km = eddyPath.ye(1:n)/1000;
plot(eddy_x_km, eddy_y_km, 'k-', 'LineWidth', 2);
plot(eddy_x_km(end), eddy_y_km(end), 'ko', 'linewidth',1.5,'MarkerSize', 8, 'MarkerFaceColor', 'none');

subplot(3,1,2);hold on
alongtrack=currentMission;
time_window=n;
t0=min(alongtrack.t);
timeo = [floor(t0):max(floor(alongtrack.t))+1];
if ~isempty(time_window)
    valid_time=logical(sum(floor(alongtrack.t)==timeo(time_window),2));
    alongtrack.t=alongtrack.t(valid_time);
    alongtrack.ssh=alongtrack.ssh(valid_time);
    if isfield(alongtrack, 'lat') && isfield(alongtrack, 'lon')
    alongtrack.lat=alongtrack.lat(valid_time);
    alongtrack.lon=alongtrack.lon(valid_time);
    end
    if isfield(alongtrack, 'x') && isfield(alongtrack, 'y')
    alongtrack.x=alongtrack.x(valid_time);
    alongtrack.y=alongtrack.y(valid_time);
    end
end

p=jpcolor(eddy_field.x/1e3, eddy_field.y/1e3, eddy_field.ssh(:,:,n)'*1e2);
set(p,'facealpha',0.4)
shading flat
axis equal tight
% xlabel('$x$ (km)', 'FontName', 'times', 'Interpreter', 'latex')
% ylabel('$y$ (km)', 'FontName', 'times', 'Interpreter', 'latex')
set(gca, 'fontname', 'times','fontsize',12)
colormap(cmap)
c = colorbar('EastOutside');
% c.Label.String = 'Time-averaged SSH (cm)';
c.Label.FontSize=16;
c.Label.FontName='times';
clim([-2.5, 12])


% Create scatter plot directly using SSH values for colors
s = scatter(alongtrack.x / 1e3, alongtrack.y / 1e3, 7, alongtrack.ssh*1e2, 'filled', ...
    'markerEdgeColor', 'none');
set(gca, 'fontname', 'times', 'fontsize', 12)
colormap(cmap)
c = colorbar('EastOutside');
% c.Label.String = 'SSH (cm)';
c.Label.FontSize = 16;
c.Label.FontName = 'times';

% Set color limits to match actual SSH range
caxis([min(alongtrack.ssh(:)*1e2), max(alongtrack.ssh(:)*1e2)])
% xlabel('$x$ (km)', 'FontName', 'times', 'Interpreter', 'latex')
% ylabel('$y$ (km)', 'FontName', 'times', 'Interpreter', 'latex')
% zlabel('SSH (cm)', 'FontName', 'times', 'Interpreter', 'latex')
axis equal tight
 xlim([-1000,1000]);ylim([-500,500])
clim([-2.5, 12])
set(gca,'xticklabel',[],'xlabel',[])

subplot(3,1,3)
alongtrack=currentMission;
time_window=n-1+[1:10];
t0=min(alongtrack.t);
timeo = [floor(t0):max(floor(alongtrack.t))+1];
if ~isempty(time_window)
    valid_time=logical(sum(floor(alongtrack.t)==timeo(time_window),2));
    alongtrack.t=alongtrack.t(valid_time);
    alongtrack.ssh=alongtrack.ssh(valid_time);
    if isfield(alongtrack, 'lat') && isfield(alongtrack, 'lon')
    alongtrack.lat=alongtrack.lat(valid_time);
    alongtrack.lon=alongtrack.lon(valid_time);
    end
    if isfield(alongtrack, 'x') && isfield(alongtrack, 'y')
    alongtrack.x=alongtrack.x(valid_time);
    alongtrack.y=alongtrack.y(valid_time);
    end
end

% Create scatter plot directly using SSH values for colors
s = scatter(alongtrack.x / 1e3, alongtrack.y / 1e3, 7, alongtrack.ssh*1e2, 'filled', ...
    'markerEdgeColor', 'none');
set(gca, 'fontname', 'times', 'fontsize', 12)
colormap(cmap)
c = colorbar('EastOutside');
% c.Label.String = 'SSH (cm)';
c.Label.FontSize = 16;
c.Label.FontName = 'times';

% Set color limits to match actual SSH range
caxis([min(alongtrack.ssh(:)*1e2), max(alongtrack.ssh(:)*1e2)])
xlabel('$x$ (km)', 'FontName', 'times', 'Interpreter', 'latex')
ylabel('$y$ (km)', 'FontName', 'times', 'Interpreter', 'latex')
% zlabel('SSH (cm)', 'FontName', 'times', 'Interpreter', 'latex')
axis equal tight
 xlim([-1000,1000]);ylim([-500,500])
clim([-2.5, 12])