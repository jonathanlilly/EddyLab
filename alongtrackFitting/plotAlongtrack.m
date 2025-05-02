function plotAlongtrack(alongtrack, eddyPath,options)
% plotAlongtrack Visualizes along-track data and eddy paths
%   This function creates plots of alongtrack data and eddy paths.
%   Can handle both lat/lon and x/y coordinate systems.
%
% Inputs:
%   alongtrack - Struct containing along-track data with fields:
%                Either (lon, lat, t, ssh) or (x, y, t, ssh)
%   eddyPath - Struct containing eddy path:
%              Either (lon, lat) or (xe, ye) as arrays or function handles
%
arguments
    alongtrack struct
    eddyPath struct
    options.time_window (1,:) {mustBeNumeric} = []
    % options.showLegend (1,1) logical = true
    % options.viewAngle (1,2) double = [-30, 35]
    % options.markerSize (1,1) double = 2
end

t0=min(alongtrack.t);
timeo = [floor(t0):max(floor(alongtrack.t))+1];
if isfield(options, 'time_window')
    valid_time=logical(sum(floor(alongtrack.t)==timeo(options.time_window),2));
    alongtrack.t=alongtrack.t(valid_time);
    alongtrack.x=alongtrack.x(valid_time);
    alongtrack.y=alongtrack.y(valid_time);
    alongtrack.ssh=alongtrack.ssh(valid_time);
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
        lono=(min(alongtrack.lon(:))+max(alongtrack.lon(:)))/2;
        lato=(min(alongtrack.lat(:))+max(alongtrack.lat(:)))/2;
    
        [eddytrack.latitude, eddytrack.longitude] = xy2latlon(eddytrack.x/1e3,eddytrack.y/1e3,lato,lono);
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

% First figure: Track visualization
figure; hold on
h(1) = plot(eddy_x, eddy_y, 'LineWidth', 3, 'Color', 0*[1 1 1]);
plot(eddy_x(1), eddy_y(1), 'ko', 'markersize', 8, 'markerfacecolor', 'k')
plot(eddy_x(1), eddy_y(1), 'wo', 'markersize', 5, 'markerfacecolor', 'w')
plot(eddy_x(end), eddy_y(end), 'wx', 'markersize', 10, 'linewidth', 2)
plot(eddy_x(end), eddy_y(end), 'kx', 'markersize', 8, 'linewidth', 1.5)
h(2) = scatter(plot_x, plot_y, 2, 'MarkerFaceColor', 'flat');
xlabel(xlabel_str, 'FontName', 'times', 'fontsize', 16,'Interpreter','latex')
ylabel(ylabel_str, 'FontName', 'times', 'fontsize', 14,'Interpreter','latex')
box on
axis tight

% Set aspect ratio correctly
if isLatLon
    set(gca, 'dataaspectratio', [1 cosd(mean(get(gca, 'ylim'))) 1])
    topoplot % continent (only for lat/lon plotting)
    latratio(30)
else
    set(gca, 'dataaspectratio', [1 1 1]) % Equal aspect ratio for x/y in km
end

set(gca, 'fontname', 'times', 'fontsize', 16)

% For optional legend
if true % set to false if no legend is desired
    clear str
    str{1} = ['Eddy path'];
    str{2} = ['Along-track'];
    legend(h, str, 'location', 'northeast', 'interpreter', 'tex')
end

% Second figure: Alongtrack plot with SSH
figure, hold on;
cmap = brewermap(256, '-Spectral');

% Create scatter plot directly using SSH values for colors
s = scatter3(plot_x, plot_y, alongtrack.ssh*1e2, 7, alongtrack.ssh*1e2, 'filled', ...
    'markerEdgeColor', 'none');
set(gca, 'fontname', 'times', 'fontsize', 14)
colormap(cmap)
c = colorbar('EastOutside');
c.Label.String = 'SSH (cm)';
c.Label.FontSize = 16;
c.Label.FontName = 'times';

% Set color limits to match actual SSH range
caxis([min(alongtrack.ssh(:)*1e2), max(alongtrack.ssh(:)*1e2)])
xlabel(xlabel_str, 'FontName', 'times', 'Interpreter', 'latex')
ylabel(ylabel_str, 'FontName', 'times', 'Interpreter', 'latex')
zlabel('SSH (cm)', 'FontName', 'times', 'Interpreter', 'latex')
grid on
xlim([min(plot_x), max(plot_x)]);
ylim([min(plot_y), max(plot_y)]);
view(-30, 35)
end
% function plotAlongtrack(alongtrack,eddyPath)
% % plotAlongtrack Visualizes along-track data and eddy paths
% %   This function creates plots of alongtrack data and eddy paths.
% %
% % Inputs:
% %   alongtrack - Struct containing along-track data with fields:
% %                x, y, t, lon, lat, ssh
% %   eddyPath_fun_t - Struct containing eddy path function handles:
% %                    xe, ye
% %   options - Optional parameters
% %
% arguments
%     alongtrack struct
%     eddyPath struct
%     % options.showLegend (1,1) logical = true
%     % options.viewAngle (1,2) double = [-30, 35]
%     % options.markerSize (1,1) double = 2
% end
% 
% % Calculate eddy positions at each alongtrack time directly
% % depending on if eddyPath is a function or an array
% elapsed_time = alongtrack.t - min(alongtrack.t);
% if isfield(eddyPath,'lat')
%     if isa(eddyPath.lat, 'function_handle')
%         eddytrack.longitude=eddyPath.lon(elapsed_time);
%         eddytrack.latitude=eddyPath.lat(elapsed_time);
%     else
%         eddytrack.longitude=eddyPath.lon;
%         eddytrack.latitude=eddyPath.lat;
%     end
% else
%     if isa(eddyPath.xe, 'function_handle')
%         eddytrack.x = eddyPath.xe(elapsed_time);
%         eddytrack.y = eddyPath.ye(elapsed_time);
%     else
%         eddytrack.x = eddyPath.xe;
%         eddytrack.y = eddyPath.ye;
%     end
%     % lato and lono are the center of the alongtrack domain
%     lono=(min(alongtrack.lon(:))+max(alongtrack.lon(:)))/2;
%     lato=(min(alongtrack.lat(:))+max(alongtrack.lat(:)))/2;
% 
%     [eddytrack.latitude, eddytrack.longitude] = xy2latlon(eddytrack.x/1e3,eddytrack.y/1e3,lato,lono);
%     eddytrack.longitude = deg180(eddytrack.longitude-lono) + lono; %avoids unwrapping issues
% end
% 
% 
% figure;hold on
% h(1)=plot(eddytrack.longitude,eddytrack.latitude,'LineWidth',3,'Color',0*[1 1 1]);
% plot(eddytrack.longitude(1),eddytrack.latitude(1),'ko','markersize',8,'markerfacecolor','k')
% plot(eddytrack.longitude(1),eddytrack.latitude(1),'wo','markersize',5,'markerfacecolor','w')
% plot(eddytrack.longitude(end),eddytrack.latitude(end),'wx','markersize',10,'linewidth',2)
% plot(eddytrack.longitude(end),eddytrack.latitude(end),'kx','markersize',8,'linewidth',1.5)
% h(2)=scatter(alongtrack.lon, alongtrack.lat,2,'MarkerFaceColor','flat');
% xlabel('Longitude','FontName', 'times','fontsize',16), ylabel('Latitude','FontName', 'times','fontsize',14), box on
% % title(['Altimeter Eddy ',eddy_id],'interpreter','tex')
% 
% %set aspect ratio correctly for midpoint of y-axes limits
% axis tight
% set(gca,'dataaspectratio',[1 cosd(mean(get(gca,'ylim'))) 1])
% set(gca, 'fontname', 'times','fontsize',16)
% topoplot %continent
% latratio(30)
% % xlim([min(alongtrack.lon),max(alongtrack.lon)]), ylim([min(alongtrack.lat),max(alongtrack.lat)]);
% %for optional legend
% if true %set to false is no legend is desired
%     clear str
%     str{1}=['Eddy path'];
%     str{2}=['Along-track'];
%     % str{3}=['continent'];
%     legend(h,str,'location','northeast','interpreter','tex')
% end
% 
% %% Alongtrack plot with SSH
% figure, hold on;
% cmap = brewermap(256, '-Spectral');
% 
% % Create scatter plot directly using SSH values for colors
% s = scatter3(alongtrack.lon, alongtrack.lat, alongtrack.ssh*1e2, 7, alongtrack.ssh*1e2, 'filled', ...
%     'markerEdgeColor', 'none');
% 
% set(gca, 'fontname', 'times','fontsize',14)
% colormap(cmap)
% c = colorbar('EastOutside');
% c.Label.String = 'SSH (cm)';
% c.Label.FontSize=16;
% c.Label.FontName='times';
% % Set color limits to match actual SSH range
% caxis([min(alongtrack.ssh(:)*1e2), max(alongtrack.ssh(:)*1e2)])
% xlabel('Longitude', 'FontName', 'times','Interpreter','latex')
% ylabel('Latitude', 'FontName', 'times','Interpreter','latex')
% zlabel('SSH (cm)', 'FontName', 'times','Interpreter','latex')
% % daspect([1,1,1.5])
% grid on
% xlim([min(alongtrack.lon),max(alongtrack.lon)]);ylim([min(alongtrack.lat),max(alongtrack.lat)])
% view(-30,35)
