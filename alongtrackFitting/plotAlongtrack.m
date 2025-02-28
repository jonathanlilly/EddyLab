function plotAlongtrack(alongtrack,eddyPath)
% plotAlongtrack Visualizes along-track data and eddy paths
%   This function creates plots of alongtrack data and eddy paths.
%
% Inputs:
%   alongtrack - Struct containing along-track data with fields:
%                x, y, t, lon, lat, ssh
%   eddyPath_fun_t - Struct containing eddy path function handles:
%                    xe, ye
%   options - Optional parameters
%
arguments
    alongtrack struct
    eddyPath struct
    % options.showLegend (1,1) logical = true
    % options.viewAngle (1,2) double = [-30, 35]
    % options.markerSize (1,1) double = 2
end

% Calculate eddy positions at each alongtrack time directly
% depending on if eddyPath is a function or an array
if isa(eddyPath.xe, 'function_handle')
    eddytrack.x = eddyPath.xe(alongtrack.t-alongtrack.t(1));
    eddytrack.y = eddyPath.ye(alongtrack.t-alongtrack.t(1));
else
    eddytrack.x = eddyPath.xe;
    eddytrack.y = eddyPath.ye;
end


% lato and lono are the center of the alongtrack domain
lono=(min(alongtrack.lon(:))+max(alongtrack.lon(:)))/2;
lato=(min(alongtrack.lat(:))+max(alongtrack.lat(:)))/2;

[eddytrack.latitude, eddytrack.longitude] = xy2latlon(eddytrack.x/1e3,eddytrack.y/1e3,lato,lono);
eddytrack.longitude = deg180(eddytrack.longitude-lono) + lono; %avoids unwrapping issues

%eddytrack
lono_eddy=eddytrack.longitude;
lato_eddy=eddytrack.latitude;

figure;hold on
h(1)=plot(eddytrack.longitude,eddytrack.latitude,'LineWidth',3,'Color',0*[1 1 1]);
plot(eddytrack.longitude(1),eddytrack.latitude(1),'ko','markersize',8,'markerfacecolor','k')
plot(eddytrack.longitude(1),eddytrack.latitude(1),'wo','markersize',5,'markerfacecolor','w')
plot(eddytrack.longitude(end),eddytrack.latitude(end),'wx','markersize',10,'linewidth',2)
plot(eddytrack.longitude(end),eddytrack.latitude(end),'kx','markersize',8,'linewidth',1.5)
h(2)=scatter(alongtrack.lon, alongtrack.lat,2,'MarkerFaceColor','flat');
xlabel('Longitude','FontName', 'times','fontsize',16), ylabel('Latitude','FontName', 'times','fontsize',14), box on
% title(['Altimeter Eddy ',eddy_id],'interpreter','tex')

%set aspect ratio correctly for midpoint of y-axes limits
axis tight
set(gca,'dataaspectratio',[1 cosd(mean(get(gca,'ylim'))) 1])
set(gca, 'fontname', 'times','fontsize',16)
topoplot %continent
latratio(30)
% xlim([min(alongtrack.lon),max(alongtrack.lon)]), ylim([min(alongtrack.lat),max(alongtrack.lat)]);
%for optional legend
if true %set to false is no legend is desired
    clear str
    str{1}=['Eddy-tracking'];
    str{2}=['along-track'];
    % str{3}=['continent'];
    legend(h,str,'location','northeast','interpreter','tex')
end

%% Alongtrack plot with SSH
figure, hold on;
cmap = brewermap(256, '-Spectral');

% Create scatter plot directly using SSH values for colors
s = scatter3(alongtrack.x/1e3, alongtrack.y/1e3, alongtrack.ssh*1e2, 7, alongtrack.ssh*1e2, 'filled', ...
    'markerEdgeColor', 'none');

set(gca, 'fontname', 'times','fontsize',14)
colormap(cmap)
c = colorbar('EastOutside');
c.Label.String = 'SSH (cm)';
c.Label.FontSize=16;
c.Label.FontName='times';
% Set color limits to match actual SSH range
caxis([min(alongtrack.ssh(:)*1e2), max(alongtrack.ssh(:)*1e2)])
xlabel('$x$ (km)', 'FontName', 'times','Interpreter','latex')
ylabel('$y$ (km)', 'FontName', 'times','Interpreter','latex')
zlabel('SSH (cm)', 'FontName', 'times','Interpreter','latex')
daspect([1,1,0.02])
grid on
xlim([min(alongtrack.x/1e3),max(alongtrack.x/1e3)]);ylim([min(alongtrack.y/1e3),max(alongtrack.y/1e3)])
view(-30,35)

%% 2D composite
% 2D statistisc on defined bins
% binsize=12.5;
% max_r=round(max(abs(xEt))/binsize)*binsize;
% [mz, xmid, ymid, numz, stdz] = twodstats(xEt, yEt, ssht, -max_r:binsize:max_r, -max_r:binsize:max_r);
