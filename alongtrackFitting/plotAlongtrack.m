function plotAlongtrack(alongtrackLonLat,eddytrack)
%alongtrackLonLat includes Lat, Lon, Time, SSH information from the
%relevant spatial window of  of along-track
%eddytrack contains eddy center information from Mason, 2014 tracking algorithm.

arguments (Input)
    alongtrackLonLat
    eddytrack 
end
% Project {lat,lon} -> {x,y}
[alongtrackXY.x, alongtrackXY.y] = latlon2xy(alongtrackLatLon.lat, alongtrackLatLon.lon, latc, lonc);

%eddytrack
lono_eddy=eddytrack.longitude;
lato_eddy=eddytrack.latitude;
track=ncread(filename,'alongtrack/track');

figure;hold on
h(1)=plot(eddytrack.longitude,eddytrack.latitude,'LineWidth',3,'Color',0*[1 1 1]);
plot(eddytrack.longitude(1),eddytrack.latitude(1),'ko','markersize',8,'markerfacecolor','k')
plot(eddytrack.longitude(1),eddytrack.latitude(1),'wo','markersize',5,'markerfacecolor','w')
plot(eddytrack.longitude(end),eddytrack.latitude(end),'wx','markersize',10,'linewidth',2)
plot(eddytrack.longitude(end),eddytrack.latitude(end),'kx','markersize',8,'linewidth',1.5)
h(2)=scatter(alongtrackLonLat.longitude, alongtrackLonLat.latitude,2,'MarkerFaceColor','flat');
xlabel('Longitude'), ylabel('Latitude'), box on
% title(['Altimeter Eddy ',eddy_id],'interpreter','tex')

%set aspect ratio correctly for midpoint of y-axes limits
set(gca,'dataaspectratio',[1 cosd(mean(get(gca,'ylim'))) 1])
topoplot %continent
latratio(30)
xlim([-40,20]), ylim([-40,-23.5]);
%for optional legend
if true %set to false is no legend is desired
    clear str
    str{1}=['Eddy-tracking'];
    str{2}=['along-track'];
    % str{3}=['continent'];
    legend(h,str,'location','best','interpreter','tex')
end

%% Alongtrack plot with SSH
figure, hold on;
cmap = brewermap(256, '-Spectral');

% Create scatter plot directly using SSH values for colors
s = scatter3(alongtrackXY.x, alongtrackXY.y, alongtrackLatLon.ssh*100, 7, alongtrackLatLon.ssh, 'o', 'filled', ...
    'markerEdgeColor', 'none');

colormap(cmap)
c = colorbar('EastOutside');
ylabel(c, 'SSH (cm)')  % Add label to colorbar
% Set color limits to match actual SSH range
caxis([min(ssht_Gauss(:)), max(ssht_Gauss(:))])
xlabel('x (km)'), ylabel('y (km)'), zlabel('ssh (cm)')
daspect([1,1,0.02])
xlim([min(alongtrackXY.x),max(alongtrackXY.x)]);ylim([min(alongtrackXY.y),max(alongtrackXY.y)])
view(-30,35)