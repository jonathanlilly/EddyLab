eddy_id='+527413';
filename = strcat(['E:/My Drive/AlongTrack/MyCode/eddy_',eddy_id,'.nc']);
timeo = ncread(filename,'eddy/time')+datenum(1950,01,01); %Initial time is also a free parameter
eddy_time = timeo-min(timeo);

clearvars eddyPath
eddyPath.lat = ncread(filename,'eddy/latitude');
eddyPath.lon = ncread(filename,'eddy/longitude');

%if you want to change the eddyPath to a function handle
eddyPath_fun_t.lon = @(t) interp1(eddy_time, eddyPath.lon, t, 'linear', 'extrap');
eddyPath_fun_t.lat = @(t) interp1(eddy_time, eddyPath.lat, t, 'linear', 'extrap');
% Function to convert any lat/lon to x/y relative to the eddy center at time t
get_eddyPath_latlon = @(x, y, t) xy2latlon(x, y, eddyPath_fun_t.lat(t),eddyPath_fun_t.lon(t));

%% Extract tracks from available Missions
radius_km=300;
alongtrackLatLon = extractAlongtrackLatLonEddyCenter_AllMissions(AlongTrackSimulator(),eddyPath,timeo,radius=radius_km,getSSH=false);
%%
plotAlongtrack(alongtrackLatLon,eddyPath);
ylim([-40,-20])
xlim([-40,20])
topoplot

%% With Mapped Product
% Set up mapped_field structure
mapped_field.filename = 'D:/UW/MyCode/aviso_madt.nc';
mapped_field.time = ncread(mapped_field.filename, 'time') + datenum(1950, 01, 01);
mapped_field.lat = ncread(mapped_field.filename, 'latitude');
mapped_field.lon = deg180(ncread(mapped_field.filename, 'longitude'));

%% Call the function (no SSH field needed in alongtrack)
plotTrackandAVISO(alongtrackLatLon, eddyPath, mapped_field, radius_km,...
    'snapshot_index', 150, ...
    'lon_bounds', [-40, 30], ...
    'lat_bounds', [-40, -20], ...
    'eddy_id', eddy_id);