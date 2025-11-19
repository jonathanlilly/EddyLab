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
alongtrackLatLon = extractAlongtrackLatLonEddyCenter_AllMissions(AlongTrackSimulator(),eddyPath,timeo,radius=300,getSSH=false);
%%
plotAlongtrack(alongtrackLatLon,eddyPath);
ylim([-40,-20])
xlim([-40,20])
topoplot
% plotTrackandAVISO(alongtrackLatLon,eddyPath);