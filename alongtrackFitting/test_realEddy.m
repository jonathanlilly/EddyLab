%% Alongtrack from OceanDB (array structure)
eddy_id='+527413';
filename = strcat(['G:/My Drive/AlongTrack/MyCode/eddy_',eddy_id,'.nc']);

eddy_time = days(ncread(filename,'eddy/time'))+datetime(1950,01,01);
eddy_lat = ncread(filename,'eddy/latitude');
eddy_lon = ncread(filename,'eddy/longitude');
eddy_amp = ncread(filename,'eddy/amplitude')*100;

at_time = ncread(filename,'alongtrack/time')+datetime(1950,01,01);
at_lat = ncread(filename,'alongtrack/latitude');
at_lon = ncread(filename,'alongtrack/longitude');
at_sla = ncread(filename,'alongtrack/sla_filtered')*100;
at_tracknumber = ncread(filename,'alongtrack/track');

lono_eddy=eddy_lon;%eddy.longitude{1};
lato_eddy=eddy_lat;%eddy.latitude{1};
date_eddy=datenum(eddy_time);%datenum(eddy.date{1});
sla_filtered=at_sla*10;%*1000;

%array of all tracks (no need to repeat)
alongtrack.time=at_time;
alongtrack.latitude=at_lat;
alongtrack.longitude=at_lon;

%% Real Eddy
%alongtrackLonLat from OceanDB
%eddyPath_track from Mason's tracking algorithm
plotAlongtrack(alongtrackLonLat,eddyPath_track);

% Gaussian fit
eddyFit_function = @(x,y,t,A,L,x0,y0,cx,cy) A.*exp(-((x-x0-cx*t).^2 + (y-y0-cy*t).^2)/L^2);
initialParams = [];
params = FitAlongTrackLatLonEddyModel(alongtrackLatLon, eddyFit_function, initialParams);