%% Alongtrack from Jonathan (3D matrix [atd, tracknumber, cycle])
%extract latitude and longitude for this region; t is for track
[latt, lont, timet] = trackextract(lat, lon, time, region);

%define a start date for the simulation
%timeo=datenum(1992,9,25)+(0:totalDays-1)';
timeo = datenum(2000, 1, 1) + (0:totalDays - 1)' + 90; %Initial time is also a free parameter
%extract time window
timet(timet < timeo(1) | timet > timeo(end)) = NaN;
timet = timet(:,:,any(~isnan(timet(:))));

lont = deg180(lont-lono) + lono; %avoids unwrapping issues

alongtrackLatLon.lat=latt;
alongtrackLatLon.lon=lont;
alongtrackLatLon.time=timet;

%% Convert alongtrack structure from 3D to array or vice versa
% maybe I need a function for this. Does Jonathan already have this?

%I think in eddy stats, 3D structure is better..?no it doesn't matter
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

lono_eddy=eddy_lon;%eddy.longitude{1};
lato_eddy=eddy_lat;%eddy.latitude{1};
date_eddy=datenum(eddy_time);%datenum(eddy.date{1});
sla_filtered=at_sla*10;%*1000;
time=at_time;
latitude=at_lat;
longitude=at_lon;
track=ncread(filename,'alongtrack/track');
%% Gaussian Eddy
% eddy shape function, and path(xe(t),ye(t)) to get SSH function (=eddy_model) with
% chosen set of parameters (x,y,t)
eddy_function = @(x,y,t,A,L,xe,ye) A.*exp(-((x-xe(t)).^2 + (y-ye(t)).^2)/L^2);
x0 = 1500e3; y0 = 1500e3; vx = -2.0e-2; vy = -0.3e-2;
eddyPath_function.xe = @(t) x0+vx*t;
eddyPath_function.ye = @(t) y0+vy*t;
params.A = 0.15; %meter
params.L = 80e3; %meter
eddy_model = analyticalEddyModel(eddyPath_function,params,eddy_function);

% apply OSSE on an eddy_model function (x,y,t)
alongtrackLatLon = analyticalEddyModelOSSE(tracks,eddy_model);
eddyPath_track = eddyPath_function(tracks.t);
plotAlongtrack(alongtrackLonLat,eddyPath_track);

% Gaussian fit
eddyFit_function = @(x,y,t,A,L,x0,y0,cx,cy) A.*exp(-((x-x0-cx*t).^2 + (y-y0-cy*t).^2)/L^2);
initialParams = [];
params = FitAlongTrackLatLonEddyModel(alongtrackLatLon, eddyFit_function, initialParams);
%% WV Eddy
%alongtrackLonLat from OceanDB
%eddyPath_track from ssh max closed contour centroid
%eddy_field is an array of eddy field in (x,y,t,ssh)
alongtrackLatLon = subsampleOSSE(tracks,eddy_field);
plotAlongtrack(alongtrackLonLat,eddyPath_track);

% Gaussian fit
eddyFit_function = @(x,y,t,A,L,x0,y0,cx,cy) A.*exp(-((x-x0-cx*t).^2 + (y-y0-cy*t).^2)/L^2);
initialParams = [];
params = FitAlongTrackLatLonEddyModel(alongtrackLatLon, eddyFit_function, initialParams);

%% Real Eddy
%alongtrackLonLat from OceanDB
%eddyPath_track from Mason's tracking algorithm
plotAlongtrack(alongtrackLonLat,eddyPath_track);

% Gaussian fit
eddyFit_function = @(x,y,t,A,L,x0,y0,cx,cy) A.*exp(-((x-x0-cx*t).^2 + (y-y0-cy*t).^2)/L^2);
initialParams = [];
params = FitAlongTrackLatLonEddyModel(alongtrackLatLon, eddyFit_function, initialParams);