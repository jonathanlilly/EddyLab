function alongtrack = analyticalEddyModelOSSE(alongtrackLatLon,eddy_model)
% tracks should be always in lat lon
arguments (Input)
    alongtrackLatLon struct
    eddy_model function_handle = []
end

% lato and lono are the center of the alongtrack domain
lono=(min(alongtrackLatLon.lon(:))+max(alongtrackLatLon.lon(:)))/2;
lato=(min(alongtrackLatLon.lat(:))+max(alongtrackLatLon.lat(:)))/2;

% Project {lat,lon} -> {x,y}
[x_km, y_km] = latlon2xy(alongtrackLatLon.lat, alongtrackLatLon.lon, lato, lono);

alongtrack.x=x_km*1000;
alongtrack.y=y_km*1000;

% Now apply the OSSE! Shift time to start from t=0;
alongtrack.ssh = eddy_model(alongtrack.x,alongtrack.y,alongtrackLatLon.time-alongtrackLatLon.time(1));

alongtrack.t=alongtrackLatLon.time;
alongtrack.lon=alongtrackLatLon.lon;
alongtrack.lat=alongtrackLatLon.lat;

