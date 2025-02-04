function params = FitAlongTrackLatLonEddyModel(alongtrackLatLon, eddy_model, initialParams)
% latc and lonc are the center of the alongtrack domain
lonc=(min(alongtrackLatLon.lon(:))+max(alongtrackLatLon.lon(:)))/2;
latc=(min(alongtrackLatLon.lat(:))+max(alongtrackLatLon.lat(:)))/2;


% Project {lat,lon} -> {x,y}
[alongtrackXY.x, alongtrackXY.y] = latlon2xy(alongtrackLatLon.lat, alongtrackLatLon.lon, latc, lonc);

% Call eddy model fit in XY
params = FitAlongTrackXYEddyModel(alongtrackXY, eddy_model, initialParams);