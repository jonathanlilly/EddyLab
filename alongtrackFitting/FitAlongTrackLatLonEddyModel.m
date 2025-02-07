function params = FitAlongTrackLatLonEddyModel(alongtrackLatLon, eddy_model, initialParams, options)
arguments (Input)
    alongtrackLatLon struct
    eddy_model function_handle
    initialParams {mustBeNumeric,mustBeReal,mustBeFinite}
    % options.windowLength="full"
end
arguments (Output)
    params {mustBeNumeric,mustBeReal,mustHaveSameSize(initalParams,params)}
end
% latc and lonc are the center of the alongtrack domain
lonc=(min(alongtrackLatLon.longitude(:))+max(alongtrackLatLon.longitude(:)))/2;
latc=(min(alongtrackLatLon.latitude(:))+max(alongtrackLatLon.latitude(:)))/2;

% Project {lat,lon} -> {x,y}
[alongtrackXY.x, alongtrackXY.y] = latlon2xy(alongtrackLatLon.latitude, alongtrackLatLon.longitude, latc, lonc);

% Call eddy model fit in XY
params = FitAlongTrackXYEddyModel(alongtrackXY, eddy_model, initialParams);