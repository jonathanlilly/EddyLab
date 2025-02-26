function params = FitAlongTrackLatLonEddyModel(alongtrack, eddyFit_fun, initialParams, options)
arguments
    alongtrack struct
    eddyFit_fun function_handle
    initialParams struct
    options.LB (1,6) double = [0, 50, -1000e3, -500e3, -10, -10]
    options.UB (1,6) double = [30, 150, 1000e3, 500e3, 0, 0]
    % options.windowLength (1,1) string = "full"
end
% lato and lono are the center of the alongtrack domain
lono=(min(alongtrackLatLon.lon(:))+max(alongtrackLatLon.lon(:)))/2;
lato=(min(alongtrackLatLon.lat(:))+max(alongtrackLatLon.lat(:)))/2;

% Project {lat,lon} -> {x,y}
[x_km, y_km] = latlon2xy(alongtrackLatLon.lat, alongtrackLatLon.lon, lato, lono);

alongtrack.x=x_km*1000;
alongtrack.y=y_km*1000;

% Call eddy model fit in XY
params = FitAlongTrackXYToEddyModel(alongtrack, eddyFit_fun, initialParams, options);