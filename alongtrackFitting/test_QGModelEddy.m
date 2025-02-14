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