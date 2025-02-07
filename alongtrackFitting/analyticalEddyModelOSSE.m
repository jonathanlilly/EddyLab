function alongtrackLatLon = analyticalEddyModelOSSE(tracks,eddy_model)
% tracks should be always in lat lon
arguments (Input)
    tracks {mustBeNumeric,mustBeReal,mustBeFinite}
    eddy_model function_handle 
end

% latc and lonc are the center of the alongtrack domain
lonc=(min(tracks.longitude(:))+max(tracks.longitude(:)))/2;
latc=(min(tracks.latitude(:))+max(tracks.latitude(:)))/2;

% Project {lat,lon} -> {x,y}
[tracks.x, tracks.y] = latlon2xy(tracks.latitude, tracks.longitude, latc, lonc);

% 
nRepeats = 5;
obs.t = zeros(nRepeats*length(tracks.x),1);
obs.x = zeros(nRepeats*length(tracks.x),1);
obs.y = zeros(nRepeats*length(tracks.x),1);
for i=1:nRepeats
    obs.t(((i-1)*length(tracks.x)+1):(i*length(tracks.x))) = tracks.t + (i-1)*nRepeats*tracks.repeatTime;
    obs.x(((i-1)*length(tracks.x)+1):(i*length(tracks.x))) = tracks.x;
    obs.y(((i-1)*length(tracks.x)+1):(i*length(tracks.x))) = tracks.y;
end

% Now apply the OSSE!
obs.ssh = eddy_model(obs.x,obs.y,obs.t);

alongtrackLatLon.longitude=tracks.longitude;
alongtrackLatLon.latitude=tracks.latitude;
alongtrackLatLon.x=tracks.x;
alongtrackLatLon.y=tracks.y;
alongtrackLatLon.ssh=tracks.ssh;

