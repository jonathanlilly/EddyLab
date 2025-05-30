function alongtrackXY = latlon2xy_centered(alongtrackLatLon)
% latlon2xy_centered Applies an Observing System Simulation Experiment
%
% This function applies an OSSE (Observing System Simulation Experiment)
% by sampling the analytical eddy model at along-track data points.
%
% Inputs:
%   track_data - Struct containing along-track data:
%               lat, lon, time - Coordinates and times of along-track points
%   eddy_model_func - Function handle @(x,y,t) representing the eddy model
%
% Output:
%   alongtrack - Struct containing sampled eddy data with fields:
%               x, y, ssh, t, lon, lat
%
arguments
    alongtrackLatLon struct
end

% lato and lono are the center of the alongtrack domain
lono=(min(alongtrackLatLon.lon(:))+max(alongtrackLatLon.lon(:)))/2;
lato=(min(alongtrackLatLon.lat(:))+max(alongtrackLatLon.lat(:)))/2;

% Project {lat,lon} -> {x,y}
[x_km, y_km] = latlon2xy(alongtrackLatLon.lat, alongtrackLatLon.lon, lato, lono);

% alongtrackXY=alongtrackLatLon;
alongtrackXY.x=x_km*1000;
alongtrackXY.y=y_km*1000;
alongtrackXY.t=alongtrackLatLon.t;



