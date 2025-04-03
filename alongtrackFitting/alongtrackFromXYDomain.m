function alongtrackLatLon = alongtrackFromXYDomain(JasonAlongTrack,x,y,timeo,options)
% alongtrackFromXYDomain Extracts along-track satellite data for a given domain
%
% This function extracts along-track satellite altimetry data within a specified
% spatial domain defined by x and y coordinates.
%
% Inputs:
%   x - Array of x coordinates defining the domain (meters)
%   y - Array of y coordinates defining the domain (meters)
%   totalDays - Number of days to include in the simulation
%   options - Optional parameters:
%     options.lono - Reference longitude (default: 305)
%     options.lato - Reference latitude (default: 24)
%
% Output:
%   alongtrackLatLon - Struct containing extracted along-track data with fields:
%                     time, lon, lat
%
arguments
    JasonAlongTrack struct
    x (1,:) {mustBeNumeric}
    y (1,:) {mustBeNumeric}
    timeo (:,:) {mustBeNumeric, mustBePositive}
    options.lono (1,1) {mustBeNumeric} = 308
    options.lato (1,1) {mustBeNumeric} = 24
    options.getSSH logical = 'false'
end
use options
% if isempty(varargin)
% lato = 24;
% tn_0 = 84; %track number
% [~, lato_i] = min(abs(lat(:, tn_0)-lato));
% lono = lon(lato_i, tn_0) - 5; %initial longitude is a free parameter
% end
%%
%first, grid x and y
xc = x - mean(x);
yc = y - mean(y);
[xg, yg] = ndgrid(xc, yc);

%covert these to longitude and latitude using a tangent plane
[latg, long] = xy2latlon(xg/1000, yg/1000, lato, lono); %spehrical geometry

%create a region enclosing minima and maxima lon and lat
region = [min(long(:)), max(long(:)), min(latg(:)), max(latg(:))];

alongtrackLatLon = alongtrackFromLatLonDomain(JasonAlongTrack,region,timeo,options);