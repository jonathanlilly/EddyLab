function alongtrackLatLon = alongtrackFromXYDomain(x,y,totalDays,options)
%Input:
% x,y arrays in meters
% totalDays constant
%Output:
%alongtrackLatLon is a struct of track information containing lat,lon,time
% arguments
%     options.lono {mustBeScalarOrEmpty}
%     options.lato {mustBeScalarOrEmpty}
% end
%% Alongtrack from Jonathan (3D matrix [atd, tracknumber, cycle])
readdir = 'G:\My Drive\AlongTrack\';
writedir = 'G:\My Drive\AlongTrack\MyCode\';
JasonAlongTrack.filename = [readdir, 'JasonAlongTrack.nc'];
lat = ncread(JasonAlongTrack.filename, 'lat');
lon = ncread(JasonAlongTrack.filename, 'lon');
atd = ncread(JasonAlongTrack.filename, 'atd');
%JML convert time to Matlab's datenum format
time = ncread(JasonAlongTrack.filename, 'time') + datenum(1950, 1, 1);
%Time is defined as beginning at 4:05 AM on Sept 23, 1992,

if isempty(options.lono)
lato = 24;
tn_0 = 84; %track number
[~, lato_i] = min(abs(lat(:, tn_0)-lato));
lono = lon(lato_i, tn_0) - 5; %initial longitude is a free parameter
end
%%
%first, grid x and y
xc = x - mean(x);
yc = y - mean(y);
[xg, yg] = ndgrid(xc, yc);

%covert these to longitude and latitude using a tangent plane
[latg, long] = xy2latlon(xg/1000, yg/1000, lato, lono); %spehrical geometry

%create a region enclosing minima and maxima lon and lat
region = [min(long(:)), max(long(:)), min(latg(:)), max(latg(:))];

%extract latitude and longitude for this region; t is for track
[latt, lont, timet] = trackextract(lat, lon, time, region);
lont = deg180(lont-lono) + lono; %avoids unwrapping issues

%define a start date for the simulation
%timeo=datenum(1992,9,25)+(0:totalDays-1)';
timeo = datenum(2000, 1, 1) + (0:totalDays - 1)' + 90; %Initial time is also a free parameter
a = find(squeeze(min(timet, [], [1, 2])) > timeo(1), 1, 'first');
b = find(squeeze(max(timet, [], [1, 2])) < timeo(end), 1, 'last');

nRepeats = length(a:b);
obs.t = zeros(nRepeats*length(lont(:)),1);
obs.lon = zeros(nRepeats*length(lont(:)),1);
obs.lat = zeros(nRepeats*length(lont(:)),1);
for i=1:nRepeats
    obs.lon(((i-1)*length(lont(:))+1):(i*length(lont(:)))) = lont(:);
    obs.lat(((i-1)*length(lont(:))+1):(i*length(lont(:)))) = latt(:);
end
obs.t = reshape(timet(:,:,a:b),[],1);
%extract time window
obs.t(obs.t < timeo(1) | obs.t > timeo(end)) = NaN;
alongtrackLatLon.time = obs.t(~isnan(obs.t(:)));
alongtrackLatLon.lon = obs.lon(~isnan(obs.t(:)));
alongtrackLatLon.lat = obs.lat(~isnan(obs.t(:)));

