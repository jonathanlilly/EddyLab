function alongtrackLatLon = alongtrackFromLatLonDomain(alongtrack,region,timeo,options)
arguments
    alongtrack struct
    region (:,:) {mustBeNumeric}
    timeo (:,:) {mustBeNumeric, mustBePositive}
    options.lono (1,1) {mustBeNumeric} = 308
    options.lato (1,1) {mustBeNumeric} = 24
    options.readdir string = 'G:\My Drive\AlongTrack\'
    options.getSSH logical = 'false'
end
use options
use alongtrack

%extract latitude and longitude for this region; t is for track
[latt, lont, timet] = trackextract(lat, lon, time, region);
lont = deg180(lont-lono) + lono; %avoids unwrapping issues

%define a start date for the simulation
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
alongtrackLatLon.t = obs.t(~isnan(obs.t(:)));
alongtrackLatLon.lon = obs.lon(~isnan(obs.t(:)));
alongtrackLatLon.lat = obs.lat(~isnan(obs.t(:)));

if getSSH
ssh = ncread(JasonAlongTrack.filename, 'sla');
alongtrackLatLon.ssh = ssh(~isnan(obs.t(:)));
end
