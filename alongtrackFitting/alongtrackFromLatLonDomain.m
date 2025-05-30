function alongtrackLatLon = alongtrackFromLatLonDomain(JasonAlongTrack,region,timeo,options)
arguments
    JasonAlongTrack struct
    region (:,:) {mustBeNumeric}
    timeo (:,:) {mustBeNumeric, mustBePositive}
    options.lono (1,1) {mustBeNumeric} = 308
    options.lato (1,1) {mustBeNumeric} = 24
    options.getSSH (1,1) logical = false
end
use options
use JasonAlongTrack

%extract latitude and longitude for this region; t is for track
[latt, lont, timet] = trackextract(lat, lon, time, region);
% [latt, lont, timet] = trackExtractCircle(lat, lon, time, region);
lont = deg180(lont-lono) + lono; %avoids unwrapping issues

%define a start date for the simulation
if numel(timeo) == 1
    % If timeo is a single value, find cycles that are closest to this time
    time_diffs = abs(squeeze(mean(timet, [1, 2], 'omitnan')) - timeo);
    [~, closest_idx] = min(time_diffs);
    a = closest_idx;
    b = closest_idx;
else %If timeo is an array
a = find(squeeze(min(timet, [], [1, 2])) > timeo(1), 1, 'first');
b = find(squeeze(max(timet, [], [1, 2])) < timeo(end), 1, 'last');
b=b+1; %allow next page
end

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
obs.t(round(obs.t) < timeo(1) | round(obs.t) > timeo(end)) = NaN;
alongtrackLatLon.t = obs.t(~isnan(obs.t(:)));
alongtrackLatLon.lon = obs.lon(~isnan(obs.t(:)));
alongtrackLatLon.lat = obs.lat(~isnan(obs.t(:)));

if getSSH
[~, ~, ssh] = trackextract(lat, lon, ssh, region);
obs.ssh = reshape(ssh(:,:,a:b),[],1);
alongtrackLatLon.ssh = obs.ssh(~isnan(obs.t(:)));
end
