function alongtrackLatLon = extractAlongtrackLatLonEddyCenter(eddyPath,timeo,options)
arguments
    eddyPath struct
    timeo (:,:) {mustBeNumeric, mustBePositive}
    options.readdir string = 'G:\My Drive\AlongTrack\'
end
use options
%% Alongtrack from Jonathan (3D matrix [atd, tracknumber, cycle])
JasonAlongTrack.filename = strcat(readdir, 'JasonAlongTrack.nc');
JasonAlongTrack.lat = ncread(JasonAlongTrack.filename, 'lat');
JasonAlongTrack.lon = ncread(JasonAlongTrack.filename, 'lon');
%JML convert time to Matlab's datenum format
JasonAlongTrack.time = ncread(JasonAlongTrack.filename, 'time') + datenum(1950, 1, 1);
%Time is defined as beginning at 4:05 AM on Sept 23, 1992,
JasonAlongTrack.ssh = ncread(JasonAlongTrack.filename, 'sla');

alongtrackLatLon.t = [];
alongtrackLatLon.lon = [];
alongtrackLatLon.lat = [];
alongtrackLatLon.ssh = [];
for n=1:length(timeo)
    latg=eddyPath.lat(n)-0.75:0.25:eddyPath.lat(n)+0.75;
    long=eddyPath.lon(n)-0.75:0.25:eddyPath.lon(n)+0.75;
    local_at = alongtrackFromLatLonDomain(JasonAlongTrack,timeo(n),lato=eddyPath.lat(n),lono=eddyPath.lon(n),readdir=readdir,getSSH='true');
    
    % Concatenate results

    alongtrackLatLon.t = [alongtrackLatLon.t; local_at.t];
    alongtrackLatLon.lon = [alongtrackLatLon.lon; local_at.lon];
    alongtrackLatLon.lat = [alongtrackLatLon.lat; local_at.lat];
    alongtrackLatLon.ssh = [alongtrackLatLon.ssh; local_at.ssh];
end