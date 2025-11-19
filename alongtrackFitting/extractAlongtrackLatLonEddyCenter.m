function alongtrackLatLon = extractAlongtrackLatLonEddyCenter(JasonAlongTrack,eddyPath,timeo,options)
arguments
    JasonAlongTrack struct
    eddyPath struct
    timeo (:,:) {mustBeNumeric, mustBePositive}
    options.radius (1,1) {mustBeNumeric} = 300
    options.getSSH = true
end
use options

%%
latc=(max(eddyPath.lat)+min(eddyPath.lat))/2;
lonc=(max(eddyPath.lon)+min(eddyPath.lon))/2;
[lat_rad,lon_rad]=xy2latlon(radius,radius,latc,lonc);
lat_rad=lat_rad-latc;
lon_rad=lon_rad-lonc;

alongtrackLatLon.t = [];
alongtrackLatLon.lon = [];
alongtrackLatLon.lat = [];
alongtrackLatLon.ssh = [];
for n=1:length(timeo)
    latg=eddyPath.lat(n)-round(lat_rad):0.25:eddyPath.lat(n)+round(lat_rad);
    long=eddyPath.lon(n)-round(lon_rad):0.25:eddyPath.lon(n)+round(lon_rad);
    %create a region enclosing minima and maxima lon and lat
    region = [min(long(:)), max(long(:)), min(latg(:)), max(latg(:))];

    local_at = alongtrackFromLatLonDomain(JasonAlongTrack,region,timeo(n),lato=eddyPath.lat(n),lono=eddyPath.lon(n),getSSH=getSSH);
    
    % Concatenate results

    alongtrackLatLon.t = [alongtrackLatLon.t; local_at.t];
    alongtrackLatLon.lon = [alongtrackLatLon.lon; local_at.lon];
    alongtrackLatLon.lat = [alongtrackLatLon.lat; local_at.lat];
    if getSSH
    alongtrackLatLon.ssh = [alongtrackLatLon.ssh; local_at.ssh];
    end
end