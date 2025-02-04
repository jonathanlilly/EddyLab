function paramsCell = FitAlongTrackLatLonToEddyModelWindowed(alongtrackLatLon, eddy_model, initialParams, windowLength)
%alongtrackLatLon: lat,lon,time
% Deduce time window length from alongtrackLatLon
at_time = datenum(alongtrackLatLon.time); %check if this is sorted in an ascending order.
totalDays = max(at_time)-min(at_time);
min_days_per_window = 45;  % Need at least 4.5 cycles for the eddy core coverage
min_total_windows = 3;      % Need at least 3 windows to see evolution
overlap = 0.5; %50% overlap between time_window
window_size = max(min_days_per_window, floor(totalDays/(min_total_windows/overlap)));

time_step=floor(window_size*(1-overlap));
totalTimeWindows=floor((totalDays - window_size) / time_step) + 1;

for i=1:totalTimeWindows
    
time_window=1+(i-1)*time_step:(i-1)*time_step+window_size;
%extract time window
at_time(at_time < time_window(1) | at_time > time_window(end)) = NaN;
at_time = at_time(any(~isnan(at_time)));

% latc and lonc are the center of the alongtrack domain
lonc=(min(alongtrackLatLon.lon(:))+max(alongtrackLatLon.lon(:)))/2;
latc=(min(alongtrackLatLon.lat(:))+max(alongtrackLatLon.lat(:)))/2;


% Project {lat,lon,date,ssh} -> {x,y,t,ssh}
[alongtrackXY.x, alongtrackXY.y] = latlon2xy(alongtrackLatLon.lat, alongtrackLatLon.lon, latc, lonc);


% Call eddy model fit in XY
params = FitAlongTrackXYEddyModel(alongtrackXY, eddy_model, initialParams);

paramsCell{i,1} = params;

end