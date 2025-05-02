function [alongtrack_timewindow,eddyPath_window]=timeWindowedAlongTrack(alongtrack,eddyPath_fun_t,time_window)
t0=min(alongtrack.t);
timeo = [floor(t0):max(floor(alongtrack.t))+1];
valid_time=logical(sum(floor(alongtrack.t)==timeo(time_window),2));
alongtrack_timewindow.t=alongtrack.t(valid_time);
alongtrack_timewindow.x=alongtrack.x(valid_time);
alongtrack_timewindow.y=alongtrack.y(valid_time);
alongtrack_timewindow.ssh=alongtrack.ssh(valid_time);

% Calculate eddy positions at each alongtrack time
elapsed_time = alongtrack_timewindow.t - t0;

% Compute the 2D composite for the current time window
eddyPath_window.xe=eddyPath_fun_t.xe(elapsed_time);
eddyPath_window.ye=eddyPath_fun_t.ye(elapsed_time);
end