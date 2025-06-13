function [alongtrack_window, window_indices] = extractAlongtrackWindow(alongtrack,window_start_day,window_end_day)

% ensure that the arrays are in ascending order in time before windowing    
alongtrack.t=floor(alongtrack.t);
[alongtrack.t,sort_idx]=sort(alongtrack.t,'ascend');
alongtrack.x=alongtrack.x(sort_idx);
alongtrack.y=alongtrack.y(sort_idx);
alongtrack.ssh=alongtrack.ssh(sort_idx);

% Find indices that correspond to times within this window
window_indices = find(alongtrack.t >= floor(window_start_day) & alongtrack.t <= floor(window_end_day));

alongtrack_window.x = alongtrack.x(window_indices);
alongtrack_window.y = alongtrack.y(window_indices);
alongtrack_window.t = alongtrack.t(window_indices);
alongtrack_window.ssh = alongtrack.ssh(window_indices);
end
