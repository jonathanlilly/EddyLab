function [window_start_day, window_end_day, totalTimeWindows] = timeWindowBounds(time, window_size, overlap)
t0=min(time);
totalDays=max(time)-t0+1;

% Ensure time_step is at least 1
time_step = max(1, floor(window_size * (1 - overlap)));

if window_size >= totalDays
    % Special case: window size equals or exceeds total days
    totalTimeWindows = 1;
else
    % Standard case
    totalTimeWindows = ceil((totalDays - window_size) / time_step);
end
for i = 1:totalTimeWindows    
window_start_day(i) = floor(t0 + (i-1)*time_step);
window_end_day(i) = floor(window_start_day(i) + window_size - 1);

if window_end_day(i) > max(time)
    window_end_day(i) = floor(max(time));
    window_start_day(i) = floor(max(t0, window_end_day(i) - window_size + 1));
end
end

end