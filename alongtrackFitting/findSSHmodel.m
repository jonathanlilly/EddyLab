function ssh_values = findSSHmodel(x, y, xmid, ymid, mz, spatial_window)

    % % Get indices for points within spatial window
    % x_valid = (xmid >= min(spatial_window) & xmid <= max(spatial_window));
    % y_valid = (ymid >= min(spatial_window) & ymid <= max(spatial_window));
    % 
    % % Filter the coordinate arrays
    % x_sWindow = xmid(x_valid);
    % y_sWindow = ymid(y_valid);
    % 
    % % Extract the corresponding subset of mz_true
    % ssh_composites = mz(x_valid, y_valid);
    % 
    % Initialize output array
    ssh_values = nan(size(x));
    
    % Process each input point
    for i = 1:numel(x)

        % Find the nearest indices in the filtered arrays
        [~, x_idx] = min(abs(xmid - x(i)));
        [~, y_idx] = min(abs(ymid - y(i)));
        
        % Get the SSH value at this nearest point
        ssh_values(i) = mz(y_idx,x_idx);
    end
    
    % Reshape to match input dimensions
    ssh_values = reshape(ssh_values, size(x));
end