function [r_core,r_shield]=zeroZetaCrossing(mz_zeta,rmid,options)
arguments
    mz_zeta (:,:) {mustBeNumeric}
    rmid (:,:) {mustBeNumeric}
    options.epsilon (1,1) {mustBeNumeric} = 0.01 %1e-4/(4e3).^2*g/fo/fo %Tolerance for shield zero-crossing, in case shield plateaus for zero-crossing zeta value. dSSH ~ 1e-2 m, dx ~ 1e4 m -> dSSH/dx*2=1e-11
end

%% zero-crossing
% For first crossing (negative to positive with possible near-zero values)

% Find first crossing (negative to positive)
r_core = NaN; % Initialize to NaN in case no crossing is found
for i = 1:length(mz_zeta)-1
    if (mz_zeta(i) < 0 && mz_zeta(i+1) >= 0)
        % Found the first crossing
        % Interpolate to find the x value
        r_core = rmid(i) + (rmid(i+1) - rmid(i)) * (-mz_zeta(i) / (mz_zeta(i+1) - mz_zeta(i)));
        r_core_idx = i;
        break;
    end
end

% Find second crossing (starts after first crossing)
peak_idx = NaN;
r_shield = NaN; % Initialize to NaN in case no crossing is found
if ~isnan(r_core)
    % Find the peak after the first crossing (maybe use findpeaks)
    [~, maxIdx] = max(mz_zeta(r_core_idx+1:end));
    peak_idx = r_core_idx + maxIdx;

    for i = peak_idx+1:length(mz_zeta)-1
        % Check if y gets close enough to zero (either crossing or just getting close)
        if (abs(mz_zeta(i)) <= options.epsilon)|| (mz_zeta(i-1) * mz_zeta(i) <= 0)
            r_shield = rmid(i);
            break;
        end
    end
end

% Display results
fprintf('First zero crossing at r = %.2f km\n', r_core/1e3);
fprintf('Second zero crossing at r = %.2f km\n', r_shield/1e3);

