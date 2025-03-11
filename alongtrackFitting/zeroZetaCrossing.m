function [r_core,r_shield]=zeroZetaCrossing(mz_zeta,rmid,options)
arguments
    mz_zeta (:,:) {mustBeNumeric}
    rmid (:,:) {mustBeNumeric}
    options.threshold = 1e-3 %threshold for zero-crossing zeta value in meter 1e-3 m = 1mm
end
%% zero-crossing
% For first crossing (negative to positive with possible near-zero values)
epsilon = 1e-6; % Tolerance for shield zero-crossing, in case shield plateaus

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
r_shield = NaN; % Initialize to NaN in case no crossing is found
if ~isnan(r_core)
    for i = r_core_idx+1:length(mz_zeta)-1
        % Check if y gets close enough to zero (either crossing or just getting close)
        if (abs(mz_zeta(i)) <= epsilon || (mz_zeta(i-1) * mz_zeta(i) <= 0))
            r_shield = rmid(i);
            break;
        end
    end
end

% Display results
fprintf('First zero crossing at r = %.2f km\n', r_core/1e3);
fprintf('Second zero crossing at r = %.2f km\n', r_shield/1e3);

