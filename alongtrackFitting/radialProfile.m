function [mzxy, rmid, numz, stdz] = radialProfile(alongtrack,eddyPath,options)
arguments
    alongtrack struct
    eddyPath struct
    options.rbin_size (1,1) {mustBeNumeric} = 12.5*1e3 %in meter
    options.showplot (1,1) logical = true %control whether to display plots
    options.data_type string = "auto" % "grid", "alongtrack", or "auto"
    options.max_r (1,1) {mustBeNumeric} = 400*1e3 %in meter
end

use alongtrack
use options

% Determine data type if set to auto
if options.data_type == "auto"
    % Automatically detect data type based on ssh structure
    if isvector(ssh) && length(ssh) == length(x) && length(ssh) == length(y) && length(ssh) == length(t)
        % If ssh is a vector with same length as x, y, t, it's alongtrack data
        data_type = "alongtrack";
    elseif ndims(ssh) == 3 && size(ssh,1) == length(x) && size(ssh,2) == length(y) && size(ssh,3) == length(t)
        % If ssh is a 3D array with dimensions matching x, y, t lengths, it's grid data
        data_type = "grid";
    elseif issparse(ssh) || isscalar(size(ssh))
        % If ssh is sparse or has a single dimension, it's alongtrack data
        data_type = "alongtrack";
    else
        warning('Could not automatically determine data type from ssh structure. Defaulting to alongtrack.');
        data_type = "alongtrack";
    end
else
    data_type = options.data_type;
end

% Calculate eddy positions at each alongtrack time directly
% depending on if eddyPath is a function or an array
elapsed_time = t - min(t);
if isa(eddyPath.xe, 'function_handle')
    xo = eddyPath.xe(elapsed_time);
    yo = eddyPath.ye(elapsed_time);
else
    xo = eddyPath.xe;
    yo = eddyPath.ye;
end

if data_type == "grid"
    % Grid data processing
    [ssh_interp, XGrid, YGrid] = interpEddyCentric(x, y, t, xo, yo, ssh, rbin_size,max_r);

 % Use interpolated coordinates and data
    xE = repmat(XGrid, [1, 1, length(t)]);
    yE = repmat(YGrid, [1, 1, length(t)]);
    ssh_data = ssh_interp;
    t_data= permute(repmat(t,[1 size(xE,1) size(xE,2)]),[2,3,1]);
    
else
    % Alongtrack array data processing
    
    % Calculate eddy-relative coordinates directly
    xE = x - xo;
    yE = y - yo;
    ssh_data = ssh;
    t_data= t;
end

% radial statistics on defined bins
max_r=round(max(sqrt(xE(:).^2 + yE(:).^2))/rbin_size)*rbin_size;
tbin = [min(t, [], "all") - .5:1:max(t, [], "all") + 0.5]'; %tbin for radialStatisticsFromScatter has to be per days for 'azimuthal'
rbin = [-rbin_size / 2:rbin_size:max_r]';

%% Compute statistics
[mzxy, rmid, numz, stdAziAvgTemp, avgAziStdTemp] = radialStatisticsFromScatter(xE, yE, t_data, ssh_data, rbin, tbin, firstAverage = 'temporal');
[mzrt, ~, ~, stdTempAvgAzi, avgTempStdAzi] = radialStatisticsFromScatter(xE, yE, t_data, ssh_data, rbin, tbin, firstAverage = 'azimuthal');
% Std Total
stdTotalxy = sqrt(stdAziAvgTemp.^2+avgAziStdTemp.^2);
stdTotalrt = sqrt(stdTempAvgAzi.^2+avgTempStdAzi.^2);

%% Plot a profile
if options.showplot %default is showplot=true
%2k : linewidth 2, k black.
%captial A-Z in rainbow order
%T:blue, U:orange, V:yellow, X:green,W:purple

figure; hold on
h = plot(rmid/1e3, [mzxy, stdTotalxy, avgAziStdTemp, stdAziAvgTemp, avgTempStdAzi, stdTempAvgAzi]*1e2);
linestyle 2k 2W 2T-- 2U-- 2V-. 2X:
hlines(0, 'k:')

xlabel('Radial distance (km)', 'FontName', 'times')
ylabel('SSH Variance (cm)', 'FontName', 'times')
set(gca, 'fontname', 'times')
lg = legend(h,'$\overline{\eta_{xy}}^\theta$', '$\Sigma_{\eta_{xy}}$', '$\overline{\sigma_\eta}^\theta$', ...
    '$\varsigma_{\overline{\eta}^t}$', ...
    '$\overline{\varsigma_\eta}^t$', '$\sigma_{\overline{\eta}^\theta}$');
set(lg, 'interpreter', 'latex', 'fontsize', 16, 'orientation', 'vertical', 'NumColumns', 2)
set(gca, 'fontsize', 16)
xlim([0, 250]); %ylim([-2, 30])
end

%save std values in a struct
stdz.stdAziAvgTemp = stdAziAvgTemp;
stdz.avgAziStdTemp = avgAziStdTemp;
stdz.stdTempAvgAzi = stdTempAvgAzi;
stdz.avgTempStdAzi = avgTempStdAzi;
stdz.stdTotalxy = stdTotalxy;
stdz.stdTotalrt = stdTotalrt;
