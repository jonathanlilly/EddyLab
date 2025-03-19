function [mzxy, rmid, numz, stdz] = radialProfile(alongtrack,eddyPath,options)
arguments
    alongtrack struct
    eddyPath struct
    options.bin_size (1,1) {mustBeNumeric} = 12.5*1e3 %in meter
    options.showplot (1,1) logical = true %control whether to display plots
end

use alongtrack
use options

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

% Calculate eddy-relative coordinates directly
xE = x - xo;
yE = y - yo;

% radial statistics on defined bins
max_r=round(max(abs(xE))/bin_size)*bin_size;
tbin = [min(t, [], "all") - .5:1:max(t, [], "all") + 0.5]'; %tbin for radialStatisticsFromScatter has to be per days for 'azimuthal'
rbin = [-bin_size / 2:bin_size:max_r]';

%% Compute statistics
[mzxy, rmid, numz, stdAziAvgTemp, avgAziStdTemp] = radialStatisticsFromScatter(xE, yE, t, ssh, rbin, tbin, firstAverage = 'temporal');
[mzrt, ~, ~, stdTempAvgAzi, avgTempStdAzi] = radialStatisticsFromScatter(xE, yE, t, ssh, rbin, tbin, firstAverage = 'azimuthal');
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
