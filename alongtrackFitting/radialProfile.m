function [mzxy, rmid, numz, stdz] = radialProfile(alongtrack,eddyPath_fun_t,options)
arguments
    alongtrack struct
    eddyPath_fun_t struct
    options.bin_size (1,1) {mustBeNumeric} = 12.5*1e3 %in meter
end

use alongtrack
use options
% eddyPath is in the increments of a day, so round the alongtrack to day
% Calculate eddy positions at each alongtrack time directly
xo = eddyPath_fun_t.xe(t-t(1));
yo = eddyPath_fun_t.ye(t-t(1));

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

%plot profile
%2k : linewidth 2, k black.
%captial A-Z in rainbow order
%T:blue, U:orange, V:yellow, X:green,W:purple
figure; hold on
h = plot(rmid/1e3, [mzxy, avgAziStdTemp, stdAziAvgTemp, stdTotalxy, mzrt, avgTempStdAzi, stdTempAvgAzi, stdTotalrt]*1e2);
linestyle 2k 2T 2U-- 2W 2k-- 2V: 2X-. 2W--
hlines(0, 'k:')

xlabel('Radial distance (km)', 'FontName', 'times')
ylabel('SSH Variance (cm)', 'FontName', 'times')
set(gca, 'fontname', 'times')
lg = legend(h, '$\overline{\eta_{xy}}^\theta$', '$\overline{\sigma_\eta}^\theta$', '$\varsigma_{\overline{\eta}^t}$', '$\Sigma_{\eta_{xy}}$', '$\overline{\eta_{rt}}^\theta$', '$\overline{\varsigma_\eta}^t$', '$\sigma_{\overline{\eta}^\theta}$', '$\Sigma_{\eta_{rt}}$');
set(lg, 'interpreter', 'latex', 'fontsize', 16, 'orientation', 'vertical', 'NumColumns', 2)
set(gca, 'fontsize', 16)
xlim([0, 250]); %ylim([-2, 30])

stdz.stdAziAvgTemp = stdAziAvgTemp;
stdz.avgAziStdTemp = avgAziStdTemp;
stdz.stdTempAvgAzi = stdTempAvgAzi;
stdz.avgTempStdAzi = avgTempStdAzi;
stdz.stdTotalxy = stdTotalxy;
stdz.stdTotalrt = stdTotalrt;