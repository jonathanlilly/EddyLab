%% Compute statistics
[mzxy, rmid, numz, stdAziAvgTemp, AvgAziStdTemp] = radialStatisticsFromScatter(xEt', yEt', timet, ssht, rbin, tbin, firstAverage = 'temporal');
[mzrt, rmid, numz, stdTempAvgAzi, AvgTempStdAzi] = radialStatisticsFromScatter(xEt', yEt', timet, ssht, rbin, tbin, firstAverage = 'azimuthal');
% Std Total
stdTotalxy = sqrt(stdAziAvgTemp.^2+AvgAziStdTemp.^2);
stdTotalrt = sqrt(stdTempAvgAzi.^2+AvgTempStdAzi.^2);

%plot profile
%2k : linewidth 2, k black.
%captial A-Z in rainbow order
%T:blue, U:orange, V:yellow, X:green,W:purple
figure; hold on
h = plot(rmid, [mzxy, AvgAziStdTemp, stdAziAvgTemp, stdTotalxy, mzrt, AvgTempStdAzi, stdTempAvgAzi, stdTotalrt]);
linestyle 2k 2T 2U-- 2W 2k-- 2V: 2X-. 2W--
hlines(0, 'k:')

xlabel('Radial distance (km)', 'FontName', 'times')
ylabel('SSH Variability (cm)', 'FontName', 'times')
set(gca, 'fontname', 'times')
lg = legend(h, '$\overline{\eta_{xy}}^\theta$', '$\overline{\sigma_\eta}^\theta$', '$\varsigma_{\overline{\eta}^t}$', '$\Sigma_{\eta_{xy}}$', '$\overline{\eta_{rt}}^\theta$', '$\overline{\varsigma_\eta}^t$', '$\sigma_{\overline{\eta}^\theta}$', '$\Sigma_{\eta_{rt}}$');
set(lg, 'interpreter', 'latex', 'fontsize', 16, 'orientation', 'vertical', 'NumColumns', 2)
set(gca, 'fontsize', 16)
xlim([0, 250]); ylim([-2, 30])