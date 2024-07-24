function [AvgVOIE, AvgAzi, stdAziAvgTemp, stdTemp2D, AvgAziStdTemp, stdTotal, bins, err] = AvgProf(xE, yE, VOIE)

%% Average
AvgVOIE = mean(VOIE, 3, 'omitnan'); %AvgVOIE(isnan(AvgVOIE))=0;
totaldays = size(VOIE, 3);

%% Azimuthal Radial profile
dbin = abs(xE(2)-xE(1));
binmax = xE(end); %250
AvgAzi = aziAvg(AvgVOIE, xE, yE, dbin, binmax);

%% temporal std
% Sum_stdTemp2D=zeros(size(AvgVOIE));
% for n=1:totaldays
%     Sum_stdTemp2D=Sum_stdTemp2D+(VOIE(:,:,n)-AvgVOIE).^2;
% end
% stdTemp2D=sqrt(Sum_stdTemp2D/totaldays);
% AvgAziStdTemp=aziAvg(stdTemp2D,xE,yE,dbin,binmax);

Sum_stdTemp2D = zeros([size(AvgVOIE), totaldays]);
for n = 1:totaldays
    Sum_stdTemp2D(:, :, n) = (VOIE(:, :, n) - AvgVOIE).^2;
end
stdTemp2D = sqrt(mean(Sum_stdTemp2D, 3, 'omitnan'));
%take out variability of zeros since there isn't enough samples and I don't
%want it to contribute to low variability.
stdTemp2D(stdTemp2D == 0) = nan;
AvgAziStdTemp = aziAvg(stdTemp2D, xE, yE, dbin, binmax);

%% azimuthal std
[X, Y] = meshgrid(xE, yE);
r = sqrt(X.^2+Y.^2);
bins = 0:dbin:binmax;%linspace(0, binmax, binnum);
bins = [bins, binmax + (bins(2) - bins(1))];
binnum=length(bins)-1;% binnum = 51;
% Initialize arrays to store radial profile and count of elements in each bin
radial_profile = zeros(1, binnum);
bin_count = zeros(1, binnum);
Sum_stdAziAvgTemp = zeros(1, binnum);
stdAziAvgTemp = zeros(1, binnum);
for i = 1:binnum
    % Identify elements in the current bin
    binIdx = (r >= bins(i)) & (r < bins(i+1));

    % Exclude NaN values from accumulation
    VOIbin = (AvgVOIE(binIdx) - AvgAzi(i)).^2;
    nonNanVOIbin = VOIbin(~isnan(VOIbin));

    % Accumulate values and count in the current bin
    Sum_stdAziAvgTemp(i) = sum(nonNanVOIbin);
    bin_count(i) = length(nonNanVOIbin); %sum(binIdx(~isnan(AvgVOIE)));
end
stdAziAvgTemp(bin_count > 0) = sqrt(Sum_stdAziAvgTemp(bin_count > 0)./bin_count(bin_count > 0));
stdAziAvgTemp(bin_count==0)=nan;
%% Total std (temporal std +azimuthal std)
stdTotal = sqrt(stdAziAvgTemp.^2+AvgAziStdTemp.^2);

%% Temporal variance of the azimuthal average
for n = 1:totaldays
    AvgAzin(:, n) = aziAvg(VOIE(:, :, n), xE, yE, dbin, binmax);
    Sum_stdTempAvgAzi(:, :, n) = (AvgAzin(:, n) - AvgAzi).^2;
end
AvgTempStdAzi = mean(std(AvgAzin), 2, 'omitnan');
stdTempAvgAzi = sqrt(mean(Sum_stdTempAvgAzi, 3, 'omitnan'));

%% Weighted Error
aziMean2D = zeros(size(VOIE, [1, 2]));
for i = 1:binnum
    % Identify elements in the current bin
    bin_indices = (r >= bins(i)) & (r < bins(i+1));
    aziMean2D(bin_indices) = AvgAzi(i);
    % Accumulate values and count in the current bin
end

%%
Sum_VOIE = zeros(size(AvgVOIE));
for n = 1:totaldays
    Sum_VOIE(:, :, n) = (VOIE(:, :, n) - aziMean2D).^2;
end
Err2D = sqrt(mean(Sum_VOIE, 3, 'omitnan'));
Err2D(isnan(Err2D)) = 0;
TotalErr = rms(Err2D,'all');
thErr = rms(stdAziAvgTemp);
tErr = rms(AvgAziStdTemp);
oErr = TotalErr^2 - tErr^2 - thErr^2;
err=[TotalErr.^2; tErr.^2; thErr.^2; oErr];
fprintf('Total error: %.3g\n Temporal error: %.3g\n Azimuthal error: %.3g\n  Model deviation error: %.3g\n', TotalErr, tErr, thErr, oErr);

%% output
bins = bins(1:end-1);
AvgAzi = AvgAzi(1:end);
stdAziAvgTemp = stdAziAvgTemp(1:end);
AvgAziStdTemp = AvgAziStdTemp(1:end);
stdTotal = stdTotal(1:end);

% %% plot profile
% prettyblue = [0, 0.4470, 0.7410];
% prettyorange = [0.8500, 0.3250, 0.0980];
% prettyyellow = [0.9290, 0.6940, 0.1250];
% prettypurple = [0.4940, 0.1840, 0.5560];
% figure(3); hold on
% h1 = plot(bins, AvgAzi, 'LineWidth', 2, 'Color', 'k');
% h2 = plot(bins, AvgAziStdTemp, 'LineWidth', 2, 'Color', prettyblue);
% h3 = plot(bins, stdAziAvgTemp, 'LineWidth', 2, 'Color', prettyorange);
% h4 = plot(bins, stdTotal, 'LineWidth', 2, 'Color', prettypurple, 'linestyle', '--');
% xlabel('Radial distance (km)', 'FontName', 'times')
% ylabel('SSH Variability (cm)', 'FontName', 'times')
% set(gca, 'fontname', 'times')
% % xlim([0,250]);ylim([-2,14])
% % zero horizontal line
% plot([bins(1), bins(end-1)], [0, 0], 'k:')
% % lg=legend([h1 h2 h3 h4],'$\overline{\psi}^\theta$','$\varsigma_{\overline{\psi}^t}$','$\overline{\sigma_\psi}^\theta$','$\Sigma_{\psi}$');
% lg = legend([h1, h2, h3, h4], '$\overline{\eta}^\theta$', '$\overline{\sigma_\eta}^\theta$','$\varsigma_{\overline{\eta}^t}$', '$\Sigma_{\eta}$');
% set(lg, 'interpreter', 'latex', 'fontsize', 15, 'orientation', 'horizontal')
% set(gca, 'fontsize', 12)
% 
% %% Composite contour
% figure(4);
% jpcolor(xE, yE, stdTemp2D); shading flat
% % hold on
% % contour(xE,yE,stdTemp2D,[0:0.1:2],'k')
% axis equal
% xlabel('Distance East (km)', 'FontName', 'times')
% ylabel('Distance North (km)', 'FontName', 'times')
% set(gca, 'fontname', 'times')
% colormap(brewermap([], '-Spectral'))
% colorbar('EastOutside');
% text(220, 220, '$\overline{\sigma_\eta}^t$ (cm)', 'fontname', 'times', 'FontSize', 12, 'interpreter', 'latex')
% clim([0,2])
% xlim([min(xE), max(xE)]); ylim([min(yE), max(yE)])
% xlim([-200,200]);ylim([-200,200])
