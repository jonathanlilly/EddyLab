%% Initialize parameters
readdir = 'G:\My Drive\AlongTrack\';
writedir = 'G:\My Drive\AlongTrack\MyCode\';
JasonAlongTrack.filename = [readdir, 'JasonAlongTrack.nc'];
lat = ncread(JasonAlongTrack.filename, 'lat');
lon = ncread(JasonAlongTrack.filename, 'lon');
atd = ncread(JasonAlongTrack.filename, 'atd');
tn_0 = 84; %track number
lato = 24;
[~, lato_i] = min(abs(lat(:, tn_0)-lato));
lono = lon(lato_i, tn_0);
time = ncread(JasonAlongTrack.filename, 'time');

%% Extract only the model region
filename = 'BetaEddyOne.nc';
x = ncread([readdir, filename], 'x') / 1000; %km
y = ncread([readdir, filename], 'y') / 1000;
ssh = permute(squeeze(ncread([readdir, filename], 'ssh'))*100, [2, 1, 3]);
totalDays = size(ssh, 3);
%redefine x and y to be centered on the domain midpoint
xc = x - mean(x);
yc = y - mean(y);

%extract a subset of tracks

%first, grid x and y
[xg, yg] = ndgrid(xc, yc);

%covert these to longitude and latitude using a tangent plane
[latg, long] = xy2latlon(xg, yg, lato, lono); %spehrical geometry
%mycode%small angle approx
% [long,latg] = xy2lonlat(x,y, mean(x), mean(y), lono, lato);
%cylindrical coordinate with axis on the equator
% [latg,long] = TransverseMercatorToLatitudeLongitude(xg,yg,lon0=lono);
% % Normalize longitude to [-180, 180] range
% long = mod(long + 180, 360) - 180;
% latg = mod(latg + 180, 360) - 180;

%create a region enclosing minima and maxima lon and lat
region = [min(long(:)), max(long(:)), min(latg(:)), max(latg(:))];

%extract latitude and longitude for this region; t is for track
[latt, lont, timet] = trackextract(lat, lon, time, region);
lont = deg180(lont-lono) + lono; %avoids unwrapping issues

figure, plot(lont, latt)

%finally, convert these latitudes and longitudes back to x and y
[xt, yt] = latlon2xy(latt, lont, lato, lono);

figure, jpcolor(xc, yc, ssh(:, :, end)), hold on
plot(xt, yt, linewidth = 2), axis tight, latratio(lato)

%% Eddy-centered x,y track matrix
% xo = 0.5 * (y(144) + y(145)); %origin defined based on eddy center
% yo = 0.5 * (y(80) + y(81));
% [center_xy, core_xy, coreE_xy]=findCenterZeroZeta(x,y,zeta2D,ssh);
% [center_xy, ~, ~]=findSSHMaxContour(xc,yc,ssh);
load('E:\Research\myCode\center_SSHMax.mat', 'center_xy')
center_xy = center_xy - [mean(x), mean(y)];
% center_xy=center_xy-center_xy(1,1);
% interp2(lonmMat,latmMat,center_lonlat(1,1),center_lonlat(1,1))

%%
%find eddy center on model lon,lat matrix
% idx_centerm=[find(x==center_xy(1, 1)),find(y==center_xy(1, 2))];
% lonm(idx_centerm(1,1)),latm(idx_centerm(1,2))
%center_lonlat has to be same as lonm, latm
for n = 1:totalDays
    % xEt(:, :, n) = xt - (center_xy(n, 1));
    % yEt(:, :, n) = yt - (center_xy(n, 2));
    lonEt(:, :, n) = lont - (center_lonlat(n, 1) - center_lonlat(1, 1));
    latEt(:, :, n) = latt - (center_lonlat(n, 2) - center_lonlat(1, 2));
    % [xEt(:, :, n),yEt(:, :, n)]=latlon2xy(latEt(:, :, n),lonEt(:, :, n),lato,lono);
    [xEt(:, :, n), yEt(:, :, n)] = lonlat2xy(lonEt(:, :, n), latEt(:, :, n), xo, yo, lono, lato);
end
xEt = xEt - (center_xy(1, 1));
yEt = yEt - (center_xy(1, 2));

%%
timesortMat = [];
sortindMat = [];
%track where there's min nan
nan_counts = sum(isnan(timet(:, :, 1)), 1);
[~, minNanCol] = min(nan_counts);
for n = 1:totalDays
    [timesort, sortind] = sort(timet(minNanCol, :, n), 2);
    timesortMat = [timesortMat; timesort];
    sortindMat = [sortindMat; sortind];
end
trackday = floor(timesortMat-min(timesortMat(:))) + 1;
ssht = NaN([size(latt), size(ssh, 3)]);
for n = 1:totalDays
    [row, col] = find(trackday == n);
    if isempty(row)
        continue
    end
    rr = unique(row);
    colkeep = sortindMat(rr(1), col);
    if size(rr) > 1
        colkeep = [colkeep, sortindMat(rr(2), col)];
    end
    lontii = NaN(size(lont));
    lontii(:, colkeep) = lont(:, colkeep);
    lattii = NaN(size(latt));
    lattii(:, colkeep) = latt(:, colkeep);
    ssht(:, :, n) = interp2(lonmMat, latmMat, ssh(:, :, n), lontii, lattii, 'linear', 0);
end



%% make video
% clearvars outputVideo3 outputVideo4 outputVideo6
% % clf(3)
% outputVideo3 = VideoWriter(strcat('E:\Research\myCode\aziProfileBin.avi'));
% fr_rate = 365 / 60 * 5; %5 days per sec
% outputVideo3.FrameRate = fr_rate;
% open(outputVideo3)
% outputVideo4 = VideoWriter(strcat('E:\Research\myCode\tempVariabilityBin.avi'));
% outputVideo4.FrameRate = fr_rate;
% open(outputVideo4)
% outputVideo6 = VideoWriter(strcat('E:\Research\myCode\AvgsshBin.avi'));
% outputVideo6.FrameRate = fr_rate;
% open(outputVideo6)

%% define bin
dbindeg = 10; %degree
dbin = round(dbindeg/(x(2) - x(1))); %index
dbinX = dbin * (x(2) - x(1)); %km
xEbound = min(abs([min(xEt(:)), max(xEt(:)), 250]));
rbin = [fliplr([0:-dbinX:-xEbound]), dbinX:dbinX:xEbound];
xEbin = rbin;
yEbin = xEbin;

%% Assign ssht to bins
sshBin = zeros(length(xEbin), length(yEbin), totalDays);
count_binAll = zeros(length(xEbin), length(yEbin), totalDays);
sum_bin = zeros(length(xEbin), length(yEbin));
count_bin = zeros(length(xEbin), length(yEbin));
for nn = 1:totalDays
    for i = 1:length(xEbin) - 1
        for j = 1:length(yEbin) - 1
            % Track elements that falls in the bin ranges
            [row, col] = find(xEt(:, :, nn) >= xEbin(i) & xEt(:, :, nn) <= xEbin(i+1) & yEt(:, :, nn) >= yEbin(j) & yEt(:, :, nn) <= yEbin(j+1));
            % if ~isempty(row) || ~isempty(col)
            %      % No data in this bin, continue to the next bin
            %      keyboard;
            %  end
            bin_data = ssht(row, col, nn); % this is where I need to edit that each snapshot takes different xEt,yEt location.
            sum_bin(j, i) = sum(bin_data, 'all', 'omitmissing');
            count_bin(j, i) = sum(~isnan(bin_data), 'all');
        end
    end
    sshBin(:, :, nn) = sum_bin ./ count_bin .* (count_bin ~= 0);
    count_binAll(:, :, nn) = count_bin;
end

%%
AvgsshAccum = zeros(size(xEt));
clearvars AvgAziAll stdAziAll stdTempAziAll stdTotalAll stdTemp2DAll
AvgsshAccumBin = zeros(length(yEbin), length(xEbin), totalDays);
countAll_lastday = sum(count_binAll(:, :, :), 3);
for n = 1:totalDays
    sshAccumBin = sshBin(:, :, 1:n);
    countAccumBin = sum(count_binAll(:, :, 1:n), 3);
    countAccumBin(countAccumBin == 0) = nan;
    % Output how many counts there are to see if it's enough sample numbers
    % and suggest alternative bin space.
    AvgsshAccumBin(:, :, n) = mean(sshAccumBin, 3, 'omitnan');
    [~, AvgAzi, stdAzi, stdTemp2D, stdTempAzi, stdTotal, bins, err] = AvgProf(xEbin, yEbin, sshAccumBin);

    AvgAziAll(:, n) = AvgAzi(1:end-1);
    stdAziAll(:, n) = stdAzi(1:end-1);
    stdTempAziAll(:, n) = stdTempAzi(1:end-1);
    stdTotalAll(:, n) = stdTotal(1:end-1);
    stdTemp2DAll(:, :, n) = stdTemp2D;
    errAll(:, n) = err;
    figure(3)
    xlim([0, 250]); ylim([-2, 14])
    text(125, 14.5, ['Day ', num2str(n)], 'FontSize', 10, 'HorizontalAlignment', 'center')
    set(gcf, 'color', 'white')
    figure(4)
    xlim([-200, 200]), ylim([-200, 200]); clim([0, 2])
    text(0, 215, ['Day ', num2str(n)], 'FontSize', 10, 'HorizontalAlignment', 'center')
    set(gcf, 'color', 'white')
    figure(5)
    jpcolor(xEbin, yEbin, countAccumBin) %log10(mat)
    text(0, 215, ['Day ', num2str(n)], 'FontSize', 10, 'HorizontalAlignment', 'center')
    shading flat
    axis equal
    xlabel('Distance East (km)', 'FontName', 'times')
    ylabel('Distance North (km)', 'FontName', 'times')
    set(gca, 'fontname', 'times')
    colormap(brewermap([], '-Spectral'))
    colorbar('EastOutside');
    xlim([-200, 200]), ylim([-200, 200]);
    clim(round([min(countAll_lastday(:)), max(countAll_lastday(:))]))
    set(gcf, 'color', 'white')
    figure(6)
    jpcolor(xEbin, yEbin, AvgsshAccumBin(:, :, n))
    shading flat
    axis equal
    xlabel('Distance East (km)', 'FontName', 'times')
    ylabel('Distance North (km)', 'FontName', 'times')
    set(gca, 'fontname', 'times')
    colormap(brewermap([], '-Spectral'))
    colorbar('EastOutside');
    xlim([-200, 200]), ylim([-200, 200]);
    clim([-2.5, 12])
    text(0, 215, ['Day ', num2str(n)], 'FontSize', 10, 'HorizontalAlignment', 'center')
    set(gcf, 'color', 'white')

    % F3 = getframe(3);
    % writeVideo(outputVideo3, F3)
    % F4 = getframe(4);
    % writeVideo(outputVideo4, F4)
    % F6 = getframe(6);
    % writeVideo(outputVideo4, F6)
    % clf(3);
    % clf(4);
    % clf(5);
    % clf(6)
end

%% Jonathan's code: radialStatisticsFromScattered
for n = 1:totalDays
    t(:, :, n) = ones(size(xEt, 1), size(xEt, 2)) * n;
end
tbin = [1:totalDays];
dbindeg = 10; %degree
dbin = round(dbindeg/(x(2) - x(1))); %index
dbinX = dbin * (x(2) - x(1)); %km
xEbound = min(abs([min(xEt(:)), max(xEt(:)), 250]));
rbin = [fliplr([0:-dbinX:-xEbound]), dbinX:dbinX:xEbound];
% Compute statistics
[mz, rmid, numz, stdAziAvgTemp, AvgAziStdTemp] = radialStatisticsFromScattered(xEt, yEt, t, ssht, rbin, tbin, firstAverage = 'temporal');
[mz, rmid, numz, stdTempAvgAzi, AvgTempStdAzi] = radialStatisticsFromScattered(xEt, yEt, t, ssht, rbin, tbin, firstAverage = 'azimuthal');
% Std Total
stdTotal1 = sqrt(stdAziAvgTemp.^2+AvgAziStdTemp.^2);
stdTotal2 = sqrt(stdTempAvgAzi.^2+AvgTempStdAzi.^2);

%%

prettyblue = [0, 0.4470, 0.7410];
prettyorange = [0.8500, 0.3250, 0.0980];
prettyyellow = [0.9290, 0.6940, 0.1250];
prettypurple = [0.4940, 0.1840, 0.5560];
prettygreen = [0.4660, 0.6740, 0.1880];
figure(2); hold on
h1 = plot(rmid, mz, 'LineWidth', 2, 'Color', 'k', 'linestyle', '-');
h2 = plot(rmid, AvgAziStdTemp, 'LineWidth', 2, 'Color', prettyblue, 'linestyle', '-');
h3 = plot(rmid, stdAziAvgTemp, 'LineWidth', 2, 'Color', prettyorange, 'linestyle', '--');
h4 = plot(rmid, AvgTempStdAzi, 'LineWidth', 2, 'Color', prettyyellow, 'linestyle', ':');
h5 = plot(rmid, stdTempAvgAzi, 'LineWidth', 2, 'Color', prettygreen, 'linestyle', '-.');
h6 = plot(rmid, stdTotal1, 'LineWidth', 2, 'Color', prettypurple, 'linestyle', '--');
h7 = plot(rmid, stdTotal2, 'LineWidth', 2, 'Color', prettypurple, 'linestyle', '-');

xlabel('Radial distance', 'FontName', 'times')
ylabel('SSH Variability', 'FontName', 'times')
set(gca, 'fontname', 'times')
% % zero horizontal line
% plot([rmid(1), rmid(end-1)], [0, 0], 'k:')
lg = legend([h1, h2, h3, h4, h5, h6, h7], '$\overline{\mathrm{z}}^\theta$', '$\overline{\sigma_\mathrm{z}}^\theta$', '$\varsigma_{\overline{\mathrm{z}}^t}$', '$\overline{\varsigma_\mathrm{z}}^t$', '$\sigma_{\overline{\mathrm{z}}^\theta}$', '$\Sigma_{z}$', '$\Sigma_{z,2}$');
set(lg, 'interpreter', 'latex', 'fontsize', 15, 'orientation', 'horizontal')
set(gca, 'fontsize', 12)