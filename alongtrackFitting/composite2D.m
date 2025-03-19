function [mz, xmid, ymid, numz, stdz, varargout] = composite2D(alongtrack,eddyPath,options)
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

% 2D statistics on defined bins
% bin_size=12.5;
max_r=round(max(abs(xE))/bin_size)*bin_size;
[mz, xmid, ymid, numz, stdz] = twodstats(xE, yE, ssh, -max_r:bin_size:max_r, -max_r:bin_size:max_r);

% sometimes there's numical precision artifacts,
% which gives stdz as a very small complex number. Do a correction:
stdz = abs(stdz);

%% Plots
if options.showplot %default is showplot=true
%% mean ssh 2D
figure;
jpcolor(xmid/1e3, ymid/1e3, mz*1e2);
% r=mean(eddy.speed_radius{1});
% th = 0:pi/50:2*pi;
% plot(r * cos(th),r*sin(th))
% legend('','mean radius')
% jpcolor(xEbin, yEbin, AvgsshAccumBin)
shading flat
axis equal
xlabel('Distance East (km)', 'FontName', 'times')
ylabel('Distance North (km)', 'FontName', 'times')
set(gca, 'fontname', 'times','fontsize',16)
if exist('brewermap')
    colormap(brewermap([], '-Spectral'))
end
c = colorbar('EastOutside');
c.Label.String = 'Time-averaged SSH (cm)';
c.Label.FontSize=16;
c.Label.FontName='times';
xlim([-250, 250]), ylim([-250, 250]);
% clim([-2.5, 25])

%% plot 2D temporal std
figure;
jpcolor(xmid/1e3, ymid/1e3, stdz*1e2);
% jpcolor(xEbin, yEbin, stdTempD);
shading flat
% hold on
% contour(xE,yE,stdTemp2D,[0:0.1:2],'k')
axis equal
xlabel('Distance East (km)', 'FontName', 'times')
ylabel('Distance North (km)', 'FontName', 'times')
set(gca, 'fontname', 'times','fontsize',16)
if exist('brewermap')
    colormap(brewermap([], '-Spectral'))
end
c = colorbar('EastOutside');
c.Label.String = 'Temporal variance (cm)';
c.Label.FontSize=16;
c.Label.FontName='times';
xlim([-250, 250]), ylim([-250, 250]); %clim([0, 8])

%% 2D track counts histogram
figure;
numz(numz == 0) = nan;
jpcolor(xmid/1e3, ymid/1e3, numz); %
% clim(round([min(numz(:)), max(numz(:))]))
shading flat
axis equal
xlabel('Distance East (km)', 'FontName', 'times')
ylabel('Distance North (km)', 'FontName', 'times')
set(gca, 'fontname', 'times','FontSize',16)
if exist('brewermap')
    colormap(brewermap([], '-Spectral'))
end
c = colorbar('EastOutside');
c.Label.String = 'Histogram (counts)';
c.Label.FontSize=16;
c.Label.FontName='times';
xlim([-250, 250]), ylim([-250, 250]);
end

