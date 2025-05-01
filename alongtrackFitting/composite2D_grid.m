function [mz, xmid, ymid, numz, stdz, varargout] = composite2D_grid(eddy_field,eddyPath,options)
arguments
    eddy_field struct
    eddyPath struct
    options.bin_size (1,1) {mustBeNumeric} = 12.5*1e3 %in meter
    options.showplot (1,1) logical = true %control whether to display plots
end

use eddy_field
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

% 2D statistics on defined bins
% bin_size=12.5;
max_r=(400*1e3/bin_size)*bin_size;
xbin=-max_r:bin_size:max_r;
xmid=(xbin+vshift(xbin,1,1))./2;
xmid=xmid(1:end-1);

% Loop over time
for i=1:length(t)

% Calculate eddy-relative coordinates directly
xE = x - xo(i);
yE = y - yo(i);

% Create meshgrid of original xE, yE positions
[XE, YE] = ndgrid(xE, yE);  % meshgrid in (x, y) order
[XGrid, YGrid] = ndgrid(xmid, xmid);
ssh_interp(:,:,i)=interpn(XE, YE, ssh(:,:,i), XGrid, YGrid,'linear',0);

    % % Flatten for interpolation
    % F = scatteredInterpolant(XE(:), YE(:), ssh(:, :, t)(:), 'linear', 'none');
    % 
    % % Interpolate onto bin center grid
    % ssh_interp(:,:,t) = F(Xgrid, Ygrid);
end

[mz, xmid, ymid, numz, stdz] = twodstats(repmat(XGrid,[1 1 length(t)]), repmat(YGrid,[1 1 length(t)]), ssh_interp, -max_r:bin_size:max_r, -max_r:bin_size:max_r);
% sometimes there's numerical precision artifacts,
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

