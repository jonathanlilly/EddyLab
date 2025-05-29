function [mz, xmid, ymid, numz, stdz, varargout] = composite2D(alongtrack,eddyPath,options)
arguments
    alongtrack struct
    eddyPath struct
    options.bin_size (1,1) {mustBeNumeric} = 12.5*1e3 %in meter
    options.showplot (1,1) logical = true %control whether to display plots
    options.data_type string = "auto" % "grid", "alongtrack", or "auto"
    options.max_r (1,1) {mustBeNumeric} = 400*1e3 %in meter
    options.interp (1,1) logical = false %interpolate scattered data onto regular grid
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
    [ssh_interp, XGrid, YGrid] = interpEddyCentric(x, y, t, xo, yo, ssh, bin_size, max_r);
    
    % Create 3D coordinate matrices varying in time
    xE = repmat(XGrid, [1, 1, length(t)]);
    yE = repmat(YGrid, [1, 1, length(t)]);
    ssh_data = ssh_interp;
    
else
    % Alongtrack array data processing
    if options.interp
    [ssh_interp, xE, yE] = interpEddyCentric(x, y, t, xo, yo, ssh, bin_size, max_r);
    ssh_data = ssh_interp;
    else
    % Calculate eddy-relative coordinates directly
    xE = x - xo;
    yE = y - yo;
    ssh_data = ssh;
    end
end
% 2D statistics on defined bins
max_r = round(max(abs([xE(:); yE(:)]))/bin_size)*bin_size;
[mz, xmid, ymid, numz, stdz] = twodstats(xE, yE, ssh_data, -max_r:bin_size:max_r, -max_r:bin_size:max_r);
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

