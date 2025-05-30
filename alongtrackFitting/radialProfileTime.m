function [mz_rt, rmid, tmid, numz_rt, stdz_rt] = radialProfileTime(alongtrack,eddyPath,options)
arguments
    alongtrack struct
    eddyPath struct
    options.rbin_size (1,1) {mustBeNumeric} = 12.5*1e3 %in meter
    options.tbin_size (1,1) {mustBeNumeric} = 9.92 %in day
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
    [ssh_interp, XGrid, YGrid] = interpEddyCentric(x, y, t, xo, yo, ssh, rbin_size, max_r);
    
    % Create 3D coordinate matrices varying in time (following your example)
    xE = repmat(XGrid, [1, 1, length(t)]);
    yE = repmat(YGrid, [1, 1, length(t)]);
    ssh_data = ssh_interp;
    tmat = repmat(permute(elapsed_time, [3, 2, 1]), [size(XGrid, 1), size(XGrid, 2)]);
    
else
    % Alongtrack array data processing
    
    % Calculate eddy-relative coordinates directly
    xE = x - xo;
    yE = y - yo;
    ssh_data = ssh;
    tmat = elapsed_time;
end

% radial statistics on defined bins

max_r = round(max(sqrt(xE(:).^2 + yE(:).^2))/rbin_size)*rbin_size;
tbin = [0.5+min(elapsed_time):tbin_size:max(elapsed_time)+0.5];
rbin = [-rbin_size / 2:rbin_size:max_r]'; %(0:rbin_size:max_r)
%% radial profile vs time
[mz_rt, rmid, tmid, numz_rt, stdz_rt] = twodstats(sqrt(xE.^2+yE.^2), tmat, ssh_data, rbin , tbin);

%% Plots
if options.showplot %default is showplot=true
%% mean radial profile vs time
figure;
hold on; 
jpcolor(rmid/1e3, [1:length(tmid)], mz_rt*1e2), xlim([0, 250]) %mz_rt'./vmean(mz_rt,2)'
xlabel('Radial distance (km)', 'FontName', 'times')
ylabel('Time (Cycle)', 'FontName', 'times')
colormap(brewermap([], '-Spectral'))
c = colorbar('EastOutside');
c.Label.Interpreter='latex';
c.Label.String = '$\overline{\mathrm{SSH}}^{\theta} (\mathrm{cm})$';
c.Label.FontSize=16;
c.Label.FontName='times';

set(gca, 'fontname', 'times','FontSize',16)
% clim([-2.5, 25])

%% std profile over time
figure;
hold on; 
jpcolor(rmid/1e3, [1:length(tmid)], abs(stdz_rt)*1e2), xlim([0, 250]) %mz_rt'./vmean(mz_rt,2)'
xlabel('Radial distance (km)', 'FontName', 'times')
ylabel('Time (Cycle)', 'FontName', 'times')
colormap(brewermap([], '-Spectral'))
c = colorbar('EastOutside');
c.Label.String = 'SSH Variance (cm)';
c.Label.FontSize=16;
c.Label.FontName='times';
set(gca, 'fontname', 'times','FontSize',16)

%% histogram profile over time
figure;
hold on; 
jpcolor(rmid/1e3, [1:length(tmid)], numz_rt), xlim([0, 250]) %mz_rt'./vmean(mz_rt,2)'
xlabel('Radial distance (km)', 'FontName', 'times')
ylabel('Time (Cycle)', 'FontName', 'times')
colormap(brewermap([], '-Spectral'))
c = colorbar('EastOutside');
c.Label.String = 'Histogram (count)';
c.Label.FontSize=16;
c.Label.FontName='times';
set(gca, 'fontname', 'times','FontSize',16)
end
