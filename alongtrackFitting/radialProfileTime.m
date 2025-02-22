function [mz_rt, rmid, tmid, numz_rt, stdz_rt] = radialProfileTime(alongtrack,eddyPath_fun_t,options)
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
tmat = t-t(1);
max_r = round(max(abs(xE))/bin_size)*bin_size;
tbin = [0.5+min(tmat):9.92:max(tmat)+0.5];
rbin = [-bin_size / 2:bin_size:max_r]'; %(0:bin_size:max_r)
%% radial profile vs time
[mz_rt, rmid, tmid, numz_rt, stdz_rt] = twodstats(sqrt(xE.^2+yE.^2), tmat, ssh, rbin , tbin);

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

% std profile over time
figure;
hold on; 
jpcolor(rmid/1e3, [1:length(tmid)], stdz_rt*1e2), xlim([0, 250]) %mz_rt'./vmean(mz_rt,2)'
xlabel('Radial distance (km)', 'FontName', 'times')
ylabel('Time (Cycle)', 'FontName', 'times')
colormap(brewermap([], '-Spectral'))
c = colorbar('EastOutside');
c.Label.String = 'SSH Variance (cm)';
c.Label.FontSize=16;
c.Label.FontName='times';
set(gca, 'fontname', 'times','FontSize',16)

% histogram profile over time
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
