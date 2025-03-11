function [mz_zeta, rmid, r_core, r_shield] = zetaProfile(fullfield,eddyPath,options)
arguments
    fullfield struct %3D matrix in x,y,t
    eddyPath struct
    options.bin_size (1,1) {mustBeNumeric} = 12.5*1e3 %in meter
    options.lato (1,1) {mustBeNumeric} = 24
    options.threshold = 1e-3 %threshold for zero-crossing zeta value in meter 1e-3 m = 1mm
end

use fullfield
use options

fo=abs(corfreq(lato))/3600; %rad/s
g=9.80665;% m/s2

% Calculate eddy positions at each alongtrack time directly
% depending on if eddyPath is a function or an array
clearvars xE yE
for n=1:length(t)
elapsed_time = t(n) - min(t(:));
if isa(eddyPath.xe, 'function_handle')
    xo = eddyPath.xe(elapsed_time);
    yo = eddyPath.ye(elapsed_time);
else
    xo = eddyPath.xe;
    yo = eddyPath.ye;
end

% Calculate eddy-relative coordinates directly
[xE(:,:,n),yE(:,:,n)] = ndgrid(x' - xo, y' - yo);

% Calculate vorticity (zeta=g/f nabla^2 eta) m/s^2 s m/m^2
[Sx, Sy] = gradient(ssh(:,:,n), x' - xo, y' - yo);
[Sxx, ~] = gradient(Sx, x' - xo, y' - yo);
[~, Syy] = gradient(Sy, x' - xo, y' - yo);
zeta(:,:,n) = (Sxx + Syy); % without the g/f
end
% radial statistics on defined bins
max_r=round(max(abs(xE(:)))/bin_size)*bin_size;
tbin = [min(t, [], "all") - .5:1:max(t, [], "all") + 0.5]'; %tbin for radialStatisticsFromScatter has to be per days for 'azimuthal'
rbin = [-bin_size / 2:bin_size:max_r]';

%% Compute statistics
tmat = permute(repmat(t',[1,length(x),length(y)]),[2,3,1]);
[mz_zeta, rmid, ~, ~, ~] = radialStatisticsFromScatter(xE, yE, tmat, zeta*g/fo, rbin, tbin, firstAverage = 'temporal');

%plot profile
%2k : linewidth 2, k black.
%captial A-Z in rainbow order
%T:blue, U:orange, V:yellow, X:green,W:purple
figure; hold on
h = plot(rmid/1e3, mz_zeta);
linestyle 2k 2T 2U-- 2W 2k-- 2V: 2X-. 2W--
hlines(0, 'k:')
xlabel('Radial distance (km)', 'FontName', 'times')
ylabel('\zeta', 'FontName', 'times')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 16)
xlim([0, 250]); %ylim([-2, 30])

%% core and shield radius
[r_core,r_shield]=zeroZetaCrossing(mz_zeta,rmid,threshold=threshold);
vlines(r_core, 'k:')
vlines(r_shield, 'k:')

