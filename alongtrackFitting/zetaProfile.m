function [mz_zeta, rmid, r_core, r_shield] = zetaProfile(fullfield,eddyPath,options)
arguments
    fullfield struct %x,y,t are arrays ssh has 3D matrix in x,y,t
    eddyPath struct
    options.bin_size (1,1) {mustBeNumeric} = 12.5*1e3 %in meter
    options.lato (1,1) {mustBeNumeric} = 24
    options.epsilon (1,1) {mustBeNumeric} = 0.01 %1e-2/(1e4).^2*g/fo %Tolerance for shield zero-crossing, in case shield plateaus for zero-crossing zeta value.  %dSSH ~ 1e-4 m, QGmodel resolution: dx ~ 4e3 m, dzeta:dSSH/dx*2*g/f/f
end

use fullfield
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

% radial statistics on defined bins
max_r=round(max(abs(xE))/bin_size)*bin_size;
tbin = [min(t, [], "all") - .5:1:max(t, [], "all") + 0.5]'; %tbin for radialStatisticsFromScatter has to be per days for 'azimuthal'
rbin = [-bin_size / 2:bin_size:max_r]';

%% Compute statistics
% tmat = permute(repmat(t',[1,length(x),length(y)]),[2,3,1]);
[mz_zeta, rmid, ~, ~, ~] = radialStatisticsFromScatter(xE, yE, t, zeta, rbin, tbin, firstAverage = 'temporal');

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
[r_core,r_shield]=zeroZetaCrossing(mz_zeta,rmid,epsilon=epsilon);
vlines(r_core/1e3, 'k:')
vlines(r_shield/1e3, 'k:')

