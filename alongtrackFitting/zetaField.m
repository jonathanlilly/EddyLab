function [zeta] = zetaField(dx,dy,ssh,options)
arguments
    dx (1,1) {mustBeNumeric}
    dy (1,1) {mustBeNumeric}
    ssh (:,:,:) {mustBeNumeric}
    options.lato (1,1) {mustBeNumeric} = 24
end

use options

fo = abs(corfreq(lato))/3600; %rad/s
g=9.80665;% m/s2

for n=1:size(ssh,3)
% Calculate vorticity (zeta=g/f nabla^2 eta) s^-1
[Sx, Sy] = gradient(ssh(:,:,n)', dx, dy);
[Sxx, ~] = gradient(Sx, dx , dy);
[~, Syy] = gradient(Sy, dx , dy);
zeta(:,:,n) = (Sxx + Syy)'*g/fo/fo; % dimensionless zeta: zeta/fo
end