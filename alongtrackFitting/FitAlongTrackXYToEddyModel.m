function params = FitAlongTrackXYToEddyModel(alongtrack, eddyFit_fun, initialParams, it_options,options)
arguments
    alongtrack struct
    eddyFit_fun function_handle
    initialParams struct
    it_options struct
    % options struct
    options.LB (1,6) double %= [0, 50, -1000e3, -500e3, -10, -10]
    options.UB (1,6) double %= [30, 150, 1000e3, 500e3, 0, 0]
end
use alongtrack
use options

p0(1)=initialParams.A;
p0(2)=initialParams.L;
p0(3)=initialParams.x0;
p0(4)=initialParams.y0;
p0(5)=initialParams.cx;
p0(6)=initialParams.cy;

% penalty_function = @(p) sum((ssh - eddyFit_fun(x,y,t,p(1),p(2),p(3),p(4),p(5),p(6))).^2);

% Change to scaled parameters
scale_factors = [1e-2, 1e3, 1e3, 1e3, 1e3, 1e3]; % Adjust based on your parameter scales
p0_scaled = p0 ./ scale_factors;

% Scaled penalty function
penalty_function_scaled = @(p_scaled) sum((ssh/1e-2 - eddyFit_fun(x/1e3, y/1e3, t, p_scaled(1)*scale_factors(1), p_scaled(2)*scale_factors(2), p_scaled(3)*scale_factors(3), p_scaled(4)*scale_factors(4), p_scaled(5)*scale_factors(5), p_scaled(6)*scale_factors(6))).^2);

if isfield(options,'LB')
    % pmin=fminsearchbnd(penalty_function, p0, LB, UB, it_options);
    pmin_scaled = fminsearchbnd(penalty_function_scaled, p0_scaled, LB, UB, it_options);

else
    % pmin=fminsearch(penalty_function, p0, it_options);
    pmin_scaled = fminsearch(penalty_function_scaled, p0_scaled, it_options);

end

% Unscale results
pmin = pmin_scaled .* scale_factors;

params.A=pmin(1);
params.L=pmin(2);
params.x0=pmin(3);
params.y0=pmin(4);
params.cx=pmin(5);
params.cy=pmin(6);

th = 0:pi/50:2*pi;
tmat = t-t(1);

% True positions
xo_true = p0(3) + p0(5)*tmat(end);
yo_true = p0(4) + p0(6)*tmat(end);

% Fit position
xo_fit = pmin(3) + pmin(5)*tmat(end);
yo_fit = pmin(4) + pmin(6)*tmat(end);

figure;hold on
plot(p0(2)*sin(th)+xo_true, p0(2)*cos(th)+yo_true,'r--');
plot(pmin(2)*sin(th)+xo_fit, pmin(2)*cos(th)+yo_fit,'b');
plot(xo_true, yo_true,'r*');
plot(xo_fit, yo_fit,'b*');
axis equal
xlim([min(x),max(x)]);ylim([min(y),max(y)])
box on
legend('True position', 'Fit position')
% Calculate error
position_error = sqrt((xo_fit - xo_true).^2 + (yo_fit - yo_true).^2);