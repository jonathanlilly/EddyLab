function params = FitAlongTrackXYToEddyModel(alongtrack, eddyFit_fun, initParams, it_options,options)
arguments
    alongtrack struct
    eddyFit_fun function_handle
    initParams struct
    it_options struct
    % options struct
    options.LB (1,6) double %= [0, 50, -1000e3, -500e3, -10, -10]
    options.UB (1,6) double %= [30, 150, 1000e3, 500e3, 0, 0]
end
use alongtrack
use options

% Add model choice section based on satellites

% Set reference time t0 and calculate elapsed time
t0 = min(t);
elapsed_time = t-t0;

p0(1)=initParams.A;
p0(2)=initParams.L;
p0(3)=initParams.x0;
p0(4)=initParams.y0;
p0(5)=initParams.cx;
p0(6)=initParams.cy;

% penalty_function = @(p) sum((ssh - eddyFit_fun(x,y,t,p(1),p(2),p(3),p(4),p(5),p(6))).^2);

% Change to s1caled parameters (A, L, x0,y0,cx,cy)
scale_factors = [1e-2, 1e3, 1e3, 1e3, 1e3, 1e3]; % Adjust based on your parameter scales
p0_scaled = p0 ./ scale_factors;

% Scaled penalty function
penalty_function_scaled = @(p_scaled) sum((ssh - eddyFit_fun(x, y, elapsed_time, p_scaled(1)*scale_factors(1), p_scaled(2)*scale_factors(2), p_scaled(3)*scale_factors(3), p_scaled(4)*scale_factors(4), p_scaled(5)*scale_factors(5), p_scaled(6)*scale_factors(6))).^2);

if isfield(options,'LB')
    % pmin=fminsearchbnd(penalty_function, p0, LB, UB, it_options);
    pmin_scaled = fminsearchbnd(penalty_function_scaled, p0_scaled, LB, UB, it_options);

else
    % pmin=fminsearch(penalty_function, p0, it_options);
    pmin_scaled = fminsearch(penalty_function_scaled, p0_scaled, it_options);

end

% Unscale pmin
pmin = pmin_scaled .* scale_factors;

params.A=pmin(1);
params.L=pmin(2);
params.x0=pmin(3);
params.y0=pmin(4);
params.cx=pmin(5);
params.cy=pmin(6);
params.t0 = t0;  % Store reference time with results
params.elapsed_time = elapsed_time;