function params = FitAlongTrackXYToEddyModel(alongtrackXY, eddy_model, initialParams)
obs=alongtrackXY;
model_eddy=eddy_model;

p(1)=initialParams.A;
p(2)=initialParams.L;
p(3)=initialParams.x0;
p(4)=initialParams.y0;
p(5)=initialParams.cx;
p(6)=initialParams.cy;
penalty_function = @(p) mean((obs.ssh - model_eddy(obs.x,obs.y,obs.t,p(1),p(2),p(3),p(4),p(5),p(6))).^2);

% feed in slightly worse values. These values are terribly scaled, but it
% seems to return exactly the right value anyway.
pmin=fminsearch(penalty_function,[0.13, 85e3, 1450e3, 1450e3, -2e-2, -3e-3]);

params.A=pmin(1);
params.L=pmin(2);
params.x0=pmin(3);
params.y0=pmin(4);
params.cx=pmin(5);
params.cy=pmin(6);