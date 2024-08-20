% Define the eddy
A = 0.10; % amplitude in m
L = 80e3; % length in m
eddy = @(r) A*exp(-r.*r/L/L);

% Create some noise
sigma = 0.01;
noiseDistribution = StudentTDistribution(sigma,3.0);
noiseDistribution = NormalDistribution(sigma);

% Now create some observations of the eddy with noise
% rng(1)
r = linspace(0,4*L,100).';
ssh = eddy(r) + noiseDistribution.rand(size(r));

% create a spline basis
% M indicates how many splines to create
K = 4; % *order* of the spline (K=4 is a cubic spline)
tKnot = BSpline.knotPointsForDataPoints(r,K=K,M=4);

% Do a simple least-squares fit to the data with these splines
constraints = struct('t',[],'D',[]);
spline = ConstrainedSpline(r,ssh,K,tKnot,noiseDistribution,constraints);

figure
tl = tiledlayout(1,2,TileSpacing="tight");
nexttile

scatter(r/1e3,ssh,'k', 'filled'), hold on, plot(tq/1e3,spline(tq))
xlabel('radial distance (km)')
ylabel('height (m)')
title('Unconstrained fit')

% The constraints let you set the value of the spline, or a derivative, to
% zero. So this example forces the fit to be zero at r_max.
constraints = struct('t',[max(r)],'D',0);
spline_constrained = ConstrainedSpline(r,ssh,K,tKnot,noiseDistribution,constraints);

nexttile
scatter(r/1e3,ssh,'k', 'filled'), hold on, plot(tq/1e3,spline_constrained(tq))
xlabel('radial distance (km)')
title('Constrained fit')