% Test analytical eddies
%% Eddy Path & extract tracks
% define eddy path: xe(t),ye(t) function handles as a function of time
x0 = 500e3; y0 = 0e3; vx = -2.0e3; vy = -0.5e3;
eddyPath_fun_t.xe = @(t) x0+vx*t;
eddyPath_fun_t.ye = @(t) y0+vy*t;

% define domain
x = linspace(-1e6,1e6,200); %meter
y = linspace(-500e3,500e3,200); %meter
totalDays = 300;

% extract track - takes somes time bc it loads the entire track matrix
alongtrackLatLon = alongtrackFromXYDomain(x,y,totalDays); %options: lato=24, lono=305
%% Test 1 - Gaussian Eddy
% define eddy shape function handle
params.A = 0.15; %meter
params.L = 80e3; %meter

% eddy_model is a function handle with a chosen set of parameters (x,y,t)
eddyShapeString = 'Gaussian';
eddy_model = analyticalEddyModel(eddyPath_fun_t,params,eddyShapeString);

% apply OSSE on an eddy_model function (x,y,t)
alongtrack = analyticalEddyModelOSSE(alongtrackLatLon,eddy_model);

%% Test 2 - Steady Elliptical Eddy:
params.A = 0.15;
L = 80e3;
params.La = 0.4*2*L;
params.Lb = 0.2*2*L;

% eddy_model is a function handle with a chosen set of parameters (x,y,t)
eddyShapeString = 'Ellipse';
eddy_model = analyticalEddyModel(eddyPath_fun_t,params,eddyShapeString);

% apply OSSE on an eddy_model function (x,y,t)
alongtrack = analyticalEddyModelOSSE(alongtrackLatLon,eddy_model);

%% Test 3 - Precessing Elliptical Eddy:
params.A = 0.15;
L = 80e3;
params.La = 0.4*2*L;
params.Lb = 0.2*2*L;
params.thetaDot= pi/totalDays;

% eddy_model is a function handle with a chosen set of parameters (x,y,t)
eddyShapeString = 'Ellipse';
eddy_model = analyticalEddyModel(eddyPath_fun_t,params,eddyShapeString);

% apply OSSE on an eddy_model function (x,y,t)
alongtrack = analyticalEddyModelOSSE(alongtrackLatLon,eddy_model);

%% plot eddy path with alongtrack, and OSSE
plotAlongtrack(alongtrack,eddyPath_fun_t);

%% make video of the propagating eddy
video_name = 'eddy_field_precessing_ellipse';
makePropagatingVideo(x,y,totalDays,eddy_model,name)

%% Gaussian fit
eddyFit_fun = @(x,y,t,A,L,x0,y0,cx,cy) A.*exp(-((x-x0-cx*t).^2 + (y-y0-cy*t).^2)/L^2);
initialParams = [];
params = FitAlongTrackLatLonEddyModel(alongtrackLatLon, eddyFit_fun, initialParams);
