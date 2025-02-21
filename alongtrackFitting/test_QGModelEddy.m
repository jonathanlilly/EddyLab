%% WV Eddy
%alongtrackLonLat from OceanDB
%eddyPath_track from ssh max closed contour centroid
%eddy_field is an array of eddy field in (x,y,t,ssh)

%% 1. Extract tracks along the eddy path
%find eddy center
[center_xoyo,amplitude,radius,core_xy] = findEddyCentroid(x, y, ssh,'thresholdratio',0.9,'GetBoundary', true);
eddyPath.xe = center_xoyo(:,1);
eddyPath.ye = center_xoyo(:,2);

% define domain
readdir = 'G:\My Drive\AlongTrack\';
filename = 'BetaEddyOne.nc';
% Load Model
x_QG = ncread([readdir, filename], 'x'); %meters
y_QG = ncread([readdir, filename], 'y'); %meters
eddy_field.ssh = squeeze(ncread([readdir, filename], 'ssh')) * 100; %cm;% matrix order in x,y,z
totalDays = size(ssh, 3);
x = x_QG - mean(x_QG);
y = y_QG - mean(y_QG);

eddy_field.x = x;
eddy_field.y = y;
eddy_field.t = [1:totalDays];

% extract track - takes somes time bc it loads the entire track matrix
alongtrackLatLon = alongtrackFromXYDomain(x,y,totalDays); %options: lato=24, lono=305

%% 2. Apply OSSE on your choice of an analytical eddy shape
%Output: alongtrack - contains arrays of lon,lat,x,y,t,ssh of OSSE 
% apply OSSE on an eddy_field from QG model in (x,y,t)
alongtrack = subsampleOSSE(alongtrackLatLon,eddy_field);

%% 3. Plots and Videos of OSSE
%% plot eddy path with alongtrack, and OSSE
plotAlongtrack(alongtrack,eddyPath_fun_t);

%% make video of the propagating eddy
video_name = 'eddy_field_precessing_ellipse';
makePropagatingVideo(x,y,totalDays,eddy_model,video_name)

%% 4. Eddy-centric coordinate system
%eddy center from eddyPath


%% 5. Eddy composites
[mz, xmid, ymid, numz, stdz] = composite2D(binsize);
radialProfile

%% 6. Gaussian fit
eddyFit_fun = @(x,y,t,A,L,x0,y0,cx,cy) A.*exp(-((x-x0-cx*t).^2 + (y-y0-cy*t).^2)/L^2);
initialParams = [];
params = FitAlongTrackLatLonEddyModel(alongtrackLatLon, eddyFit_fun, initialParams);
% paramsCell = FitAlongTrackLatLonToEddyModelWindowed(alongtrackLatLon, eddyFit_fun, initialParams, windowLength);