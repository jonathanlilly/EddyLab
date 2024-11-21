% (xe,ye) are both function handles, that take time as a parameter. So
% although the shape of the eddy is limited to a two parameter gaussian,
% it's propagation path is entirely unspecified.
test_eddy = @(x,y,t,A,L,xe,ye) A.*exp(-((x-xe(t)).^2 + (y-ye(t)).^2)/L^2);

% Let's define a path analytically... but we could use an existing path
% using interp
x0 = 1500e3; y0 = 1500e3; vx = -2.0e-2; vy = -0.3e-2;
x_lin = @(t) x0+vx*t;
y_lin = @(t) y0+vy*t;

% Lets actually make an eddy with a chosen set of parameters.
my_eddy = @(x,y,t) test_eddy(x,y,t,0.15,80e3,x_lin,y_lin);

% Now lets sample the eddy at the appropriate along-track spots.
tracks = load('tracks-lat-30-Lx-2000-Ly-2000.mat');

% the track file only contains one repeat cycle, so we need to repeat it.
% I'm sure there's a better way to do this, but it works.
nRepeats = 5;
obs.t = zeros(nRepeats*length(tracks.x),1);
obs.x = zeros(nRepeats*length(tracks.x),1);
obs.y = zeros(nRepeats*length(tracks.x),1);
for i=1:nRepeats
    obs.t(((i-1)*length(tracks.x)+1):(i*length(tracks.x))) = tracks.t + (i-1)*nRepeats*tracks.repeatTime;
    obs.x(((i-1)*length(tracks.x)+1):(i*length(tracks.x))) = tracks.x;
    obs.y(((i-1)*length(tracks.x)+1):(i*length(tracks.x))) = tracks.y;
end

% Now apply the OSSE!
obs.ssh = my_eddy(obs.x,obs.y,obs.t);

figure, scatter3(obs.x/1e3,obs.y/1e3,obs.ssh*1e2,[],obs.ssh*1e2,'filled'), colorbar('eastoutside')

% Now fit this to a simple model eddy...
model_eddy = @(x,y,t,A,L,x0,y0,cx,cy) A.*exp(-((x-x0-cx*t).^2 + (y-y0-cy*t).^2)/L^2);
penalty_function = @(p) mean((obs.ssh - model_eddy(obs.x,obs.y,obs.t,p(1),p(2),p(3),p(4),p(5),p(6))).^2);

% feed in slightly worse values. These values are terribly scaled, but it
% seems to return exactly the right value anyway.
pmin=fminsearch(penalty_function,[0.13, 85e3, 1450e3, 1450e3, -2e-2, -3e-3]);