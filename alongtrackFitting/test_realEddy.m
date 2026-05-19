%% Alongtrack from OceanDB (array structure)
%eddyPath_track from eddy tracking algorithm
% no eddy field
clear all

1. Load extracted tracks along the eddy path
eddy_id='+527413';
filename = strcat(['E:/My Drive/AlongTrack/MyCode/eddy_',eddy_id,'.nc']);
timeo = ncread(filename,'eddy/time')+datenum(1950,01,01); %Initial time is also a free parameter
eddy_time = timeo-min(timeo);

clearvars eddyPath
eddyPath.lat = ncread(filename,'eddy/latitude');
eddyPath.lon = ncread(filename,'eddy/longitude');

%if you want to change the eddyPath to a function handle
eddyPath_fun_t.lon = @(t) interp1(eddy_time, eddyPath.lon, t, 'linear', 'extrap');
eddyPath_fun_t.lat = @(t) interp1(eddy_time, eddyPath.lat, t, 'linear', 'extrap');
% Function to convert any lat/lon to x/y relative to the eddy center at time t
get_eddyPath_xy = @(lat, lon, t) latlon2xy(lat, lon, eddyPath_fun_t.lat(t),eddyPath_fun_t.lon(t));


% %already extracted alongtrack
% at_time = ncread(filename,'alongtrack/time')+datenum(1950,01,01);
% at_lat = ncread(filename,'alongtrack/latitude');
% at_lon = ncread(filename,'alongtrack/longitude');
% at_sla = ncread(filename,'alongtrack/sla_filtered');
% at_tracknumber = ncread(filename,'alongtrack/track');
% 
% %array of all tracks (no need to repeat)
% alongtrackLatLon.t=at_time;
% alongtrackLatLon.lat=at_lat;
% alongtrackLatLon.lon=at_lon;
% alongtrackLatLon.ssh=at_sla;


% % other option to get it from Jonathan's nc 3D file.
% readdir = 'E:\My Drive\AlongTrack\';
% %% Alongtrack from Jonathan (3D matrix [atd, tracknumber, cycle])
% JasonAlongTrack.filename = strcat(readdir, 'JasonAlongTrack.nc');
% JasonAlongTrack.lat = ncread(JasonAlongTrack.filename, 'lat');
% JasonAlongTrack.lon = ncread(JasonAlongTrack.filename, 'lon');
% %JML convert time to Matlab's datenum format
% JasonAlongTrack.time = ncread(JasonAlongTrack.filename, 'time') + datenum(1950, 1, 1);
% %Time is defined as beginning at 4:05 AM on Sept 23, 1992,
% JasonAlongTrack.ssh = ncread(JasonAlongTrack.filename, 'sla');
% alongtrackLatLon = extractAlongtrackLatLonEddyCenter(JasonAlongTrack,eddyPath,timeo,radius=300); %radius in km

% % if loading jasonalongtrack above takes too long, load the mat file 
current_file_path = matlab.desktop.editor.getActiveFilename;
[current_dir, ~, ~] = fileparts(current_file_path);
load(strcat(current_dir,'\alongtrackLatLon_real.mat'))


Alongtrack in x,y projection about the eddy center
[x_km, y_km]=get_eddyPath_xy(alongtrackLatLon.lat,alongtrackLatLon.lon,alongtrackLatLon.t-min(alongtrackLatLon.t));
alongtrackXY.x=x_km*1e3;
alongtrackXY.y=y_km*1e3;
%if you want to use eddyPath_fun_t.xe,ye as a function handle
%alongtrack are already adjusted to be eddy centered so xo,yo are zeros
eddyPath_fun_t.xe=@(t) 0*t;
eddyPath_fun_t.ye=@(t) 0*t;

3. Plots and Videos of OSSE
plotAlongtrack(alongtrackLatLon,eddyPath);
% clim([-6,25])

%% make video of the propagating eddy
video_name = 'eddy_field_real';
makePropagatingVideo(x,y,eddy_time,eddy_field,video_name)


4. Eddy composites from along-track
%eddy center from eddyPath_track, which is from ssh max closed contour centroid
alongtrack.t=alongtrackLatLon.t;
alongtrack.x=alongtrackXY.x;
alongtrack.y=alongtrackXY.y;
alongtrack.ssh=alongtrackLatLon.ssh;

%If you want to narrow the time window
time_window=200:500;%1:length(timeo)-1;
window_idx = find(alongtrack.t-min(alongtrack.t) >= time_window(1) & alongtrackLatLon.t-min(alongtrack.t) <= time_window(end));
alongtrack_window.t=alongtrack.t(window_idx);
alongtrack_window.x=alongtrack.x(window_idx);
alongtrack_window.y=alongtrack.y(window_idx);
alongtrack_window.ssh=alongtrack.ssh(window_idx);
% time-averaged eddy composite
[mz_xy, xmid_xy, ymid_xy, numz_xy, stdz_xy] = composite2D(alongtrack_window,eddyPath_fun_t);% options: bin_size=12.5*1e3
% time-averaged radial profile
[mz_r, rmid_r, numz_r, stdz_r] = radialProfile(alongtrack_window,eddyPath_fun_t);
% radialProfile over time
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(alongtrack,eddyPath_fun_t);

4. Eddy composites from mapped product
filename = 'D:/UW/MyCode/aviso_madt.nc';
% filename = 'eddy_-41.nc';
ni=ncinfo(filename);
mapped_field.time = ncread(filename,'time')+datenum(1950,01,01);
mapped_field.lat = ncread(filename,'latitude'); %ascending
mapped_field.lon = deg180(ncread(filename,'longitude')); % ascending


Plot mapped product first
%% Call the function (no SSH field needed in alongtrack)
plotTrackandAVISO(alongtrackLatLon, eddyPath, mapped_field, radius_km,...
    'snapshot_index', 150, ...
    'lon_bounds', [-40, 30], ...
    'lat_bounds', [-40, -20], ...
    'eddy_id', eddy_id);

% Set up eddy composite bins
radius_km=300;%km
bin_size=12.5*1e3;
max_r=(radius_km*1e3/bin_size)*bin_size;
xbin=-max_r:bin_size:max_r;
xmid=(xbin+vshift(xbin,1,1))./2;
xmid=xmid(1:end-1);
ymid=xmid;
[XGrid, YGrid] = ndgrid(xmid, ymid);
clearvars ssh_interp

for n = 1:length(timeo)
    % Find time index closest to current eddy time
    [~, time_idx] = min(abs(mapped_field.time - timeo(n)));

    % Find indices for region around current eddy position
    [~, center_lat] = min(abs(mapped_field.lat - eddyPath.lat(n)));
    [~, center_lon] = min(abs(mapped_field.lon - eddyPath.lon(n)));

    % Calculate the lat-lon box that corresponds to our desired km radius
    [lat_box, lon_box] = xy2latlon([-radius_km, radius_km], [-radius_km, radius_km], eddyPath.lat(n), eddyPath.lon(n));
     
    % Find indices in the grid that match these boundaries
    [~, start_lat] = min(abs(mapped_field.lat - lat_box(1)));
    [~, last_lat] = min(abs(mapped_field.lat - lat_box(2)));
    [~, start_lon] = min(abs(mapped_field.lon - lon_box(1)));
    [~, last_lon] = min(abs(mapped_field.lon - lon_box(2)));
    
    % Handle longitude wrap-around if needed
    if start_lon > last_lon % Crossing the longitude boundary
        % First part: end of longitude array to the end
        start1 = [start_lon, start_lat, time_idx];
        count1 = [length(mapped_field.lon) - start_lon + 1, last_lat - start_lat + 1, 1];
        ssh_part1 = ncread(filename, 'sla', start1, count1); %m
        lon_part1 = mapped_field.lon(start_lon:end);
        
        % Second part: beginning of longitude array to last_lon
        start2 = [1, start_lat, time_idx];
        count2 = [last_lon, last_lat - start_lat + 1, 1];
        ssh_part2 = ncread(filename, 'sla', start2, count2);%m
        lon_part2 = mapped_field.lon(1:last_lon);
        
        % Combine the two parts
        current_ssh = [ssh_part1; ssh_part2];
        current_lon = [lon_part1; lon_part2];
    else
        % Normal case - extraction window within longitude bounds
        start = [start_lon, start_lat, time_idx];
        count = [last_lon - start_lon + 1, last_lat - start_lat + 1, 1];
        current_ssh = ncread(filename, 'sla', start, count);%m
        current_lon = mapped_field.lon(start_lon:last_lon);
    end
    
    current_lat = mapped_field.lat(start_lat:last_lat); %lat doesn't have a if statement bc we don't do eddies that cross the equator
    % Create meshgrid for extracted region
    [lon_grid, lat_grid] = ndgrid(current_lon, current_lat);
    [x_km, y_km]=get_eddyPath_xy(lat_grid,lon_grid,timeo(n)-min(timeo));

    % Interpolate onto a grid
    ssh_interp(:,:,n) = griddata(double(x_km * 1e3), double(y_km * 1e3), current_ssh, XGrid, YGrid, 'linear');
    % ssh_interp(:,:,n)=interpn(x_km * 1e3, y_km * 1e3, current_ssh, XGrid, YGrid,'linear',0);   
end
    % Store mapped grid data instead of an array
    mapped_extracted.ssh = ssh_interp;
    mapped_extracted.x = xmid;
    mapped_extracted.y = ymid;
    mapped_extracted.t = timeo;%permute(repmat(timeo,[1,length(xmid),length(xmid)]),[2,3,1]);

% time_window=1:length(timeo)-1;
mapped_window.t=mapped_extracted.t(time_window);
mapped_window.x=mapped_extracted.x;
mapped_window.y=mapped_extracted.y;
mapped_window.ssh=mapped_extracted.ssh(:,:,time_window);

% time-averaged eddy composite
% [mz_xy, xmid_xy, ymid_xy, numz_xy, stdz_xy] = composite2D_grid(mapped_window,eddyPath_fun_t);
[mz_xy, xmid_xy, ymid_xy, numz_xy, stdz_xy] = composite2D(mapped_window,eddyPath_fun_t);% options: bin_size=12.5*1e3
% time-averaged radial profile
[mz_r, rmid_r, numz_r, stdz_r] = radialProfile(mapped_window,eddyPath_fun_t);

time_window=1:length(timeo)-1;
mapped_window.t=mapped_extracted.t(time_window);
mapped_window.x=mapped_extracted.x;
mapped_window.y=mapped_extracted.y;
mapped_window.ssh=mapped_extracted.ssh(:,:,time_window);

% radialProfile over time
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(mapped_window,eddyPath_fun_t,"tbin_size",10); %mapped product are in the increments of days

5. Parameter extraction from mapped product
%%
% [center_zeta,amplitude,radius,core_zeta,shield_zeta] = findEddyCentroid(x, y, ssh,'BoundaryType','zero-zeta','lato',24,'GetBoundary', true,'epsilon',0.01);
amplitude = ncread(filename,'eddy/amplitude')*100;
radius = ncread(filename,'eddy/speed_radius')/1e3;

% Calculate grid spacing
dx = mapped_extracted.x(2) - mapped_extracted.x(1);
dy = mapped_extracted.y(2) - mapped_extracted.y(1);
% calculate zeta in (x,y,t) matrix
[zeta] = zetaField(dx,dy,mapped_extracted.ssh,lato=mean(eddyPath.lat));

[center_zeta,amplitude,radius,core_zeta,shield_zeta] = findEddyCentroid(mapped_window.x, mapped_window.y, mapped_window.ssh,'BoundaryType','zero-zeta','lato',24,'GetBoundary', true,'epsilon',0.01);

n=25;
zetaContours.core=core_zeta{n};
% zetaContours.shield=shield_zeta{n};
plotComposite(mapped_window.x,mapped_window.y,mapped_window.t(n),zeta(:,:,n),eddyPath_fun_t,contour=zetaContours)
hold on
text(-245,270,['Day ',num2str(n)],'FontSize',14,'FontName','times')

mapped_window.zeta = reshape(zeta,[],1);

[mz_zeta, rmid] = zetaProfile(mapped_window,eddyPath_fun_t);

5. Convergence Rate
% how many cycles are needed to reach convergence to time-averaged profile
[convergence] = convergenceRate(mz_rt, numz_rt);


6. Gaussian fit
    % mapped_field_extracted.ssh = ssh_interp;
    % mapped_field_extracted.x = xmid;
    % mapped_field_extracted.y = xmid;
    % mapped_field_extracted.t = timeo;
eddyFit_fun = @(x,y,t,A,L,x0,y0,cx,cy) A.*exp(-((x-x0-cx*t).^2 + (y-y0-cy*t).^2)/L^2);
[~,amplitude,radius] = findEddyCentroid(mapped_extracted.x, mapped_extracted.y, mapped_extracted.ssh,'thresholdratio',0.9,'GetBoundary', true);
eddyParams.A=amplitude;
eddyParams.L=radius;
eddyParams.xe=eddyPath_fun_t.xe;
eddyParams.ye=eddyPath_fun_t.ye;

trueParams.A=mean(eddyParams.A);
trueParams.L=mean(eddyParams.L);
trueParams.x0=eddyParams.xe(0);
trueParams.y0=eddyParams.ye(0);
trueParams.cx = mean(vdiff(eddyParams.xe(eddy_time),1));
trueParams.cy = mean(vdiff(eddyParams.ye(eddy_time),1));

% new inital parameters for this window
initParams.A = trueParams.A+2*(rand-0.5)*1e-2; %random uncertainty +/- 1e-2
initParams.L = trueParams.L+2*(rand-0.5)*1e3; %random uncertainty +/- 1e3
%assuming you roughly know the eddy center from eddy-tracking algorithm
%beginning of this particular window minus t0 offset of the entire eddy lifetime
% initParams_window.x0 = initParams.x0+i*initParams.cx*max(elapsed_time_window);
% initParams_window.y0 = initParams.y0+i*initParams.cy*max(elapsed_time_window);
initParams.x0 = trueParams.x0+2*(rand-0.5)*10e3; %random uncertainty +/- 10e3
initParams.y0 = trueParams.y0+2*(rand-0.5)*10e3;
initParams.cx = trueParams.cx+2*(rand-0.5)*1e2;  %random uncertainty +/- 1e2
initParams.cy = trueParams.cy+2*(rand-0.5)*1e2;


% initParams.A = true_params.A-0.
% initParams.L = true_params.L+5e3;
% initParams.x0 = true_params.x0-50e3;
% initParams.y0 = true_params.y0+50e3;
% initParams.cx = true_params.cx-0.2e3;
% initParams.cy = true_params.cy+0.2e3;

plot_func = @(xtrans, optimValues, state) plotParamsIteration_ref(xtrans, optimValues, state, [], trueParams);
it_options = optimset('OutputFcn', plot_func,'TolX',1e-3,'TolFun',1e-3);
paramsFit = FitAlongTrackXYToEddyModel(alongtrack, eddyFit_fun, initParams, it_options);

plotFitPosition(eddyPath_fun_t,paramsFit,trueParams)
xlim([min(x_km(:)),max(x_km(:))]);ylim([min(y_km(:)),max(y_km(:))])

%% time windowed fit
it_options = optimset('TolX',1e-3,'TolFun',1e-3);
[paramsCell, trueParamsCell, window_center] = FitAlongTrackXYToEddyModelWindowed(alongtrack, eddyFit_fun, eddyParams, it_options,"window",100);

plotFitPosition(eddyPath_fun_t,paramsCell,trueParamsCell)
xlim([min(x_km(:)),max(x_km(:))]);ylim([min(y_km(:)),max(y_km(:))])
plotFitWindowed(paramsCell,trueParamsCell, eddyPath_fun_t)

% % To analyze a specific window in detail (e.g., window 3):
 plotSingleWindowFit(alongtrack, eddyPath_fun_t, paramsCell, trueParamsCell, 3);



%% 7. Publication figures (D3 + D5.2) ----------------------------------------
% Depends on sections 4 (along-track) and 4b (mapped product) having run.
%
% Struct reference after mlx sync:
%   alongtrack        — full-life along-track (t,x,y,ssh as 1D arrays, t in datenum)
%   alongtrack_window — pre-windowed subset (time_window=200:500 days), used for 2D composites above
%   mapped_extracted  — grid truth: .ssh [nx,ny,nt] in m, .x/.y [1D] in m, .t [nt] datenum
%   mapped_window     — last assignment is full-life (time_window=1:end-1); used for r-t above
%
% NOTE: mz_rt/stdz_rt in workspace = mapped full-life (overwritten from along-track version).
%       Section 7 recomputes both with distinct _at / _mp suffixes to avoid confusion.

% --- 7a. MSE sweep: optimal time window (truth = mapped_extracted grid) ---
window_days_arr = [10,20,50,100,300,1000];               % candidate window lengths in cycles
rmse_all = nan(size(window_days_arr));
for k = 1:numel(window_days_arr)
    % compute_model_error(alongtrack, eddy_field_truth, eddyPath, window_days, model)
    % mapped_extracted is the grid truth (auto-detected as grid by composite2D inside)
    [~, rmse_all(k)] = compute_model_error(alongtrack, mapped_extracted, ...
        eddyPath_fun_t, window_days_arr(k), 'composite', showplot=true);
end
figure;
plot(window_days_arr, rmse_all, 'k-o', 'LineWidth',2, 'MarkerFaceColor','k');
xlabel('Time window (cycles)','FontName','times');
ylabel('RMSE vs mapped product','FontName','times');
title('MSE sweep','FontName','times');
set(gca,'FontName','times','FontSize',14); grid on;
[~, k_opt]   = min(rmse_all);
W_opt_days   = window_days_arr(k_opt);
fprintf('MSE argmin: W* = %.1f days — inspect curve and override if elbow differs\n', ...
    W_opt_days);
% W_opt_cycles = ??;   % <-- uncomment and set manually if elbow != argmin
% W_opt_days   = W_opt_cycles * cycle_days;
%%
% --- 7b. Subset along-track and mapped_extracted to optimal window ---
W_opt_days=100; W_opt_cycles=W_opt_days/10;
% Along-track: subset full-life `alongtrack` by elapsed days
t0_at  = min(alongtrack.t);
in_W_at = (alongtrack.t - t0_at) <= W_opt_days;
alongtrack_W.t   = alongtrack.t(in_W_at);
alongtrack_W.x   = alongtrack.x(in_W_at);
alongtrack_W.y   = alongtrack.y(in_W_at);
alongtrack_W.ssh = alongtrack.ssh(in_W_at);

% Mapped: subset mapped_extracted grid by same elapsed days
t0_mp   = min(mapped_extracted.t);
in_W_mp = (mapped_extracted.t - t0_mp) <= W_opt_days;
mapped_W_grid.x   = mapped_extracted.x;
mapped_W_grid.y   = mapped_extracted.y;
mapped_W_grid.t   = mapped_extracted.t(in_W_mp);
mapped_W_grid.ssh = mapped_extracted.ssh(:,:,in_W_mp);   % 3D grid subset

% --- 7c. 2D composites at optimal window (panels a,b,d,e) ---
[mz_xy_at_W, xmid_pub, ymid_pub, ~, stdz_xy_at_W] = composite2D(alongtrack_W,  eddyPath_fun_t, showplot=false);
[mz_xy_mp_W, ~,         ~,        ~, stdz_xy_mp_W] = composite2D(mapped_W_grid, eddyPath_fun_t, showplot=false);

% --- 7d. Full-life r-t composites (panels c,f) ---
% Recompute with distinct variable names (_at / _mp) to avoid overwrite confusion.
[~, rmid_pub, tmid_at, ~, stdz_rt_at] = radialProfileTime(alongtrack,      eddyPath_fun_t, showplot=false);
[~, ~,        tmid_mp, ~, stdz_rt_mp] = radialProfileTime(mapped_extracted, eddyPath_fun_t, showplot=false, tbin_size=10);
tmid_cyc_at = (tmid_at - min(tmid_at)) / cycle_days;   % days -> cycles
tmid_cyc_mp = (tmid_mp - min(tmid_mp)) / cycle_days;

% --- 7e. Combined 2x3 publication figure (new Fig 5) ---
% CHECK clim values once figure is generated — adjust to actual data range.
clim_ssh = [-5 20];   % cm; adjust to data
clim_std = [0  5];    % cm; adjust to data

fig5 = figure('Position',[100 100 1500 800]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile; imagesc(xmid_pub/1e3, ymid_pub/1e3, mz_xy_at_W');  axis xy image; clim(clim_ssh);
colorbar; xlabel('x (km)'); ylabel('y (km)'); title('(a) along-track \overline\eta^t')

nexttile; imagesc(xmid_pub/1e3, ymid_pub/1e3, stdz_xy_at_W'); axis xy image; clim(clim_std);
colorbar; xlabel('x (km)'); title('(b) along-track \sigma_\eta')

nexttile; imagesc(tmid_cyc_at, rmid_pub/1e3, stdz_rt_at); axis xy;
colorbar; xlabel('Cycle'); ylabel('r (km)'); title('(c) along-track \varsigma_{\eta''} (r,t)')
hold on; xline(W_opt_cycles,'w--','LineWidth',2,'Label',sprintf('W^*=%d',W_opt_cycles)); hold off

nexttile; imagesc(xmid_pub/1e3, ymid_pub/1e3, mz_xy_mp_W');  axis xy image; clim(clim_ssh);
colorbar; xlabel('x (km)'); ylabel('y (km)'); title('(d) mapped \overline\eta^t')

nexttile; imagesc(xmid_pub/1e3, ymid_pub/1e3, stdz_xy_mp_W'); axis xy image; clim(clim_std);
colorbar; xlabel('x (km)'); title('(e) mapped \sigma_\eta')

nexttile; imagesc(tmid_cyc_mp, rmid_pub/1e3, stdz_rt_mp); axis xy;
colorbar; xlabel('Cycle'); ylabel('r (km)'); title('(f) mapped \varsigma_{\eta''} (r,t)')
hold on; xline(W_opt_cycles,'w--','LineWidth',2,'Label',sprintf('W^*=%d',W_opt_cycles)); hold off

sgtitle(sprintf('Agulhas Ring +527413  |  W^* = %d cycles (%.0f days)', W_opt_cycles, W_opt_days), ...
    'FontName','times','FontSize',14)
% exportgraphics(fig5, fullfile(fileparts(mfilename('fullpath')), ...
%     '../../AlongTrack-GRL/fig/real_eddy_composite.png'), 'Resolution',300);
% fprintf('Saved real_eddy_composite.png\n');

%% 8. D5.2 — Velocity & KE underestimation numbers --------------------------
% Compares along-track vs mapped radial profiles at the optimal window.
% Truth for the mapped product = mapped_W_grid (same window as along-track composite).
% Units: ssh in m (mapped_extracted is raw ncread in m); check along-track units match.

[mz_r_at, rmid_r_at] = radialProfile(alongtrack_W,  eddyPath_fun_t, showplot=false);
[mz_r_mp, rmid_r_mp] = radialProfile(mapped_W_grid, eddyPath_fun_t, showplot=false);

% Raw finite-difference d(eta)/dr — no B-spline (D4 decision)
dhdr_at = diff(mz_r_at) ./ diff(rmid_r_at);  r_c_at = (rmid_r_at(1:end-1)+rmid_r_at(2:end))/2;
dhdr_mp = diff(mz_r_mp) ./ diff(rmid_r_mp);  r_c_mp = (rmid_r_mp(1:end-1)+rmid_r_mp(2:end))/2;

% Headline: peak |d(eta)/dr| — physically: peak rotational velocity location
v_peak_at = max(abs(dhdr_at),[],'omitnan');
v_peak_mp = max(abs(dhdr_mp),[],'omitnan');
pct_v_peak  = 100*(v_peak_at - v_peak_mp) / v_peak_mp;
pct_KE_peak = 100*((v_peak_at/v_peak_mp)^2 - 1);

% Robustness check: core-mean within r <= Le
Le      = 80e3;   % eddy core radius (m) — matches QG model init in eq:Gaussian
core_at = abs(r_c_at) <= Le;
core_mp = abs(r_c_mp) <= Le;
pct_v_core = 100*(mean(abs(dhdr_at(core_at)),'omitnan') - ...
                  mean(abs(dhdr_mp(core_mp)),'omitnan')) / ...
                  mean(abs(dhdr_mp(core_mp)),'omitnan');

fprintf('\n--- D5.2 Velocity & KE (single-mission Agulhas, W*=%d cycles) ---\n', W_opt_cycles)
fprintf('Peak |d(eta)/dr| gain (along-track vs mapped): %.1f%%\n', pct_v_peak)
fprintf('Implied peak-KE underestimation in mapped:     %.1f%%\n', pct_KE_peak)
fprintf('Core-mean |d(eta)/dr| gain (r<=Le=%.0fkm):     %.1f%%\n', Le/1e3, pct_v_core)
fprintf('-> discussion.tex:32 and :42-43: XX=%%.0f, KE-YY=%%.0f (revisit with multi-mission data)\n')