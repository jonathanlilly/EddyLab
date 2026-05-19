%% realEddy_composites_radtime.m
% Generates the combined Fig 5 (real_eddy_composite.png) for the publication.
% Style matches QGEddy_composites_radtime.m exactly.
%
% Run AFTER test_realEddy.mlx (section 7) so the workspace contains:
%   alongtrack        — full-life along-track struct (t,x,y,ssh)
%   mapped_extracted  — grid truth struct (ssh [nx,ny,nt], x/y [1D] m, t datenum)
%   eddyPath_fun_t    — eddy path function handles
%   W_opt_cycles      — MSE-optimal window in cycles (set in Section 7 of test_realEddy)
%
% Saves: D:\UW\AlongTrack-GRL\fig\real_eddy_composite.png  (.eps)

%% --- Check for required workspace variables --------------------------------
W_opt_days = 200;   % fallback — set manually if MSE sweep not yet run

cycle_days   = 10;
W_opt_cycles   = W_opt_days * cycle_days;
time_window=200+[1:W_opt_days];%1:length(timeo)-1;
%% --- Build optimal-window subsets ------------------------------------------
% Along-track: subset full-life struct to first W* days
%If you want to narrow the time window

in_W_at = find(alongtrack.t-min(alongtrack.t) >= time_window(1) & alongtrackLatLon.t-min(alongtrack.t) <= time_window(end));

alongtrack_W.t   = alongtrack.t(in_W_at);
alongtrack_W.x   = alongtrack.x(in_W_at);
alongtrack_W.y   = alongtrack.y(in_W_at);
alongtrack_W.ssh = alongtrack.ssh(in_W_at);

% Mapped: subset mapped_extracted grid to same elapsed days
mapped_W_grid.x   = mapped_extracted.x;
mapped_W_grid.y   = mapped_extracted.y;
mapped_W_grid.t   = mapped_extracted.t(time_window);
mapped_W_grid.ssh = mapped_extracted.ssh(:,:,time_window);

%% --- Compute composites ----------------------------------------------------
% Along-track at W* (panels a,b — 2D spatial)
[mz, xmid, ymid, ~, stdz] = composite2D(alongtrack_W, eddyPath_fun_t, showplot=0);
[mzxy, rmid, ~, stdz_rt]  = radialProfile(alongtrack_W, eddyPath_fun_t, showplot=0);
data{1} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

% Mapped at W* (panels d,e — 2D spatial)
[mz, xmid, ymid, ~, stdz] = composite2D(mapped_W_grid, eddyPath_fun_t, showplot=0);
[mzxy, rmid, ~, stdz_rt]  = radialProfile(mapped_W_grid, eddyPath_fun_t, showplot=0);
data{2} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

% Along-track full-life r-t (panel c — azimuthal variance)
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(alongtrack, eddyPath_fun_t, showplot=0);
data2{1} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

% Mapped full-life r-t (panel f — azimuthal variance)
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(mapped_extracted, eddyPath_fun_t, showplot=0, tbin_size=10);
data2{2} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

%% --- Figure setup ----------------------------------------------------------
row      = 2;
figname  = 'real_eddy_composite_radtime';
plot_data  = data;
plot_data2 = data2;

main_fig = figure('color','w');
set(gcf,'Units','centimeters');
set(gcf,'Position',[0 0 18 12]);

labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};
figure(main_fig);
clf(main_fig)

%% --- Shared colorbars at top -----------------------------------------------
cbar1_pos = [1.5,  10.7, 4.3, 0.4];
cbar2_pos = [6.8,  10.7, 4.3, 0.4];
cbar3_pos = [12.5, 10.7, 4.3, 0.4];

% Colorbar 1: time-averaged SSH
cbar1_ax = axes('Units','centimeters','Position', cbar1_pos);
colormap(cbar1_ax, brewermap([], '-Spectral'));
cbar1 = colorbar('peer', cbar1_ax, 'Location', 'North');
cbar1.FontName = 'times';  cbar1.FontSize = 12;
cbar1.Position(1) = 0.095;
cbar1.Position(2) = cbar1.Position(2) - 0.06;
cbar1.Position(4) = 0.02;
clim_cbar1 = [-10, 30];   % cm — adjust if Agulhas range differs
caxis(cbar1_ax, clim_cbar1);
axis(cbar1_ax, 'off');
title(cbar1_ax, 'Time-averaged SSH (cm)', 'fontsize',12,'fontname','times');

% Colorbar 2: temporal variance
cbar2_ax = axes('Units','centimeters','Position', cbar2_pos);
colormap(cbar2_ax, brewermap([], '-Spectral'));
cbar2 = colorbar('peer', cbar2_ax, 'Location', 'North');
cbar2.FontName = 'times';  cbar2.FontSize = 12;
cbar2.Position(1) = 0.38;
cbar2.Position(2) = cbar2.Position(2) - 0.06;
cbar2.Position(4) = 0.02;
clim_cbar2 = [0, 10];     % cm — adjust to data
caxis(cbar2_ax, clim_cbar2);
axis(cbar2_ax, 'off');
title(cbar2_ax, 'Temporal variance (cm)', 'fontsize',12,'fontname','times');

% Colorbar 3: azimuthal variance (r-t panels)
cbar3_ax = axes('Units','centimeters','Position', cbar3_pos);
colormap(cbar3_ax, brewermap([], '-Spectral'));
cbar3 = colorbar('peer', cbar3_ax, 'Location', 'North');
cbar3.FontName = 'times';  cbar3.FontSize = 12;
cbar3.Position(1) = 0.7;
cbar3.Position(2) = cbar3.Position(2) - 0.06;
cbar3.Position(4) = 0.02;
clim_cbar3 = [0, 10];     % cm — adjust to data
caxis(cbar3_ax, clim_cbar3);
axis(cbar3_ax, 'off');
title(cbar3_ax, 'Azimuthal variance (cm)', 'fontsize',12,'fontname','times');

%% --- Subplot layout (matches QGEddy_composites_radtime.m) ------------------
subplot_width  = 6;
subplot_height = 4;
vgap       = 0.3;
col1_start = 0.7;
col2_start = 5.8;
col3_start = 12.5;

row_starts = 1.5 + (subplot_height + vgap) .* [1, 0];

for i = 1:row
    % Unpack 2D composite data
    mz       = plot_data{i}{1};
    xmid     = plot_data{i}{2};
    ymid     = plot_data{i}{3};
    stdz     = plot_data{i}{4};

    % Unpack full-life r-t data
    rmid_rt  = plot_data2{i}{2};
    tmid_rt  = plot_data2{i}{3};
    stdz_rt2 = plot_data2{i}{5};   % azimuthal variance (matches QG col-3 variable name)
    n_cycles = length(tmid_rt);

    % --- Column 1: time-averaged SSH (2D x-y) ---
    axes('Units','centimeters','Position',[col1_start, row_starts(i), subplot_width, subplot_height]);
    hold on;
    jpcolor(xmid/1e3, ymid/1e3, mz*1e2);
    shading flat;
    axis equal;
    if i ~= row
        set(gca,'XTickLabel',[])
    else
        xlabel('Distance East (km)', 'FontName','times');
    end
    ylabel('Distance North (km)', 'FontName','times');
    set(gca, 'fontname','times', 'fontsize',12);
    colormap(brewermap([], '-Spectral'));
    clim(clim_cbar1);
    xlim([-250, 250]);  ylim([-250, 250]);
    text(-230, 210, labels{3*(i-1)+1}, 'fontsize',14,'fontname','times','color','k');

    % --- Column 2: temporal variance (2D x-y) ---
    axes('Units','centimeters','Position',[col2_start, row_starts(i), subplot_width, subplot_height]);
    hold on;
    jpcolor(xmid/1e3, ymid/1e3, stdz*1e2);
    shading flat;
    axis equal;
    if i ~= row
        set(gca,'XTickLabel',[])
    else
        xlabel('Distance East (km)', 'FontName','times');
    end
    set(gca,'YTickLabel',[]);
    set(gca, 'fontname','times', 'fontsize',12);
    colormap(brewermap([], '-Spectral'));
    clim(clim_cbar2);
    xlim([-250, 250]);  ylim([-250, 250]);
    text(-230, 210, labels{3*(i-1)+2}, 'fontsize',14,'fontname','times','color','k');

    % --- Column 3: azimuthal variance r-t (full lifetime) ---
    axes('Units','centimeters','Position',[col3_start, row_starts(i), subplot_width-1.7, subplot_height]);
    hold on;
    jpcolor(rmid_rt/1e3, 1:n_cycles, abs(stdz_rt2)*1e2);
    shading flat;
    % Mark the optimal window with a dashed white box spanning the full radial extent.
    % time_window is in days; divide by cycle_days to convert to cycle-axis units.
    box_y0 = time_window(1)  / cycle_days;   % start cycle of window
    box_y1 = time_window(end)/ cycle_days;   % end cycle of window
    rectangle('Position', [0, box_y0, 250, box_y1 - box_y0], ...
        'EdgeColor', [0.7,0.4,0.8], 'LineWidth', 1.5, 'LineStyle', '-');
    if i ~= row
        set(gca,'XTickLabel',[])
    else
        xlabel('Radial distance (km)', 'FontName','times');
    end
    ylabel('Time (Cycle)', 'FontName','times');
    set(gca, 'fontname','times', 'fontsize',12);
    colormap(brewermap([], '-Spectral'));
    clim(clim_cbar3);
    xlim([0, 250]);  ylim([1, n_cycles]-1);
    text(10, n_cycles*0.9, labels{3*(i-1)+3}, 'fontsize',14,'fontname','times','color','k');
end

%% --- Save ------------------------------------------------------------------
folder_name = 'D:\UW\AlongTrack-GRL\fig';
set(gcf,'Renderer','opengl');
set(gca, 'Color', 'w');
set(gcf,'inverthardcopy','off')
set(1,'paperpositionmode','auto')
print('-dpng','-r600', strcat(folder_name,'\',figname,'.png'))
print('-depsc','-r1200', strcat(folder_name,'\',figname,'.eps'))
fprintf('Saved %s\n', fullfile(folder_name, figname))
