%% realEddy_mappedComparison.m
% Generates Fig 7 (mappedComparison.png) for the publication.
% Style matches plotFitProfile.m and radialProfile.m internal plot.
%
% 2×2 layout:
%   Row 1 (top):    QG model multi-mission OSSE
%   Row 2 (bottom): Agulhas Ring +527413 (real, single-mission)
%   Col 1 (left):   Time-averaged radial SSH profile  \bar{\eta}^{t,\theta}(r)
%   Col 2 (right):  Raw radial SSH gradient  d\eta/dr
%
% Each row overlays non-overlapping fixed-width windows sliding through the
% eddy lifetime.  Alpha increases with time (earlier = lighter, later = darker),
% showing whether the composite stabilises regardless of when the window falls.
%   QG model  : W* = 10 days  → windows [1–10, 11–20, 21–30, …]
%   Agulhas   : W* = 200 days → windows [1–200, 201–400, …]
% Truth (full QG field / CMEMS mapped) is plotted in black.
%
% Run AFTER test_multiMissions.mlx AND test_realEddy.mlx so workspace contains:
%   QG row   — fullfield_mm      : QG truth grid struct (ssh[nx,ny,nt], x,y m, t)
%              alongtrack_mm     : multi-mission along-track struct (t,x,y,ssh)
%              eddyPath_fun_t_mm : eddy path function handles (.xe, .ye)
%              (falls back to loading alongtrack_QG.mat if not found)
%   Agulhas  — mapped_extracted  : CMEMS grid struct (ssh[nx,ny,nt], x,y m, t)
%              alongtrack        : single-mission along-track struct (t,x,y,ssh)
%              eddyPath_fun_t    : eddy path function handles (.xe, .ye)
%
% Saves: D:\UW\AlongTrack-GRL\fig\mappedComparison.png  (.eps)

%% --- Parameters -----------------------------------------------------------
W_qg = 10;    % QG optimal window width (days) — 1 reference cycle
W_re = 200;   % Agulhas optimal window width (days) — W* from MSE sweep

figname     = 'QGmodel_realEddy_radialProfile';
folder_name = 'D:\UW\AlongTrack-GRL\fig';

%% --- Load / assign QG data ------------------------------------------------
% CHECK: verify alongtrack_mm / fullfield_mm are the multi-mission workspace
%        variables from test_multiMissions.mlx.
if exist('alongtrack_mm','var') && exist('fullfield_mm','var')
    qg_at_full  = alongtrack_mm;
    qg_truth    = fullfield_mm;
    eddyPath_qg = eddyPath_fun_t_mm;
else
    warning('alongtrack_mm / fullfield_mm not found — loading QG mat files (single-mission fallback)');
    qg_vars     = load('D:\UW\EddyLab\alongtrackFitting\alongtrack_QG.mat');
    qg_at_full  = qg_vars.alongtrack;
    qg_truth    = qg_vars.fullfield;
    eddyPath_qg = qg_vars.eddyPath_fun_t;
end

%% --- QG: build non-overlapping W_qg-day windows --------------------------
t0_qg   = min(qg_at_full.t);
tmax_qg = max(qg_at_full.t);
edges_qg = t0_qg : W_qg : tmax_qg;           % window-start edges (datenum)
n_win_qg = numel(edges_qg) - 1;

[mz_r_qg_truth, rmid_qg] = radialProfile(qg_truth, eddyPath_qg, showplot=false);

mz_r_qg = cell(n_win_qg, 1);
for k = 1:n_win_qg
    in_w = qg_at_full.t >= edges_qg(k) & qg_at_full.t < edges_qg(k+1);
    if sum(in_w) < 2, continue; end
    sub = struct('t', qg_at_full.t(in_w), 'x', qg_at_full.x(in_w), ...
                 'y', qg_at_full.y(in_w), 'ssh', qg_at_full.ssh(in_w));
    [mz_r_qg{k}, ~] = radialProfile(sub, eddyPath_qg, showplot=false);
end

%% --- Agulhas: build non-overlapping W_re-day windows ---------------------
% CHECK: mapped_extracted.ssh and alongtrack.ssh must be in the same units (m).
t0_re   = min(alongtrack.t);
tmax_re = max(alongtrack.t);
edges_re = t0_re : W_re : tmax_re;
n_win_re = numel(edges_re) - 1;

[mz_r_re_truth, rmid_re] = radialProfile(mapped_extracted, eddyPath_fun_t, showplot=false);

mz_r_re = cell(n_win_re, 1);
for k = 1:n_win_re
    in_w = alongtrack.t >= edges_re(k) & alongtrack.t < edges_re(k+1);
    if sum(in_w) < 2, continue; end
    sub = struct('t', alongtrack.t(in_w), 'x', alongtrack.x(in_w), ...
                 'y', alongtrack.y(in_w), 'ssh', alongtrack.ssh(in_w));
    [mz_r_re{k}, ~] = radialProfile(sub, eddyPath_fun_t, showplot=false);
end

%% --- Figure setup ---------------------------------------------------------
main_fig = figure('color','w');
set(gcf,'Units','centimeters');
set(gcf,'Position',[0 0 18 12]);

subplot_width  = 7.5;
subplot_height = 4.3;
col1_start     = 0.8;
col2_start     = 9.8;
vgap           = 0.8;
row_starts     = 1.2 + (subplot_height + vgap) .* [1, 0];

labels     = {'(a)', '(b)', '(c)', '(d)'};
row_labels = {'QG model', 'Agulhas Ring'};

%% --- Plot rows ------------------------------------------------------------
for i = 1:2   % i=1 QG, i=2 Agulhas

    if i == 1
        mz_truth = mz_r_qg_truth;
        rmid     = rmid_qg;
        mz_win   = mz_r_qg;
        n_win    = n_win_qg;
        W_days   = W_qg;
    else
        mz_truth = mz_r_re_truth;
        rmid     = rmid_re;
        mz_win   = mz_r_re;
        n_win    = n_win_re;
        W_days   = W_re;
    end

    % Drop any empty cells (windows with < 2 points)
    valid = ~cellfun(@isempty, mz_win);
    mz_win   = mz_win(valid);
    n_valid  = sum(valid);
    % Alpha: earlier windows lighter, later darker
    alpha_vals = linspace(0.25, 1.0, n_valid);
    % Red ramp: light salmon → brick red
    red_ramp = [linspace(0.95, 0.65, n_valid)', ...
                linspace(0.30, 0.00, n_valid)', ...
                linspace(0.20, 0.00, n_valid)'];

    % ---- Column 1: time-averaged SSH profile ----
    axes('Units','centimeters','Position',[col1_start, row_starts(i), subplot_width, subplot_height]);
    hold on;

    for k = 1:n_valid
        plot(rmid/1e3, mz_win{k}*1e2, '-', 'LineWidth', 1.5, ...
             'Color', [red_ramp(k,:), alpha_vals(k)]);
    end
    h_truth = plot(rmid/1e3, mz_truth*1e2, 'k-', 'LineWidth', 2.5);
    hlines(0, 'k:');

    if i == 2
        xlabel('Radial distance (km)', 'FontName','times');
    else
        set(gca,'XTickLabel',[]);
    end
    ylabel('SSH (cm)', 'FontName','times');
    set(gca,'FontName','times','FontSize',16);
    xlim([0, 250]);
    yl = ylim;
    text(8, yl(1) + 0.92*(yl(2)-yl(1)), labels{2*(i-1)+1}, ...
         'FontSize',14,'FontName','times','FontWeight','bold');
    text(8, yl(1) + 0.80*(yl(2)-yl(1)), row_labels{i}, ...
         'FontSize',11,'FontName','times','Color',[0.45 0.45 0.45]);

    % Legend on top-left (QG) panel only
    if i == 1
        % Proxy handles for first, last window, and truth
        h_first = plot(nan,nan,'-','LineWidth',1.5,'Color',[red_ramp(1,:), alpha_vals(1)]);
        h_last  = plot(nan,nan,'-','LineWidth',1.5,'Color',[red_ramp(end,:), alpha_vals(end)]);
        lg = legend([h_truth, h_first, h_last], ...
            {'$\overline{\eta}^{t,\theta}$ truth', ...
             sprintf('$t_1$–$%g$ d', W_days), ...
             sprintf('later windows')}, ...
            'Location','northeast');
        set(lg,'interpreter','latex','FontName','times','FontSize',13, ...
               'orientation','vertical','NumColumns',1);
    end

    % ---- Column 2: raw d(SSH)/dr ----
    axes('Units','centimeters','Position',[col2_start, row_starts(i), subplot_width, subplot_height]);
    hold on;

    for k = 1:n_valid
        plot(rmid/1e3, vdiff(mz_win{k}*1e2, 1), '-', 'LineWidth', 1.5, ...
             'Color', [red_ramp(k,:), alpha_vals(k)]);
    end
    plot(rmid/1e3, vdiff(mz_truth*1e2, 1), 'k-', 'LineWidth', 2.5);
    hlines(0, 'k:');   % zero reference — matches radialProfile.m style

    if i == 2
        xlabel('Radial distance (km)', 'FontName','times');
    else
        set(gca,'XTickLabel',[]);
    end
    ylabel('$\partial\eta/\partial r$ (cm\,km$^{-1}$)', ...
           'interpreter','latex','FontName','times');
    set(gca,'FontName','times','FontSize',16);
    xlim([0, 250]);
    yl = ylim;
    text(8, yl(1) + 0.92*(yl(2)-yl(1)), labels{2*(i-1)+2}, ...
         'FontSize',14,'FontName','times','FontWeight','bold');
end

%% --- Save -----------------------------------------------------------------
set(gcf,'Renderer','opengl');
set(gcf,'inverthardcopy','off');
set(1,'paperpositionmode','auto');
print('-dpng','-r600',  strcat(folder_name, '\', figname, '.png'));
print('-depsc','-r1200', strcat(folder_name, '\', figname, '.eps'));
fprintf('Saved %s\n', fullfile(folder_name, figname));
