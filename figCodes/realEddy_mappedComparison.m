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
W_qg      = 10;             % QG optimal window width (days) — 1 reference cycle
W_re      = 200;            % Agulhas optimal window width (days) — W* from MSE sweep
% Window start days (1-indexed; window k covers days [starts(k), starts(k)+W-1]).
% Hand-picked to sample early / mid / late life of each eddy.
starts_qg = [1, 101, 251];  % QG : windows [1,10], [100,109], [250,259]
starts_re = [1, 201, 801];  % Agulhas: windows [1,200], [200,399], [800,999]
n_show    = numel(starts_qg);

figname     = 'QGmodel_realEddy_radialProfile';
folder_name = 'D:\UW\AlongTrack-GRL\fig';

%% --- Load / assign QG data ------------------------------------------------
% CHECK: verify alongtrack_mm / fullfield_mm are the multi-mission workspace
%        variables from test_multiMissions.mlx.

    qg_at_full  = currentMission;
    qg_truth    = eddy_field;
    eddyPath_qg = eddyPath_fun_t_mm;


%% --- QG: 3 paired (truth, along-track) windows ---------------------------
% For each window we build a *shifted* eddyPath so that radialProfile's
% internal  elapsed_time = t - min(t)  evaluates the path at the correct
% point in the eddy's life, not at day 0.  Both the along-track subset and
% the full-field grid are subset to the same window so that the black curves
% (truth) and red curves (along-track) are paired at the same life stage.
t0_qg = min(qg_at_full.t);

mz_r_qg       = cell(n_show, 1);   % along-track radial profile per window
rmid_r_qg     = cell(n_show, 1);
mz_r_qg_tr    = cell(n_show, 1);   % truth radial profile per window
rmid_r_qg_tr  = cell(n_show, 1);
labels_qg     = cell(n_show, 1);

for k = 1:n_show
    s = starts_qg(k);           % 1-indexed start day
    t_lo = t0_qg + (s - 1);     % datenum of first day in window
    t_hi = t_lo + W_qg;         % datenum of one-past-last day in window
    labels_qg{k} = sprintf('$[%d,%d]$ day', s, s + W_qg - 1);

    % Shifted eddyPath — offset in elapsed days = (s-1)
    ep_k.xe = @(dt) eddyPath_qg.xe(dt + (s-1));
    ep_k.ye = @(dt) eddyPath_qg.ye(dt + (s-1));

    % Along-track subset
    in_at = qg_at_full.t >= t_lo & qg_at_full.t < t_hi;
    if sum(in_at) >= 2
        sub_at = struct('t', qg_at_full.t(in_at), 'x', qg_at_full.x(in_at), ...
                        'y', qg_at_full.y(in_at), 'ssh', qg_at_full.ssh(in_at));
        [mz_r_qg{k}, rmid_r_qg{k}] = radialProfile(sub_at, ep_k, showplot=false);
    end

    % Truth (grid) subset — same time window
    in_tr = qg_truth.t >= t_lo-t0_qg & qg_truth.t < t_hi-t0_qg;
    if sum(in_tr) >= 1
        sub_tr = struct('x', qg_truth.x, 'y', qg_truth.y, ...
                        't', qg_truth.t(in_tr), ...
                        'ssh', qg_truth.ssh(:, :, in_tr));
        [mz_r_qg_tr{k}, rmid_r_qg_tr{k}] = radialProfile(sub_tr, ep_k, showplot=false);
    end
end

%% --- Agulhas: 3 paired (truth, along-track) windows ----------------------
% CHECK: mapped_extracted.ssh and alongtrack.ssh must be in the same units (m).
t0_re = min(alongtrack.t);

mz_r_re       = cell(n_show, 1);
rmid_r_re     = cell(n_show, 1);
mz_r_re_tr    = cell(n_show, 1);
rmid_r_re_tr  = cell(n_show, 1);
labels_re     = cell(n_show, 1);

for k = 1:n_show
    s = starts_re(k);
    t_lo = t0_re + (s - 1);
    t_hi = t_lo + W_re;
    labels_re{k} = sprintf('$[%d,%d]$ day', s, s + W_re - 1);

    ep_k.xe = @(dt) eddyPath_fun_t.xe(dt + (s-1));
    ep_k.ye = @(dt) eddyPath_fun_t.ye(dt + (s-1));

    in_at = alongtrack.t >= t_lo & alongtrack.t < t_hi;
    if sum(in_at) >= 2
        sub_at = struct('t', alongtrack.t(in_at), 'x', alongtrack.x(in_at), ...
                        'y', alongtrack.y(in_at), 'ssh', alongtrack.ssh(in_at));
        [mz_r_re{k}, rmid_r_re{k}] = radialProfile(sub_at, ep_k, showplot=false);
    end

    in_tr = mapped_extracted.t >= t_lo & mapped_extracted.t < t_hi;
    if sum(in_tr) >= 1
        sub_tr = struct('x', mapped_extracted.x, 'y', mapped_extracted.y, ...
                        't', mapped_extracted.t(in_tr), ...
                        'ssh', mapped_extracted.ssh(:, :, in_tr));
        [mz_r_re_tr{k}, rmid_r_re_tr{k}] = radialProfile(sub_tr, ep_k, showplot=false);
    end
end

%% --- Figure setup ---------------------------------------------------------
% Positions baked in from the MATLAB-generated createfigure_realEddy_mappedComparison
% (i.e. the layout you tweaked by hand in the figure editor and exported).
main_fig = figure('color','w');
set(gcf,'Units','centimeters');
set(gcf,'Position',[0 0 19 12.2]);   % wider so col2 + width fits

subplot_width  = 7.50;
subplot_height = 4.30;
col1_start     = 1.55;   % matches axes1 / axes3 from createfigure_*
col2_start     = 11.23;  % matches axes2 / axes4
vgap           = 1.38;   % 7.34 - 1.66 - 4.30
row_starts     = 1.66 + (subplot_height + vgap) .* [1, 0];  % [7.34, 1.66]

labels     = {'(a)', '(b)', '(c)', '(d)'};
row_labels = {'QG model', 'Agulhas Ring'};

%% --- Plot rows ------------------------------------------------------------
for i = 1:2   % i=1 QG, i=2 Agulhas

    if i == 1
        mz_win      = mz_r_qg;
        rmid_win    = rmid_r_qg;
        mz_tr_win   = mz_r_qg_tr;
        rmid_tr_win = rmid_r_qg_tr;
        win_labels  = labels_qg;
        truth_kind  = 'Full field';
    else
        mz_win      = mz_r_re;
        rmid_win    = rmid_r_re;
        mz_tr_win   = mz_r_re_tr;
        rmid_tr_win = rmid_r_re_tr;
        win_labels  = labels_re;
        truth_kind  = 'CMEMS';
    end

    % Drop any empty cells (windows that contained too few points)
    valid       = ~cellfun(@isempty, mz_win) & ~cellfun(@isempty, mz_tr_win);
    mz_win      = mz_win(valid);
    rmid_win    = rmid_win(valid);
    mz_tr_win   = mz_tr_win(valid);
    rmid_tr_win = rmid_tr_win(valid);
    win_labels  = win_labels(valid);
    n_valid     = sum(valid);
    % alpha: very faint early → fully opaque late
    alpha_vals = linspace(0.35, 1.0, n_valid);
    % Red ramp (along-track): pinkish red (faint) → pure deep red.
    % R held high, G = B descend together so the hue stays red (no orange).
    red_ramp = [linspace(1.00, 0.70, n_valid)', ...   % R
                linspace(0.70, 0.00, n_valid)', ...   % G
                linspace(0.70, 0.00, n_valid)'];      % B  (= G keeps it red, not orange)
    % Black ramp (truth): light grey → pure black
    blk_ramp = repmat(linspace(0.65, 0.00, n_valid)', 1, 3);

    % ---- Column 1: time-averaged SSH profile ----
    axes('Units','centimeters','Position',[col1_start, row_starts(i), subplot_width, subplot_height]);
    hold on;

    h_tr_arr = gobjects(n_valid,1);
    h_at_arr = gobjects(n_valid,1);
    for k = 1:n_valid
        % Truth window: black gradient, thicker
        h_tr_arr(k) = plot(rmid_tr_win{k}/1e3, mz_tr_win{k}*1e2, '-', ...
            'LineWidth', 2.0, 'Color', [blk_ramp(k,:), alpha_vals(k)]);
        % Along-track window: red gradient
        h_at_arr(k) = plot(rmid_win{k}/1e3, mz_win{k}*1e2, '-', ...
            'LineWidth', 1.5, 'Color', [red_ramp(k,:), alpha_vals(k)]);
    end
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

    % Legend: pair each window with both a black (truth) and red (along-track) swatch.
    % Use the along-track handles as legend entries with the window-day label,
    % plus a single black proxy at the top labelling the truth source.
    h_truth_proxy = plot(nan, nan, '-', 'LineWidth', 2.0, 'Color', [0 0 0]);
    h_at_proxy    = plot(nan, nan, '-', 'LineWidth', 1.5, 'Color', [0.75 0.1 0.05]);
    leg_handles = [h_truth_proxy, h_at_proxy, h_tr_arr(:)'];
    leg_strs    = [{['Truth (' truth_kind ')']}, {'Along-track'}, win_labels(:)'];
    lg = legend(leg_handles, leg_strs, 'Location','northeast');
    set(lg,'interpreter','latex','FontName','times','FontSize',10, ...
           'orientation','vertical','NumColumns',1);

    % ---- Column 2: raw d(SSH)/dr ----
    axes('Units','centimeters','Position',[col2_start, row_starts(i), subplot_width, subplot_height]);
    hold on;

    for k = 1:n_valid
        % Truth window gradient (black)
        plot(rmid_tr_win{k}/1e3, vdiff(mz_tr_win{k}*1e2, 1), '-', ...
            'LineWidth', 2.0, 'Color', [blk_ramp(k,:), alpha_vals(k)]);
        % Along-track window gradient (red)
        plot(rmid_win{k}/1e3, vdiff(mz_win{k}*1e2, 1), '-', ...
            'LineWidth', 1.5, 'Color', [red_ramp(k,:), alpha_vals(k)]);
    end
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
print('-dpng','-r1200',  strcat(folder_name, '\', figname, '.png'));
print('-depsc','-r1200', strcat(folder_name, '\', figname, '.eps'));
fprintf('Saved %s\n', fullfile(folder_name, figname));
