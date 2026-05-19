%% pub_numbers_multiMissions.m
% Run AFTER test_multiMissions.mlx has executed so workspace variables are live.
%
% Computes:
%   D2  — sample-density number to replace histogram in method-result_multi.tex:55
%   D5.1 — convergence cycles for T/P-only vs full constellation (discussion.tex:12-15, 40)
%
% Assumes workspace contains (from test_multiMissions.mlx):
%   numz_xy          — 2D bin count from composite2D on multi-mission OSSE
%   xmid_xy, ymid_xy — bin-centre vectors (m) matching numz_xy
%   alongtrack_TPonly — along-track struct for T/P-only mission
%   alongtrack_multi  — along-track struct for full constellation
%   eddyPath_fun_t    — eddy-path function handle struct

%% D2 — Sample density inside eddy core (replaces histogram, tex:55) ---------
Le = 80e3;   % effective eddy radius used in QG init (m)
[Xb, Yb] = meshgrid(xmid_xy, ymid_xy);
core_mask = sqrt(Xb.^2 + Yb.^2) <= Le;
samples_core = numz_xy(core_mask);
samples_core = samples_core(~isnan(samples_core) & samples_core > 0);

med_samples = median(samples_core);
q25 = quantile(samples_core, 0.25);
q75 = quantile(samples_core, 0.75);

fprintf('\n--- D2 Sample Density (r <= Le = %.0f km) ---\n', Le/1e3)
fprintf('Median: %.1f samples per (1/8 deg)^2 bin\n', med_samples)
fprintf('IQR:    [%.1f, %.1f]\n', q25, q75)
fprintf('--> Insert into tex:55: "median ~%.0f observations per bin (IQR [%.0f, %.0f])"\n', ...
        med_samples, q25, q75)

%% D5.1 — Convergence curves: T/P-only vs full constellation (tex:12-15) ----
[mz_rt_tpo, ~, ~, numz_rt_tpo] = radialProfileTime(alongtrack_TPonly, eddyPath_fun_t, showplot=false);
[mz_rt_mm,  ~, ~, numz_rt_mm ] = radialProfileTime(alongtrack_multi,  eddyPath_fun_t, showplot=false);

conv_tpo = convergenceRate(mz_rt_tpo, numz_rt_tpo);
conv_mm  = convergenceRate(mz_rt_mm,  numz_rt_mm);

% --- Plot side-by-side so elbow can be read visually ---
figure('Position',[100 100 800 400]);
semilogy(1:numel(conv_tpo), conv_tpo, 'k-', 'LineWidth',2, 'DisplayName','T/P-only'); hold on
semilogy(1:numel(conv_mm),  conv_mm,  'r-', 'LineWidth',2, 'DisplayName','Multi-mission');
yline(0.10,'k--','10% threshold','LabelVerticalAlignment','bottom');
xlabel('Cycle','FontName','times');
ylabel('RMS relative convergence error','FontName','times');
legend('Location','northeast','FontName','times'); grid on;
set(gca,'FontName','times','FontSize',14);
title('Convergence: T/P-only vs full constellation','FontName','times')
% --> Inspect curve. The "elbow" is where the rate of improvement slows sharply.
%     Set n_tpo and n_mm below after inspection.

% --- Fill in after inspecting the plot ---
n_tpo = NaN;   % <-- set to elbow cycle for T/P-only
n_mm  = NaN;   % <-- set to elbow cycle for multi-mission

if ~isnan(n_tpo) && ~isnan(n_mm)
    pct_reduction = 100*(1 - n_mm/n_tpo);
    fprintf('\n--- D5.1 Convergence ---\n')
    fprintf('T/P-only:    %d cycles (%.0f days) at elbow\n', n_tpo, n_tpo*9.92)
    fprintf('Multi-miss:  %d cycles (%.0f days) at elbow\n', n_mm,  n_mm*9.92)
    fprintf('Time reduction with full constellation: %.0f%%\n', pct_reduction)
    fprintf('\n--> discussion.tex:12-13: "approximately %d cycles (~%.0f days)"\n', n_tpo, n_tpo*9.92)
    fprintf('--> discussion.tex:14-15: "reduces required observation time by ~%.0f%%, achieving comparable accuracy in %d cycles"\n', pct_reduction, n_mm)
else
    fprintf('\n[Set n_tpo and n_mm after inspecting the convergence plot]\n')
end
