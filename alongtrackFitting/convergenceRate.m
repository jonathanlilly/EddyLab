function [convergence] = convergenceRate(mz_rt, numz_rt)
% mz_all is time-averaged profile
% mz_rt is radial profile as a function of time
% Plot convergence rate to see how many tracks is needed
% Accumulated tracks at each cycle
cumAvg_track = cumsum(mz_rt.*numz_rt, 'omitnan') ./ cumsum(numz_rt);
Avg_all = cumAvg_track(end,:);
% Calculate convergence rate as relative difference from final value, which
% is the time-averaged radial profile
% Relative Convergence Error of Time-Azimuthal Average $m_z(r,t)$
%$\overline{h}^{\theta}(r,t)-\overline{h}(r)$ (\%)
convergence = rms(abs(cumAvg_track-Avg_all)./abs(Avg_all),2,'omitmissing');

% cumAvg_all = cumsum(mz_all.*numz_all, 'omitnan') ./ cumsum(numz_all);
% dmz = rms((cumAvg_track - cumAvg_all)./repmat(cumAvg_all(end, :), [size(cumAvg_all, 1), 1]), 2, 'omitmissing');
% convergence_rate = dmz/rms(cumAvg_all(end, :));
figure;
plot(1:size(mz_rt, 1), convergence, 'LineWidth', 2, 'Color', 'k');
xlabel('Track cycles', 'FontName', 'times')
ylabel('Relative convergence error of $\overline{h}^{\theta}(r,t)$', 'FontName', 'times','Interpreter','latex')
set(gca, 'fontname', 'times','FontSize',16)