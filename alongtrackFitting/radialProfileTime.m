
%% radial profile vs time
tmat = timet-timet(1);
[mz_rt, rmid, tmid, numz_rt, stdz_rt] = twodstats(sqrt(xEt.^2+yEt.^2), tmat', ssht, (0:binsize:max_r), 0.5+min(tmat):9.92:max(tmat)+0.5);

figure;
hold on; 
jpcolor(rmid, [1:length(tmid)], mz_rt), xlim([0, 250]) %mz_rt'./vmean(mz_rt,2)'
xlabel('Radial distance (km)', 'FontName', 'times')
ylabel('Time (Cycle)', 'FontName', 'times')
colormap(brewermap([], '-Spectral'))
c = colorbar('EastOutside');
c.Label.String = 'Azi-Avg SSH (cm)';
c.Label.FontSize=16;
c.Label.FontName='times';
set(gca, 'fontname', 'times','FontSize',16)
clim([-2.5, 25])

% std profile over time
figure;
hold on; 
jpcolor(rmid, [1:length(tmid)], stdz_rt), xlim([0, 250]) %mz_rt'./vmean(mz_rt,2)'
xlabel('Radial distance (km)', 'FontName', 'times')
ylabel('Time (Cycle)', 'FontName', 'times')
colormap(brewermap([], '-Spectral'))
c = colorbar('EastOutside');
c.Label.String = 'Azi-Avg Variance (cm)';
c.Label.FontSize=16;
c.Label.FontName='times';
set(gca, 'fontname', 'times','FontSize',16)

