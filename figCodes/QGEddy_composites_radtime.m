load('D:\UW\EddyLab\alongtrackFitting\alongtrack_QG.mat')
load('D:\UW\EddyLab\alongtrackFitting\alongtrackLatLon_QG.mat')

%%
%Full
[mz, xmid, ymid, ~, stdz] = composite2D(fullfield, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(fullfield, eddyPath_fun_t,showplot=0);
data{1} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

%OSSE
[mz, xmid, ymid, ~, stdz] = composite2D(alongtrack, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(alongtrack, eddyPath_fun_t,showplot=0);
data{2} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

%Full
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(fullfield,eddyPath_fun_t,showplot=0);
data2{1} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

%OSSE
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(alongtrack,eddyPath_fun_t,showplot=0);
data2{2} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

row=2;
%%
figname='QG_eddy_composite_radtime';
plot_data = data;
plot_data2 = data2;
%% subplots
main_fig = figure('color','w');
set(gcf,'Units','centimeters');  %Ensure all dimensions are in cm
set(gcf,'Position',[0 0 18 12]);  %Size of the figure. The first two numbers indicate the location of the lower left-hand corner.
                                    %The second two numbers indicate size

% Subplot labels (a) through (i)
labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'};
figure(main_fig); % Make sure we're on the main figure before creating subplots
clf(main_fig)
% Create axes for the colorbars at the top
cbar1_pos = [1.5, 10.7, 4.3, 0.4]; % Position for the first colorbar
cbar2_pos = [6.8, 10.7, 4.3, 0.4];  % Position for the second colorbar
cbar3_pos = [12, 10.7, 4.3, 0.4]; % Position for the legend at the top

% First colorbar for SSH
cbar1_ax = axes('Units','centimeters','Position', cbar1_pos);
colormap(cbar1_ax, brewermap([], '-Spectral'));
cbar1 = colorbar('peer', cbar1_ax, 'Location', 'North');
cbar1.FontName='times';
cbar1.FontSize=12;
cbar1.Position(1)=0.095;
cbar1.Position(2)=cbar1.Position(2)-0.06;
cbar1.Position(4)=0.02;
clim_cbar1=[-2, 13];
caxis(cbar1_ax, clim_cbar1);
axis(cbar1_ax, 'off');
title(cbar1_ax, 'Time-averaged SSH (cm)', 'fontsize', 12, 'fontname', 'times');

% Second colorbar for variance
cbar2_ax = axes('Units','centimeters','Position', cbar2_pos);
colormap(cbar2_ax, brewermap([], '-Spectral'));
cbar2 = colorbar('peer', cbar2_ax, 'Location', 'North');
cbar2.FontName='times';
cbar2.FontSize=12;
cbar2.Position(1)=0.38;
cbar2.Position(2)=cbar2.Position(2)-0.06;
cbar2.Position(4)=0.02;
% clim setting the min as max(std_rt_test 1) and max as max(std_rt_test 3) 
clim_cbar2 = [0,1.5];%round([max([min(OSSE_data{1}{4},[],'all'),min(fullfield_data{1}{4},[],'all')]), max([max(OSSE_data{1}{4},[],'all'),max(fullfield_data{1}{4},[],'all')])]*1e2,1);
caxis(cbar2_ax, clim_cbar2);
axis(cbar2_ax, 'off');
title(cbar2_ax, 'Temporal variance (cm)', 'fontsize', 12, 'fontname', 'times');

% Third colorbar for variance
% Create an invisible axes for the legend at the top
cbar3_ax = axes('Units','centimeters','Position', cbar3_pos);
colormap(cbar3_ax, brewermap([], '-Spectral'));
cbar3 = colorbar('peer', cbar3_ax, 'Location', 'North');
cbar3.Position(1)=0.68;
cbar3.Position(4)=0.02;
cbar3.Position(2)=cbar3.Position(2)-0.06;
cbar3.FontName = 'times';
cbar3.FontSize = 12;
clim_cbar3 = [0,2];
caxis(cbar3_ax, clim_cbar3);
axis(cbar3_ax, 'off');
title(cbar3_ax, 'Azimuthal variance (cm)', 'fontsize', 12, 'fontname', 'times');

% Create tighter subplots
% Define new positions for the subplots
subplot_width = 6;
subplot_height = 4;
hgap = 0.01;
vgap = 0.3;
col1_start = 0.7;
col2_start = 5.8;
col3_start = 12;

row_starts = 1.5+(subplot_height+vgap).*[1,0];

for i = 1:row
    % Unpack data for this test
    mz = plot_data{i}{1};
    xmid = plot_data{i}{2};
    ymid = plot_data{i}{3};
    stdz = plot_data{i}{4};
    mzxy = plot_data{i}{5};
    rmid = plot_data{i}{6};
    stdz_rt = plot_data{i}{7};

    mz_rt = plot_data2{i}{1};
    rmid_rt = plot_data2{i}{2};
    tmid_rt = plot_data2{i}{3};
    numz_rt = plot_data2{i}{4};
    stdz_rt2 = plot_data2{i}{5};

    % Column 1: Plot mz (time-averaged SSH)
    axes('Units','centimeters','Position', [col1_start, row_starts(i), subplot_width, subplot_height]);
    %subplot(3, 3, 3*(i-1)+1);
    hold on;
    jpcolor(xmid/1e3, ymid/1e3, mz*1e2);
    shading flat;
    axis equal;
    if(i~=row)
        set(gca,'XTickLabel',[])
    else
        xlabel('Distance East (km)', 'FontName', 'times');
    end
    ylabel('Distance North (km)', 'FontName', 'times');
    set(gca, 'fontname', 'times', 'fontsize', 12);
    colormap(brewermap([], '-Spectral'));
    % colorbar('EastOutside');
    clim(clim_cbar1)
    xlim([-250, 250]);
    ylim([-250, 250]);
    text(-230, 210, labels{3*(i-1)+1}, 'fontsize', 14, 'fontname', 'times','color','k');
    
    % Column 2: Plot stdz (temporal std)
    axes('Units','centimeters','Position', [col2_start, row_starts(i), subplot_width, subplot_height]);
    % subplot(3, 3, 3*(i-1)+2);
    hold on;
    jpcolor(xmid/1e3, ymid/1e3, stdz*1e2);
    shading flat;
    axis equal;
    if(i~=row)
        set(gca,'XTickLabel',[])
    else
        xlabel('Distance East (km)', 'FontName', 'times');
    end
    % ylabel('Distance North (km)', 'FontName', 'times');
    set(gca,'YTickLabel',[])
    
    set(gca, 'fontname', 'times', 'fontsize', 12);
    colormap(brewermap([], '-Spectral'));
    % colorbar('EastOutside');
    clim(clim_cbar2)
    xlim([-250, 250]);
    ylim([-250, 250]);
    text(-230, 210, labels{3*(i-1)+2}, 'fontsize', 14, 'fontname', 'times','color','k');

    % Column 3: Plot stdz in time (azi std)
    axes('Units','centimeters','Position', [col3_start, row_starts(i), subplot_width-1.7, subplot_height]);
    % subplot(3, 3, 3*(i-1)+3);
    hold on;
    
    % Plot the radial profiles (assuming avgAziStdTemp, etc. are available)
    jpcolor(rmid_rt/1e3, [1:length(tmid_rt)], abs(stdz_rt2)*1e2) %mz_rt'./vmean(mz_rt,2)'
    shading flat;
    % axis equal;
    if(i~=row)
        set(gca,'XTickLabel',[])
    else
        xlabel('Radial distance (km)', 'FontName', 'times')
    end
    ylabel('Time (Cycle)', 'FontName', 'times')
    set(gca,'YTickLabel',[])
    set(gca, 'fontname', 'times', 'fontsize', 12);
    colormap(brewermap([], '-Spectral'));
    % colorbar('EastOutside');
    clim(clim_cbar3)
    xlim([0, 250]);ylim([1,length(tmid_rt)]-1)
    text(10, 32, labels{3*(i-1)+3}, 'fontsize', 14, 'fontname', 'times','color','k');
end

%%
% Save
folder_name='D:\UW\AlongTrack-GRL\fig';%'E:\Research\AlongTrack-GRL\fig';%'D:\UW\AlongTrack-GRL\fig';
set(gcf,'Renderer','opengl');         % somehow painters or completely vectorize creates white lines
set(gca, 'Color', 'w'); % Sets axis background to white
set(gcf,'inverthardcopy','off') %to save white text as white
set(1,'paperpositionmode','auto')
print('-dpng','-r1200',strcat(folder_name,'\',figname,'.png'))
print('-depsc','-r1200',strcat(folder_name,'\',figname,'.eps'))
% exportgraphics(gcf, fullfile(folder_name, strcat(figname, '.eps')), 'ContentType', 'vector');