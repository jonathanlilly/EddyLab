load('D:\UW\EddyLab\alongtrack_QG.mat')
load('D:\UW\EddyLab\alongtrackLatLon_QG.mat')

%%
%Full
[mz, xmid, ymid, ~, stdz] = composite2D(fullfield, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(fullfield, eddyPath_fun_t,showplot=0);
data{1} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

%OSSE
[mz, xmid, ymid, ~, stdz] = composite2D(alongtrack, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(alongtrack, eddyPath_fun_t,showplot=0);
data{2} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

row=2;
%%
figname='QG_eddy_composite';
plot_data = data;
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
legend_pos = [12.8, 10.7, 4.3, 0.4]; % Position for the legend at the top

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

% Third legend for variance
% Create an invisible axes for the legend at the top
legend_ax = axes('Units','centimeters','Position', legend_pos);

% Create dummy lines with the same formatting as will be used in the plots
% h_dummy = cell(7,1);
h = plot(NaN, nan(1,6)*1e2);
linestyle 2k 2W 2T-- 2U-- 2V-. 2X:
% Create the legend at the top
legend_top = legend(legend_ax, ...
    '$\overline{\eta_{xy}}^\theta$', '$\Sigma_{\eta_{xy}}$', '$\overline{\sigma_\eta}^\theta$', ...
    '$\varsigma_{\overline{\eta}^t}$', ...
    '$\overline{\varsigma_\eta}^t$', '$\sigma_{\overline{\eta}^\theta}$',...
    'interpreter', 'latex', 'fontsize', 10, 'orientation', 'vertical', 'NumColumns', 3);
% lg = legend(h_dummy, '$\overline{\eta_{xy}}^\theta$', '$\overline{\sigma_\eta}^\theta$', '$\varsigma_{\overline{\eta}^t}$', '$\Sigma_{\eta_{xy}}$', '$\overline{\varsigma_\eta}^t$', '$\sigma_{\overline{\eta}^\theta}$', '$\Sigma_{\eta_{rt}}$');
%         set(lg, 'interpreter', 'latex', 'fontsize', 12, 'orientation', 'vertical', 'NumColumns', 2);
% set(legend_top,'Position', legend_pos)
        % 
        % Position the legend box
 legend_top.Position=[0.67,0.84,0.32,0.055];
title(legend_ax, 'SSH variance (cm)', 'fontsize', 12, 'fontname', 'times');
axis(legend_ax, 'off');

% Create tighter subplots
% Define new positions for the subplots
subplot_width = 6;
subplot_height = 4;
hgap = 0.01;
vgap = 0.3;
col1_start = 0.7;
col2_start = 5.8;
col3_start = 12.5;

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

    % Column 3: Plot radial profiles
    axes('Units','centimeters','Position', [col3_start, row_starts(i), subplot_width-1.5, subplot_height-0.2]);
    % subplot(3, 3, 3*(i-1)+3);
    hold on;
    
    % Plot the radial profiles (assuming avgAziStdTemp, etc. are available)
    % Note: You might need to adjust this part depending on how these variables are defined in your script
    use stdz_rt
    % Apply line styles
    h = plot(rmid/1e3, [mzxy,  stdTotalxy, avgAziStdTemp, stdAziAvgTemp, avgTempStdAzi, stdTempAvgAzi]*1e2);
    linestyle 2k 2W 2T-- 2U-- 2V-. 2X:
    hlines(0, 'k:')

    if(i~=row)
        set(gca,'XTickLabel',[])
    else
        xlabel('Radial distance (km)', 'FontName', 'times');
    end
    % ylabel('SSH Variance (cm)', 'FontName', 'times');
    set(gca, 'ticklength', [0.03, 0.02], 'fontname', 'times', 'fontsize', 12);
    
    % if i == 1  % Only add legend to the fi`rst radial profile plot
    %     lg = legend(h, '$\overline{\eta_{xy}}^\theta$', '$\overline{\sigma_\eta}^\theta$', '$\varsigma_{\overline{\eta}^t}$', '$\Sigma_{\eta_{xy}}$', '$\overline{\varsigma_\eta}^t$', '$\sigma_{\overline{\eta}^\theta}$', '$\Sigma_{\eta_{rt}}$');
    %     set(lg, 'interpreter', 'latex', 'fontsize', 12, 'orientation', 'vertical', 'NumColumns', 2);
    % end
    
    xlim([0, 250]);ylim([-2,13])
    text(-60, max(get(gca, 'YLim')), labels{3*(i-1)+3}, 'fontsize', 14, 'fontname', 'times');
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