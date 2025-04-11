n=[0,100,234,250,320]+1;
for i=1:length(n)
t_i=eddy_field.t(n(i));
    xo(i) = eddyPath_fun_t.xe(t_i);
    yo(i) = eddyPath_fun_t.ye(t_i);
end
% figure;
% hold on
% jpcolor(xE/1e3, yE/1e3, zeta');
% if isfield(options,'contour') %core
%     plot((options.contour.core(1,:)-xo)/1e3,(options.contour.core(2,:)-yo)/1e3,'k','LineWidth',1)
% 
%     if isfield(options.contour,'shield') % shield
%     plot((options.contour.shield(1,:)-xo)/1e3,(options.contour.shield(2,:)-yo)/1e3,'k','LineWidth',1)
%     end
% end
% % r=mean(eddy.speed_radius{1});
% % th = 0:pi/50:2*pi;
% % plot(r * cos(th),r*sin(th))
% % legend('','mean radius')
% % jpcolor(xEbin, yEbin, AvgsshAccumBin)
% shading flat
% axis equal
% xlabel('Distance East (km)', 'FontName', 'times')
% ylabel('Distance North (km)', 'FontName', 'times')
% set(gca, 'fontname', 'times','fontsize',16)
% if exist('brewermap')
%     colormap(brewermap([], '-Spectral'))
% end
% c = colorbar('EastOutside');
% c.Label.String = '\zeta/f_0';
% c.Label.FontSize=16;
% c.Label.FontName='times';
% xlim([-250, 250]), ylim([-250, 250]);
row=1;
%%
figname='QG_eddy_zeta_snapshot';
%% subplots
main_fig = figure('color','w');
set(gcf,'Units','centimeters');  %Ensure all dimensions are in cm
set(gcf,'Position',[0 0 24 8]);  %Size of the figure. The first two numbers indicate the location of the lower left-hand corner.
                                    %The second two numbers indicate size

% Subplot labels (a) through (i)
labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'};
figure(main_fig); % Make sure we're on the main figure before creating subplots
clf(main_fig)
% Create axes for the colorbars at the top
cbar2_pos = [6.8, 5.5, 4.3, 0.4];  % Position for the second colorbar


% Second colorbar for variance
cbar2_ax = axes('Units','centimeters','Position', cbar2_pos);
colormap(cbar2_ax, brewermap([], '-Spectral'));
cbar2 = colorbar('peer', cbar2_ax, 'Location', 'North');
cbar2.FontName='times';
cbar2.FontSize=12;
cbar2.Position(1)=0.38;
cbar2.Position(2)=cbar2.Position(2)+0.06;
cbar2.Position(4)=0.04;
% clim setting the min as max(std_rt_test 1) and max as max(std_rt_test 3) 
clim_cbar2 = [-0.25,0.05];%round([max([min(OSSE_data{1}{4},[],'all'),min(fullfield_data{1}{4},[],'all')]), max([max(OSSE_data{1}{4},[],'all'),max(fullfield_data{1}{4},[],'all')])]*1e2,1);
caxis(cbar2_ax, clim_cbar2);
axis(cbar2_ax, 'off');
title(cbar2_ax, '\zeta/f_0', 'fontsize', 12, 'fontname', 'times');

% Create tighter subplots
% Define new positions for the subplots
subplot_width = 6;
subplot_height = 4;
hgap = 0.01;
vgap = 0.3;
col1_start = 0.7;
col2_start = col1_start+5.1;
col3_start = col2_start+5.1;
col4_start = col3_start+5.1;

row_starts = 1.5+(subplot_height+vgap).*[0];
i = row;
j=1;
zetaContours.core=core_zeta{n(j)};
zetaContours.shield=shield_zeta{n(j)};
    % Column 1: Plot mz (time-averaged SSH)
    axes('Units','centimeters','Position', [col1_start, row_starts(i), subplot_width, subplot_height]);
    %subplot(3, 3, 3*(i-1)+1);
    hold on;
    jpcolor((eddy_field.x - xo(j))/1e3, (eddy_field.y - yo(j))/1e3, zeta(:,:,n(j))');
    plot((zetaContours.core(1,:)-xo(j))/1e3,(zetaContours.core(2,:)-yo(j))/1e3,'k','LineWidth',1)
    plot((zetaContours.shield(1,:)-xo(j))/1e3,(zetaContours.shield(2,:)-yo(j))/1e3,'k','LineWidth',1)
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
    clim(clim_cbar2)
    xlim([-250, 250]);
    ylim([-250, 250]);
    text(-230, 210, labels{j}, 'fontsize', 14, 'fontname', 'times','color','k');
    text(-245,270,['Day ',num2str(n(j)-1)],'FontSize',14,'FontName','times')
    
    % Column 2: Plot stdz (temporal std)
    j=2;
    zetaContours.core=core_zeta{n(j)};
    zetaContours.shield=shield_zeta{n(j)};
    axes('Units','centimeters','Position', [col2_start, row_starts(i), subplot_width, subplot_height]);
    % subplot(3, 3, 3*(i-1)+2);
    hold on;
    jpcolor((eddy_field.x - xo(j))/1e3, (eddy_field.y - yo(j))/1e3, zeta(:,:,n(j))');
    plot((zetaContours.core(1,:)-xo(j))/1e3,(zetaContours.core(2,:)-yo(j))/1e3,'k','LineWidth',1)
    plot((zetaContours.shield(1,:)-xo(j))/1e3,(zetaContours.shield(2,:)-yo(j))/1e3,'k','LineWidth',1)
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
    text(-230, 210, labels{j}, 'fontsize', 14, 'fontname', 'times','color','k');
    text(-245,270,['Day ',num2str(n(j)-1)],'FontSize',14,'FontName','times')

    % Column 3: Plot stdz (temporal std)
    j=3;
    zetaContours.core=core_zeta{n(j)};
    zetaContours.shield=shield_zeta{n(j)};
    axes('Units','centimeters','Position', [col3_start, row_starts(i), subplot_width, subplot_height]);
    % subplot(3, 3, 3*(i-1)+2);
    hold on;
    jpcolor((eddy_field.x - xo(j))/1e3, (eddy_field.y - yo(j))/1e3, zeta(:,:,n(j))');
    plot((zetaContours.core(1,:)-xo(j))/1e3,(zetaContours.core(2,:)-yo(j))/1e3,'k','LineWidth',1)
    plot((zetaContours.shield(1,:)-xo(j))/1e3,(zetaContours.shield(2,:)-yo(j))/1e3,'k','LineWidth',1)
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
    text(-230, 210, labels{j}, 'fontsize', 14, 'fontname', 'times','color','k');
    text(-245,270,['Day ',num2str(n(j)-1)],'FontSize',14,'FontName','times')

    % Column 4: Plot stdz (temporal std)
    j=4;
    zetaContours.core=core_zeta{n(j)};
    zetaContours.shield=shield_zeta{n(j)};
    axes('Units','centimeters','Position', [col4_start, row_starts(i), subplot_width, subplot_height]);
    % subplot(3, 3, 3*(i-1)+2);
    hold on;
    jpcolor((eddy_field.x - xo(j))/1e3, (eddy_field.y - yo(j))/1e3, zeta(:,:,n(j))');
    plot((zetaContours.core(1,:)-xo(j))/1e3,(zetaContours.core(2,:)-yo(j))/1e3,'k','LineWidth',1)
    plot((zetaContours.shield(1,:)-xo(j))/1e3,(zetaContours.shield(2,:)-yo(j))/1e3,'k','LineWidth',1)
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
    text(-230, 210, labels{j}, 'fontsize', 14, 'fontname', 'times','color','k');
    text(-245,270,['Day ',num2str(n(j)-1)],'FontSize',14,'FontName','times')
    
    % if i == 1  % Only add legend to the fi`rst radial profile plot
    %     lg = legend(h, '$\overline{\eta_{xy}}^\theta$', '$\overline{\sigma_\eta}^\theta$', '$\varsigma_{\overline{\eta}^t}$', '$\Sigma_{\eta_{xy}}$', '$\overline{\varsigma_\eta}^t$', '$\sigma_{\overline{\eta}^\theta}$', '$\Sigma_{\eta_{rt}}$');
    %     set(lg, 'interpreter', 'latex', 'fontsize', 12, 'orientation', 'vertical', 'NumColumns', 2);
    % end
   
%%
% Save
folder_name='E:\Research\AlongTrack-GRL\fig';%'D:\UW\AlongTrack-GRL\fig';%;%'D:\UW\AlongTrack-GRL\fig';
set(gcf,'Renderer','opengl');         % somehow painters or completely vectorize creates white lines
set(gca, 'Color', 'w'); % Sets axis background to white
set(gcf,'inverthardcopy','off') %to save white text as white
set(1,'paperpositionmode','auto')
print('-dpng','-r1200',strcat(folder_name,'\',figname,'.png'))
print('-depsc','-r1200',strcat(folder_name,'\',figname,'.eps'))
% exportgraphics(gcf, fullfile(folder_name, strcat(figname, '.eps')), 'ContentType', 'vector');