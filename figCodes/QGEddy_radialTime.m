%% WV Eddy
%alongtrackLonLat from OceanDB
%eddyPath_track from ssh max closed contour centroid
%eddy_field is an array of eddy field in (x,y,t,ssh)
clear all

%% 1. Extract tracks along the eddy path
% define domain
current_file_path = matlab.desktop.editor.getActiveFilename;
[current_dir, ~, ~] = fileparts(current_file_path);
readdir = 'G:\My Drive\AlongTrack\';
filename = 'BetaEddyOne.nc';
% Load Model
x_QG = ncread([readdir, filename], 'x'); %meters
y_QG = ncread([readdir, filename], 'y'); %meters
ssh = squeeze(ncread([readdir, filename], 'ssh')); %meter;% matrix order in x,y,z
lato = 24;

totalDays = size(ssh, 3);
eddy_time = ([1:totalDays]-1)';
x = x_QG - mean(x_QG);
y = y_QG - mean(y_QG);

eddy_field.x = x;
eddy_field.y = y;
eddy_field.t = eddy_time;
eddy_field.ssh = ssh;

%find eddy center
[center_xoyo,amplitude] = findEddyCentroid(x, y, ssh,'thresholdratio',0.9,'GetBoundary', false);
eddyPath.xe = center_xoyo(:,1);
eddyPath.ye = center_xoyo(:,2);

%if you want to change the eddyPath to a function handle
eddyPath_fun_t.xe = @(t) interp1(eddy_field.t, eddyPath.xe, t, 'linear', 'extrap');
eddyPath_fun_t.ye = @(t) interp1(eddy_field.t, eddyPath.ye, t, 'linear', 'extrap');

% % extract track - takes somes time bc it loads the entire track matrix

%timeo=datenum(1992,9,25)+(0:totalDays-1)';
timeo = datenum(2000, 1, 1) + eddy_time' + 90; %Initial time is also a free parameter

%% Alongtrack from Jonathan (3D matrix [atd, tracknumber, cycle])
% JasonAlongTrack.filename = strcat(readdir, 'JasonAlongTrack.nc');
% JasonAlongTrack.lat = ncread(JasonAlongTrack.filename, 'lat');
% JasonAlongTrack.lon = ncread(JasonAlongTrack.filename, 'lon');
% %JML convert time to Matlab's datenum format
% JasonAlongTrack.time = ncread(JasonAlongTrack.filename, 'time') + datenum(1950, 1, 1);
% 
% alongtrackLatLon = alongtrackFromXYDomain(JasonAlongTrack,x,y,timeo,lato=lato,lono=-40); %options: lato=24, lono=305
% save(strcat(current_dir,'\alongtrackLatLon_QG.mat'),'alongtrackLatLon')

load(strcat(current_dir,'\alongtrackLatLon_QG.mat'))

alongtrackXY = latlon2xy_centered(alongtrackLatLon);

%% 2. Apply OSSE on your choice of an analytical eddy shape
%Output: alongtrack - contains arrays of lon,lat,x,y,t,ssh of OSSE 
% apply OSSE on an eddy_field from QG model in (x,y,t)
alongtrack = subsampleOSSE(alongtrackXY,eddy_field);

% load('D:\UW\EddyLab\alongtrack_QG.mat')

%% eddy center from eddyPath
% Full field struct
[X,Y,T] = ndgrid(1:length(eddy_field.x), 1:length(eddy_field.y), 1:length(eddy_field.t));
fullfield.x = eddy_field.x(reshape(X,[],1));
fullfield.y = eddy_field.y(reshape(Y,[],1));
fullfield.t = eddy_field.t(reshape(T,[],1)); % Subtract 1 to match your indexing
fullfield.ssh = reshape(eddy_field.ssh,[],1);

%%
%Full
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(fullfield,eddyPath_fun_t,showplot=0);
data{1} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

%OSSE
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(alongtrack,eddyPath_fun_t,showplot=0);
data{2} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

row=2;
%%
figname='QG_eddy_radtime';
plot_data = data;
%% subplots
main_fig = figure('color','w');
set(gcf,'Units','centimeters');  %Ensure all dimensions are in cm
set(gcf,'Position',[0 0 13 12]);  %Size of the figure. The first two numbers indicate the location of the lower left-hand corner.
                                    %The second two numbers indicate size
num_cols=2;
% Subplot labels (a) through (i)
labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'};
figure(main_fig); % Make sure we're on the main figure before creating subplots
clf(main_fig)

% Compute subplot width dynamically
total_margin_fraction = 0.225;  % Total margin space as fraction of figure width
subplot_width = 6;%(1 - total_margin_fraction) / num_cols;  % Adjust width based on column count
subplot_height = 4;  % Keep height same

% Compute dynamic margins
margins = (1 - num_cols * subplot_width) / (num_cols + 1);

% Compute column positions dynamically
col_starts = 0.03+ margins + (0:num_cols-1) * (subplot_width + margins);

% Compute center positions of each column for colorbars
cbar_width = subplot_width * 0.9;  % Slightly smaller than subplot width
cbar_height = 0.4;  % Fixed small height
cbar_y = 10.7;  % Fixed height near the top

cbar_positions = [col_starts + subplot_width / 2 - cbar_width / 2; 
                  repmat(cbar_y, 1, num_cols);  
                  repmat(cbar_width, 1, num_cols);
                  repmat(cbar_height, 1, num_cols)];

% Create first colorbar
cbar1_ax = axes('Units','centimeters','Position', cbar_positions(:,1)');
colormap(cbar1_ax, brewermap([], '-Spectral'));
cbar1 = colorbar('peer', cbar1_ax, 'Location', 'North');
cbar1.Position(4)=0.02;
cbar1.Position(2)=cbar1.Position(2)-0.06;
cbar1.FontName = 'times';
cbar1.FontSize = 12;
clim_cbar1 = [-2, 13];
caxis(cbar1_ax, clim_cbar1);
axis(cbar1_ax, 'off');
title(cbar1_ax, 'Time-averaged SSH (cm)', 'fontsize', 12, 'fontname', 'times');

% Create second colorbar
cbar2_ax = axes('Units','centimeters','Position', cbar_positions(:,2)');
colormap(cbar2_ax, brewermap([], '-Spectral'));
cbar2 = colorbar('peer', cbar2_ax, 'Location', 'North');
cbar2.Position(4)=0.02;
cbar2.Position(2)=cbar2.Position(2)-0.06;
cbar2.FontName = 'times';
cbar2.FontSize = 12;
clim_cbar2 = [0,2];
caxis(cbar2_ax, clim_cbar2);
axis(cbar2_ax, 'off');
title(cbar2_ax, 'Azimuthal variance (cm)', 'fontsize', 12, 'fontname', 'times');

% Set legend position to match third column (or last column if fewer than 3)
% legend_pos = cbar_positions(:, min(num_cols, 3))'; 
% cbar3_ax = axes('Position', legend_pos);
% axis(cbar3_ax, 'off');  % Replace with actual legend command if needed
% title(cbar3_ax, 'Legend', 'fontsize', 12, 'fontname', 'times');



% col3_start = col2_start + subplot_width + margin;
row_starts = 0.08+(subplot_height+0.02).*[1,0];


for i = 1:row
    % Unpack data for this test
    mz_rt = plot_data{i}{1};
    rmid_rt = plot_data{i}{2};
    tmid_rt = plot_data{i}{3};
    numz_rt = plot_data{i}{4};
    stdz_rt = plot_data{i}{5};
    
    % Column 1: Plot mz (time-averaged SSH)
    axes('Units','centimeters','Position', [col_starts(1), row_starts(i), subplot_width, subplot_height]);
    %subplot(3, 3, 3*(i-1)+1);
    hold on;
    jpcolor(rmid_rt/1e3, [1:length(tmid_rt)], mz_rt*1e2)%mz_rt'./vmean(mz_rt,2)'
    shading flat;
    % axis equal;
    if(i~=row)
        set(gca,'XTickLabel',[])
    else
        xlabel('Radial distance (km)', 'FontName', 'times')
    end
    ylabel('Time (Cycle)', 'FontName', 'times')
    set(gca, 'fontname', 'times', 'fontsize', 12);
    colormap(brewermap([], '-Spectral'));
    % colorbar('EastOutside');
    clim(clim_cbar1)
    xlim([0, 250]);ylim([1,length(tmid_rt)]-1)
    text(215, 31, labels{3*(i-1)+1}, 'fontsize', 14, 'fontname', 'times','color','w');
    
    % Column 2: Plot stdz (temporal std)
    axes('Units','centimeters', [col_starts(2), row_starts(i), subplot_width, subplot_height]);
    % subplot(3, 3, 3*(i-1)+2);
    hold on;
    jpcolor(rmid_rt/1e3, [1:length(tmid_rt)], stdz_rt*1e2) %mz_rt'./vmean(mz_rt,2)'
    shading flat;
    % axis equal;
    if(i~=row)
        set(gca,'XTickLabel',[])
    else
        xlabel('Radial distance (km)', 'FontName', 'times')
    end
    % ylabel('Time (Cycle)', 'FontName', 'times')
    set(gca,'YTickLabel',[])
    set(gca, 'fontname', 'times', 'fontsize', 12);
    colormap(brewermap([], '-Spectral'));
    % colorbar('EastOutside');
    clim(clim_cbar2)
    xlim([0, 250]);ylim([1,length(tmid_rt)]-1)
    text(215, 31, labels{3*(i-1)+2}, 'fontsize', 14, 'fontname', 'times','color','w');
    
    % % Column 3: Plot numz (histogram)
    % subplot('Position', [col_starts(3), row_starts(i), subplot_width, subplot_height]);
    % % subplot(3, 3, 3*(i-1)+2);
    % hold on;
    % jpcolor(rmid_rt/1e3, [1:length(tmid_rt)], numz_rt) %mz_rt'./vmean(mz_rt,2)'
    % shading flat;
    % axis equal;
    % if(i==1||i==2)
    %     set(gca,'XTickLabel',[])
    % else
    %     xlabel('Radial distance (km)', 'FontName', 'times')
    % end
    % % ylabel('Time (Cycle)', 'FontName', 'times')
    % set(gca,'YTickLabel',[])
    % 
    % set(gca, 'fontname', 'times', 'fontsize', 12);
    % colormap(brewermap([], '-Spectral'));
    % % colorbar('EastOutside');
    % clim(clim_cbar3)
    % xlim([0, 250])
    % text(-230, 210, labels{3*(i-1)+2}, 'fontsize', 14, 'fontname', 'times','color','w');
end

%  ZDATA1:  surface zdata
%  YDATA1:  surface ydata
%  XDATA1:  surface xdata
%  CDATA1:  surface cdata
%  CDATA2:  surface cdata
%  ZDATA2:  surface zdata
%  YDATA2:  surface ydata
%  XDATA2:  surface xdata
%  CDATA3:  surface cdata
%  CDATA4:  surface cdata
%%
% Get all subplot axes
ax = findobj(gcf, 'Type', 'axes');

% Select the subplot (e.g., first subplot)
ax1 = ax(4);  % Adjust index based on your layout
ax2 = ax(3);
ax3 = ax(2);
ax4 = ax(1);
% Find the surface object within the selected subplot
hSurf = findobj(ax1, 'Type', 'surface');
hSurf2 = findobj(ax2, 'Type', 'surface');
hSurf3 = findobj(ax3, 'Type', 'surface');
hSurf4 = findobj(ax4, 'Type', 'surface');



% Extract the required data
ZData1 = get(hSurf, 'ZData');
YData1 = get(hSurf, 'YData');
XData1 = get(hSurf, 'XData');
CData1 = get(hSurf, 'CData');  % First color data
CData2 = get(hSurf2, 'CData');
ZData2 = get(hSurf3, 'ZData');
YData2 = get(hSurf3, 'YData');
XData2 = get(hSurf3, 'XData');
CData3 = get(hSurf3, 'CData');
CData4 = get(hSurf4, 'CData');
%%
createfigure_radialProfileTime_QG(ZData1, YData1, XData1, CData1, CData2, ZData2, YData2, XData2, CData3, CData4)
%%
% Save
folder_name='E:\Research\AlongTrack-GRL\fig';%'D:\UW\AlongTrack-GRL\fig';
set(gcf,'Renderer','opengl');         %Instead of painter,opengl\
set(gca, 'Color', 'w'); % Sets axis background to white
set(gcf,'inverthardcopy','off') %to save white text as white
set(1,'paperpositionmode','auto')
print('-dpng','-r1200',strcat(folder_name,'\',figname,'.png'))
print('-depsc','-r1200',strcat(folder_name,'\',figname,'.eps'))
% exportgraphics(gcf, fullfile(folder_name, strcat(figname, '.eps')), 'ContentType', 'vector');

% epsclean(strcat(folder_name,'\',figname,'.eps'),'closeGaps',true)