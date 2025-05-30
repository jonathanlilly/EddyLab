%% Load eddy center tracking
eddy_id='+527413';
readdir = 'G:\My Drive\AlongTrack\';
filename = strcat([readdir,'MyCode\eddy_',eddy_id,'.nc']);
timeo = ncread(filename,'eddy/time')+datenum(1950,01,01); %Initial time is also a free parameter
eddy_time = timeo-min(timeo);

clearvars eddyPath
eddyPath.lat = ncread(filename,'eddy/latitude');
eddyPath.lon = ncread(filename,'eddy/longitude');

%if you want to change the eddyPath to a function handle
eddyPath_fun_t.lon = @(t) interp1(eddy_time, eddyPath.lon, t, 'linear', 'extrap');
eddyPath_fun_t.lat = @(t) interp1(eddy_time, eddyPath.lat, t, 'linear', 'extrap');
% Function to convert any lat/lon to x/y relative to the eddy center at time t
get_eddyPath_xy = @(lat, lon, t) latlon2xy(lat, lon, eddyPath_fun_t.lat(t),eddyPath_fun_t.lon(t));

%if you want to use eddyPath_fun_t.xe,ye as a function handle
%alongtrack are already adjusted to be eddy centered so xo,yo are zeros
eddyPath_fun_t.xe=@(t) 0*t;
eddyPath_fun_t.ye=@(t) 0*t;
%% Load alongtrack
JasonAlongTrack.filename = strcat(readdir, 'JasonAlongTrack.nc');
JasonAlongTrack.lat = ncread(JasonAlongTrack.filename, 'lat');
JasonAlongTrack.lon = ncread(JasonAlongTrack.filename, 'lon');
%JML convert time to Matlab's datenum format
JasonAlongTrack.time = ncread(JasonAlongTrack.filename, 'time') + datenum(1950, 1, 1);
%Time is defined as beginning at 4:05 AM on Sept 23, 1992,
JasonAlongTrack.ssh = ncread(JasonAlongTrack.filename, 'sla');
alongtrackLatLon = extractAlongtrackLatLonEddyCenter(JasonAlongTrack,eddyPath,timeo,radius=300); %radius in km

% % if loading jasonalongtrack above takes too long, load the mat file 
% current_file_path = matlab.desktop.editor.getActiveFilename;
% [current_dir, ~, ~] = fileparts(current_file_path);
% load(strcat(current_dir,'\alongtrackLatLon_real.mat'))

[x_km, y_km]=get_eddyPath_xy(alongtrackLatLon.lat,alongtrackLatLon.lon,alongtrackLatLon.t-min(alongtrackLatLon.t));
alongtrackXY.x=x_km*1e3;
alongtrackXY.y=y_km*1e3;

%% Load mapped field
filename = 'E:\Research\myCode\aviso_madt.nc';
ni=ncinfo(filename);
mapped_field.time = ncread(filename,'time')+datenum(1950,01,01);
mapped_field.lat = ncread(filename,'latitude'); %ascending
mapped_field.lon = deg180(ncread(filename,'longitude')); % ascending

radius_km=300;%km
bin_size=12.5*1e3;
max_r=(radius_km*1e3/bin_size)*bin_size;
xbin=-max_r:bin_size:max_r;
xmid=(xbin+vshift(xbin,1,1))./2;
xmid=xmid(1:end-1);
ymid=xmid;
[XGrid, YGrid] = ndgrid(xmid, ymid);
clearvars ssh_interp

for n = 1:length(timeo)
    % Find time index closest to current eddy time
    [~, time_idx] = min(abs(mapped_field.time - timeo(n)));

    % Find indices for region around current eddy position
    [~, center_lat] = min(abs(mapped_field.lat - eddyPath.lat(n)));
    [~, center_lon] = min(abs(mapped_field.lon - eddyPath.lon(n)));

    % Calculate the lat-lon box that corresponds to our desired km radius
    [lat_box, lon_box] = xy2latlon([-radius_km, radius_km], [-radius_km, radius_km], eddyPath.lat(n), eddyPath.lon(n));
     
    % Find indices in the grid that match these boundaries
    [~, start_lat] = min(abs(mapped_field.lat - lat_box(1)));
    [~, last_lat] = min(abs(mapped_field.lat - lat_box(2)));
    [~, start_lon] = min(abs(mapped_field.lon - lon_box(1)));
    [~, last_lon] = min(abs(mapped_field.lon - lon_box(2)));
    
    % Handle longitude wrap-around if needed
    if start_lon > last_lon % Crossing the longitude boundary
        % First part: end of longitude array to the end
        start1 = [start_lon, start_lat, time_idx];
        count1 = [length(mapped_field.lon) - start_lon + 1, last_lat - start_lat + 1, 1];
        ssh_part1 = ncread(filename, 'sla', start1, count1); %m
        lon_part1 = mapped_field.lon(start_lon:end);
        
        % Second part: beginning of longitude array to last_lon
        start2 = [1, start_lat, time_idx];
        count2 = [last_lon, last_lat - start_lat + 1, 1];
        ssh_part2 = ncread(filename, 'sla', start2, count2);%m
        lon_part2 = mapped_field.lon(1:last_lon);
        
        % Combine the two parts
        current_ssh = [ssh_part1; ssh_part2];
        current_lon = [lon_part1; lon_part2];
    else
        % Normal case - extraction window within longitude bounds
        start = [start_lon, start_lat, time_idx];
        count = [last_lon - start_lon + 1, last_lat - start_lat + 1, 1];
        current_ssh = ncread(filename, 'sla', start, count);%m
        current_lon = mapped_field.lon(start_lon:last_lon);
    end
    
    current_lat = mapped_field.lat(start_lat:last_lat); %lat doesn't have a if statement bc we don't do eddies that cross the equator
    % Create meshgrid for extracted region
    [lon_grid, lat_grid] = ndgrid(current_lon, current_lat);
    [x_km, y_km]=get_eddyPath_xy(lat_grid,lon_grid,timeo(n)-min(timeo));

    % Interpolate onto a grid
    ssh_interp(:,:,n) = griddata(double(x_km * 1e3), double(y_km * 1e3), current_ssh, XGrid, YGrid, 'linear');
end
    % Store mapped grid data instead of an array
    mapped_extracted.ssh = ssh_interp;
    mapped_extracted.x = xmid;
    mapped_extracted.y = ymid;
    mapped_extracted.t = timeo;%permute(repmat(timeo,[1,length(xmid),length(xmid)]),[2,3,1]);
%%
time_window=1:length(timeo)-1;
window_idx = find(alongtrackLatLon.t-min(alongtrackLatLon.t) >= time_window(1) & alongtrackLatLon.t-min(alongtrackLatLon.t) <= time_window(end));
alongtrack.t=alongtrackLatLon.t(window_idx);
alongtrack.x=alongtrackXY.x(window_idx);
alongtrack.y=alongtrackXY.y(window_idx);
alongtrack.ssh=alongtrackLatLon.ssh(window_idx);

mapped_window.t=mapped_extracted.t(time_window);
mapped_window.x=mapped_extracted.x;
mapped_window.y=mapped_extracted.y;
mapped_window.ssh=mapped_extracted.ssh(:,:,time_window);

%%
%Full
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(alongtrack,eddyPath_fun_t,showplot=0);
data{1} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

%OSSE
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(mapped_window,eddyPath_fun_t,showplot=0);
data{2} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

row=2;
%%
figname='real_eddy_radtime';
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
clim_cbar1 = [-10,40];
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
clim_cbar2 = [0,15];
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
    axes('Units','centimeters','Position', [col_starts(2), row_starts(i), subplot_width, subplot_height]);
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
createfigure_radialProfileTime_real(ZData1, YData1, XData1, CData1, CData2, ZData2, YData2, XData2, CData3, CData4)
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