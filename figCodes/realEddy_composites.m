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
time_window=100:300;%1:length(timeo)-1;
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
%Alongtrack
[mz, xmid, ymid, ~, stdz] = composite2D(alongtrack, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(alongtrack, eddyPath_fun_t,showplot=0);
data{1} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

%Mapped
[mz, xmid, ymid, ~, stdz] = composite2D(mapped_window, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(mapped_window, eddyPath_fun_t,showplot=0);
data{2} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

row=2;
%%
figname='real_eddy_composite';
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
clim_cbar1=[-10, 30];
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
clim_cbar2 = [0,10];%round([max([min(OSSE_data{1}{4},[],'all'),min(fullfield_data{1}{4},[],'all')]), max([max(OSSE_data{1}{4},[],'all'),max(fullfield_data{1}{4},[],'all')])]*1e2,1);
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
    
    xlim([0, 250]);ylim([-5,35])
    text(-60, max(get(gca, 'YLim')), labels{3*(i-1)+3}, 'fontsize', 14, 'fontname', 'times');
end

%%
% Save
folder_name='E:\Research\AlongTrack-GRL\fig';%'D:\UW\AlongTrack-GRL\fig';
set(gcf,'Renderer','opengl');         % somehow painters or completely vectorize creates white lines
set(gca, 'Color', 'w'); % Sets axis background to white
set(gcf,'inverthardcopy','off') %to save white text as white
set(1,'paperpositionmode','auto')
print('-dpng','-r1200',strcat(folder_name,'\',figname,'.png'))
print('-depsc','-r1200',strcat(folder_name,'\',figname,'.eps'))
% exportgraphics(gcf, fullfile(folder_name, strcat(figname, '.eps')), 'ContentType', 'vector');