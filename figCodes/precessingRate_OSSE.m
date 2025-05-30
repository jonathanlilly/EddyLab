%% precessingRate_OSSE 
%% Extract tracks along the defined eddy path
% define eddy path: xe(t),ye(t) function handles as a function of time
x0 = 500e3; y0 = 0e3; vx = -2.0e3; vy = -0.5e3;
eddyPath_fun_t.xe = @(t) x0+vx*t;
eddyPath_fun_t.ye = @(t) y0+vy*t;

% define domain
x = linspace(-1e6,1e6,200); %meter
y = linspace(-500e3,500e3,200); %meter
totalDays = 366;

% extract track - takes somes time bc it loads the entire track matrix
% my_readdir = 'G:\My Drive\AlongTrack\';
% alongtrackLatLon = alongtrackFromXYDomain(x,y,totalDays,lono=-40,readdir=my_readdir); %options: lato=24, lono=308

% If alongtrackFromXYDomain takes too long or you don't have
% JasonAlongTrack.nc file,
% To ensure that you're getting the correct directory:
current_file_path = matlab.desktop.editor.getActiveFilename;
[file_dir, ~, ~] = fileparts(current_file_path);
% Move one level up from the current directory
current_dir = fileparts(file_dir);
load(strcat(current_dir,'\alongtrackFitting\alongtrackLatLon.mat'))

% convert lat,lon to x,y to generate an eddy
alongtrackXY = latlon2xy_centered(alongtrackLatLon);
alongtrack = alongtrackXY; %pass on the ground track lat,lon,t arrays

%% loop precessing rate with Test 3 case
thetaDots= [-10*pi/365,-5*pi/365,-pi/365,-0.5*pi/365];%linspace(-10*pi/365,-0.5*pi/365,4);%-10*pi/365; %similar to QG model ~5 rotations per year. 0.01-0.2 deg per day 
% Test 3 - Precessing Elliptical Eddy:
clearvars params
params.A = 0.15;
L = 80e3;
lambda = 0.5;
La=L/sqrt(lambda);
params.La = L/sqrt(lambda);%0.4*2*L;
params.Lb = lambda*params.La;%0.2*2*L;

OSSE_data = cell(length(thetaDots), 1);
for i=1:length(thetaDots)
params.thetaDot= thetaDots(i);
eddy_model = analyticalEddyModel(eddyPath_fun_t,params);
alongtrack.ssh = eddy_model(alongtrackXY.x,alongtrackXY.y,alongtrackXY.t-min(alongtrackXY.t));
[mz, xmid, ymid, ~, stdz] = composite2D(alongtrack, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(alongtrack, eddyPath_fun_t,showplot=0);
OSSE_data{i} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};
end
%%
figname='test3_OSSE_thetaDot';
plot_data = OSSE_data;

%% subplots in multiple columns
main_fig = figure('color','w');
set(gcf,'Units','centimeters'); %Ensure all dimensions are in cm
set(gcf,'Position',[0 0 15.74 6]); %Size of the figure. Modified to be wider for single-row layout

% Subplot labels (a) through (f)
labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
figure(main_fig); % Make sure we're on the main figure before creating subplots
clf(main_fig)

% Create axes for the colorbar at the top
cbar_pos = [0.4, 0.88, 0.2, 0.02]; % Position for the colorbar

% Multiple column layout
num_rows = 1; % 1 row
num_cols = length(thetaDots); % 6 columns

% Compute subplot dimensions and margins dynamically
total_margin_fraction = 0.2; % Adjust based on figure size
subplot_width = (1 - total_margin_fraction) / num_cols;  
h_margin = (1 - num_cols * subplot_width) / (num_cols + 1); 
% Compute column positions dynamically
col_starts = h_margin + (0:num_cols-1) * (subplot_width + h_margin);

subplot_height = 0.5; % Height of subplot (increased for single row)
v_margin = 0.15; % Fixed positive vertical margin

% Create the colorbar
cbar_ax = axes('Position', cbar_pos);
colormap(cbar_ax, brewermap([], '-Spectral'));
cbar = colorbar('peer', cbar_ax, 'Location', 'North');
cbar.Position(4) = 0.015;
cbar.Position(2) = cbar.Position(2) - 0.01;
cbar.FontName = 'times';
cbar.FontSize = 12;
clim_cbar = round([0.9,max(OSSE_data{1}{4},[],'all')*1e2],1);
caxis(cbar_ax, clim_cbar);
axis(cbar_ax, 'off');
title(cbar_ax, 'Temporal SSH variance (cm)', 'fontsize', 12, 'fontname', 'times');

% Plot all subplots in grid layout
for i = 1:length(thetaDots)
    % Calculate row and column indices
    row = ceil(i/num_cols);
    col = mod(i-1, num_cols) + 1;
    
    % Calculate position for this subplot
    x_pos = col_starts(i);
    y_pos = 1 - (row * (subplot_height + v_margin));
    
    % Unpack data for this test - we only need stdz for the plots
    xmid = plot_data{i}{2};
    ymid = plot_data{i}{3};
    stdz = plot_data{i}{4};
    
    % Plot stdz (temporal std) in grid layout
    subplot('Position', [x_pos, y_pos, subplot_width, subplot_height]);
    hold on;
    jpcolor(xmid/1e3, ymid/1e3, stdz*1e2);
    shading flat;
    axis equal;
    
    % Only show y-axis label for the leftmost column
    if (col == 1)
        ylabel('Distance North (km)', 'FontName', 'times');
    else
        set(gca,'YTickLabel',[])
    end
    
    % Show x-axis labels for all plots in the single row
    xlabel('Distance East (km)', 'FontName', 'times');
    
    set(gca, 'fontname', 'times', 'fontsize', 12);
    colormap(brewermap([], '-Spectral'));
    clim(clim_cbar)
    xlim([-250, 250]);
    ylim([-250, 250]);
    
    % Add subplot label in top-left corner
    text(-230, 210, labels{i}, 'fontsize', 14, 'fontname', 'times', 'color', 'w');
end

% % Adjust figure to ensure proper spacing
% set(gcf,'Units','centimeters');
% set(gcf,'Position',[0 0 16 12]);
%%
% Get all subplot axes
ax = findobj(gcf, 'Type', 'axes');

% Select the subplot (e.g., first subplot)
ax1 = ax(1);  % Adjust index based on your layout
ax2 = ax(2);
ax3 = ax(3);
ax4 = ax(4); 
% Find the surface object within the selected subplot
hSurf = findobj(ax1, 'Type', 'surface');
hSurf2 = findobj(ax2, 'Type', 'surface');
hSurf3 = findobj(ax3, 'Type', 'surface');
hSurf4 = findobj(ax4, 'Type', 'surface');

% Extract the required data
ZData1 = get(hSurf, 'ZData');
YData1 = get(hSurf, 'YData');
CData1 = get(hSurf, 'CData');  % First color data
CData2 = get(hSurf2, 'CData');
CData3 = get(hSurf3, 'CData');
CData4 = get(hSurf4, 'CData');

%%
createfigure_precessingRate_OSSE(ZData1, YData1, CData4, CData3, CData2, CData1);
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