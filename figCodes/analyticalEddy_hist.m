% print histogram from analytical eddy
main_fig = figure('color','w');
set(gcf,'Units','centimeters');  %Ensure all dimensions are in cm
% set(gcf,'Position',[0 0 18 16]);  %Size of the figure. The first two numbers indicate the location of the lower left-hand corner.
                                    %The second two numbers indicate size
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
%% Test data
% Test 1 - Gaussian Eddy
clearvars params
params.A = 0.15; %meter
params.L = 80e3; %meter

eddy_model = analyticalEddyModel(eddyPath_fun_t,params);
alongtrack.ssh = eddy_model(alongtrackXY.x,alongtrackXY.y,alongtrackXY.t-min(alongtrackXY.t));
[mz, xmid, ymid, numz, stdz] = composite2D(alongtrack, eddyPath_fun_t,showplot=0);
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(alongtrack,eddyPath_fun_t,showplot=0);
%% 2D track counts histogram
%% subplots
main_fig = figure('color','w');
set(gcf,'Units','centimeters');  %Ensure all dimensions are in cm
% set(gcf,'Position',[0 0 13 16]);  %Size of the figure. The first two numbers indicate the location of the lower left-hand corner.
subplot(1,2,1)
numz(numz == 0) = nan;
jpcolor(xmid/1e3, ymid/1e3, numz); %
% clim(round([min(numz(:)), max(numz(:))]))
shading flat
axis equal
xlabel('Distance East (km)', 'FontName', 'times')
ylabel('Distance North (km)', 'FontName', 'times')
set(gca, 'fontname', 'times','FontSize',16)
if exist('brewermap')
    colormap(brewermap([], '-Spectral'))
end
c = colorbar('EastOutside');
c.Label.String = 'Histogram (counts)';
c.Label.FontSize=16;
c.Label.FontName='times';
xlim([-250, 250]), ylim([-250, 250]);

text(-230, 210, '(a)', 'fontsize', 14, 'fontname', 'times','color','w');

%%
% Rad-Time Hist
subplot(1,2,2)
% subplot(3, 3, 3*(i-1)+2);
hold on;
jpcolor(rmid_rt/1e3, [1:length(tmid_rt)], numz_rt) %mz_rt'./vmean(mz_rt,2)'
shading flat;
% axis equal;
ylabel('Time (Cycle)', 'FontName', 'times')
xlabel('Radial distance (km)', 'FontName', 'times')
% ylabel('Time (Cycle)', 'FontName', 'times')
set(gca, 'fontname', 'times', 'fontsize', 16);
colormap(brewermap([], '-Spectral'));
colorbar('EastOutside');
% clim([1,30])
xlim([0, 250]);ylim([1,length(tmid_rt)]-0.5)
text(10, 31.5, '(b)', 'fontsize', 14, 'fontname', 'times','color','w');

%%
% Get all subplot axes
ax = findobj(gcf, 'Type', 'axes');

% Select the subplot (e.g., first subplot)
ax2 = ax(1);  % Adjust index based on your layout
ax1 = ax(2);
% Find the surface object within the selected subplot
hSurf = findobj(ax1, 'Type', 'surface');
hSurf2 = findobj(ax2, 'Type', 'surface');


% Extract the required data
ZData1 = get(hSurf, 'ZData');ZData2 = get(hSurf2, 'ZData');
YData1 = get(hSurf, 'YData');XData2 = get(hSurf2, 'XData');
YData2 = get(hSurf2, 'YData');
CData1 = get(hSurf, 'CData');  % First color data
CData2 = get(hSurf2, 'CData');
%%
createfigure_hist(ZData1, YData1, CData1, ZData2, YData2, XData2, CData2)
%% Save
folder_name='E:\Research\AlongTrack-GRL\fig';%'D:\UW\AlongTrack-GRL\fig';
figname = 'composite_hist';
set(gcf,'Renderer','opengl');         %Instead of painter,opengl
set(gca, 'Color', 'w'); % Sets axis background to white
set(gcf,'inverthardcopy','off') %to save white text as white
set(1,'paperpositionmode','auto')
print('-dpng','-r600',strcat(folder_name,'\',figname,'.png'))
print('-depsc','-r600',strcat(folder_name,'\',figname,'.eps'))