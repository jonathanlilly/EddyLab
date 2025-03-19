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

%% Test data
% Test 1 - Gaussian Eddy
clearvars params
params.A = 0.15; %meter
params.L = 80e3; %meter

eddy_model = analyticalEddyModel(eddyPath_fun_t,params);
alongtrack = analyticalEddyModelOSSE(alongtrackLatLon,eddy_model);
[mz, xmid, ymid, numz, stdz] = composite2D(alongtrack, eddyPath_fun_t,showplot=0);

%% 2D track counts histogram
figure(main_fig)
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
%% Save
folder_name='D:\UW\AlongTrack-GRL\fig';
figname = 'composite_hist';
set(gcf,'Renderer','painter');         %Instead of painter,opengl
set(gcf,'inverthardcopy','off') %to save white text as white
set(1,'paperpositionmode','auto')
print('-dpng','-r600',strcat(folder_name,'\',figname,'.png'))
print('-depsc','-r600',strcat(folder_name,'\',figname,'.eps'))