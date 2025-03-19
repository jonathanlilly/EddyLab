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

% full field
fullfield.x = reshape(repmat(x',[1,length(y),totalDays]),[],1);
fullfield.y = reshape(repmat(y,[length(x),1,totalDays]),[],1);
fullfield.t = reshape(permute(repmat([1:totalDays]'-1,[1,length(x),length(y)]),[2,3,1]),[],1);

fullfield_data = cell(3, 1);
OSSE_data = cell(3, 1);
%% Test data
% Test 1 - Gaussian Eddy
clearvars params
params.A = 0.15; %meter
params.L = 80e3; %meter

eddy_model = analyticalEddyModel(eddyPath_fun_t,params);
alongtrack = analyticalEddyModelOSSE(alongtrackLatLon,eddy_model);
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(alongtrack,eddyPath_fun_t,showplot=0);
OSSE_data{1} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

fullfield.ssh = eddy_model(fullfield.x,fullfield.y,fullfield.t);
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(fullfield,eddyPath_fun_t,showplot=0);
fullfield_data{1} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

% Test 2 - Steady Elliptical Eddy:
clearvars params
params.A = 0.15;
L = 80e3;
lambda = 0.5;
La=L/sqrt(lambda);
params.La = L/sqrt(lambda);%0.4*2*L;
params.Lb = lambda*params.La;%0.2*2*L;

eddy_model = analyticalEddyModel(eddyPath_fun_t,params);
alongtrack = analyticalEddyModelOSSE(alongtrackLatLon,eddy_model);
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(alongtrack,eddyPath_fun_t,showplot=0);
OSSE_data{2} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

fullfield.ssh = eddy_model(fullfield.x,fullfield.y,fullfield.t);
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(fullfield,eddyPath_fun_t,showplot=0);
fullfield_data{2} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

% Test 3 - Precessing Elliptical Eddy:
clearvars params
params.A = 0.15;
L = 80e3;
lambda = 0.5;
La=L/sqrt(lambda);
params.La = L/sqrt(lambda);%0.4*2*L;
params.Lb = lambda*params.La;%0.2*2*L;
params.thetaDot= -10*pi/365; %similar to QG model ~5 rotations per year
eddy_model = analyticalEddyModel(eddyPath_fun_t,params);
alongtrack = analyticalEddyModelOSSE(alongtrackLatLon,eddy_model);
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(alongtrack,eddyPath_fun_t,showplot=0);
OSSE_data{3} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};

fullfield.ssh = eddy_model(fullfield.x,fullfield.y,fullfield.t);
[mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt] = radialProfileTime(fullfield,eddyPath_fun_t,showplot=0);
fullfield_data{3} = {mz_rt, rmid_rt, tmid_rt, numz_rt, stdz_rt};
%%
% figname='analytical_radtime_full';
% plot_data = fullfield_data;
% 
figname='analytical_radtime_OSSE';
plot_data = OSSE_data;
%% subplots
main_fig = figure('color','w');
set(gcf,'Units','centimeters');  %Ensure all dimensions are in cm
set(gcf,'Position',[0 0 13 16]);  %Size of the figure. The first two numbers indicate the location of the lower left-hand corner.
                                    %The second two numbers indicate size
num_cols=2;
% Subplot labels (a) through (i)
labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'};
figure(main_fig); % Make sure we're on the main figure before creating subplots
clf(main_fig)
% % Create axes for the colorbars at the top
% cbar1_pos = [0.1, 0.93, 0.2, 0.02]; % Position for the first colorbar
% cbar2_pos = [0.38, 0.93, 0.2, 0.02];  % Position for the second colorbar
% legend_pos = [0.71, 0.93, 0.2, 0.02]; % Position for the legend at the top
% cbar3_pos=legend_pos;
% % First colorbar for SSH
% cbar1_ax = axes('Position', cbar1_pos);
% colormap(cbar1_ax, brewermap([], '-Spectral'));
% cbar1 = colorbar('peer', cbar1_ax, 'Location', 'North');
% cbar1.FontName='times';
% cbar1.FontSize=12;
% cbar1.Position(1)=0.115;
% cbar1.Position(2)=cbar1.Position(2)-0.04;
% cbar1.Position(4)=0.015;
% clim_cbar1=[-2.5, 15];
% caxis(cbar1_ax, clim_cbar1);
% axis(cbar1_ax, 'off');
% title(cbar1_ax, 'Time-averaged SSH (cm)', 'fontsize', 12, 'fontname', 'times');
% 
% % Second colorbar for variance
% cbar2_ax = axes('Position', cbar2_pos);
% colormap(cbar2_ax, brewermap([], '-Spectral'));
% cbar2 = colorbar('peer', cbar2_ax, 'Location', 'North');
% cbar2.FontName='times';
% cbar2.FontSize=12;
% cbar2.Position(1)=0.387;
% cbar2.Position(2)=cbar2.Position(2)-0.04;
% cbar2.Position(4)=0.015;
% % clim_cbar2=[0.5,2.5];
% clim_cbar2=round([max(fullfield_data{1}{5},[],'all'),max(fullfield_data{3}{5},[],'all')]*1e2,1);
% caxis(cbar2_ax, clim_cbar2);
% % set(cbar2, 'Ticks', [clim_cbar2(1), round((clim_cbar2(1)+clim_cbar2(2))/2,1),clim_cbar2(2)]);
% axis(cbar2_ax, 'off');
% title(cbar2_ax, 'Azimuthal variance (cm)', 'fontsize', 12, 'fontname', 'times');
% 
% % % Third legend for variance
% % cbar3_ax = axes('Position', cbar3_pos);
% % colormap(cbar3_ax, brewermap([], '-Spectral'));
% % cbar3 = colorbar('peer', cbar3_ax, 'Location', 'North');
% % cbar3.FontName='times';
% % cbar3.FontSize=12;
% % cbar3.Position(1)=0.115;
% % cbar3.Position(2)=cbar3.Position(2)-0.04;
% % cbar3.Position(4)=0.015;
% % clim_cbar3=[-2.5, 15];
% % caxis(cbar3_ax, clim_cbar3);
% % axis(cbar3_ax, 'off');
% % title(cbar3_ax, 'Histogram (counts)', 'fontsize', 12, 'fontname', 'times');

% Create tighter subplots
% Define new positions for the subplots
% subplot_width = 0.22;
% subplot_height = 0.25;
% col1_start = 0.1;
% col2_start = 0.37;
% col3_start = 0.69;
% row_starts = 0.08+(subplot_height+0.02).*[2,1,0];

% Compute center positions of each column for colorbars
cbar_width = subplot_width * 0.9;  % Slightly smaller than subplot width
cbar_height = 0.02;  % Fixed small height
cbar_y = 0.93;  % Fixed height near the top

cbar_positions = [col_starts + subplot_width / 2 - cbar_width / 2; 
                  repmat(cbar_y, 1, num_cols);  
                  repmat(cbar_width, 1, num_cols);
                  repmat(cbar_height, 1, num_cols)];

% Create first colorbar
cbar1_ax = axes('Position', cbar_positions(:,1)');
colormap(cbar1_ax, brewermap([], '-Spectral'));
cbar1 = colorbar('peer', cbar1_ax, 'Location', 'North');
cbar1.Position(4)=0.015;
cbar1.Position(2)=cbar1.Position(2)-0.04;
cbar1.FontName = 'times';
cbar1.FontSize = 12;
caxis(cbar1_ax, [-2.5, 15]);
axis(cbar1_ax, 'off');
title(cbar1_ax, 'Time-averaged SSH (cm)', 'fontsize', 12, 'fontname', 'times');

% Create second colorbar
cbar2_ax = axes('Position', cbar_positions(:,2)');
colormap(cbar2_ax, brewermap([], '-Spectral'));
cbar2 = colorbar('peer', cbar2_ax, 'Location', 'North');
cbar2.Position(4)=0.015;
cbar2.Position(2)=cbar2.Position(2)-0.04;
cbar2.FontName = 'times';
cbar2.FontSize = 12;
clim_cbar2 = round([max([max(OSSE_data{1}{5},[],'all'),max(fullfield_data{1}{5},[],'all')]), max([max(OSSE_data{3}{5},[],'all'),max(fullfield_data{3}{5},[],'all')])]*1e2,1);
caxis(cbar2_ax, clim_cbar2);
axis(cbar2_ax, 'off');
title(cbar2_ax, 'Azimuthal variance (cm)', 'fontsize', 12, 'fontname', 'times');

% Set legend position to match third column (or last column if fewer than 3)
% legend_pos = cbar_positions(:, min(num_cols, 3))'; 
% cbar3_ax = axes('Position', legend_pos);
% axis(cbar3_ax, 'off');  % Replace with actual legend command if needed
% title(cbar3_ax, 'Legend', 'fontsize', 12, 'fontname', 'times');

% Compute subplot width dynamically
total_margin_fraction = 0.225;  % Total margin space as fraction of figure width
subplot_width = (1 - total_margin_fraction) / num_cols;  % Adjust width based on column count
subplot_height = 0.25;  % Keep height same

% Compute dynamic margins
margin = (1 - num_cols * subplot_width) / (num_cols + 1);

% Compute column positions dynamically
col_starts = 0.03+ margin + (0:num_cols-1) * (subplot_width + margin);
% col3_start = col2_start + subplot_width + margin;
% row_starts = 0.08+(subplot_height+0.02).*[2,1,0];


for i = 1:3
    % Unpack data for this test
    mz_rt = plot_data{i}{1};
    rmid_rt = plot_data{i}{2};
    tmid_rt = plot_data{i}{3};
    numz_rt = plot_data{i}{4};
    stdz_rt = plot_data{i}{5};
    
    % Column 1: Plot mz (time-averaged SSH)
    subplot('Position', [col_starts(1), row_starts(i), subplot_width, subplot_height]);
    %subplot(3, 3, 3*(i-1)+1);
    hold on;
    jpcolor(rmid_rt/1e3, [1:length(tmid_rt)], mz_rt*1e2)%mz_rt'./vmean(mz_rt,2)'
    shading flat;
    % axis equal;
    if(i==1||i==2)
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
    text(10, 32, labels{3*(i-1)+1}, 'fontsize', 14, 'fontname', 'times','color','w');
    
    % Column 2: Plot stdz (temporal std)
    subplot('Position', [col_starts(2), row_starts(i), subplot_width, subplot_height]);
    % subplot(3, 3, 3*(i-1)+2);
    hold on;
    jpcolor(rmid_rt/1e3, [1:length(tmid_rt)], stdz_rt*1e2) %mz_rt'./vmean(mz_rt,2)'
    shading flat;
    % axis equal;
    if(i==1||i==2)
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
    text(10, 32, labels{3*(i-1)+2}, 'fontsize', 14, 'fontname', 'times','color','w');
    
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

%%
% Save
folder_name='D:\UW\AlongTrack-GRL\fig';
set(gcf,'Renderer','painter');         %Instead of painter,opengl
set(gcf,'inverthardcopy','off') %to save white text as white
set(1,'paperpositionmode','auto')
print('-dpng','-r1200',strcat(folder_name,'\',figname,'.png'))
print('-depsc','-r1200',strcat(folder_name,'\',figname,'.eps'))