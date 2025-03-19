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
[mz, xmid, ymid, ~, stdz] = composite2D(alongtrack, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(alongtrack, eddyPath_fun_t,showplot=0);
OSSE_data{1} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

fullfield.ssh = eddy_model(fullfield.x,fullfield.y,fullfield.t);
[mz, xmid, ymid, ~, stdz] = composite2D(fullfield, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(fullfield, eddyPath_fun_t,showplot=0);
fullfield_data{1} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

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
[mz, xmid, ymid, ~, stdz] = composite2D(alongtrack, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(alongtrack, eddyPath_fun_t,showplot=0);
OSSE_data{2} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

fullfield.ssh = eddy_model(fullfield.x,fullfield.y,fullfield.t);
[mz, xmid, ymid, ~, stdz] = composite2D(fullfield, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(fullfield, eddyPath_fun_t,showplot=0);
fullfield_data{2} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

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
[mz, xmid, ymid, ~, stdz] = composite2D(alongtrack, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(alongtrack, eddyPath_fun_t,showplot=0);
OSSE_data{3} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};

fullfield.ssh = eddy_model(fullfield.x,fullfield.y,fullfield.t);
[mz, xmid, ymid, ~, stdz] = composite2D(fullfield, eddyPath_fun_t,showplot=0);
[mzxy, rmid, ~, stdz_rt] = radialProfile(fullfield, eddyPath_fun_t,showplot=0);
fullfield_data{3} = {mz, xmid, ymid, stdz, mzxy, rmid, stdz_rt};
%%
figname='analytical_eddies_full';
plot_data = fullfield_data;
% 
% figname='analytical_eddies_OSSE';
% plot_data = OSSE_data;
%% subplots
main_fig = figure('color','w');
set(gcf,'Units','centimeters');  %Ensure all dimensions are in cm
set(gcf,'Position',[0 0 18 16]);  %Size of the figure. The first two numbers indicate the location of the lower left-hand corner.
                                    %The second two numbers indicate size

% Subplot labels (a) through (i)
labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'};
figure(main_fig); % Make sure we're on the main figure before creating subplots
clf(main_fig)
% Create axes for the colorbars at the top
cbar1_pos = [0.1, 0.93, 0.2, 0.02]; % Position for the first colorbar
cbar2_pos = [0.38, 0.93, 0.2, 0.02];  % Position for the second colorbar
legend_pos = [0.71, 0.93, 0.2, 0.02]; % Position for the legend at the top

% First colorbar for SSH
cbar1_ax = axes('Position', cbar1_pos);
colormap(cbar1_ax, brewermap([], '-Spectral'));
cbar1 = colorbar('peer', cbar1_ax, 'Location', 'North');
cbar1.FontName='times';
cbar1.FontSize=12;
cbar1.Position(1)=0.115;
cbar1.Position(2)=cbar1.Position(2)-0.04;
cbar1.Position(4)=0.015;
clim_cbar1=[-2.5, 15];
caxis(cbar1_ax, clim_cbar1);
axis(cbar1_ax, 'off');
title(cbar1_ax, 'Time-averaged SSH (cm)', 'fontsize', 12, 'fontname', 'times');

% Second colorbar for variance
cbar2_ax = axes('Position', cbar2_pos);
colormap(cbar2_ax, brewermap([], '-Spectral'));
cbar2 = colorbar('peer', cbar2_ax, 'Location', 'North');
cbar2.FontName='times';
cbar2.FontSize=12;
cbar2.Position(1)=0.387;
cbar2.Position(2)=cbar2.Position(2)-0.04;
cbar2.Position(4)=0.015;
clim_cbar2=[0, 2];
caxis(cbar2_ax, clim_cbar2);
axis(cbar2_ax, 'off');
title(cbar2_ax, 'Temporal variance (cm)', 'fontsize', 12, 'fontname', 'times');

% Third legend for variance
% Create an invisible axes for the legend at the top
legend_ax = axes('Position', legend_pos);

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
 legend_top.Position=[0.64,0.88,0.32,0.055];
title(legend_ax, 'SSH variance (cm)', 'fontsize', 12, 'fontname', 'times');
axis(legend_ax, 'off');

% Create tighter subplots
% Define new positions for the subplots
subplot_width = 0.22;
subplot_height = 0.25;
col1_start = 0.1;
col2_start = 0.37;
col3_start = 0.69;
row_starts = 0.08+(subplot_height+0.02).*[2,1,0];

for i = 1:3
    % Unpack data for this test
    mz = plot_data{i}{1};
    xmid = plot_data{i}{2};
    ymid = plot_data{i}{3};
    stdz = plot_data{i}{4};
    mzxy = plot_data{i}{5};
    rmid = plot_data{i}{6};
    stdz_rt = plot_data{i}{7};
    
    % Column 1: Plot mz (time-averaged SSH)
    subplot('Position', [col1_start, row_starts(i), subplot_width, subplot_height]);
    %subplot(3, 3, 3*(i-1)+1);
    hold on;
    jpcolor(xmid/1e3, ymid/1e3, mz*1e2);
    shading flat;
    axis equal;
    if(i==1||i==2)
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
    text(-230, 210, labels{3*(i-1)+1}, 'fontsize', 14, 'fontname', 'times','color','w');
    
    % Column 2: Plot stdz (temporal std)
    subplot('Position', [col2_start, row_starts(i), subplot_width, subplot_height]);
    % subplot(3, 3, 3*(i-1)+2);
    hold on;
    jpcolor(xmid/1e3, ymid/1e3, stdz*1e2);
    shading flat;
    axis equal;
    if(i==1||i==2)
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
    text(-230, 210, labels{3*(i-1)+2}, 'fontsize', 14, 'fontname', 'times','color','w');

    % Column 3: Plot radial profiles
    subplot('Position', [col3_start, row_starts(i), subplot_width+0.02, subplot_height-0.02]);
    % subplot(3, 3, 3*(i-1)+3);
    hold on;
    
    % Plot the radial profiles (assuming avgAziStdTemp, etc. are available)
    % Note: You might need to adjust this part depending on how these variables are defined in your script
    use stdz_rt
    % Apply line styles
    h = plot(rmid/1e3, [mzxy,  stdTotalxy, avgAziStdTemp, stdAziAvgTemp, avgTempStdAzi, stdTempAvgAzi]*1e2);
    linestyle 2k 2W 2T-- 2U-- 2V-. 2X:
    hlines(0, 'k:')

    if(i==1||i==2)
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
    
    xlim([0, 250]);ylim([-2.5,15])
    text(-70, max(get(gca, 'YLim')), labels{3*(i-1)+3}, 'fontsize', 14, 'fontname', 'times');
end

%%
% Save
folder_name='D:\UW\AlongTrack-GRL\fig';
set(gcf,'Renderer','painter');         %Instead of painter,opengl
set(gcf,'inverthardcopy','off') %to save white text as white
set(1,'paperpositionmode','auto')
print('-dpng','-r1200',strcat(folder_name,'\',figname,'.png'))
print('-depsc','-r1200',strcat(folder_name,'\',figname,'.eps'))