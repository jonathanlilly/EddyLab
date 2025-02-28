%% Test 1 - Gaussian Eddy
% define eddy shape function handle
clearvars params
params.A = 0.15; %meter
params.L = 80e3; %meter

% eddy_model is a function handle with a chosen set of parameters (x,y,t)
% eddyShapeString = 'Gaussian';
eddy_model = analyticalEddyModel(eddyPath_fun_t,params);

% apply OSSE on an eddy_model function (x,y,t)
alongtrack = analyticalEddyModelOSSE(alongtrackLatLon,eddy_model);


%% Test 2 - Steady Elliptical Eddy:
clearvars params
params.A = 0.15;
L = 80e3;
params.La = 0.4*2*L;
params.Lb = 0.2*2*L;

% eddy_model is a function handle with a chosen set of parameters (x,y,t)
% eddyShapeString = 'Ellipse';
eddy_model = analyticalEddyModel(eddyPath_fun_t,params);

% apply OSSE on an eddy_model function (x,y,t)
alongtrack = analyticalEddyModelOSSE(alongtrackLatLon,eddy_model);

%% Test 3 - Precessing Elliptical Eddy:
clearvars params
params.A = 0.15;
L = 80e3;
params.La = 0.4*2*L;
params.Lb = 0.2*2*L;
params.thetaDot= -2*pi/totalDays; %per day or same as eddy path function time increment (i.e. relative to propation speed)

% eddy_model is a function handle with a chosen set of parameters (x,y,t)
% eddyShapeString = 'Ellipse';
eddy_model = analyticalEddyModel(eddyPath_fun_t,params);

% apply OSSE on an eddy_model function (x,y,t)
alongtrack = analyticalEddyModelOSSE(alongtrackLatLon,eddy_model);
%% subplots
[mz, xmid, ymid, ~, stdz] = composite2D(alongtrack,eddyPath_fun_t);
[mzxy, rmid, ~, stdz_rt] = radialProfile(alongtrack,eddyPath_fun_t);
use stdz_rt

figure(1);hold on
figure('color','w')
figure(1)

set(1,'Units','centimeters');  %Ensure all dimensions are in cm
set(1,'Position',[0 0 17 24.5]);  %Size of the figure. The first two numbers indicate the location of the lower left-hand corner.
                                    %The second two numbers indicate size
set(1,'Renderer','opengl');         %Instead of painter

h=subplot(1,3,1);hold on 
set(h,'units','centimeters')
jpcolor(xmid/1e3, ymid/1e3, mz*1e2);
% r=mean(eddy.speed_radius{1});
% th = 0:pi/50:2*pi;
% plot(r * cos(th),r*sin(th))
% legend('','mean radius')
text(-270, 270, '(a)', 'fontsize', 12, 'fontname', 'times')

shading flat
axis equal
xlabel('Distance East (km)', 'FontName', 'times')
ylabel('Distance North (km)', 'FontName', 'times')
set(gca, 'fontname', 'times','fontsize',16)
colormap(brewermap([], '-Spectral'))
c = colorbar('EastOutside');
c.Label.String = 'Time-averaged SSH (cm)';
c.Label.FontSize=16;
c.Label.FontName='times';
xlim([-250, 250]), ylim([-250, 250]);
set(h,'position',[13 19.5-4.45*(i-1) 3 3.5])

% plot 2D temporal std
figure(1)
subplot(1,3,2);hold on 
jpcolor(xmid/1e3, ymid/1e3, stdz*1e2);
text(-270, 270, '(a)', 'fontsize', 12, 'fontname', 'times')

shading flat
% hold on
% contour(xE,yE,stdTemp2D,[0:0.1:2],'k')
axis equal
xlabel('Distance East (km)', 'FontName', 'times')
ylabel('Distance North (km)', 'FontName', 'times')
set(gca, 'fontname', 'times','fontsize',16)
colormap(brewermap([], '-Spectral'))
c = colorbar('EastOutside');
c.Label.String = 'Temporal variance (cm)';
c.Label.FontSize=16;
c.Label.FontName='times';
xlim([-250, 250]), ylim([-250, 250]); %clim([0, 8])

figure(1);
subplot(1,3,3);hold on 
h = plot(rmid/1e3, [mzxy, avgAziStdTemp, stdAziAvgTemp, stdTotalxy, avgTempStdAzi, stdTempAvgAzi, stdTotalrt]*1e2);
linestyle 2k 2T 2U-- 2W 2k-- 2V: 2X-. 2W--
hlines(0, 'k:')

xlabel('Radial distance (km)', 'FontName', 'times')
ylabel('SSH Variance (cm)', 'FontName', 'times')
set(gca, 'ticklength', [0.03,0.02],'fontname', 'times','fontsize', 16)
lg = legend(h, '$\overline{\eta_{xy}}^\theta$', '$\overline{\sigma_\eta}^\theta$', '$\varsigma_{\overline{\eta}^t}$', '$\Sigma_{\eta_{xy}}$', '$\overline{\varsigma_\eta}^t$', '$\sigma_{\overline{\eta}^\theta}$', '$\Sigma_{\eta_{rt}}$');
set(lg, 'interpreter', 'latex', 'fontsize', 16, 'orientation', 'vertical', 'NumColumns', 2)
xlim([0, 250]); %ylim([-2, 30])
%%
% Save
folder_name='E:\Research\AlongTrack-GRL\fig';

set(gcf,'inverthardcopy','off') %to save white text as white
set(1,'paperpositionmode','auto')
figname='analytical_Gaussian';
print('-dpng','-r600',strcat(folder_name,'\',figname,'.png'))
% print('-deps','-r600',strcat(folder_name,'\',figname,'.eps'))