function plotComposite(x,y,t_i,zeta,eddyPath,options)
arguments
    x (:,1) {mustBeNumeric}
    y (:,1) {mustBeNumeric}
    t_i (1,1) {mustBeNumeric}
    zeta (:,:) {mustBeNumeric}
    eddyPath struct
    options.contour struct
end

if isa(eddyPath.xe, 'function_handle')
    xo = eddyPath.xe(t_i);
    yo = eddyPath.ye(t_i);
else
    xo = eddyPath.xe;
    yo = eddyPath.ye;
end

% Calculate eddy-relative coordinates directly
xE = x - xo;
yE = y - yo;
figure;hold on
jpcolor(xE/1e3, yE/1e3, zeta');
if isfield(options,'contour') %core
    plot((options.contour.core(1,:)-xo)/1e3,(options.contour.core(2,:)-yo)/1e3,'k','LineWidth',1)
    
    if isfield(options.contour,'shield') % shield
    plot((options.contour.shield(1,:)-xo)/1e3,(options.contour.shield(2,:)-yo)/1e3,'k','LineWidth',1)
    end
end
% r=mean(eddy.speed_radius{1});
% th = 0:pi/50:2*pi;
% plot(r * cos(th),r*sin(th))
% legend('','mean radius')
% jpcolor(xEbin, yEbin, AvgsshAccumBin)
shading flat
axis equal
xlabel('Distance East (km)', 'FontName', 'times')
ylabel('Distance North (km)', 'FontName', 'times')
set(gca, 'fontname', 'times','fontsize',16)
if exist('brewermap')
    colormap(brewermap([], '-Spectral'))
end
c = colorbar('EastOutside');
c.Label.String = '\zeta/f_0';
c.Label.FontSize=16;
c.Label.FontName='times';
xlim([-250, 250]), ylim([-250, 250]);