function plotFitWindowed(paramsCell,true_params, eddyPath_fun_t)
totalTimeWindows=length(paramsCell);
window_size=max(floor(paramsCell{1}.elapsed_time))+1;
scale_factors = [1e-2, 1e3, 1e3, 1e3, 1e3, 1e3];
if totalTimeWindows>1
param_label = {'A','L','x_o','y_o','v_x','v_y'};
param_var = {'A','L','x0','y0','cx','cy'};

if length(initParams.A)==1
    initParams.A=true_params.A;
    initParams.L=true_params.L;
else
    initParams.A=true_params.A(window_start_day(i)-t0+1);
    initParams.L=true_params.L(window_start_day(i)-t0+1);
end

% new inital parameters for this window
initParams.A = A+2*(rand-0.5)*1e-2; %random uncertainty +/- 1e-2
initParams.L = L+2*(rand-0.5)*1e3; %random uncertainty +/- 1e3
%assuming you roughly know the eddy center from eddy-tracking algorithm
%beginning of this particular window minus t0 offset of the entire eddy lifetime
initParams.x0 = eddyPath_fun_t.xe(t0_window-t0)+2*(rand-0.5)*10e3; %random uncertainty +/- 10e3
initParams.y0 = eddyPath_fun_t.ye(t0_window-t0)+2*(rand-0.5)*10e3;
initParams.cx = cx(1)+2*(rand-0.5)*1e2;  %random uncertainty +/- 1e2
initParams.cy = cy(1)+2*(rand-0.5)*1e2;
initParamsCell{i,1} = initParams;

xlimits = [true_params.A+[-0.05,0.05];true_params.L+[-10e3,10e3];true_params.x0+[-200e3,200e3];true_params.y0+[-200e3,200e3];true_params.cx+[-0.2e3,0.2e3];true_params.cy+[-0.2e3,0.2e3]];
for j = 1:totalTimeWindows
    %redefine true params positions
    true_params.x0=eddyPath_fun_t.xe(paramsCell{j}.t0-paramsCell{1}.t0);
    true_params.y0=eddyPath_fun_t.ye(paramsCell{j}.t0-paramsCell{1}.t0);
    for i = 1:6
        paramValues(j, i) = paramsCell{j}.(param_var{i});
        initParamValues(j, i) = initParamsCell{i,1}.(param_var{i});
    end
end
figure;hold on
for i = 1:6
    subplot(6,1,i)  % 6 rows, 1 column
    hold on
    plot([1:totalTimeWindows],paramValues(:, i)./scale_factors(i), 'k.', 'MarkerSize', 15);
    plot([1:totalTimeWindows],initParamValues(:, i)./scale_factors(i), 'r:', 'LineWidth', 1.5)  % Initial value
    set(gca, 'fontname', 'times','fontsize',11)
    hold off
    ylabel(param_label{i},'fontname', 'times','fontsize',14)
    if i~=3 && i~=4
    %skip x0 and y0
    % ylim([xlimits(i,1:2)]./scale_factors(i))
    end
     if i == 6
        xlabel(strcat('Time (window=',num2str(window_size),' Days)'), 'fontname', 'times','fontsize',14)
     else
         set(gca,'xticklabel',[])
     end
     xlim([1,totalTimeWindows])
    h = get(gca, 'Position');
    set(gca, 'Position', [h(1) h(2)+(6-i)*0.01 h(3) h(4)+ 0.01])
end
% 
end
%  % Spatial window figure
%     figure(n),hold on
%     if i==1 && i~=totalTimeWindows
%         % jpcolor(xc, yc, ssh(:,:,idx)')
%         % plot(xt, yt, linewidth = 2,Color=[0.7,0.7,0.7]), axis tight, axis equal,
%         plot([xilim(1) xilim(1) xilim(2) xilim(2) xilim(1)],[yilim(1) yilim(2) yilim(2) yilim(1) yilim(1)],linewidth=2)
%         set(gca, 'fontname', 'times')
%         xlim([-1000,1000]);ylim([-500,500])
%         % colormap(brewermap([], '-Spectral'))
%         % c = colorbar('EastOutside');
%         % c.Label.String = 'ssh (cm)';
%         % clim([-2.5, 12])
%         xlabel('x (km)')
%         ylabel('y (km)')
% 
%     elseif i==totalTimeWindows
%         jpcolor(xc, yc, ssh(:,:,idx)')
%         plot(xt, yt, linewidth = 2,Color=[0.7,0.7,0.7]), axis tight, axis equal,
%         plot([xilim(1) xilim(1) xilim(2) xilim(2) xilim(1)],[yilim(1) yilim(2) yilim(2) yilim(1) yilim(1)],linewidth=2)
%         set(gca, 'fontname', 'times')
%         xlim([-1000,1000]);ylim([-500,500])
%         colormap(brewermap([], '-Spectral'))
%         c = colorbar('EastOutside');
%         c.Label.String = 'ssh (cm)';
%         clim([-2.5, 12])
%         xlabel('x (km)')
%         ylabel('y (km)')
%         title(strcat('Time window = [',num2str(time_window(1)),', ',num2str(time_window(end)),'] cycles'))  
%     else
%         plot([xilim(1) xilim(1) xilim(2) xilim(2) xilim(1)],[yilim(1) yilim(2) yilim(2) yilim(1) yilim(1)],linewidth=2)
%     end
% 
%     p0=[A,L,xoi,vxi,yoi,vyi];%[max(ssht),max(r),r(1),0];
%     p0Series=[p0Series;p0];
% 
%     [p]=gaussianFreeCenter(x1(~isnan(ssht_range)),y1(~isnan(ssht_range)),t(~isnan(ssht_range)),ssht_range(~isnan(ssht_range)),p0,LB,UB,options);
%     pSeries(i,:)=p;
% 
%     xc_fit=x_fit2-(p(3)+p(4).*t_fit2);
%     yc_fit=y_fit2-(p(5)+p(6).*t_fit2);
%     ssh_fit=p(1)*exp(-(xc_fit.^2+yc_fit.^2)/p(2)^2);
%     ro_fit2=sqrt(xc_fit.^2+yc_fit.^2);
%     ro_series(:,i)=ro_fit2(:);%reshape(ro_fit2,[],size(ro_fit2,3));
%     ssh_series(:,i)=ssh_fit(:);
% end

