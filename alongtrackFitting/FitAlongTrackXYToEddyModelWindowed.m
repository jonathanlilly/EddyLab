function [paramsCell, initParamsCell] = FitAlongTrackXYToEddyModelWindowed(alongtrack, eddyFit_fun, initParams, eddyPath_fun_t, it_options, options)
arguments
    alongtrack struct
    eddyFit_fun function_handle 
    initParams struct
    eddyPath_fun_t struct
    it_options struct
    options.LB (1,6) double %= [0, 50, -1000e3, -500e3, -10, -10]
    options.UB (1,6) double %= [30, 150, 1000e3, 500e3, 0, 0]
end
%alongtrackLatLon: lat,lon,time
% Deduce time window length from alongtrackLatLon
totalDays = max(alongtrack.t)-min(alongtrack.t);

min_days_per_window = 45;  % Need at least 4.5 cycles for the eddy core coverage
min_total_windows = 3;      % Need at least 3 windows to see evolution
overlap = 0.5; %50% overlap between time_window
window_size = max(min_days_per_window, floor(totalDays/(min_total_windows/overlap)));

time_step=floor(window_size*(1-overlap));
totalTimeWindows=floor((totalDays - window_size) / time_step) + 1;

% ensure that the arrays are in ascending order in time before windowing
[alongtrack.t,sort_idx]=sort(alongtrack.t,'ascend');
alongtrack.x=alongtrack.x(sort_idx);
alongtrack.y=alongtrack.y(sort_idx);
alongtrack.lon=alongtrack.lon(sort_idx);
alongtrack.lat=alongtrack.lat(sort_idx);
alongtrack.ssh=alongtrack.ssh(sort_idx);

% pre-allocate cells for number of windows
paramsCell = cell(totalTimeWindows, 1);

for i=1:totalTimeWindows
    % Calculate the start and end times for this window in days
    window_start_day = min(alongtrack.t) + (i-1)*time_step;
    window_end_day = window_start_day + window_size;

    % Find indices that correspond to times within this window
    window_indices = find(alongtrack.t >= window_start_day & alongtrack.t <= window_end_day);

% Extract time window
at_window.x = alongtrack.x(window_indices);
at_window.y = alongtrack.y(window_indices);
at_window.t = alongtrack.t(window_indices);
at_window.ssh = alongtrack.ssh(window_indices);

% Define t0 for this specific window
t0_window = min(at_window.t);
elapsed_time_window = at_window.t-t0_window;

% % latc and lonc are the center of the alongtrack domain
% lonc=(min(at_window.lon(:))+max(at_window.lon(:)))/2;
% latc=(min(at_window.lat(:))+max(at_window.lat(:)))/2;
% 
% % Project {lat,lon} -> {x,y}
% [at_window.x, at_window.y] = latlon2xy(at_window.latitude, at_window.longitude, latc, lonc);

% new inital parameters for this window
initParams_window.A = initParams.A;
initParams_window.L = initParams.L;
% initParams_window.x0 = initParams.x0+i*initParams.cx*max(elapsed_time_window);
% initParams_window.y0 = initParams.y0+i*initParams.cy*max(elapsed_time_window);
%assuming you roughly know the eddy center from eddy-tracking algorithm
%beginning of this particular window minus t0 offset of the entire eddy lifetime
initParams_window.x0 = eddyPath_fun_t.xe(t0_window-min(alongtrack.t))+rand*20e3-10e3; %random uncertainty +/- 10e3
initParams_window.y0 = eddyPath_fun_t.ye(t0_window-min(alongtrack.t))+rand*20e3-10e3;
initParams_window.cx = initParams.cx;
initParams_window.cy = initParams.cy;
initParamsCell{i,1} = initParams_window;

% bound


% Call eddy model fit in XY
params = FitAlongTrackXYToEddyModel(at_window, eddyFit_fun, initParams_window, it_options);

paramsCell{i,1} = params;
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

if totalTimeWindows>1
param_label = {'A','L','x_o','y_o','v_x','v_y'};
param_var = {'A','L','x0','y0','cx','cy'};
% % xlimits = [8,18;50,100;min(center_xoyo(:,1)),max(center_xoyo(:,1));min(vxo),max(vxo);min(center_xoyo(:,2)),max(center_xoyo(:,2));min(vyo),max(vyo)];
for j = 1:totalTimeWindows
    for i = 1:6
        paramValues(j, i) = paramsCell{j}.(param_var{i});
        initParamValues(j, i) = initParamsCell{j}.(param_var{i});
    end
end
figure;hold on
for i = 1:6
    subplot(6,1,i)  % 6 rows, 1 column
    hold on
    plot([1:totalTimeWindows],paramValues(:, i), 'k.', 'MarkerSize', 15);
    plot([1:totalTimeWindows],initParamValues(:, i), 'r:', 'LineWidth', 1.5)  % Initial value
    set(gca, 'fontname', 'times','fontsize',11)
    hold off
    ylabel(param_label{i},'fontname', 'times','fontsize',14)
    % ylim([xlimits(i,1:2)])
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