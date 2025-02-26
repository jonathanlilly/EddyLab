function paramsCell = FitAlongTrackXYToEddyModelWindowed(alongtrack, eddyFit_fun, initialParams)
arguments
    alongtrack struct
    eddyFit_fun function_handle 
    initialParams struct
    % options.windowLength (1,1) = 
end
%alongtrackLatLon: lat,lon,time
% Deduce time window length from alongtrackLatLon
at_time = datenum(alongtrack.t); %check if this is sorted in an ascending order.
totalDays = max(at_time)-min(at_time);

min_days_per_window = 45;  % Need at least 4.5 cycles for the eddy core coverage
min_total_windows = 3;      % Need at least 3 windows to see evolution
overlap = 0.5; %50% overlap between time_window
window_size = max(min_days_per_window, floor(totalDays/(min_total_windows/overlap)));

time_step=floor(window_size*(1-overlap));
totalTimeWindows=floor((totalDays - window_size) / time_step) + 1;

for i=1:totalTimeWindows
time_window=1+(i-1)*time_step:(i-1)*time_step+window_size;

% Extract time window
at_window.x = alongtrack.x(time_window);
at_window.y = alongtrack.y(time_window);
at_window.t = alongtrack.t(time_window);
at_window.ssh = alongtrack.ssh(time_window);

% % latc and lonc are the center of the alongtrack domain
% lonc=(min(at_window.lon(:))+max(at_window.lon(:)))/2;
% latc=(min(at_window.lat(:))+max(at_window.lat(:)))/2;
% 
% % Project {lat,lon} -> {x,y}
% [at_window.x, at_window.y] = latlon2xy(at_window.latitude, at_window.longitude, latc, lonc);

% new inital parameters for this window
initialParams_window.A = initialParams.A;
initialParams_window.L = initialParams.L;
initialParams_window.x0 = initialParams.x0+initialParams.cx*time_window(end);
initialParams_window.y0 = initialParams.y0+initialParams.cy*time_window(end);
initialParams_window.cx = initialParams.cx;
initialParams_window.cy = initialParams.cy;

% bound


% Call eddy model fit in XY
params = FitAlongTrackXYToEddyModel(at_window, eddyFit_fun, initialParams_window);

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
param_names = {'A','L','x_o','v_x','y_o','v_y'};
xlimits = [8,18;50,100;min(center_xoyo(:,1)),max(center_xoyo(:,1));min(vxo),max(vxo);min(center_xoyo(:,2)),max(center_xoyo(:,2));min(vyo),max(vyo)];
figure;hold on
for i = 1:6
    subplot(6,1,i)  % 6 rows, 1 column
    plot([1:totalTimeWindows],pSeries(:,i), 'k.', 'MarkerSize', 15);
    hold on
    plot([1:totalTimeWindows], p0Series(:,i), 'r:', 'LineWidth', 1.5)  % Initial value
    hold off
    ylabel(param_names{i},'fontsize',12)
    ylim([xlimits(i,1:2)])
     if i == 6
        xlabel('Time (2 cycles)')
     else
         set(gca,'xticklabel',[])
     end
     xlim([1,totalTimeWindows])
    % h = get(gca, 'Position');
    % set(gca, 'Position', [h(1) h(2) h(3) 0.1])
end

end