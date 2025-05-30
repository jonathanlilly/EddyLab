function plotFitPosition(eddyPath_fun_t, paramsCell, trueParamsCell)
%%
totalTimeWindows=length(paramsCell);
th = 0:pi/50:2*pi;
figure;hold on
for j = 1:totalTimeWindows
    if iscell(paramsCell)
        trueParams=trueParamsCell{j};
        params=paramsCell{j};
        t0=paramsCell{1}.t0;
        current_window_t0=paramsCell{j}.t0;
        window_start_day=current_window_t0-t0;
    else
        trueParams=trueParamsCell;
        params=paramsCell;
        window_start_day=0;
    end
    
    % True positions on the last day of the window
    xo_true = mean(eddyPath_fun_t.xe(unique(params.elapsed_time)+window_start_day));
    yo_true = mean(eddyPath_fun_t.ye(unique(params.elapsed_time)+window_start_day));
    A_true = mean(trueParams.A);
    L_true = mean(trueParams.L);

    % Fit position on the last day of the window
    xo_fit = mean(params.x0 + params.cx*unique(params.elapsed_time));
    yo_fit = mean(params.y0 + params.cy*unique(params.elapsed_time));
    A_fit = params.A;
    L_fit = params.L;
    
    plot((L_true*sin(th)+xo_true)/1e3, (L_true*cos(th)+yo_true)/1e3,'r');
    plot((L_fit*sin(th)+xo_fit)/1e3, (L_fit*cos(th)+yo_fit)/1e3,'b--');

    plot(xo_true/1e3, yo_true/1e3,'rx');
    plot(xo_fit/1e3, yo_fit/1e3,'b+');

    % Calculate error
    position_error(j) = sqrt((xo_fit - xo_true).^2 + (yo_fit - yo_true).^2);
end

axis equal

% xlim([min(x),max(x)]);ylim([min(y),max(y)])
box on
legend('True position', 'Fit position')
set(gca, 'fontname', 'times','fontsize',14)
xlabel('$x$ (km)', 'FontName', 'times','Interpreter','latex')
ylabel('$y$ (km)', 'FontName', 'times','Interpreter','latex')