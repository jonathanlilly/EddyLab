function plotFitPosition(paramsCell,initParamsCell)
totalTimeWindows=length(paramsCell);
th = 0:pi/50:2*pi;
figure;hold on
for j = 1:totalTimeWindows
    if totalTimeWindows==1
        initParams=initParamsCell;
        params=paramsCell;
    else
        initParams=initParamsCell{j};
        params=paramsCell{j};
    end
    % True positions
    xo_true = initParams.x0 + initParams.cx*max(params.elapsed_time);
    yo_true = initParams.y0 + initParams.cy*max(params.elapsed_time);
    L_true = initParams.L;

    % Fit position
    xo_fit = params.x0 + params.cx*max(params.elapsed_time);
    yo_fit = params.y0 + params.cy*max(params.elapsed_time);
    L_fit = params.L;

    plot(L_true*sin(th)+xo_true, L_true*cos(th)+yo_true,'r--');
    plot(L_fit*sin(th)+xo_fit, L_fit*cos(th)+yo_fit,'b');

    plot(xo_true, yo_true,'r*');
    plot(xo_fit, yo_fit,'b*');

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