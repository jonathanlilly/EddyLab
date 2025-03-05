% Create a custom plot function to visualize parameter values
function stop = plotParamsIteraction_ref(xtrans, optimValues, state, x, true_params,flag)
persistent data p0_vals ref_params
    if strcmp(state,'init')
        data = [];
        p0_vals=xtrans;
        ref_params = true_params; % Store the reference parameters
    end
    if strcmp(state,'iter')
        data = [data, xtrans(:)];
        param_label = {'A','L','x_o','y_o','v_x','v_y'};
        param_var = {'A','L','x0','y0','cx','cy'};
        scale_factors = [1e-2, 1e3, 1e3, 1e3, 1e3, 1e3];
        for i = 1:6
            subplot(6,1,i)  % 6 rows, 1 column
            plot(data(i,:), '.', 'MarkerSize', 15);
            hold on
            plot(xlim, [ref_params.(param_var{i}); ref_params.(param_var{i})]./scale_factors(i), 'r:', 'LineWidth', 1.5)  % Initial value
            hold off
            ylabel(param_label{i},'fontsize',12)
            grid on

            % % Make plots taller by adjusting height
            % h = get(gca, 'Position');
            % bottom = 0.07 + (6-i) * ((h(4)+0.1) + 0.02); %bottomMargin +(6-i)*(subplotHeight + btwnPlotSpace)
            % set(gca, 'Position', [h(1) h(2)+0.12*i h(3) (h(4)+0.1)])
            % Only put xlabel on bottom plot
            if i == 6
                xlabel('Iteration')
            end
            
        end
        
        % plot(data.', '.', 'MarkerSize', 15);  % Just dots, no lines

        % legend(arrayfun(@(i) sprintf('p%d', i), 1:length(xtrans), 'UniformOutput', false));
    end
    stop = 0;
end