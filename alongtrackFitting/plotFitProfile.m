%% Compute MSE per day
window_size = 50;
model = 'Gaussian';
rbin_size=12.5*1e3;
max_r=400*1e3;
rbin = [-rbin_size / 2:rbin_size:max_r]';
[rmse_daily,rmse,getModelForDay,window_center]=compute_model_error(alongtrack,eddy_field, eddyPath_fun_t, window_size, model);

% Create a function that returns an interpolated model for any day
% getModelForDay = @(day) createInterpolatedModel(day, window_center, ssh_model);

% Initialize MSE array
mse_t = nan(totalDays, 1);

for n = 1:totalDays
    % Get model for this day
    % if totalTimeWindows == 1
    %     ssh_model_interp = ssh_model{1};
    % else
        ssh_model_interp = getModelForDay(t0 + n - 1);
    % end

    % eddy-centered coordinates
    x_rel = eddy_field.x - eddyPath_fun_t.xe(n-1);
    y_rel = eddy_field.y - eddyPath_fun_t.ye(n-1);
    
    % Find points within spatial window
    in_window_x = x_rel >= min(spatial_window) & x_rel <= max(spatial_window);
    in_window_y = y_rel >= min(spatial_window) & y_rel <= max(spatial_window);
    
    ssh_true_n = eddy_field.ssh(in_window_x, in_window_y, n)';
    % Handle coordinates based on model type
    switch model
        case 'composite'
            % Composite models expect eddy-centered coordinates
            [x_grid, y_grid] = ndgrid(x_rel(in_window_x), y_rel(in_window_y));
            ssh_model_n = ssh_model_interp(x_grid', y_grid', n-1);
            
        case 'Gaussian'
            % Gaussian models expect earth fixed coordinates
            [x_grid, y_grid] = ndgrid(eddy_field.x, eddy_field.y);
            ssh_model_fulldomain = ssh_model_interp(x_grid', y_grid', n-1);
            ssh_model_n = ssh_model_fulldomain(in_window_y, in_window_x);
    end
    [x_eddycenter, y_eddycenter] = ndgrid(x_rel(in_window_x), y_rel(in_window_y));
    [mz, ~, ~, ~, ~] = radialStatisticsFromScatter(x_eddycenter', y_eddycenter', ones(size(x_eddycenter')), ssh_model_n, rbin, [0,1], firstAverage = 'temporal');
    mz_model{1,n}=mz;
    [mz, rmid, ~, ~, ~] = radialStatisticsFromScatter(x_eddycenter', y_eddycenter', ones(size(x_eddycenter')), ssh_true_n, rbin, [0,1], firstAverage = 'temporal');
    mz_true{1,n}=mz;
    rmid_all{1,n}=rmid;
end
%% 
f1 = figure;
f2 = figure;

% Define the time indices to plot
time_indices = [1,256+1];
n_plots = length(time_indices);

for i = 1:n_plots
    n = time_indices(i);
    
    % Calculate transparency (alpha) - starts low, increases with time
    % Alpha ranges from 0.2 to 1.0
    alpha_val = round(0.2 + 0.8 * (i-1) / (n_plots-1),2);
    
figure(f1);
    plot(rmid_all{1,n}/1e3, mz_true{1,n}*1e2, 'k-', 'LineWidth', 2, 'Color', [0, 0, 0, alpha_val]);
    hold on;
    plot(rmid_all{1,n}/1e3, mz_model{1,n}*1e2, 'b-', 'LineWidth', 2, 'Color', [0.07,0.44,0.75, alpha_val]);
    % 
    figure(f2);
    plot(rmid_all{1,n}/1e3, vdiff(mz_true{1,n}*1e2, 1), 'k-', 'LineWidth', 2, 'Color', [0, 0, 0, alpha_val]);
    hold on;
    plot(rmid_all{1,n}/1e3, vdiff(mz_model{1,n}*1e2, 1), 'b-', 'LineWidth', 2, 'Color', [0, 0, 1, alpha_val]);
end

% Format Figure 1
figure(f1);
xlabel('Radial Distance (km)', 'FontName', 'times');
ylabel('SSH (cm)', 'FontName', 'times');
% title('SSH Profiles: True vs Model', 'FontName', 'times', 'FontSize', 14);
legend('True$_{t=0}$', 'Composite$_{t=0}$','True$_{t=256}$', 'Composite$_{t=256}$','Gaussian$_{t=256}$', 'interpreter','latex','Location', 'northeast','orientation','horizontal','NumColumns',2);
xlim([0,250])
set(gca,'FontName', 'times', 'FontSize', 16);

% Format Figure 2
figure(f2);
xlabel('Radial Distance (km)', 'FontName', 'times');
ylabel('SSH Gradient $\partial \eta / \partial r$','interpreter','latex','FontName', 'times');
% title('SSH Gradients: True vs Model', 'FontName', 'times', 'FontSize', 14);
legend('True$_{t=0}$', 'Composite$_{t=0}$','True$_{t=256}$', 'Composite$_{t=256}$','Gaussian$_{t=256}$', 'interpreter','latex','Location', 'southeast','orientation','horizontal','NumColumns',2);
set(gca, 'FontName', 'times', 'FontSize', 16);
xlim([0,250])
%%
window_size_array=[366	100	50	25	10	5	1];
model='composite';
data=alongtrack;

clearvars rmse_array rmse_daily_array

f1=figure;hold on
for w=1:length(window_size_array)
    window_size=window_size_array(w);
    [rmse_daily,rmse,~,~]=compute_model_error(data,eddy_field, eddyPath_fun_t, window_size, model);
    rmse_array(w)=rmse;
    rmse_daily_array(:,w)=rmse_daily;
    figure(f1);plot([1:totalDays]-1,rmse_daily(:)*1e2,'linewidth',2)
end

xlabel('Time (day)')
ylabel('RMSE (cm)')
set(gca,'fontname','times','fontsize',16)
lg=legend(strcat(['$t_w=',num2str(window_size_array(1)),'$ days']),strcat(['$t_w=',num2str(window_size_array(2)),'$ days']),strcat(['$t_w=',num2str(window_size_array(3)),'$ days']),strcat(['$t_w=',num2str(window_size_array(4)),'$ days']),strcat(['$t_w=',num2str(window_size_array(5)),'$ days']),strcat(['$t_w=',num2str(window_size_array(6)),'$ days']));
set(lg,'interpreter','latex','fontname','times','fontsize',16)
xlim([1,totalDays]-1)

figure;plot(window_size_array,rmse_array*1e2,'k-o','MarkerFaceColor','black','Color','black','MarkerSize',6,'LineWidth',2)
xlabel('Time window size (days)')
ylabel('RMSE (cm)')
set(gca,'fontname','times','fontsize',16)

