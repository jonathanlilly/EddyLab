Tol = 1e-2; %tolerance for variance to vanish
prettyblue=[0, 0.4470, 0.7410];
prettyorange=[0.8500, 0.3250, 0.0980];
prettyyellow=[0.9290, 0.6940, 0.1250];
prettypurple=[0.4940, 0.1840, 0.5560];
prettygreen=[0.4660 0.6740 0.1880];
%% Test 1: Azimuthally symmetric eddy with radially symmetric perturbations
A = 1;
R = 1 / sqrt(-2*log(Tol)); %so that sqrt(2)/R=0. Roughly R=0.33
totalDays = 100;
x = -1:0.01:1;
y = x;
[xMat, yMat] = meshgrid(x, y);
r = sqrt(xMat.^2+yMat.^2);

% Azimuthally symmetric eddy
z0 = A * exp(-0.5*(r / R).^2);

% Add radial perturbation that changes sinusoidally each day
for t = 1:totalDays
    % Radially symmetric perturbation that varies with time
    perturb = 0.01 * sin(t*6*pi/totalDays);
    z(:, :, t) = A * exp(-0.5*(r * (1 + perturb) / R).^2);
end

% Compute statistics
[mz, rmid, numz, stdz1, stdz2] = radialStatisticsFromGridded(x, y, z, firstAverage = 'temporal');
stdAziAvgTemp = stdz1; %azimuthal variance of the temporal average
AvgAziStdTemp = stdz2; %azimuthal average of the temporal variance
[mz, rmid, numz, stdz1, stdz2] = radialStatisticsFromGridded(x, y, z, firstAverage = 'azimuthal');
stdTempAvgAzi = stdz1; %temporal variance of the azimuthal average
AvgTempStdAzi = stdz2; %temporal average of the azimuthal variance

% Verify the results
disp('Test1:')
v1 = aresame(rms(AvgAziStdTemp), 0, Tol);
v2 = aresame(rms(stdAziAvgTemp), 0, Tol);
v3 = aresame(rms(AvgTempStdAzi), 0, Tol);
v4 = aresame(rms(stdTempAvgAzi), 0, Tol);

if v1
    disp('Azimuthal average of the temporal variance vanished')
else
    disp('Azimuthal average of the temporal variance is not vanished')
end
if v2
    disp('Azimuthal variance of the temporal average vanished')
else
    disp('Azimuthal variance of the temporal average is not vanished')
end
%azimuthal variation zero
if v3
    disp('Temporal average of the azimuthal variance vanished')
else
    disp('Temporal average of the azimuthal variance is not vanished')
end
if v4
    disp('Temporal variance of the azimuthal average vanished')
else
    disp('Temporal variance of the azimuthal average is not vanished')
end

if v2 && v3
    disp('Test 1 successful. Eddy is azimuthally symmetric')
else
    disp('Test 1 failed. Eddy is not azimuthally symmetric')
end

% Std Total
stdTotal = sqrt(stdAziAvgTemp.^2+AvgAziStdTemp.^2);

% Plot Test 1 profiles
figure(1); hold on
h1 = plot(rmid, mz, 'LineWidth', 2, 'Color', 'k', 'linestyle', '-');
h2 = plot(rmid, AvgAziStdTemp, 'LineWidth', 2, 'Color', prettyblue, 'linestyle', '-');
h3 = plot(rmid, stdAziAvgTemp, 'LineWidth', 2, 'Color', prettyorange, 'linestyle', '--');
h4 = plot(rmid, AvgTempStdAzi, 'LineWidth', 2, 'Color', prettyyellow, 'linestyle', ':');
h5 = plot(rmid, stdTempAvgAzi, 'LineWidth', 2, 'Color', prettygreen, 'linestyle', '-.');
h6 = plot(rmid, stdTotal, 'LineWidth', 2, 'Color', prettypurple, 'linestyle', '-');
xlabel('Radial distance', 'FontName', 'times')
ylabel('SSH Variability', 'FontName', 'times')
set(gca, 'fontname', 'times')
% % zero horizontal line
% plot([rmid(1), rmid(end-1)], [0, 0], 'k:')
lg = legend([h1, h2, h3, h4, h5, h6], '$\overline{\mathrm{z}}^\theta$','$\overline{\sigma_\mathrm{z}}^\theta$', '$\varsigma_{\overline{\mathrm{z}}^t}$', '$\overline{\varsigma_\mathrm{z}}^t$', '$\sigma_{\overline{\mathrm{z}}^\theta}$', '$\overline{\mathrm{z}}^\theta$','$\Sigma_{z}$');
set(lg, 'interpreter', 'latex', 'fontsize', 15, 'orientation', 'horizontal')
set(gca, 'fontsize', 12)

%% Test 2: Steady elliptical eddy
a = 0.4;
b = 0.2;
z = A * exp(-0.5*(xMat / a).^2-0.5*(yMat / b).^2);

% No temporal variation. Repeat z for other days
z = repmat(z, 1, 1, totalDays);

% Compute statistics
[mz, rmid, numz, stdz1, stdz2] = radialStatisticsFromGridded(x, y, z, firstAverage = 'temporal');
stdAziAvgTemp = stdz1; %azimuthal variance of the temporal average
AvgAziStdTemp = stdz2; %azimuthal average of the temporal variance
[mz, rmid, numz, stdz1, stdz2] = radialStatisticsFromGridded(x, y, z, firstAverage = 'azimuthal');
stdTempAvgAzi = stdz1; %temporal variance of the azimuthal average
AvgTempStdAzi = stdz2; %temporal average of the azimuthal variance

% Verify the results
disp('Test2:')
v1 = aresame(rms(AvgAziStdTemp), 0, Tol);
v2 = aresame(rms(stdAziAvgTemp), 0, Tol);
v3 = aresame(rms(AvgTempStdAzi), 0, Tol);
v4 = aresame(rms(stdTempAvgAzi), 0, Tol);

if v1
    disp('Azimuthal average of the temporal variance vanished')
else
    disp('Azimuthal average of the temporal variance is not vanished')
end
if v2
    disp('Azimuthal variance of the temporal average vanished')
else
    disp('Azimuthal variance of the temporal average is not vanished')
end
%azimuthal variation zero
if v3
    disp('Temporal average of the azimuthal variance vanished')
else
    disp('Temporal average of the azimuthal variance is not vanished')
end
if v4
    disp('Temporal variance of the azimuthal average vanished')
else
    disp('Temporal variance of the azimuthal average is not vanished')
end


if v1 && v4
    disp('Test 2 successful. Eddy is an steady elipse')
else
    disp('Test 2 failed. Eddy is not an steady elipse')
end

% Std Total
stdTotal = sqrt(stdAziAvgTemp.^2+AvgAziStdTemp.^2);

% Plot Test 2 profiles
figure(2); hold on
h1 = plot(rmid, mz, 'LineWidth', 2, 'Color', 'k', 'linestyle', '-');
h2 = plot(rmid, AvgAziStdTemp, 'LineWidth', 2, 'Color', prettyblue, 'linestyle', '-');
h3 = plot(rmid, stdAziAvgTemp, 'LineWidth', 2, 'Color', prettyorange, 'linestyle', '--');
h4 = plot(rmid, AvgTempStdAzi, 'LineWidth', 2, 'Color', prettyyellow, 'linestyle', ':');
h5 = plot(rmid, stdTempAvgAzi, 'LineWidth', 2, 'Color', prettygreen, 'linestyle', '-.');
h6 = plot(rmid, stdTotal, 'LineWidth', 2, 'Color', prettypurple, 'linestyle', '-');
xlabel('Radial distance', 'FontName', 'times')
ylabel('SSH Variability', 'FontName', 'times')
set(gca, 'fontname', 'times')
% % zero horizontal line
% plot([rmid(1), rmid(end-1)], [0, 0], 'k:')
lg = legend([h1, h2, h3, h4, h5, h6], '$\overline{\mathrm{z}}^\theta$','$\overline{\sigma_\mathrm{z}}^\theta$', '$\varsigma_{\overline{\mathrm{z}}^t}$', '$\overline{\varsigma_\mathrm{z}}^t$', '$\sigma_{\overline{\mathrm{z}}^\theta}$', '$\overline{\mathrm{z}}^\theta$','$\Sigma_{z}$');
set(lg, 'interpreter', 'latex', 'fontsize', 15, 'orientation', 'horizontal')
set(gca, 'fontsize', 12)

%% Test 3: Uniformly precessing elliptical eddy
z = zeros(size(xMat, 1), size(xMat, 2), totalDays);

for t = 1:totalDays
    x_rot = xMat * cos(t*pi/totalDays) - yMat * sin(t*pi/totalDays);
    y_rot = xMat * sin(t*pi/totalDays) + yMat * cos(t*pi/totalDays);
    z(:, :, t) = A * exp(-0.5*(x_rot / a).^2-0.5*(y_rot / b).^2);
end

% Compute statistics
[mz, rmid, numz, stdz1, stdz2] = radialStatisticsFromGridded(x, y, z, firstAverage = 'temporal');
stdAziAvgTemp = stdz1; %azimuthal variance of the temporal average
AvgAziStdTemp = stdz2; %azimuthal average of the temporal variance
[mz, rmid, numz, stdz1, stdz2] = radialStatisticsFromGridded(x, y, z, firstAverage = 'azimuthal');
stdTempAvgAzi = stdz1; %temporal variance of the azimuthal average
AvgTempStdAzi = stdz2; %temporal average of the azimuthal variance

% Verify the results
disp('Test3:')
v1 = aresame(rms(AvgAziStdTemp), 0, Tol);
v2 = aresame(rms(stdAziAvgTemp), 0, Tol);
v3 = aresame(rms(AvgTempStdAzi), 0, Tol);
v4 = aresame(rms(stdTempAvgAzi), 0, Tol);

if v1
    disp('Azimuthal average of the temporal variance vanished')
else
    disp('Azimuthal average of the temporal variance is not vanished')
end
if v2
    disp('Azimuthal variance of the temporal average vanished')
else
    disp('Azimuthal variance of the temporal average is not vanished')
end
%azimuthal variation zero
if v3
    disp('Temporal average of the azimuthal variance vanished')
else
    disp('Temporal average of the azimuthal variance is not vanished')
end
if v4
    disp('Temporal variance of the azimuthal average vanished')
else
    disp('Temporal variance of the azimuthal average is not vanished')
end


if v2 && v4
    disp('Test 3 successful. Eddy is eliptical and uniformly precessing')
else
    disp('Test 3 failed. Eddy is not eliptical and uniformly precessing')
end

% Std Total
stdTotal = sqrt(stdAziAvgTemp.^2+AvgAziStdTemp.^2);

% Plot Test 3 profiles
figure(3); hold on
h1 = plot(rmid, mz, 'LineWidth', 2, 'Color', 'k', 'linestyle', '-');
h2 = plot(rmid, AvgAziStdTemp, 'LineWidth', 2, 'Color', prettyblue, 'linestyle', '-');
h3 = plot(rmid, stdAziAvgTemp, 'LineWidth', 2, 'Color', prettyorange, 'linestyle', '--');
h4 = plot(rmid, AvgTempStdAzi, 'LineWidth', 2, 'Color', prettyyellow, 'linestyle', ':');
h5 = plot(rmid, stdTempAvgAzi, 'LineWidth', 2, 'Color', prettygreen, 'linestyle', '-.');
h6 = plot(rmid, stdTotal, 'LineWidth', 2, 'Color', prettypurple, 'linestyle', '-');
xlabel('Radial distance', 'FontName', 'times')
ylabel('SSH Variability', 'FontName', 'times')
set(gca, 'fontname', 'times')
% % zero horizontal line
% plot([rmid(1), rmid(end-1)], [0, 0], 'k:')
lg = legend([h1, h2, h3, h4, h5, h6], '$\overline{\mathrm{z}}^\theta$','$\overline{\sigma_\mathrm{z}}^\theta$', '$\varsigma_{\overline{\mathrm{z}}^t}$', '$\overline{\varsigma_\mathrm{z}}^t$', '$\sigma_{\overline{\mathrm{z}}^\theta}$', '$\overline{\mathrm{z}}^\theta$','$\Sigma_{z}$');
set(lg, 'interpreter', 'latex', 'fontsize', 15, 'orientation', 'horizontal')
set(gca, 'fontsize', 12)