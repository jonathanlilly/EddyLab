function makePropagatingVideo(x,y,time,eddy_model,video_name)
% makePropagatingVideo Creates a video of a propagating eddy field
%
% Inputs:
%   x - Array of x coordinates (meters)
%   y - Array of y coordinates (meters)
%   totalDays - Total simulation days
%   eddy_model - Function handle representing the eddy model
%   video_name - Name for the output video file
%
arguments
    x (1,:) {mustBeNumeric}
    y (1,:) {mustBeNumeric}
    time (:,:) {mustBeNumeric, mustBePositive}
    eddy_model function_handle
    video_name (1,1) string = "eddy_field"
end

%name is video file name
[xMat, yMat] = ndgrid(x, y);

figure;
% Create VideoWriter object
v = VideoWriter(strcat(video_name,'.mp4'), 'MPEG-4');
v.FrameRate = 30;
open(v);

% figure loop
for n = 1:length(time)
ssh = eddy_model(xMat',yMat',time(n)-min(time));
jpcolor(x/1e3,y/1e3,ssh*1e2), shading flat
set(gca, 'fontname', 'times','fontsize',16)
colormap(brewermap([], '-Spectral'))
c = colorbar('EastOutside');
c.Label.String = 'SSH (cm)';
c.Label.FontSize=16;
xlabel('Distance East (km)', 'FontName', 'times')
ylabel('Distance North (km)', 'FontName', 'times')
axis equal,axis tight 
clim([-2.5, 12])
colorbar
% Capture the frame
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v)