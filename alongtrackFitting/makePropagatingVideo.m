function makePropagatingVideo(x,y,totalDays,eddy_model,name)
%name is video file name
[xMat, yMat] = ndgrid(x, y);

figure;
% Create VideoWriter object
v = VideoWriter(strcat(name,'.mp4'), 'MPEG-4');
v.FrameRate = 10;
open(v);

% figure loop
for n = 1:totalDays
ssh = eddy_model(xMat',yMat',n);
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