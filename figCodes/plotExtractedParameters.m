%figure for parameters in time

figure;
subplot(3,1,2);
plot([1:totalDays], radius, 'k-');
xlabel('Time (Days)');
ylabel('Radius (km)');
title('SSH e-fold contour radius');

subplot(3,1,3);
plot([1:totalDays], amplitude, 'k-');
xlabel('Time (Days)');
ylabel('Amplitude (cm)');
title('Max SSH');

subplot(3,2,1);hold on
plot([1:totalDays],at,'k','LineWidth',2)
plot([1:totalDays],bt,'k--','LineWidth',2)
xlabel('Time (Days)');
ylabel('a,b (km)');
legend('$a$','$b$','Location','northwest','interpreter','latex')

subplot(3,2,2);
plot([1:totalDays], lambdat, 'k-','LineWidth',2);
xlabel('Time (Days)');
ylabel('$\lambda$','interpreter','latex');

subplot(3,2,3);
plot([1:totalDays], degunwrap(jrad2deg(thetat)), 'k.','LineWidth',2);
xlabel('Time (Days)');
ylabel('$\theta$','interpreter','latex');

omega = convertAngularRate(jrad2deg(thetadot));
subplot(3,2,4);
plot([1:totalDays], thetadot, 'k.','LineWidth',2);
xlabel('Time (Days)');
ylabel('$\dot{\theta}$','interpreter','latex');

subplot(3,2,5);hold on
plot([1:totalDays],vdiff(center_xoyo(:,1),1),'k.','LineWidth',2)
plot([1:totalDays],vdiff(center_xoyo(:,2),1),'b.','LineWidth',2)
xlabel('Time (Days)');
ylabel('Velocity (km/day)');
legend('vx','vy','Location','northoutside','Orientation','horizontal')

subplot(3,2,6);
plot([1:totalDays],eddy_speed,'k.','LineWidth',2)
xlabel('Time (Days)');
ylabel('Velocity (km/day)');