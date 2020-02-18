% r1 = xlsread('..\data\solar+load\NY_NE_1244665_solarPV_power_density.csv');
% figure
% hold on
% x = 0:24:24*365;
% plot(1:(30*24),r1((x(140)+1):(x(140)+(30*24)),2))
% hold off
% figure
% hold on
% x = 0:24:24*365;
% plot(1:24,r1((x(26)+1):x(27),2))
% plot(1:24,r1((x(16)+1):x(17),2))
% plot(1:24,r1((x(160)+1):x(161),2))
% plot(1:24,r1((x(144)+1):x(145),2))
% xlabel('time (hours)')
% ylabel('power density')
% legend('winter-sunny','winter-cloudy','summer-sunny','summer-cloudy')
% box off
% legend boxoff
% set(gca,'FontSize',14)
% % xlim([1,24])
% hold off
r2 = xlsread('..\results\RestorationPolish_Case_redo_ex_1.csv');
figure()
hold on
plot(r2(:,1),r2(:,2),'LineWidth',1.5)
xlim([-1,300])
xlabel('minutes')
title('Load shed (MW)','FontWeight','Normal')
%ylabel('load shed (MW)')
grid on
box off
set(gca,'FontSize',14)