%% Make plots to check distribution behavior
figure
lineR = xlsread('results/checks/LineOutageRealiz2_v1-3.csv');
lineD = xlsread('results/checks/LinesDist_v3.csv');
% lineR = xlsread('results/checks/LineOutageRealiz_v2.csv');
% lineD = xlsread('results/checks/LinesDist_v2.csv');
rtimeR = xlsread('results/checks/RecTimeRealiz_v2.csv');
rtimeD = xlsread('results/checks/RecTimeDist_v2.csv');
% lineR = xlsread('results/checks/LineOutageRealiz.csv')
% lineD = xlsread('results/checks/LinesDist.csv')
% rtimeR = xlsread('results/checks/RecTimeRealiz.csv');
% rtimeD = xlsread('results/checks/RecTimeDist.csv');

LR_ao = lineR(:,1);
LR_nao = lineR(:,2);
LR_p1 = lineR(:,3);
P = lineR(:,4);
subplot(1,2,1)
% [lineRX,lineRY] = stairs(lineR(:,2),lineR(:,3));
% semilogx(lineR(:,2),lineR(:,3),'r')
% semilogx(lineRX,lineRY,'r')
semilogx(LR_p1,P,'r')
hold on
% semilogx(LR_nao,P,'g')
% semilogx(LR_p1,P,'k')
[DX,DY] = stairs(lineD(:,1),lineD(:,3));
semilogx(DX,DY,'b')
% semilogx(lineD(:,1),lineD(:,3),'b.')
legend('empirical-allow 0 plus 1','analytic')
title('No. Line Outages'); xlabel('k'); ylabel('Prob(no. lines \geq k)')
set(gca, 'fontsize',13)
hold off
subplot(1,2,2)
semilogx(rtimeR(:,2),rtimeR(:,3),'r')
hold on
semilogx(rtimeD(:,1),rtimeD(:,3),'b')
legend('empirical','analytic')
title('Restoration Times'); xlabel('T (min)'); ylabel('Prob(restoration time \geq T)')
set(gca, 'fontsize',13)
hold off


%Plot of line outages only
% figure
% plot(LR_ao,P,'r')
% hold on
% plot(LR_nao,P,'g')
% plot(LR_p1,P,'k')
% % [DX,DY] = stairs(lineD(:,1),lineD(:,3));
% plot(lineD(:,1),lineD(:,3),'b.')
% legend('empirical-allow 0 out','empirical-fair sample','empirical-plus 1','analytic')
% title('No. Line Outages'); xlabel('k'); ylabel('Prob(no. lines \geq k)')
% set(gca, 'fontsize',13)
% axis([0,15,0,1])
% hold off
