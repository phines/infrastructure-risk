%% Make plots to check distribution behavior
figure
% lineR = xlsread('results/checks/LineOutageRealiz_v2.csv');
% lineD = xlsread('results/checks/LinesDist_v2.csv');
% rtimeR = xlsread('results/checks/RecTimeRealiz_v2.csv')
% rtimeD = xlsread('results/checks/RecTimeDist_v2.csv')
lineR = xlsread('results/checks/LineOutageRealiz.csv');
lineD = xlsread('results/checks/LinesDist.csv');
rtimeR = xlsread('results/checks/RecTimeRealiz.csv')
rtimeD = xlsread('results/checks/RecTimeDist.csv')

subplot(1,2,1)
semilogx(lineR(:,2),lineR(:,3),'r')
hold on
semilogx(lineD(:,1),lineD(:,3),'b*')
hold off
subplot(1,2,2)
semilogx(rtimeR(:,2),rtimeR(:,3),'r')
hold on
semilogx(rtimeD(:,1),rtimeD(:,3),'b*')
hold off