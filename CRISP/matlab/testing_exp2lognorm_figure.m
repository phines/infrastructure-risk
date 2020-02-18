r1 = xlsread('..\results\test_exponential2lognormal_p1_0.3.csv');
r2 = xlsread('..\results\test_exponential2lognormal_p1_0.4.csv');
r3 = xlsread('..\results\test_exponential2lognormal_p1_0.5.csv');
r4 = xlsread('..\results\test_exponential2lognormal_p1_0.8.csv');
r5 = xlsread('..\results\test_exponential2lognormal_p1_0.9.csv');
r6 = xlsread('..\results\test_exponential2lognormal_p1_1.0.csv');
r7 = xlsread('..\results\test_exponential2lognormal_p1_1.1.csv');
results0 = r1(:,1);
results1 = r1(:,2);
results2 = r2(:,2);
results3 = r3(:,2);
results4 = r4(:,2);
results5 = r5(:,2);
results6 = r6(:,2);
results7 = r7(:,2);
n = length(results0);

% figure(1); clf;
% hist(results1,20);
% hold on
% hist(results1,20,'r');
% legend('lognormal(3.22,2.43)','weibull(1000,10),factor 1.5')
% hold off

figure(2); clf;
x = sort(results0);
p = (n:-1:1)/n;
loglog(x,p,'LineWidth',2);
hold on
% x = sort(results1);
% loglog(x,p);
x = sort(results2);
loglog(x,p,'LineWidth',2);
% x = sort(results3);
% loglog(x,p);
% x = sort(results4);
% loglog(x,p);
% x = sort(results5);
% loglog(x,p);
% x = sort(results6);
% loglog(x,p);
% x = sort(results7);
% loglog(x,p);
legend('lognormal(3.22,2.43)','weibull(0.4,60),factor 1.1')
% legend('lognormal(3.22,2.43)','weibull(0.3,60),factor 1.1',...
%     'weibull(0.4,60),factor 1.1','weibull(0.5,60),factor 1.1')
% legend('lognormal(3.22,2.43)','weibull(0.3,10),factor 1.1',...
%     'weibull(0.4,10),factor 1.1','weibull(0.5,10),factor 1.1',...
%     'weibull(0.8,10),factor 1.1','weibull(0.9,10),factor 1.1','weibull(1.0,10)','weibull(1.1,10)')
box off
legend boxoff
set(gca,'FontSize',14)
hold off

%% p2

r1 = xlsread('..\results\test_exponential2lognormal_p2_40.csv');
r2 = xlsread('..\results\test_exponential2lognormal_p2_45.csv');
r3 = xlsread('..\results\test_exponential2lognormal_p2_50.csv');
r4 = xlsread('..\results\test_exponential2lognormal_p2_60.csv');
r5 = xlsread('..\results\test_exponential2lognormal_p2_70.csv');
results0 = r1(:,1);
results1 = r1(:,2);
results2 = r2(:,2);
results3 = r3(:,2);
results4 = r4(:,2);
results5 = r5(:,2);

figure(3)
x = sort(results0);
p = (n:-1:1)/n;
loglog(x,p);
hold on
x = sort(results1);
loglog(x,p);
x = sort(results2);
loglog(x,p);
x = sort(results3);
loglog(x,p);
x = sort(results4);
loglog(x,p);
x = sort(results5);
loglog(x,p);
%legend('lognormal(3.22,2.43)','weibull(1,15),factor 1.1',...
    %'weibull(1,20),factor 1.1','weibull(1,25),factor 1.1',...
    %'weibull(1,30),factor 1.1','weibull(1,40),factor 1.1')

legend('lognormal(3.22,2.43)','weibull(1,40),factor 1.1',...
    'weibull(1,45),factor 1.1','weibull(1,50),factor 1.1',...
    'weibull(1,60),factor 1.1','weibull(1,70),factor 1.1')
hold off