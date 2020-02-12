r1 = xlsread('..\results\test_exponential2lognormal_p1_0.1.csv');
r2 = xlsread('..\results\test_exponential2lognormal_p1_0.2.csv');
r3 = xlsread('..\results\test_exponential2lognormal_p1_0.3.csv');
r4 = xlsread('..\results\test_exponential2lognormal_p1_0.4.csv');
r5 = xlsread('..\results\test_exponential2lognormal_p1_0.5.csv');
results0 = r1(:,1);
results1 = r1(:,2);
results2 = r2(:,2);
results3 = r3(:,2);
results4 = r4(:,2);
results5 = r5(:,2);
n = length(results);

% figure(1); clf;
% hist(results1,20);
% hold on
% hist(results,20,'r');
% legend('lognormal(3.22,2.43)','weibull(1000,10),factor 1.5')
% hold off

figure(2); clf;
x = sort(results0);
p = (n:-1:1)/n;
loglog(x,p);
hold on
x = sort(results1);
p = (n:-1:1)/n;
loglog(x,p);
x = sort(results2);
p = (n:-1:1)/n;
loglog(x,p);
x = sort(results3);
p = (n:-1:1)/n;
loglog(x,p);
x = sort(results4);
p = (n:-1:1)/n;
loglog(x,p);
x = sort(results5);
p = (n:-1:1)/n;
loglog(x,p);
legend('lognormal(3.22,2.43)','weibull(0.1,10),factor 1.1',...
    'weibull(0.2,10),factor 1.1','weibull(0.3,10),factor 1.1',...
    'weibull(0.4,10),factor 1.1','weibull(0.5,10),factor 1.1')
hold off

%% p2

r1 = xlsread('..\results\test_exponential2lognormal_p2_15.csv');
r2 = xlsread('..\results\test_exponential2lognormal_p2_20.csv');
r3 = xlsread('..\results\test_exponential2lognormal_p2_25.csv');
r4 = xlsread('..\results\test_exponential2lognormal_p2_30.csv');
r5 = xlsread('..\results\test_exponential2lognormal_p2_40.csv');
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
p = (n:-1:1)/n;
loglog(x,p);
x = sort(results2);
p = (n:-1:1)/n;
loglog(x,p);
x = sort(results3);
p = (n:-1:1)/n;
loglog(x,p);
x = sort(results4);
p = (n:-1:1)/n;
loglog(x,p);
x = sort(results5);
p = (n:-1:1)/n;
loglog(x,p);
legend('lognormal(3.22,2.43)','weibull(0.3,15),factor 1.1',...
    'weibull(0.3,20),factor 1.1','weibull(0.3,25),factor 1.1',...
    'weibull(0.3,30),factor 1.1','weibull(0.3,40),factor 1.1')
hold off