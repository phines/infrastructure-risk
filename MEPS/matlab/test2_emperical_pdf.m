%load data
% costs1 = xlsread('..\results\experiments\6\res_out_case39_p1_2.csv');
% costs2 = xlsread('..\results\experiments\6\res_out_case39_05PV_p1_2.csv');
% costs3 = xlsread('..\results\experiments\6\res_out_case39_20PV_p1_2.csv');

costs1 = xlsread('..\results\experiments\8\res_out_case39_n-1.csv');
costs2 = xlsread('..\results\experiments\8\res_out_case39_n-1_05PV.csv');
costs3 = xlsread('..\results\experiments\8\res_out_case39_n-1_20PV.csv');

%remove any NaN
costs1(isnan(costs1))=0;
costs2(isnan(costs2))=0;
costs3(isnan(costs3))=0;

%find the pdf
[x1,pdf1] = emperical_pdf(costs1,80);
[x2,pdf2] = emperical_pdf(costs2,80);
[x3,pdf3] = emperical_pdf(costs3,80);

% make plots
figure(1); clf;
plot(x1,pdf1)
hold on
plot(x2,pdf2)
plot(x3,pdf3)
legend('original', '5% load added in DG', '20% load added in DG')
title('Resilience 39 bus PDF')
xlabel('ENS'); ylabel('PDF(ENS)');
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca, 'fontsize',13)

%double check area is 1 -> outputs NaN...
% area1 = trapz(x1,pdf1)
% area2 = trapz(x2,pdf2)
% area3 = trapz(x3,pdf3)
% 
% figure(2); clf;
% [yvalues,xvalues] = histc(x,20)
% [p,x] = histcounts(x1,20,'Normalization','cdf');
% hist(p,x);
% set(gca,'xscale','log');
% set(gca,'yscale','log');