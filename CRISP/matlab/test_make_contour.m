% Resilience bar charts for MEPS conference paper
% Molly Kelly-Gorham 5/21/2019

%clear

% % Original 39 bus p1 sampling.
% data1 = xlsread('..\results\experiments\9\res_out_case39_p1_lines.csv');
% data2 = xlsread('..\results\experiments\9\res_out_case39_05PV_p1_lines.csv');
% data3 = xlsread('..\results\experiments\9\res_out_case39_20PV_p1_lines.csv');
% 
% costs1 = xlsread('..\results\experiments\9\res_out_case39_p1.csv');
% costs2 = xlsread('..\results\experiments\9\res_out_case39_05PV_p1.csv');
% costs3 = xlsread('..\results\experiments\9\res_out_case39_20PV_p1.csv');

% N-1 Secure 39 bus p1 sampling.
data1 = xlsread('..\results\experiments\9\res_out_case39_n-1_p1_lines.csv');
data2 = xlsread('..\results\experiments\9\res_out_case39_n-1_05PV_p1_lines.csv');
data3 = xlsread('..\results\experiments\9\res_out_case39_n-1_20PV_p1_lines.csv');

costs1 = xlsread('..\results\experiments\9\res_out_case39_n-1_p1.csv');
costs2 = xlsread('..\results\experiments\9\res_out_case39_n-1_05PV_p1.csv');
costs3 = xlsread('..\results\experiments\9\res_out_case39_n-1_20PV_p1.csv');

lines1 = data1(:,1);
lines2 = data2(:,1);
lines3 = data3(:,1);


time1 = data1(:,2);
time2 = data2(:,2);
time3 = data3(:,2);

ls1 = data1(:,3);
ls2 = data2(:,3);
ls3 = data3(:,3);

costs1(isnan(costs1))=0;
costs2(isnan(costs2))=0;
costs3(isnan(costs3))=0;

figure
subplot(3,1,1)
[Costs, bin_ls, bin_rt] = make_contour(ls1,time1,costs1);
subplot(3,1,2)
[Costs2, bin_ls2, bin_rt2] = make_contour(ls2,time2,costs2);
subplot(3,1,3)
[Costs3, bin_ls3, bin_rt3] = make_contour(ls3,time3,costs3);