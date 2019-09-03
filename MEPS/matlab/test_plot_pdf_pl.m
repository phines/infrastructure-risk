% test plot_pdf_pl
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

% plot_pdf_pl(X,n_bins,min_per_bin)
plot_pdf_pl(data1,50,20)