% Resilience bar charts for MEPS conference paper
% Molly Kelly-Gorham 5/21/2019

clear

% Original 39 bus p1 sampling.
% data1 = xlsread('..\results\experiments\9\res_out_case39_p1_4_lines.csv');
% data2 = xlsread('..\results\experiments\9\res_out_case39_05PV_p1_4_lines.csv');
% data3 = xlsread('..\results\experiments\9\res_out_case39_20PV_p1_4_lines.csv');

% costs1 = xlsread('..\results\experiments\9\res_out_case39_p1_4.csv');
% costs2 = xlsread('..\results\experiments\9\res_out_case39_05PV_p1_4.csv');
% costs3 = xlsread('..\results\experiments\9\res_out_case39_20PV_p1_4.csv');

% N-1 Secure 39 bus p1 sampling.
data1 = xlsread('..\results\experiments\9\res_out_case39_n-1_p1_4_lines.csv');
data2 = xlsread('..\results\experiments\9\res_out_case39_n-1_05PV_p1_4_lines.csv');
data3 = xlsread('..\results\experiments\9\res_out_case39_n-1_20PV_p1_4_lines.csv');

lines1 = data1(:,1);
lines2 = data2(:,1);
lines3 = data3(:,1);


time1 = data1(:,2);
time2 = data2(:,2);
time3 = data3(:,2);

costs1 = xlsread('..\results\experiments\9\res_out_case39_n-1_p1_4.csv');
costs2 = xlsread('..\results\experiments\9\res_out_case39_n-1_05PV_p1_4.csv');
costs3 = xlsread('..\results\experiments\9\res_out_case39_n-1_20PV_p1_4.csv');

costs1(isnan(costs1))=0;
costs2(isnan(costs2))=0;
costs3(isnan(costs3))=0;

% find eens for different sized events
small = 1:2;
medium = 3:7;
large = 7:60;

d(1,1)= mean(sort_data(small,lines1,costs1));
d(1,2)= mean(sort_data(medium,lines1,costs1));
d(1,3)= mean(sort_data(large,lines1,costs1));

d(2,1)= mean(sort_data(small,lines2,costs2));
d(2,2)= mean(sort_data(medium,lines2,costs2));
d(2,3)= mean(sort_data(large,lines2,costs2));

d(3,1)= mean(sort_data(small,lines3,costs3));
d(3,2)= mean(sort_data(medium,lines3,costs3));
d(3,3)= mean(sort_data(large,lines3,costs3));

figure
bar(d')
hold on
legend("original", "+5% load in DG", "+20% load in DG")
xlabel(['small                  ' 'medium                 ' 'large'])
set(gca,'xticklabels',[])
ylabel('EENS (MWh)')
set(gca, 'fontsize',13)
legend boxoff
box off
hold off


% % find the pdf
% [x1,pdf1] = emperical_pdf(costs1,80);
% [x2,pdf2] = emperical_pdf(costs2,80);
% [x3,pdf3] = emperical_pdf(costs3,80);
% 
% % make ccdf
% N = length(costs1);
% Pr=(N:-1:1)/N;
% sorted_costs1 = sort(costs1);
% sorted_costs2 = sort(costs2);
% sorted_costs3 = sort(costs3);
% for jj = 1:N
%     if sorted_costs1(jj)<0.001 || isnan(sorted_costs1(jj))
%     sorted_costs1(jj) = 0;
%     end
%     if sorted_costs2(jj)<0.001 || isnan(sorted_costs2(jj))
%         sorted_costs2(jj) = 0;
%     end
%     if sorted_costs3(jj)<0.001 || isnan(sorted_costs3(jj))
%         sorted_costs3(jj) = 0;
%     end
% end
% 
% figure
% semilogx(sorted_costs1,Pr)
% hold on
% % title("CCDF of Resilience of N-1 secure 39 bus cases")
% semilogx(sorted_costs2,Pr)
% semilogx(sorted_costs3,Pr)
% legend("original", "+5% load in DG", "+20% load in DG")
% ylabel("Prob(Energy lost \geq ENS)")
% set(gca, 'fontsize',13)
% xlim([1 10^7.5])
% box off
% legend boxoff
% set(gca,'xticklabels',[])
% hold off
