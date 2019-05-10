% simple script to plot histograms when in a hurry
%case 6ww:
%costs1 = xlsread('results\case6ww\test\res_out_test_001.csv');
%costs2 = xlsread('results\case6ww\test\res_out_test_002.csv');

%case39:
%figuring out weird change in distribution
%costs1 = xlsread('results\case39\resilience_case39.csv');
%costs2 = xlsread('results\case39\resilience_case39_4.csv');
%costs3 = xlsread('results\case39\resilience_case39_5.csv');

%costs1 = xlsread('results\case39\test\res_out_test_010.csv'); 
%costs2 = xlsread('results\case39\test\res_out_test_011.csv'); 

%Compare - no lines out
%costs1 = xlsread('results\case39\resilience_out_noLinesOut.csv');
%costs2 = xlsread('results\case39_05PV\resilience_noLinesOut.csv');
%costs3 = xlsread('results\case39_30PV\resilience_noLinesOut.csv');

%costs1 = xlsread('results\experiments\1\res_out_case39.csv');
%costs2 = xlsread('results\experiments\1\res_out_case39_05PV.csv');
%costs3 = xlsread('results\experiments\1\res_out_case39_30PV.csv');

% costs1 = xlsread('results\experiments\1\res_out_case6ww.csv');
% costs2 = xlsread('results\experiments\1\res_out_case6ww_05PV.csv');
% costs3 = xlsread('results\experiments\1\res_out_case6ww_05PV_RS.csv');

% costs1 = xlsread('results\experiments\2\res_out_case39_orig_bias_samp1.csv');
% costs2 = xlsread('results\experiments\2\res_out_case39_allow_0out.csv');
% costs3 = xlsread('results\experiments\2\res_out_case39_fair_sample_no0.csv');

% costs1 = xlsread('results\experiments\1\res_out_case6ww_A0O.csv');
% costs2 = xlsread('results\experiments\1\res_out_case6ww_05PV_A0O.csv');
% costs3 = xlsread('results\experiments\1\res_out_case6ww_100PV_A0O.csv');
% costs3 = xlsread('results\experiments\1\res_out_case6ww_05PV_RS_A0O.csv');

% costs1 = xlsread('results\experiments\3\res_out_case6ww_A0O_1.csv');
% costs2 = xlsread('results\experiments\3\res_out_case6ww_05PV_v2_A0O_1.csv');
% costs3 = xlsread('results\experiments\3\res_out_case6ww_20PV_A0O_1.csv');
% costs4 = xlsread('results\experiments\3\res_out_case6ww_100PV_A0O_1.csv');

% costs1 = xlsread('results\experiments\4\res_out_case39_A0O_1.csv');
% costs2 = xlsread('results\experiments\4\res_out_case39_05PV_A0O_1.csv');
% costs3 = xlsread('results\experiments\4\res_out_case39_20PV_A0O_1.csv');
% costs4 = xlsread('results\experiments\4\res_out_case39_100PV_A0O_1.csv');

% costs1 = xlsread('results\experiments\5\res_out_case39_A0O_2.csv');
% costs2 = xlsread('results\experiments\5\res_out_case39_05PV_A0O_2.csv');
% costs3 = xlsread('results\experiments\5\res_out_case39_20PV_A0O_2.csv');
% costs4 = xlsread('results\experiments\5\res_out_case39_20PV_A0O_2.csv');

costs1 = xlsread('results\experiments\6\res_out_case39_p1.csv');
costs2 = xlsread('results\experiments\6\res_out_case39_05PV_p1.csv');
costs3 = xlsread('results\experiments\6\res_out_case39_20PV_p1.csv');
costs4 = xlsread('results\experiments\6\res_out_case39_100PV_p1.csv');

% costs1 = xlsread('results\experiments\5\res_out_case39_A0O_3.csv');
% costs2 = xlsread('results\experiments\5\res_out_case39_05PV_A0O_3.csv');
% costs3 = xlsread('results\experiments\5\res_out_case39_20PV_A0O_3.csv');
% costs4 = xlsread('results\experiments\5\res_out_case39_20PV_A0O_3.csv');

% costs1 = xlsread('results\experiments\0\res_out_case39_0_out.csv');
% costs2 = xlsread('results\experiments\0\res_out_case39_fair_sample_no_0_out.csv');
% costs3 = xlsread('results\experiments\0\res_out_case39_0_out_p1.csv');

% costs1 = xlsread('results\experiments\0\res_out_case6ww_0_out2.csv');
% costs2 = xlsread('results\experiments\0\res_out_case6ww_fair_sample_no_0_out2.csv');
% costs3 = xlsread('results\experiments\0\res_out_case6ww_0_out_p1.csv');

costs1(isnan(costs1))=0;
costs2(isnan(costs2))=0;
costs3(isnan(costs3))=0;
costs4(isnan(costs4))=0;

% make ccdf plost
N = length(costs1);
Pr=(N:-1:1)/N;
sorted_costs1 = sort(costs1);
sorted_costs2 = sort(costs2);
sorted_costs3 = sort(costs3);
sorted_costs4 = sort(costs4);

for jj = 1:N
    if sorted_costs1(jj)<0.001 || isnan(sorted_costs1(jj))
    sorted_costs1(jj) = 0;
    end
    if sorted_costs2(jj)<0.001 || isnan(sorted_costs2(jj))
        sorted_costs2(jj) = 0;
    end
    if sorted_costs3(jj)<0.001 || isnan(sorted_costs3(jj))
        sorted_costs3(jj) = 0;
    end
    if sorted_costs4(jj)<0.001 || isnan(sorted_costs4(jj))
        sorted_costs4(jj) = 0;
    end
end
numbins = 10;
cost = sorted_costs3(sorted_costs3~=0);
c1 = log10(min(cost));
c2 = log10(max(costs1));
bins = logspace(c1,c2,numbins);
P1 = zeros(numbins,1);
P1(1) = sum(sorted_costs1 <= bins(1));
P2 = zeros(numbins,1);
P2(1) = sum(sorted_costs2 <= bins(1));
P3 = zeros(numbins,1);
P3(1) = sum(sorted_costs3 <= bins(1));
P4 = zeros(numbins,1);
P4(1) = sum(sorted_costs4 <= bins(1));
for n = 2:numbins
    P1(n) = sum(sorted_costs1 <= bins(n))-sum(P1(1:n-1));
    P2(n) = sum(sorted_costs2 <= bins(n))-sum(P2(1:n-1));
    P3(n) = sum(sorted_costs3 <= bins(n))-sum(P3(1:n-1));
    P4(n) = sum(sorted_costs4 <= bins(n))-sum(P4(1:n-1));
end
% for n = 1:numbins
%     P1(n) = sum(sorted_costs1 >= bins(n))-sum(P1(1:n-1));
%     P2(n) = sum(sorted_costs2 >= bins(n))-sum(P2(1:n-1));
%     P3(n) = sum(sorted_costs3 >= bins(n))-sum(P3(1:n-1));
%     P4(n) = sum(sorted_costs4 >= bins(n))-sum(P4(1:n-1));
% end
P1 = P1./10000;
P2 = P2./10000;
P3 = P3./10000;
% P4 = P4./10000;
figure
% subplot(1,2,1)
% loglog(bins(2:end),P1(2:end))
% hold on
% loglog(bins(2:end),P2(2:end))
% loglog(bins(2:end),P3(2:end))
loglog(bins,P1)
hold on
loglog(bins,P2)
loglog(bins,P3)
% loglog(bins,P4)
title("Resilience measures of 39 bus cases")
% legend("original", "original2")
legend("original", "+5% load in DG", "+20% load in DG")%, "39 bus +100% load in DG")
% legend("39 bus", "39 bus 5% PV", "39 bus 20% PV", "39 bus 100% PV")%, "39 bus 100% PV")
% legend("6 bus", "6 bus 5% PV", "6 bus 20% PV")
% legend("6 bus", "6 bus 5% PV", "6 bus 20% PV", "6 bus 100% PV")
% legend("6 bus", "6 bus 5% PV", "6 bus 5% RL");
% legend("6 bus - allow 0 outages", "6 bus - fair sample", "6 bus - allow 0 outages and add 1");
% legend("39 bus - allow 0 outages", "39 bus - fair sample", "39 bus - allow 0 outages and add 1"); %"39 bus - bias sample 1 line out", 
xlabel("ENS (MWh)"); ylabel("Prob(ENS)")
set(gca, 'fontsize',13)
hold off