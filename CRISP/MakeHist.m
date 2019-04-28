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

costs1 = xlsread('results\experiments\5\res_out_case39_A0O_2.csv');
% costs2 = xlsread('results\experiments\5\res_out_case39_A0O_3.csv');
costs2 = xlsread('results\experiments\5\res_out_case39_05PV_A0O_2.csv');
costs3 = xlsread('results\experiments\5\res_out_case39_20PV_A0O_2.csv');
costs4 = xlsread('results\experiments\5\res_out_case39_20PV_A0O_2.csv');


% costs1 = xlsread('results\experiments\5\res_out_case39_A0O_3.csv');
% costs2 = xlsread('results\experiments\5\res_out_case39_05PV_A0O_3.csv');
% costs3 = xlsread('results\experiments\5\res_out_case39_20PV_A0O_3.csv');
% costs4 = xlsread('results\experiments\5\res_out_case39_20PV_A0O_3.csv');

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

figure(1)
subplot(1,2,1)
semilogx(sorted_costs1,Pr)
hold on
semilogx(sorted_costs2,Pr)
semilogx(sorted_costs3,Pr)
% semilogx(sorted_costs4,Pr)
title("Resilience 39 bus - 1")
% legend("original", "original2")
legend("original", "+5% load in DG", "+20% load in DG")%, "39 bus +100% load in DG")
% legend("39 bus", "39 bus 5% PV", "39 bus 20% PV", "39 bus 100% PV")%, "39 bus 100% PV")
% legend("6 bus", "6 bus 5% PV", "6 bus 20% PV")
% legend("6 bus", "6 bus 5% PV", "6 bus 20% PV", "6 bus 100% PV")
% legend("6 bus", "6 bus 5% PV", "6 bus 5% RL");
%legend("39 bus - bias sample 1 line out", "39 bus - allow 0 outages", "39 bus - fair sample");
xlabel("EL (MWh)"); ylabel("Prob(Energy lost \geq EL)")
set(gca, 'fontsize',13)
hold off