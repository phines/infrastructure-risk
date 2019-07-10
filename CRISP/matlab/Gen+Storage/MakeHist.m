% simple script to plot histograms when in a hurry
%case 6ww:
costs1 = xlsread('..\..\results\experiments_gen\6\res_out_case73_storage.csv');
costs2 = xlsread('..\..\results\experiments_gen\6\res_out_case73_noPWS.csv');


costs1(isnan(costs1))=0;
costs2(isnan(costs2))=0;
% costs3(isnan(costs3))=0;
% costs4(isnan(costs4))=0;

% make ccdf plost
N = length(costs1);
Pr=(N:-1:1)/N;
sorted_costs1 = sort(costs1);
sorted_costs2 = sort(costs2);
% sorted_costs3 = sort(costs3);
% sorted_costs4 = sort(costs4);
for jj = 1:N
    if sorted_costs1(jj)<0.001 || isnan(sorted_costs1(jj))
    sorted_costs1(jj) = 0;
    end
    if sorted_costs2(jj)<0.001 || isnan(sorted_costs2(jj))
        sorted_costs2(jj) = 0;
    end
%     if sorted_costs3(jj)<0.001 || isnan(sorted_costs3(jj))
%         sorted_costs3(jj) = 0;
%     end
%     if sorted_costs4(jj)<0.001 || isnan(sorted_costs4(jj))
%         sorted_costs4(jj) = 0;
%     end
end

figure(1)
% subplot(1,2,1)
% hold on
% loglog(sorted_costs1,Pr)
% hold on
% loglog(sorted_costs2,Pr)
% loglog(sorted_costs3,Pr)
semilogx(sorted_costs1,Pr)
hold on
semilogx(sorted_costs2,Pr)
% semilogx(sorted_costs3,Pr)
% semilogx(sorted_costs4,Pr)
title("CCDF of Resilience of N-1 secure 73 bus cases")
% legend("original", "original2")
% legend("original", "+5% load in DG", "+20% load in DG")%, "39 bus +100% load in DG")
% legend("39 bus", "39 bus 5% PV", "39 bus 20% PV", "39 bus 100% PV")%, "39 bus 100% PV")
% legend("6 bus", "6 bus 5% PV", "6 bus 20% PV")
% legend("6 bus", "6 bus 5% PV", "6 bus 20% PV", "6 bus 100% PV")
% legend("6 bus", "6 bus 5% PV", "6 bus 5% RL");
% legend("6 bus - allow 0 outages", "6 bus - fair sample", "6 bus - allow 0 outages and add 1");
% legend("39 bus - allow 0 outages", "39 bus - fair sample", "39 bus - allow 0 outages and add 1"); %"39 bus - bias sample 1 line out", 
xlabel("ENS (MWh)"); ylabel("Prob(Energy lost \geq ENS")
set(gca, 'fontsize',13)
hold off