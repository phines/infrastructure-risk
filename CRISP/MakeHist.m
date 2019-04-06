% simple script to plot histograms when in a hurry

costs1 = xlsread('results\case39\test\res_out_test_008.csv'); %('results\case39\resilience_out_noLinesOut.csv');
costs2 = xlsread('results\case39\test\res_out_test_009.csv'); %('results\case39_05PV\resilience_noLinesOut.csv');
%costs3 = xlsread('results\case39_30PV\resilience_noLinesOut.csv');

%costs1 = xlsread('results\experiments\1\res_out_case39.csv'); %('results\case39\resilience_out_noLinesOut.csv');
%costs2 = xlsread('results\experiments\1\res_out_case39_05PV.csv'); %('results\case39_05PV\resilience_noLinesOut.csv');

% make ccdf plost
N = length(costs1);
Pr=(N:-1:1)/N;
sorted_costs1 = sort(costs1);
sorted_costs2 = sort(costs2);
%sorted_costs3 = sort(costs3);
for jj = 1:N
    if sorted_costs1(jj)<0.001 || sorted_costs1(jj) ==NaN
    sorted_costs1(jj) = 0;
    end
    if sorted_costs2(jj)<0.001 || sorted_costs2(jj) ==NaN
        sorted_costs2(jj) = 0;
    end
    %if sorted_costs3(jj)<0.001 || sorted_costs3(jj) ==NaN
    %    sorted_costs3(jj) = 0;
    %end
end

figure
semilogx(sorted_costs1,Pr)
hold on
semilogx(sorted_costs2,Pr)
%semilogx(sorted_costs3,Pr)
title("Resilience")
legend("39 bus-old code", "39 bus-new code")%, "39 bus 30% PV");
xlabel("C"); ylabel("Prob(Cost \geq C)")
hold off