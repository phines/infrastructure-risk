% simple script to plot histograms when in a hurry

costs = xlsread('results\case39\resilience_case39.csv')
figure
hist(costs)
hold on
title("Resilience Distribution 39 bus");
xlabel("cost"); ylabel("number of events")
hold off

% make ccdf plost
N = length(costs);
Pr=(N:-1:1)/N;
sorted_costs = sort(costs);
semilogx(sorted_costs,Pr)
hold on
title("Resilience Distribution 39 bus");
xlabel("cost"); ylabel("CCDF")
hold off