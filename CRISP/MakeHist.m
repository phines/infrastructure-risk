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
for jj = 1:N
    if sorted_costs(jj)<0.01
    sorted_costs(jj) = 0;
    end
end

semilogx(sorted_costs,Pr)
hold on
title("Resilience Distribution 39 bus");
xlabel("C"); ylabel("Prob(Cost \geq C)")
hold off