% simple script to plot histograms when in a hurry

costs1 = xlsread('results\case39\resilience_out3.csv');
costs2 = xlsread('results\case39_05PV\resilience_out1.csv');
costs3 = xlsread('results\case39_30PV\resilience_out1.csv');

% make ccdf plost
N = length(costs1);
Pr=(N:-1:1)/N;
sorted_costs1 = sort(costs1);
sorted_costs2 = sort(costs2);
sorted_costs3 = sort(costs3);
for jj = 1:N
    if sorted_costs1(jj)<0.001 || sorted_costs1(jj) ==NaN
    sorted_costs1(jj) = 0;
    end
    if sorted_costs2(jj)<0.001 || sorted_costs2(jj) ==NaN
        sorted_costs2(jj) = 0;
    end
    if sorted_costs3(jj)<0.001 || sorted_costs3(jj) ==NaN
        sorted_costs3(jj) = 0;
    end
end

figure
subplot(1,3,1)
semilogx(sorted_costs1,Pr)
hold on
title("Resilience, 39 bus");
xlabel("C"); ylabel("Prob(Cost \geq C)")
hold off
subplot(1,3,2)
semilogx(sorted_costs2,Pr)
hold on
title("Resilience, 39 bus 5% PV");
xlabel("C"); ylabel("Prob(Cost \geq C)")
hold off
subplot(1,3,3)
semilogx(sorted_costs3,Pr)
hold on
title("Resilience, 39 bus 30% PV");
xlabel("C"); ylabel("Prob(Cost \geq C)")
hold off