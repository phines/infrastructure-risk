% Resilience figures for MEPS conference paper
% Molly Kelly-Gorham 5/21/2019

% Original 39 bus p1 sampling.
% costs1 = xlsread('..\results\experiments\6\res_out_case39_p1_2.csv');
% costs2 = xlsread('..\results\experiments\6\res_out_case39_05PV_p1_2.csv');
% costs3 = xlsread('..\results\experiments\6\res_out_case39_20PV_p1_2.csv');

% N-1 Secure 39 bus p1 sampling.
costs1 = xlsread('..\results\experiments\8\res_out_case39_n-1.csv');
costs2 = xlsread('..\results\experiments\8\res_out_case39_n-1_05PV.csv');
costs3 = xlsread('..\results\experiments\8\res_out_case39_n-1_20PV.csv');

costs1(isnan(costs1))=0;
costs2(isnan(costs2))=0;
costs3(isnan(costs3))=0;

% find the pdf
[x1,pdf1] = emperical_pdf(costs1,80);
[x2,pdf2] = emperical_pdf(costs2,80);
[x3,pdf3] = emperical_pdf(costs3,80);

% make ccdf
N = length(costs1);
Pr=(N:-1:1)/N;
sorted_costs1 = sort(costs1);
sorted_costs2 = sort(costs2);
sorted_costs3 = sort(costs3);
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
end

figure(2)
subplot(3,1,1)
semilogx(sorted_costs1,Pr)
hold on
% title("CCDF of Resilience of N-1 secure 39 bus cases")
semilogx(sorted_costs2,Pr)
semilogx(sorted_costs3,Pr)
legend("original", "+5% load in DG", "+20% load in DG")
ylabel("Prob(Energy lost \geq ENS)")
set(gca, 'fontsize',13)
xlim([1 10^7.5])
box off
legend boxoff
set(gca,'xticklabels',[])
hold off
subplot(3,1,2)
loglog(sorted_costs1,Pr)
hold on
loglog(sorted_costs2,Pr)
loglog(sorted_costs3,Pr)
legend("original", "+5% load in DG", "+20% load in DG")
ylabel("Prob(Energy lost \geq ENS)")
set(gca, 'fontsize',13)
xlim([1 10^7.5])
box off
legend boxoff
set(gca,'xticklabels',[])
hold off
% make plots
subplot(3,1,3)
plot(x1,pdf1)
hold on
plot(x2,pdf2)
plot(x3,pdf3)
legend('original', '5% load added in DG', '20% load added in DG')
%title('Resilience 39 bus PDF')
xlabel("ENS (MWh)"); ylabel('PDF(ENS)');
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca, 'fontsize',13)
xlim([1 10^7.5])
box off
legend boxoff
hold off
