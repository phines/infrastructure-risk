% Resilience figures for MEPS conference paper
% Molly Kelly-Gorham 5/21/2019

% Original 39 bus p1 sampling.
% costs1 = xlsread('..\results\experiments\6\res_out_case39_p1_2.csv');
% costs2 = xlsread('..\results\experiments\6\res_out_case39_05PV_p1_2.csv');
% costs3 = xlsread('..\results\experiments\6\res_out_case39_20PV_p1_2.csv');

% N-1 Secure 73 bus original model sampling.
% costs1 = xlsread('..\results\100\case73\data\res_out_case39_n-1_p1.csv');
% costs2 = xlsread('..\results\100\case73_load1.5\res_out_case73_n-1_p1.csv');
% costs3 = xlsread('..\results\100\case73_load2\res_out_case73_n-1_p1.csv');

% N-1 Secure 73 bus case generator outages
costs1 = xlsread('..\VACC\results\experiments\mh\resilience_case73_n-1_3.csv');

costs1(isnan(costs1))=0;
costs1 = -1*costs1;
% costs2(isnan(costs2))=0;
% costs3(isnan(costs3))=0;

% find the pdf
[x1,pdf1] = emperical_pdf(costs1,80);
% [x2,pdf2] = emperical_pdf(costs2,80);
% [x3,pdf3] = emperical_pdf(costs3,80);

% make ccdf
N = length(costs1);
% Nn = length(costs2);
Pr=(N:-1:1)/N;
% Pr2=(Nn:-1:1)/Nn;
sorted_costs1 = sort(costs1);
% sorted_costs2 = sort(costs2);
% sorted_costs3 = sort(costs3);
for jj = 1:N
    if sorted_costs1(jj)<0.001
    sorted_costs1(jj) = 0;
    end
end
% for jj = 1:Nn
%     if sorted_costs2(jj)<0.001 || isnan(sorted_costs2(jj))
%         sorted_costs2(jj) = 0;
%     end
%     if sorted_costs3(jj)<0.001 || isnan(sorted_costs3(jj))
%         sorted_costs3(jj) = 0;
%     end
% end


loglog(sorted_costs1,Pr,'b--')
hold on
% loglog(sorted_costs2,Pr2,'r--')
% loglog(sorted_costs3,Pr2,'--', 'color',[1.000    0.7000    0.0000])
% legend("original", "1.5*load", "2.0*load")
% legend("original", "+5% load in DG", "+20% load in DG","n-1 secure", "n-1 +5% load in DG", "n-1 +20% load in DG")
ylabel("Prob(Energy lost \geq ENS)")
set(gca, 'fontsize',13)
%xlim([1 10^7.5])
box off
legend boxoff
% set(gca,'xticklabels',[])

% make pdf plots
% subplot(3,1,3)
% plot(x1,pdf1)
% hold on
% plot(x2,pdf2)
% plot(x3,pdf3)
% legend('original', '5% load added in DG', '20% load added in DG')
% %title('Resilience 39 bus PDF')
xlabel("ENS (MWh)"); 
% ylabel('PDF(ENS)');
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% set(gca, 'fontsize',13)
% xlim([1 10^7.5])
% box off
% legend boxoff
% hold off
hold off