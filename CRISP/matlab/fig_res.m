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

% N-1 Secure 73 bus case generator outages storage
% costs1 = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1.csv');
% costs2 = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1+S5.csv');
% costs3 = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1+S20.csv');
% costs = xlsread('..\VACC\results\experiments\mh\RiskResults_case73_noPWS_lx2_n-1.csv');
% costs1 = costs(:,1)
% costs = xlsread('..\VACC\results\experiments\mh\RiskResults_case73_noPWS_lx2_n-1+S5.csv');
% costs2 = costs(:,1)
% costs = xlsread('..\VACC\results\experiments\mh\RiskResults_case73_noPWS_lx2_n-1+S20.csv');
% costs3 = costs(:,1)

% N-1 Secure 73 bus case generator outages DG
% costs1 = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1.csv');
% costs2 = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1+PV5.csv');
% costs3 = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1+PV20.csv');
% costs = xlsread('..\VACC\results\experiments\mh\RiskResults_case73_noPWS_lx2_n-1.csv');
% costs1 = costs(:,1)
% costs = xlsread('..\VACC\results\experiments\mh\RiskResults_case73_noPWS_lx2_n-1+PV5.csv');
% costs2 = costs(:,1)
% costs = xlsread('..\VACC\results\experiments\mh\RiskResults_case73_noPWS_lx2_n-1+PV20.csv');
% costs3 = costs(:,1)

% N-1 Secure 73 bus case generator outages DG 5 and storage
% costs1 = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1+PV20.csv');
% costs2 = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1+PV20+S5.csv');
% costs3 = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1+PV20+S20.csv');

% N-1 Secure 73 bus case generator outages storage 5 and DG
% costs = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1.csv');
% costs1 = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1+S20.csv');
% costs2 = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1+PV5+S20.csv');
% costs3 = xlsread('..\VACC\results\experiments\mh\cascade_set\res_73_noPWS_lx2_n-1+PV20+S20.csv');
% costs = xlsread('..\VACC\results\experiments\mh\RiskResults_case73_noPWS_lx2_n-1+S20.csv');
% costs1 = costs(:,1)
% costs = xlsread('..\VACC\results\experiments\mh\RiskResults_case73_noPWS_lx2_n-1+PV5+S20.csv');
% costs2 = costs(:,1)
% costs = xlsread('..\VACC\results\experiments\mh\RiskResults_case73_noPWS_lx2_n-1+PV20+S20.csv');
% costs3 = costs(:,1)


% N-1 Secure 73 bus case generator outages communication system increased
% restoration times
costs = xlsread('..\VACC\results\experiments\mh\communications_interactions\res_73_noPWS_lx2_n-1.csv');
costs1 = xlsread('..\VACC\results\experiments\mh\communications_interactions\res_73_noPWS_lx2_n-1+S20.csv');
costs2 = xlsread('..\VACC\results\experiments\mh\communications_interactions\res_73_noPWS_lx2_n-1+PV5+S20.csv');
costs3 = xlsread('..\VACC\results\experiments\mh\communications_interactions\res_73_noPWS_lx2_n-1+PV20+S20.csv');

costs(isnan(costs))=0;
costs1(isnan(costs1))=0;
costs2(isnan(costs2))=0;
costs3(isnan(costs3))=0;

% find the pdf
[x1,pdf1] = emperical_pdf(costs1,80);
[x2,pdf2] = emperical_pdf(costs2,80);
[x3,pdf3] = emperical_pdf(costs3,80);

% make ccdf
N0 = length(costs);
N = length(costs1);
N2 = length(costs2);
N3 = length(costs3);
Pr0 = (N0:-1:1)/N0;
Pr=(N:-1:1)/N;
Pr2=(N2:-1:1)/N2;
Pr3=(N3:-1:1)/N3;
sorted_costs1 = sort(costs1);
sorted_costs2 = sort(costs2);
sorted_costs3 = sort(costs3);
sorted_costs0 = sort(costs);
for jj = 1:N0
    if sorted_costs0(jj)<0.001 || isnan(sorted_costs0(jj))
    sorted_costs0(jj) = 0;
    end
end
for jj = 1:N
    if sorted_costs1(jj)<0.001 || isnan(sorted_costs1(jj))
    sorted_costs1(jj) = 0;
    end
end
for jj = 1:N2
    if sorted_costs2(jj)<0.001 || isnan(sorted_costs2(jj))
        sorted_costs2(jj) = 0;
    end
end
for jj = 1:N3
    if sorted_costs3(jj)<0.001 || isnan(sorted_costs3(jj))
        sorted_costs3(jj) = 0;
    end
end
Pr0 = Pr0(sorted_costs0 ~= 0) 
sorted_costs0 = sorted_costs0(sorted_costs0 ~= 0)
Pr = Pr(sorted_costs1 ~= 0)
sorted_costs1 = sorted_costs1(sorted_costs1 ~= 0)
Pr2 = Pr2(sorted_costs2 ~= 0)
sorted_costs2 = sorted_costs2(sorted_costs2 ~= 0)
Pr3 = Pr3(sorted_costs3 ~= 0)
sorted_costs3 = sorted_costs3(sorted_costs3 ~= 0)
% loglog([10^(-1); sorted_costs1],[Pr(1) Pr],'b')
loglog([10^(-1); sorted_costs0],[Pr0(1) Pr0],'b')
hold on
loglog([10^(-1); sorted_costs1],[Pr(1) Pr],'color',[1.000    0.7000    0.0000])
loglog([10^(-1); sorted_costs2],[Pr2(1) Pr2],'color',[0.1 1.0 0.1]) %rgb(238,130,238) = violet
loglog([10^(-1); sorted_costs3],[Pr3(1) Pr3],'color',[0.5    0.50    0.8000]) %rgb(255,69,0) orange-red
% loglog([10^(-1); sorted_costs2],[Pr2(1) Pr2],'c')
% loglog([10^(-1); sorted_costs3],[Pr3(1) Pr3], 'm')
% loglog([10^(-1); sorted_costs2],[Pr2(1) Pr2],'r')
% loglog([10^(-1); sorted_costs3],[Pr3(1) Pr3], 'color',[1.000    0.7000    0.0000])
% legend("original", "+5% load in DG", "+20% load in DG","n-1 secure", "n-1 +5% load in DG", "n-1 +20% load in DG")
% legend("base case", "+5% PV", "+20% PV")
% legend("base case", "+5% storage", "+20% storage")
legend("base case","+20% storage","+20% storage +5% PV", " +20% storage +20% PV")
ylabel("Prob(ENS \geq x)")
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
xlabel("x (MWh)"); 
% ylabel('PDF(ENS)');
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca, 'fontsize',13)
xlim([1 10^5])
box off
legend boxoff
% hold off
hold off