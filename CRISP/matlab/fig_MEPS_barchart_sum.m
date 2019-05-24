% Resilience bar charts for MEPS conference paper
% Molly Kelly-Gorham 5/21/2019

%clear

% % Original 39 bus p1 sampling.
% data1 = xlsread('..\results\experiments\9\res_out_case39_p1_lines.csv');
% data2 = xlsread('..\results\experiments\9\res_out_case39_05PV_p1_lines.csv');
% data3 = xlsread('..\results\experiments\9\res_out_case39_20PV_p1_lines.csv');
% 
% costs1 = xlsread('..\results\experiments\9\res_out_case39_p1.csv');
% costs2 = xlsread('..\results\experiments\9\res_out_case39_05PV_p1.csv');
% costs3 = xlsread('..\results\experiments\9\res_out_case39_20PV_p1.csv');

% N-1 Secure 39 bus p1 sampling.
data1 = xlsread('..\results\experiments\9\res_out_case39_n-1_p1_lines.csv');
data2 = xlsread('..\results\experiments\9\res_out_case39_n-1_05PV_p1_lines.csv');
data3 = xlsread('..\results\experiments\9\res_out_case39_n-1_20PV_p1_lines.csv');

costs1 = xlsread('..\results\experiments\9\res_out_case39_n-1_p1.csv');
costs2 = xlsread('..\results\experiments\9\res_out_case39_n-1_05PV_p1.csv');
costs3 = xlsread('..\results\experiments\9\res_out_case39_n-1_20PV_p1.csv');

lines1 = data1(:,1);
lines2 = data2(:,1);
lines3 = data3(:,1);


time1 = data1(:,2);
time2 = data2(:,2);
time3 = data3(:,2);

ls1 = data1(:,3);
ls2 = data2(:,3);
ls3 = data3(:,3);

costs1(isnan(costs1))=0;
costs2(isnan(costs2))=0;
costs3(isnan(costs3))=0;

% find eens for different sized events
short = [0,1000];
long =  [1000,max([time1; time2; time3])];
smallMW = [0,100];
largeMW = [100,max([costs1; costs2; costs3])];

d1(1,1) = sort_and_sum_data(short,time1,smallMW,ls1,costs1);
d1(1,2) = sort_and_sum_data(short,time1,largeMW,ls1,costs1);
d1(1,3) = sort_and_sum_data(long,time1,smallMW,ls1,costs1);
d1(1,4) = sort_and_sum_data(long,time1,largeMW,ls1,costs1);

d1(2,1) = sort_and_sum_data(short, time2, smallMW,ls2, costs2);
d1(2,2) = sort_and_sum_data(short,time2,largeMW,ls2,costs2);
d1(2,3) = sort_and_sum_data(long,time2,smallMW,ls2,costs2);
d1(2,4) = sort_and_sum_data(long,time2,largeMW,ls2,costs2);

d1(3,1) = sort_and_sum_data(short, time3, smallMW, ls3, costs3);
d1(3,2) = sort_and_sum_data(short,time3,largeMW,ls3,costs3);
d1(3,3) = sort_and_sum_data(long,time3,smallMW,ls3,costs3);
d1(3,4) = sort_and_sum_data(long,time3,largeMW,ls3,costs3);

figure
data2D = d1;
H=bar(data2D, 'stack');
P=findobj(gca,'type','patch');
myC= [0 0 1
    1 0 0
    1 0.4 0
    0 0.8 1
    0.6 0 1
    0 1 0 ];
for n= 1 : length(P) 
    set(P(n),'facecolor',myC(n,:));
end
a = ['< ' num2str(long(1)) ' min, <'  num2str(smallMW(2)) 'MW'];
b = ['< ' num2str(long(1)) ' min, >'  num2str(smallMW(2)) 'MW'];
c = ['> ' num2str(long(1)) ' min, <' num2str(smallMW(2)) 'MW'];
e = ['> ' num2str(long(1)) ' min, >'  num2str(smallMW(2)) 'MW'];
legend(a,b,c,e,'Location','eastoutside')
set(gca,'xticklabels',["+0% "; "+5% "; "+20%"])
%xlabel(['small                  ' 'medium                 ' 'large'])
%legend('< 50 min, < 10 MW','< 50 min, > 10 MW','> 50 min, < 10 MW','> 50 min, > 10 MW','Location','eastoutside')
%AX=legend(a,b,c,e);
%LEG = findobj(AX,'type','text');
%set(LEG,'FontSize',8);
%set(gca,'yscale','log');
ylabel('1/n \Sigma ENS (MWh)')
set(gca, 'fontsize',13)
legend boxoff
box off


% % find the pdf
% [x1,pdf1] = emperical_pdf(costs1,80);
% [x2,pdf2] = emperical_pdf(costs2,80);
% [x3,pdf3] = emperical_pdf(costs3,80);
% 
% % make ccdf
% N = length(costs1);
% Pr=(N:-1:1)/N;
% sorted_costs1 = sort(costs1);
% sorted_costs2 = sort(costs2);
% sorted_costs3 = sort(costs3);
% for jj = 1:N
%     if sorted_costs1(jj)<0.001 || isnan(sorted_costs1(jj))
%     sorted_costs1(jj) = 0;
%     end
%     if sorted_costs2(jj)<0.001 || isnan(sorted_costs2(jj))
%         sorted_costs2(jj) = 0;
%     end
%     if sorted_costs3(jj)<0.001 || isnan(sorted_costs3(jj))
%         sorted_costs3(jj) = 0;
%     end
% end
% 
% figure
% semilogx(sorted_costs1,Pr)
% hold on
% % title("CCDF of Resilience of N-1 secure 39 bus cases")
% semilogx(sorted_costs2,Pr)
% semilogx(sorted_costs3,Pr)
% legend("original", "+5% load in DG", "+20% load in DG")
% ylabel("Prob(Energy lost \geq ENS)")
% set(gca, 'fontsize',13)
% xlim([1 10^7.5])
% box off
% legend boxoff
% set(gca,'xticklabels',[])
% hold off
