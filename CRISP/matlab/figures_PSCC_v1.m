% Resilience figures for PSCC conference paper
% Molly Kelly-Gorham 12/18/2019

% N-1 Secure 73 bus original model sampling.
%% Added PV
costs = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1.csv');
costs = costs(:,1)
costs1 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1.csv');
costs1 = costs1(:,1);
costs2 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+PV5.csv');
costs2 = costs2(:,1);
costs3 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+PV20.csv');
costs3 = costs3(:,1);

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
EENS_0 = sum(costs)/N0;
EENS_1 = sum(costs1)/N;
EENS_2 = sum(costs2)/N2;
EENS_3 = sum(costs3)/N3;
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
figure
loglog([10^(-1); sorted_costs1],[Pr(1) Pr],'b')
hold on
loglog([10^(-1); sorted_costs2],[Pr2(1) Pr2],'c')
loglog([10^(-1); sorted_costs3],[Pr3(1) Pr3], 'm')
legend("base case", "+5% PV", "+20% PV")
xlabel("x")
ylabel("Prob(ENS \geq x)")
set(gca, 'fontsize',13)
%xlim([1 10^7.5])
box off
legend boxoff
% set(gca,'xticklabels',[])


%% Added Storage

costs1 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1.csv');
costs1 = costs1(:,1);
costs2 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+S5.csv')
costs2 = costs2(:,1);
costs3 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+S20.csv');
costs3 = costs3(:,1);

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
EENS_0 = sum(costs)/N0;
EENS_1 = sum(costs1)/N;
EENS_2 = sum(costs2)/N2;
EENS_3 = sum(costs3)/N3;
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
figure
loglog([10^(-1); sorted_costs1],[Pr(1) Pr],'b')
hold on
loglog([10^(-1); sorted_costs2],[Pr2(1) Pr2],'r')
loglog([10^(-1); sorted_costs3],[Pr3(1) Pr3], 'color',[1.000    0.7000    0.0000])

legend("base case", "+5% storage", "+20% storage")
xlabel("x")
ylabel("Prob(ENS \geq x)")
set(gca, 'fontsize',13)
box off
legend boxoff

%% PV and Storage% Added PV
costs = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1.csv');
costs = costs(:,1)
costs1 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+S20.csv');
costs1 = costs1(:,1);
costs2 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+PV5+S20.csv');
costs2 = costs2(:,1);
costs3 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+PV20+S20.csv');
costs3 = costs3(:,1);

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
EENS_0 = sum(costs)/N0;
EENS_1 = sum(costs1)/N;
EENS_2 = sum(costs2)/N2;
EENS_3 = sum(costs3)/N3;
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
figure
loglog([10^(-1); sorted_costs0],[Pr0(1) Pr0],'b')%v-')
hold on
loglog([10^(-1); sorted_costs1],[Pr(1) Pr],'color',[1.000 0.7000 0.0000])%'s-',
loglog([10^(-1); sorted_costs2],[Pr2(1) Pr2],'color',[0.1 1.0 0.1])%'x-', %rgb(238,130,238) = violet
loglog([10^(-1); sorted_costs3],[Pr3(1) Pr3],'color',[0.5 0.50 0.8000])%'d-', %rgb(255,69,0) orange-red
legend("base case","+20% storage","+20% storage +5% PV", " +20% storage +20% PV")
xlabel("x")
ylabel("Prob(ENS \geq x)")
set(gca, 'fontsize',13)
box off
legend boxoff

%% Make Square plots
data1 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1.csv');
data11 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+S5.csv');
data21 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+S20.csv');
data2 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+PV5.csv');
data3 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+PV20.csv');
data31 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+PV5+S20.csv');
data41 = xlsread('..\VACC\results\experiments\mh\casc2_1\RiskResults_case73_noPWS_lx2_n-1+PV20+S20.csv');

costs11 = data11(:,1);
costs21 = data21(:,1);
costs31 = data31(:,1);
costs41 = data41(:,1);
costs1 = data1(:,1);
costs2 = data2(:,1);
costs3 = data3(:,1);

time1 = data1(:,3);
time2 = data2(:,3);
time3 = data3(:,3);

time11 = data11(:,3);
time21 = data21(:,3);
time31 = data31(:,3);
time41 = data41(:,3);

ls1 = data1(:,2);
ls2 = data2(:,2);
ls3 = data3(:,2);

ls11 = data11(:,2);
ls21 = data21(:,2);
ls31 = data31(:,2);
ls41 = data41(:,2);

costs1(isnan(costs1))=0;
costs2(isnan(costs2))=0;
costs3(isnan(costs3))=0;

costs11(isnan(costs11))=0;
costs21(isnan(costs21))=0;
costs31(isnan(costs31))=0;
costs41(isnan(costs41))=0;

% find eens for different sized events
short = [0,60*48];
long =  [60*48,max([time1; time2; time3])];
smallMW = [0,250];
largeMW = [250,max([costs1; costs2; costs3])];

d1(1,1) = sort_and_sum_data(short,time11,smallMW,ls11,costs11);
d1(1,2) = sort_and_sum_data(short,time11,largeMW,ls11,costs11);
d1(1,3) = sort_and_sum_data(long,time11,smallMW,ls11,costs11);
d1(1,4) = sort_and_sum_data(long,time11,largeMW,ls11,costs11);

d1(2,1) = sort_and_sum_data(short, time21, smallMW,ls21, costs21);
d1(2,2) = sort_and_sum_data(short,time21,largeMW,ls21,costs21);
d1(2,3) = sort_and_sum_data(long,time21,smallMW,ls21,costs21);
d1(2,4) = sort_and_sum_data(long,time21,largeMW,ls21,costs21);

d1(3,1) = sort_and_sum_data(short, time31, smallMW, ls31, costs31);
d1(3,2) = sort_and_sum_data(short,time31,largeMW,ls31,costs31);
d1(3,3) = sort_and_sum_data(long,time31,smallMW,ls31,costs31);
d1(3,4) = sort_and_sum_data(long,time31,largeMW,ls31,costs31);

d1(4,1) = sort_and_sum_data(short,time41,smallMW,ls41,costs41);
d1(4,2) = sort_and_sum_data(short,time41,largeMW,ls41,costs41);
d1(4,3) = sort_and_sum_data(long,time41,smallMW,ls41,costs41);
d1(4,4) = sort_and_sum_data(long,time41,largeMW,ls41,costs41);

d1(5,1) = 0;
d1(5,2) = 0;
d1(5,3) = 0;
d1(5,4) = 0;

d1(6,1) = sort_and_sum_data(short,time1,smallMW,ls1,costs1);
d1(6,2) = sort_and_sum_data(short,time1,largeMW,ls1,costs1);
d1(6,3) = sort_and_sum_data(long,time1,smallMW,ls1,costs1);
d1(6,4) = sort_and_sum_data(long,time1,largeMW,ls1,costs1);

d1(7,1) = sort_and_sum_data(short, time2, smallMW,ls2, costs2);
d1(7,2) = sort_and_sum_data(short,time2,largeMW,ls2,costs2);
d1(7,3) = sort_and_sum_data(long,time2,smallMW,ls2,costs2);
d1(7,4) = sort_and_sum_data(long,time2,largeMW,ls2,costs2);

d1(8,1) = sort_and_sum_data(short, time3, smallMW, ls3, costs3);
d1(8,2) = sort_and_sum_data(short,time3,largeMW,ls3,costs3);
d1(8,3) = sort_and_sum_data(long,time3,smallMW,ls3,costs3);
d1(8,4) = sort_and_sum_data(long,time3,largeMW,ls3,costs3);


% make a 6x4 grid of square plots, with six squares in each
data = d1; %rand(6,4)

figure(4); clf;
arrangement = [2 4]
square_plots(data,arrangement)