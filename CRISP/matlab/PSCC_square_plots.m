% Original 39 bus p1 sampling.
data1 = xlsread('..\VACC\results\experiments\mh\cascade_Set\RiskResults_case73_noPWS_lx2_n-1.csv');
data11 = xlsread('..\VACC\results\experiments\mh\cascade_Set\RiskResults_case73_noPWS_lx2_n-1+S5.csv');
data21 = xlsread('..\VACC\results\experiments\mh\cascade_Set\RiskResults_case73_noPWS_lx2_n-1+S20.csv');
% N-1 Secure 39 bus p1 sampling.
data2 = xlsread('..\VACC\results\experiments\mh\cascade_Set\RiskResults_case73_noPWS_lx2_n-1+PV5.csv');
data3 = xlsread('..\VACC\results\experiments\mh\cascade_Set\RiskResults_case73_noPWS_lx2_n-1+PV20.csv');
data31 = xlsread('..\VACC\results\experiments\mh\cascade_Set\RiskResults_case73_noPWS_lx2_n-1+PV5+S20.csv');
data41 = xlsread('..\VACC\results\experiments\mh\cascade_Set\RiskResults_case73_noPWS_lx2_n-1+PV20+S20.csv');

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


% make a 4x4 grid of square plots, with four squares in each
data = d1; %rand(6,4)

figure(1); clf;
arrangement = [2 4]
square_plots(data,arrangement)

