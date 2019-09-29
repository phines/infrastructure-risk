% Original 39 bus p1 sampling.
data11 = xlsread('..\VACC\results\experiments\mh\cascade_set\RiskResults_case73_noPWS_lx2_n-1.csv');
data21 = xlsread('..\results\experiments\9\res_out_case39_05PV_p1_lines.csv');
data31 = xlsread('..\results\experiments\9\res_out_case39_20PV_p1_lines.csv');

costs11 = xlsread('..\results\experiments\9\res_out_case39_p1.csv');
costs21 = xlsread('..\results\experiments\9\res_out_case39_05PV_p1.csv');
costs31 = xlsread('..\results\experiments\9\res_out_case39_20PV_p1.csv');

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

lines1 = data1(:,1);
lines2 = data2(:,1);
lines3 = data3(:,1);

time1 = data1(:,2);
time2 = data2(:,2);
time3 = data3(:,2);

time11 = data11(:,2);
time21 = data21(:,2);
time31 = data31(:,2);

ls1 = data1(:,3);
ls2 = data2(:,3);
ls3 = data3(:,3);

ls11 = data11(:,3);
ls21 = data21(:,3);
ls31 = data31(:,3);

costs1(isnan(costs1))=0;
costs2(isnan(costs2))=0;
costs3(isnan(costs3))=0;

costs11(isnan(costs11))=0;
costs21(isnan(costs21))=0;
costs31(isnan(costs31))=0;

% find eens for different sized events
short = [0,60*24];
long =  [60*24,max([time1; time2; time3])];
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

d1(4,1) = sort_and_sum_data(short,time11,smallMW,ls11,costs11);
d1(4,2) = sort_and_sum_data(short,time11,largeMW,ls11,costs11);
d1(4,3) = sort_and_sum_data(long,time11,smallMW,ls11,costs11);
d1(4,4) = sort_and_sum_data(long,time11,largeMW,ls11,costs11);

d1(5,1) = sort_and_sum_data(short, time21, smallMW,ls21, costs21);
d1(5,2) = sort_and_sum_data(short,time21,largeMW,ls21,costs21);
d1(5,3) = sort_and_sum_data(long,time21,smallMW,ls21,costs21);
d1(5,4) = sort_and_sum_data(long,time21,largeMW,ls21,costs21);

d1(6,1) = sort_and_sum_data(short, time31, smallMW, ls31, costs31);
d1(6,2) = sort_and_sum_data(short,time31,largeMW,ls31,costs31);
d1(6,3) = sort_and_sum_data(long,time31,smallMW,ls31,costs31);
d1(6,4) = sort_and_sum_data(long,time31,largeMW,ls31,costs31);


% make a 4x4 grid of square plots, with four squares in each
data = d1; %rand(6,4)

figure(1); clf;
arrangement = [2 3]
square_plots(data,arrangement)

