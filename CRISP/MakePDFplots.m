% not accurate! use code in Matlab folder
clear all
close all

%% Load data

%case 6ww:
%costs1 = xlsread('results\case6ww\test\res_out_test_001.csv');
%costs2 = xlsread('results\case6ww\test\res_out_test_002.csv');

%case39:
%figuring out weird change in distribution
%costs1 = xlsread('results\case39\resilience_case39.csv');
%costs2 = xlsread('results\case39\resilience_case39_4.csv');
%costs3 = xlsread('results\case39\resilience_case39_5.csv');

%costs1 = xlsread('results\case39\test\res_out_test_010.csv'); 
%costs2 = xlsread('results\case39\test\res_out_test_011.csv'); 

%Compare - no lines out
%costs1 = xlsread('results\case39\resilience_out_noLinesOut.csv');
%costs2 = xlsread('results\case39_05PV\resilience_noLinesOut.csv');
%costs3 = xlsread('results\case39_30PV\resilience_noLinesOut.csv');

%costs1 = xlsread('results\experiments\1\res_out_case39.csv');
%costs2 = xlsread('results\experiments\1\res_out_case39_05PV.csv');
%costs3 = xlsread('results\experiments\1\res_out_case39_30PV.csv');

% costs1 = xlsread('results\experiments\1\res_out_case6ww.csv');
% costs2 = xlsread('results\experiments\1\res_out_case6ww_05PV.csv');
% costs3 = xlsread('results\experiments\1\res_out_case6ww_05PV_RS.csv');

% costs1 = xlsread('results\experiments\2\res_out_case39_orig_bias_samp1.csv');
% costs2 = xlsread('results\experiments\2\res_out_case39_allow_0out.csv');
% costs3 = xlsread('results\experiments\2\res_out_case39_fair_sample_no0.csv');

% costs1 = xlsread('results\experiments\1\res_out_case6ww_A0O.csv');
% costs2 = xlsread('results\experiments\1\res_out_case6ww_05PV_A0O.csv');
% costs3 = xlsread('results\experiments\1\res_out_case6ww_100PV_A0O.csv');
% costs3 = xlsread('results\experiments\1\res_out_case6ww_05PV_RS_A0O.csv');

% costs1 = xlsread('results\experiments\3\res_out_case6ww_A0O_1.csv');
% costs2 = xlsread('results\experiments\3\res_out_case6ww_05PV_v2_A0O_1.csv');
% costs3 = xlsread('results\experiments\3\res_out_case6ww_20PV_A0O_1.csv');
% costs4 = xlsread('results\experiments\3\res_out_case6ww_100PV_A0O_1.csv');

% costs1 = xlsread('results\experiments\4\res_out_case39_A0O_1.csv');
% costs2 = xlsread('results\experiments\4\res_out_case39_05PV_A0O_1.csv');
% costs3 = xlsread('results\experiments\4\res_out_case39_20PV_A0O_1.csv');
% costs4 = xlsread('results\experiments\4\res_out_case39_100PV_A0O_1.csv');

% costs1 = xlsread('results\experiments\5\res_out_case39_A0O_2.csv');
% costs2 = xlsread('results\experiments\5\res_out_case39_05PV_A0O_2.csv');
% costs3 = xlsread('results\experiments\5\res_out_case39_20PV_A0O_2.csv');
% costs4 = xlsread('results\experiments\5\res_out_case39_20PV_A0O_2.csv');

% costs1 = xlsread('results\experiments\6\res_out_case39_fair-sample.csv');
% costs2 = xlsread('results\experiments\6\res_out_case39_05PV_fair-sample.csv');
% costs3 = xlsread('results\experiments\6\res_out_case39_20PV_fair-sample.csv');
% costs4 = xlsread('results\experiments\6\res_out_case39_100PV_fair-sample.csv');

costs1 = xlsread('results\experiments\6\res_out_case39_p1_2.csv');
costs2 = xlsread('results\experiments\6\res_out_case39_05PV_p1_2.csv');
costs3 = xlsread('results\experiments\6\res_out_case39_20PV_p1_2.csv');
costs4 = xlsread('results\experiments\6\res_out_case39_100PV_p1_2.csv');

% costs1 = xlsread('results\experiments\5\res_out_case39_A0O_3.csv');
% costs2 = xlsread('results\experiments\5\res_out_case39_05PV_A0O_3.csv');
% costs3 = xlsread('results\experiments\5\res_out_case39_20PV_A0O_3.csv');
% costs4 = xlsread('results\experiments\5\res_out_case39_20PV_A0O_3.csv');

% costs1 = xlsread('results\experiments\0\res_out_case39_0_out.csv');
% costs2 = xlsread('results\experiments\0\res_out_case39_fair_sample_no_0_out.csv');
% costs3 = xlsread('results\experiments\0\res_out_case39_0_out_p1.csv');

% costs1 = xlsread('results\experiments\0\res_out_case6ww_0_out2.csv');
% costs2 = xlsread('results\experiments\0\res_out_case6ww_fair_sample_no_0_out2.csv');
% costs3 = xlsread('results\experiments\0\res_out_case6ww_0_out_p1.csv');

%% Clean data

costs1(isnan(costs1))=0;
costs2(isnan(costs2))=0;
costs3(isnan(costs3))=0;
costs4(isnan(costs4))=0;

% make ccdf plost
N = length(costs1);
Pr=(N:-1:1)/N;
sorted_costs1 = sort(costs1);
sorted_costs2 = sort(costs2);
sorted_costs3 = sort(costs3);
sorted_costs4 = sort(costs4);

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
    if sorted_costs4(jj)<0.001 || isnan(sorted_costs4(jj))
        sorted_costs4(jj) = 0;
    end
end


%% Separating out each case resilience data into while loops to measure probability function

numbins = 3;
cost = sorted_costs3(sorted_costs3~=0);
c1 = log10(min(cost));
c2 = log10(max(costs1));
bins = logspace(c1,c2,numbins);
P1 = zeros(numbins,1);
P1(1) = sum(sorted_costs1 <= bins(1));
P2 = zeros(numbins,1);
P2(1) = sum(sorted_costs2 <= bins(1));
P3 = zeros(numbins,1);
P3(1) = sum(sorted_costs3 <= bins(1));
P4 = zeros(numbins,1);
P4(1) = sum(sorted_costs4 <= bins(1));
% 
% binsize0 = 10;
% binsize = binsize0;
% % save bins
% bin = zeros(2,1); m=1; bin(m) = binsize0;
% binstep = 10;
% bin_factor = 10;
% % minimum number of events within a bin
% mineventspbin = 20;
% % count the number of events between 0 and 100 MWh
% n=1;
% while (length(sorted_costs1) - n)>0
%     m = m+1;
%     bin(m) = binsize;
%     P1(m) = sum(sorted_costs1 <= binsize)-sum(P1(1:m-1));
%     P2(m) = sum(sorted_costs2 <= binsize)-sum(P2(1:m-1));
%     P3(m) = sum(sorted_costs3 <= binsize)-sum(P3(1:m-1));
%     if sum(sorted_costs1 <= binsize) <= sum(sorted_costs2 <= binsize) && sum(sorted_costs1 <= binsize) <= sum(sorted_costs3 <= binsize)
%         n = sum(sorted_costs1 <= binsize);
%     elseif sum(sorted_costs2 <= binsize) <= sum(sorted_costs1 <= binsize) && sum(sorted_costs2 <= binsize) <= sum(sorted_costs3 <= binsize)
%         n = sum(sorted_costs2 <= binsize);
%     else
%         n = sum(sorted_costs3 <= binsize);
%     end
%     if (P1(m) >= 20 && P2(m) >= 20 && P3(m) >= 20 && P4(m) >= 20) || length(sorted_costs1)-n <=20
%         binsize = binsize + binstep;
%     else
%         binstep = binstep*bin_factor;
%         binsize = binsize + binstep;
%         bin(m) = binsize;
%         P1(m) = sum(sorted_costs1 <= binsize)-sum(P1(1:m-1));
%         P2(m) = sum(sorted_costs2 <= binsize)-sum(P2(1:m-1));
%         P3(m) = sum(sorted_costs3 <= binsize)-sum(P3(1:m-1));
%         if sum(sorted_costs1 <= binsize) <= sum(sorted_costs2 <= binsize) && sum(sorted_costs1 <= binsize) <= sum(sorted_costs3 <= binsize)
%             n = sum(sorted_costs1 <= binsize);
%         elseif sum(sorted_costs2 <= binsize) <= sum(sorted_costs1 <= binsize) && sum(sorted_costs2 <= binsize) <= sum(sorted_costs3 <= binsize)
%             n = sum(sorted_costs2 <= binsize);
%         else
%             n = sum(sorted_costs3 <= binsize);
%         end
%     end
% end


binsize0 = 1;
binstep = 1;
bin_factor = 2;
binsize = binsize0;
% save bins
bin1 = zeros(2,1); m=0; %bin1(m) = binsize0;
% minimum number of events within a bin
mineventspbin = 20;
% count the number of events between 0 and 100 MWh
n=1;
while (length(sorted_costs1) - n)>0
    m = m+1;
    bin1(m) = binsize;
    P1(m) = sum(sorted_costs1 <= binsize)-sum(P1(1:m-1));
    n = sum(sorted_costs1 <= binsize);
    if (P1(m) >= mineventspbin) || sum(P1)==length(sorted_costs1)
        binsize = binsize + binstep;
    else
        while (P1(m) <= mineventspbin) && sum(P1)~=length(sorted_costs1)
        binstep = binstep*bin_factor;
        binsize = binsize + binstep;
        bin1(m) = binsize;
        P1(m) = sum(sorted_costs1 <= binsize)-sum(P1(1:m-1));
        n = sum(sorted_costs1 <= binsize);
        end
    end
end
% Next set
binsize = binsize0;
% save bins
bin2 = zeros(10,1); m=0; %bin2(m) = binsize0;
binstep = 1;
% count the number of events between 0 and 100 MWh
n=1;
while (length(sorted_costs2) - n)>0
    m = m+1;
    bin2(m) = binsize;
    P2(m) = sum(sorted_costs2 <= binsize)-sum(P2(1:m-1));
    n = sum(sorted_costs2 <= binsize);
    if (P2(m) >= mineventspbin) || sum(P2)==length(sorted_costs1)
        binsize = binsize + binstep; % calculate end of next bin
    else
        while (P2(m) <= mineventspbin) && sum(P2)~=length(sorted_costs1)
        binstep = binstep*bin_factor; % increase binstep size
        binsize = binsize + binstep; % calculate end of next bin
        bin2(m) = binsize; % save end of bin
        P2(m) = sum(sorted_costs2 <= binsize)-sum(P2(1:m-1));
        n = sum(sorted_costs2 <= binsize);
        end
    end
end

%Next set:
binsize = binsize0;
% save bins
bin3 = zeros(2,1); m=0; %bin3(m) = binsize0;
binstep = 1;
% minimum number of events within a bin
mineventspbin = 20;
% count the number of events between 0 and 100 MWh
n=1;
while (length(sorted_costs1) - n)>0
    m = m+1;
    bin3(m) = binsize;
    P3(m) = sum(sorted_costs3 <= binsize)-sum(P3(1:m-1));
    n = sum(sorted_costs3 <= binsize);
    if (P3(m) >= mineventspbin) || sum(P3)==length(sorted_costs1)
        binsize = binsize + binstep;
    else
        while (P3(m) <= mineventspbin) && sum(P3)~=length(sorted_costs1)
        binstep = binstep*bin_factor;
        binsize = binsize + binstep;
        bin3(m) = binsize;
        P3(m) = sum(sorted_costs3 <= binsize)-sum(P3(1:m-1));
        n = sum(sorted_costs3 <= binsize);
        end
    end
end

% Normalizing histograms

% Probability density function
[N1, BIN] = histc(P1, bin1);
pdf1 = N1./trapz(bin1,N1);
[N2, BIN] = histc(P2, bin2);
pdf2 = N2./trapz(bin2,N2);
[N3, BIN] = histc(P3, bin3);
pdf3 = N3./trapz(bin3,N3);

% % Probability distribution
% P1 = P1./length(sorted_costs1);
% P2 = P2./length(sorted_costs1);
% P3 = P3./length(sorted_costs1);

% % probability distribution
% figure
% % subplot(1,2,1)
% % loglog(bins(2:end),P1(2:end))
% % hold on
% % loglog(bins(2:end),P2(2:end))
% % loglog(bins(2:end),P3(2:end))
% % [bins1, P1_0] = stairs(bin1,P1);
% % [bins2, P2_0] = stairs(bin2,P2);
% % [bins3, P3_0] = stairs(bin3,P3);
% % loglog(bins1,P1_0)
% % hold on
% % loglog(bins2,P2_0)
% % loglog(bins3,P3_0)
% title("Resilience measures of 39 bus cases")
% % legend("original", "original2")
% legend("original", "+5% load in DG", "+20% load in DG")%, "39 bus +100% load in DG")
% % legend("39 bus", "39 bus 5% PV", "39 bus 20% PV", "39 bus 100% PV")%, "39 bus 100% PV")
% % legend("6 bus", "6 bus 5% PV", "6 bus 20% PV")
% % legend("6 bus", "6 bus 5% PV", "6 bus 20% PV", "6 bus 100% PV")
% % legend("6 bus", "6 bus 5% PV", "6 bus 5% RL");
% % legend("6 bus - allow 0 outages", "6 bus - fair sample", "6 bus - allow 0 outages and add 1");
% % legend("39 bus - allow 0 outages", "39 bus - fair sample", "39 bus - allow 0 outages and add 1"); %"39 bus - bias sample 1 line out", 
% xlabel("ENS (MWh)"); ylabel("Prob(ENS)")
% set(gca, 'fontsize',13)
% hold off

% stairs plot of PDF
[bins1, PDF1] = stairs(bin1,pdf1);
[bins2, PDF2] = stairs(bin2,pdf2);
[bins3, PDF3] = stairs(bin3,pdf3);

figure
semilogx(bins1,PDF1)
hold on
semilogx(bins2,PDF2)
semilogx(bins3,PDF3)
title("Resilience measures of 39 bus cases")
xlabel("ENS (MWh)"); ylabel("PDF(ENS)")
set(gca, 'fontsize',13)
hold off

figure
p = histcounts(P1,[0; bin1],'Normalization','pdf')
semilogx(bin1,p);