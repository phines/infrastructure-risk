% loop through markovrecover
close all
clear

% indroduce parameters
geog = 100;% size of square geography
Ns = [10 20 50 100 200]; % number of nodes allowable states for each node are 0, 1 
k = [1 2.5 4]; % number of edges fully connected directed would be N^2, unconnected N(N-1)/2
NF = 3; % number of failures ot implement
maxtimesteps = 1000; % number of timesteps to iterate probabilities
coefNeighbors = 0.9;
averages = 1000;
r = 10; % proportional to radius of hurricane
debug = 0;
costloadshed = 10;
recoveryRate = 0.5;

avgrecov = zeros(length(Ns),length(coefNeighbors));
avgcost = avgrecov;
recovstdev = avgrecov;
coststdev = avgrecov;
for n = 1:length(Ns)
    N = Ns(n);
    K = N;
    figure
    [avgrecov(n,:), avgcost(n,:),recovstdev(n,:),coststdev(n,:)] = markovRecover(geog,N,K,...
        maxtimesteps,coefNeighbors,averages,r,costloadshed,recoveryRate,debug);
end
figure
for l = 1:length(coefNeighbors)
subplot(1,length(coefNeighbors),l)
hold on
errorbar(Ns,avgrecov(:,l),recovstdev(:,l),'mo')
%errorbar(Ns,avgcost(:,l),coststdev(:,l),'co')
title('Average Recovery Times over Several Network Sizes')
xlabel('Network size, N'); ylabel('Average Recovery Time')
hold off
end

%% Ks
% loop through markovrecover
close all
clear

% indroduce parameters
geog = 100;% size of square geography
N = 50; % 20 50 100 200]; % number of nodes allowable states for each node are 0, 1 
k = [0 10 50 100 500]; % number of edges fully connected directed would be N^2, unconnected N(N-1)/2
maxtimesteps = 100; % number of timesteps to iterate probabilities
coefNeighbors = 0.9;
averages = 100;
r = 10; % proportional to radius of hurricane
debug = 0;
costloadshed = 10;
recoveryRate = 0.5;
avgrecov = zeros(length(k),length(coefNeighbors));
avgcost = avgrecov;
recovstdev = avgrecov;
coststdev = avgrecov;
for n = 1:length(k)
    K = k(n);
    figure
    [avgrecov(n,:), avgcost(n,:),recovstdev(n,:),coststdev(n,:)] = markovRecover(geog,N,K,...
        maxtimesteps,coefNeighbors,averages,r,costloadshed,recoveryRate,debug);
end
figure
for l = 1:length(coefNeighbors)
subplot(1,length(coefNeighbors),l)
hold on
lll = k/N+1;
errorbar(lll,avgrecov(:,l),recovstdev(:,l),'mo')
%errorbar(Ns,avgcost(:,l),coststdev(:,l),'co')
title('Average Recovery Times over Several Network Degrees')
xlabel('Network degree, K'); ylabel('Average Recovery Time')
hold off
end

%% new one
%housekeeping
clear
close all
% indroduce parameters
geog = 100;% size of square geography
N = 50; % number of nodes allowable states for each node are 0, 1 
k = [100 500]; % number of edges fully connected directed would be N^2, unconnected N(N-1)/2
NF = 3; % number of failures ot implement
maxtimesteps = 100; % number of timesteps to iterate probabilities
coefNeighbors = 0.9;
averages = 100;
r = 10; % proportional to radius of hurricane
debug = 0;
costloadshed = 10;
recoveryRate = 0.5;


Rrate = k; %differnt values of K
Recovery{length(Rrate),1} = [];
Cost{length(Rrate),1} = [];
Prob{length(Rrate),1} = [];
for r = 1:length(Rrate)
    K = Rrate(r);
    figure; title(['Recovery rate = ' num2str(K)])
    [avgrecov(r,:), avgcost(r,:),recovstdev(r,:),coststdev(r,:),Recovery{r},...
        Cost{r},Prob{r},Bins{r}] = markovRecover ...
    (geog,N,K,maxtimesteps,coefNeighbors,averages,r,costloadshed,recoveryRate,debug);
    Rtemp = Prob{r};
    Btemp = Bins{r};
    figure
    for n = 1:length(coefNeighbors)
        Precover = (Rtemp{n});
        bin = (Btemp{n});
        val = polyfit(log(bin(2:end-3))',log(Precover(2:end-3)),1);
        subplot(1,length(coefNeighbors),n)
        loglog(bin,Precover)
        hold on
        Kk = K/N+1;
        legend(['K = ' num2str(Kk) ', Slope = ' num2str(val(1))])
        title(['Wieght of Neighbors = ' num2str(coefNeighbors(n))]); xlabel('Recovery time'); ylabel('Number of Cases Recovered in <= x(Recovery time)')
        hold off
    end
end