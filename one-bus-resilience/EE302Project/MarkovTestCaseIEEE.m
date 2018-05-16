function [recovery, cost, Precover,Bins] = MarkovTestCaseIEEE(coefNeighbors, averages,r,costloadshed,Recoveryrate)
% attempt to recreate using markov modeling why recovery times follow a
% power law in real life systems
% the typical model for recovery times is an exponentialdecay curve


% indroduce parameters
geog = 100; % size of square geography
N = 73; % number of nodes allowable states for each node are 0, 1
K = []; % number of edges fully connected directed would be N^2, unconnected N(N-1)/2
maxtimesteps = 1000; % number of timesteps to iterate probabilities
%coefNeighbors = 0.5:0.1:0.9;
%averages = 100;
%r = 10; % proportional to radius of hurricane
debug = 0;
%costloadshed = 10;
%Recoveryrate = 1;

L = load('IEEErts96BusLocations.mat');
location = L.IEEErts96BusLocations;
edges = IEEEedges();
GenDemand = [];

% introduce stress to network, the probability of each node failing(0->1)
% is based on stress at the node and the connections of the node that have
% failed
recovery = zeros(length(coefNeighbors),averages);
Cost = zeros(length(coefNeighbors),averages);
Precover{length(coefNeighbors)} = [];
Bins{length(coefNeighbors)} = [];
for W = 1:length(coefNeighbors)
    neighborWeight = coefNeighbors(W);
    recoverytime = zeros(averages,1);
    cost = zeros(averages,1);
    for m = 1:averages
        Hurricane = hurricane2dcont1( geog,r,location,debug ); % gives a percent of the intesity of the wind at the locations of each of the nodes
        startingstate = zeros(N,1);
        if isempty(Recoveryrate)
        [statematrix,recoverytime(m)] = MarkovModel0(N,edges,Hurricane,neighborWeight,maxtimesteps,startingstate);
        else
        [statematrix,recoverytime(m)] = MarkovModel1(N,edges,Hurricane,neighborWeight,maxtimesteps,startingstate,Recoveryrate);
        end
            cost(m) = LostLoad(N,edges,GenDemand,statematrix,costloadshed);
    end
    recovery(W,:) = sort(recoverytime);
    Cost(W,:) = sort(cost);
    bins = min(recovery(W,:)):1:max(recovery(W,:));

    Pr = zeros(length(bins));
    binsbackwards = max(recovery(W,:)):-1:min(recovery(W,:));
    Pr(1) = sum(bins(end) == recovery(W,:));%/length(recovery(W,:));
    for ii = 2:length(bins)
        Pr(ii) = Pr(ii-1) + sum(binsbackwards(ii) == recovery(W,:));%/length(recovery(W,:));
        %Pr(ii) = sum(bins(ii) == recovery(W,:));%/length(recovery(W,:));
    end
    Pr = Pr(Pr~=0);
    Bins{W} = binsbackwards(Pr~=0);
    Precover{W} = Pr;

    binscost = min(Cost(W,:)):1:max(Cost(W,:));
    figure(gcf)
    subplot(2,length(coefNeighbors),W)
    hold on
    hist(recovery(W,:),bins)
    title(['Weight of Neighbors = ' num2str(neighborWeight)])
    xlabel('Recovery Time '); ylabel('Number of Cases')
    hold off
    subplot(2,length(coefNeighbors),length(coefNeighbors)+W)
    hold on
    hist(Cost(W,:),binscost)
    title(['Weight of Neighbors = ' num2str(neighborWeight)])
    xlabel('Cost') ; ylabel('Number of Cases')
    hold off
end
end
