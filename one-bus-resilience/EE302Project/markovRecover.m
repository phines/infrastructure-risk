function [avgrecov, avgcost,recovstdev,coststdev,recovery,cost,Precover,Bins] = markovRecover ...
    (geog,N,K,maxtimesteps,coefNeighbors,averages,r,costloadshed,recoveryRate,debug)
% indroduce parameters
% geog = 100; % size of square geography
% N = 1000; % number of nodes allowable states for each node are 0, 1
% K = 2500; % number of edges fully connected directed would be N^2, unconnected N(N-1)/2
% NF = 3; % number of failures ot implement
% maxtimesteps = 100; % number of timesteps to iterate probabilities
% coefNeighbors = 0:0.2:1;
% averages = 100;
% r = 10; % proportional to radius of hurricane
% debug = 0;
% costloadshed = 10;

% introduce small network of markov nodes
connections = randi(N,K,2);
edges{N,1} = [];
% Make sure each node has at least one edge
for jj = 1:N
    Rand = randi(N,1);
    if Rand == jj
        if jj == N
            Rand = 1;
        else
            Rand = jj+1;
        end
        edges{jj,1} = Rand;
        edges{Rand,1} = jj;
    else
        edges{jj,1} = Rand;
        edges{Rand,1} = jj;
    end
end

for n = 1:K
    line = connections(n,:);
    if line(1) == line(2) % discard selfloops
        if line(1) == N
            line(2) = line(2) - 1;
        else
            line(2) = line(1) + 1;
        end
        edges{line(1)} = [edges{line(1)} line(2)];
        edges{line(2)} = [edges{line(2)} line(1)];
    else
        edges{line(1)} = [edges{line(1)} line(2)];
        edges{line(2)} = [edges{line(2)} line(1)];
    end
end


%creating locations for each node in the geography
row = geog*rand((N),1);
col = geog*rand((N),1);
location = [row col];  % saving the locations of the nodes

%creating current demand and generation before the hurricane (assuming 
% these characteristics are static throughout the model)
GenDemand = [];


% introduce stress to network, the probability of each node failing(0->1)
% is based on stress at the node and the connections of the node that have
% failed
recovery = zeros(length(coefNeighbors),averages);
Cost = zeros(length(coefNeighbors),averages);
avgrecov = zeros(length(coefNeighbors),1);
avgcost = avgrecov;
recovstdev = avgrecov;
coststdev = avgrecov;
Precover{length(coefNeighbors)} = [];
Bins{length(coefNeighbors)} = [];
for W = 1:length(coefNeighbors)
    neighborWeight = coefNeighbors(W);
    recoverytime = zeros(averages,1);
    cost = zeros(averages,1);
    for m = 1:averages
        Hurricane = hurricane2dcont1( geog,r,location,debug ); % gives a percent of the intesity of the wind at the locations of each of the nodes
        startingstate = zeros(N,1);
        [statematrix,recoverytime(m)] = MarkovModel0(N,edges,Hurricane,neighborWeight,maxtimesteps,startingstate);
        cost(m) = LostLoad(N,edges,GenDemand,statematrix,costloadshed);
    end
    recovery(W,:) = sort(recoverytime);
    Cost(W,:) = sort(cost);
    avgrecov(W) = mean(recoverytime);
    avgcost(W) = mean(cost);
    recovstdev(W) = std(recoverytime);
    coststdev(W) = std(cost);
    bins = min(recovery(W,:)):1:max(recovery(W,:));
    binscost = min(Cost(W,:)):1:max(Cost(W,:));
    
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
    
    
    figure(gcf)
    subplot(2,length(coefNeighbors),W)
    hold on
    hist(recovery(W,:),bins)
    title(['Weight of Neighbors = ' num2str(neighborWeight)])
    xlabel('Recovery Time '); ylabel('Number of Nodes')
    hold off
    subplot(2,length(coefNeighbors),length(coefNeighbors)+W)
    hold on
    hist(Cost(W,:),binscost)
    title(['Weight of Neighbors = ' num2str(neighborWeight)])
    xlabel('Cost') ; ylabel('Number of Nodes')
    hold off
end


end

