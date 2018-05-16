%run MarkovTestCase

% house keeping
clear
close all
clear

% indroduce parameters
geog = 100; % size of square geography
N = 73; % number of nodes allowable states for each node are 0, 1
K = []; % number of edges fully connected directed would be N^2, unconnected N(N-1)/2
maxtimesteps = 1000; % number of timesteps to iterate probabilities
coefNeighbors = 0.9:0.01:0.97;
averages = 1000;
radius = 10; % proportional to radius of hurricane
debug = 0;
costloadshed = 10;
Rrate = 0.5;%0.1:.1:1;
Recovery{length(Rrate),1} = [];
Cost{length(Rrate),1} = [];
Prob{length(Rrate),1} = [];
for r = 1:length(Rrate)
    Recoveryrate = Rrate(r);
    figure; title(['Recovery rate = ' num2str(Recoveryrate)])
    [Recovery{r}, Cost{r},Prob{r},Bins{r}] = MarkovTestCaseIEEE(coefNeighbors, averages,radius,costloadshed,Recoveryrate);
    Rtemp = Prob{r};
    Btemp = Bins{r};
    figure
    for n = 1:length(coefNeighbors)
        Precover = (Rtemp{n});
        bin = (Btemp{n});
        val = polyfit(log(bin(2:end-13))',log(Precover(2:end-13)),1);
        subplot(1,length(coefNeighbors),n)
        loglog(bin,Precover)
        hold on
        legend(['Slope = ' num2str(val(1))])
        title(['Wieght of Neighbors = ' num2str(coefNeighbors(n))]); xlabel('Recovery time'); 
        ylabel('Number of Times Recovered in <= x(Recovery time)')
        hold off
    end
end
for r = 1:length(Rrate)
    Recoveryrate = [];
    figure; title(['Recovery rate = ' num2str(Recoveryrate)])
    [Recovery{r}, Cost{r},Prob{r},Bins{r}] = MarkovTestCaseIEEE(coefNeighbors, averages,radius,costloadshed,Recoveryrate);
    Rtemp = Prob{r};
    Btemp = Bins{r};
    figure
    for n = 1:length(coefNeighbors)
        Precover = (Rtemp{n});
        bin = (Btemp{n});
        val = polyfit(log(bin(2:end-13))',log(Precover(2:end-13)),1);
        subplot(1,length(coefNeighbors),n)
        loglog(bin,Precover)
        hold on
        legend(['Slope = ' num2str(val(1))])
        title(['Wieght of Neighbors = ' num2str(coefNeighbors(n))]); xlabel('Recovery time'); 
        ylabel('Number of Times Recovered in <= x(Recovery time)')
        hold off
    end
end


