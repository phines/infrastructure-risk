% RunningTestcases
%run FailRecoverContDriver.m


%   0, 1, or 2 for figures from main driver debug = 1; for both figures from 
%   main and figures from functions run from driver debug =2;
debug = 2;
geog = 100; % size of space
Mag = 5; %magnitude of hurricane on rictor scale
r = 10; %radius of hurricane -- should depends on space size and units
j = 10; % number of iterations
numgen = 33; % number of generators
numload = 51; % number of loads
robustness = 0.4*ones(numgen+numload,1); %information about reliability of generators and loads
recoverystats = 0.5*ones(numgen+numload,1); %information about ease of recovery in generators and loads  ---- percent resources
NormalRecovery = 0;

[ Recovery, TotFails, Hurricane, location ] = TestCaseIEEE96_76( geog, Mag, r, j, robustness, recoverystats, NormalRecovery, debug);

index = 1:length(Hurricane);
values = [index' location Hurricane' TotFails' Recovery'];
%fprintf( '    index   x and y location     Hurricane    TotFails   Recovery')
%values


figure
subplot(1,3,1)
hold on
imagesc(Recovery)
colormap(jet);
set(gca,'clim',[0 max(max(Recovery))]); % make sure the color limits down change dynamically
ch=colorbar;
set(ch,'Ytick',[0 1 max(max(Recovery))],'Yticklabel',{'No Failure','1 Timestep Recovery',[num2str(max(max(Recovery))) ' Timestep Recovery']})
title('Map of Recovery Times over each hurricane(down) and component (across)')
hold off

subplot(1,3,2)
hold on
plot(1:j, sum(TotFails'))
title('Number of failures due to hurricane2dcont')
xlabel('Hurricane number')
ylabel('Number of Failures')
hold off

cost = 1000;
Cost = findcost(Recovery, cost);
subplot(1,3,3)
hold on
plot(1:j,Cost)
title('Cost due to each hurricane')
xlabel('Hurricane number')
ylabel('Total Cost')
hold off
