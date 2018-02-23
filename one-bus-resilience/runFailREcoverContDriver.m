%run FailRecoverContDriver.m


%   0, 1, or 2 for figures from main driver debug = 1; for both figures from 
%   main and figures from functions run from driver debug =2;
debug = 2;
geog = 100; % size of space
Mag = 4; %magnitude of hurricane on rictor scale
r = 10; %radius of hurricane -- should depends on space size and units
j = 10; % number of iterations
numgen = 15; % number of generators
numload = 15; % number of loads
robustness = 0.2*ones(numgen+numload,1); %information about reliability of generators and loads
recoverystats = 0.2*ones(numgen+numload,1); %information about ease of recovery in generators and loads

[ Recovery, TotFails ] = FailRecoverContDriver( geog, Mag, r, j, numgen, numload, robustness, recoverystats, debug)

figure
imagesc(Recovery)
