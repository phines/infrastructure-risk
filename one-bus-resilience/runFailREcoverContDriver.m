%run FailRecoverContDriver.m


%   0, 1, or 2 for figures from main driver debug = 1; for both figures from 
%   main and figures from functions run from driver debug =2;
debug = 2;
geog = 100; % size of space
Mag = 4; %magnitude of hurricane on rictor scale
r = 1; %radius of hurricane -- should depends on space size and units
j = 5; % number of iterations
numgen = 10; % number of generators
numload = 10; % number of loads
stats = 0.5*ones(numgen+numload,1); %information about reliability of generators and loads

[ Recovery, TotFails ] = FailRecoverContDriver( geog, Mag, r, j, numgen, numload, stats, debug);

Recovery
imagesc(Recovery)
