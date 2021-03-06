function [ Recovery, TotFails, Hurricane, location ] = FailRecoverContDriver( geog, hurricaneMagnitude, ...
    hurricaneSize, numHurricane, numgen, numload, robustness, recoveryStats, NormalRecovery, debug)
%[ Recovery, TotFails ] = FailRecoverDriver( n, M, r, j, numgen, numload,debug)
% Generates random hurricanes, checks part failures based on those
% hurricanes using arbitrary probabilities of 70% for generators and 50%
% for lines
%   Inputs:
%   geog is the size of the square matrix that crudely represents our geography
%   hurricaneMagnitude is the magnitude of the generated hurricanes
%   numgen and num load are the number of generators represented in the
%   geography and the number of loads represented in the geography,
%   respectively
%   numHurricane is the number of hurricane events to simulate
%   debug creates plots when it is not false, 0, or not input at all
%   Outputs:
%   Recovery - recovery time for each part (column) for each hurricane
%   (row); if the part did not fail the value will be 1 otherwise the value
%   will be the timestep at which recovery occured.
%   TotFails - matrix of failures of each part (column) for each hurricane
%   (row)


if nargin == 10
    if debug >= 2
        Debug = 1;
        debugger =1;
    elseif debug >=1
        Debug = 0;
        debugger =1;
    else
        Debug = 0;
        debugger = 0;
    end
else
    Debug = 0;
    debugger = 0;
end
if debugger
    figure(1)
end

Recovery = zeros(numHurricane,numgen+numload);
TotFails = zeros(numHurricane,numgen+numload);
for ii = 1:numHurricane
    row = geog*rand((numgen+numload),1);
    col = geog*rand((numgen+numload),1);
    location = [row col];  % saving the locations of generation and loads
    Hurricane = hurricane2dcont(geog,hurricaneMagnitude,hurricaneSize,location,Debug);
    [ failures ] = fail( Hurricane, robustness, Debug );
    recoveries = zeros(numgen+numload,1); n=1;
    recoveries(failures == 0) = 1; % these parts of the system are operational and wasn't critically damaged in hurricane
    if NormalRecovery
    while ~isempty(recoveries(recoveries == 0)) || n <= 100
        n = n+1;
        [ recoveries ] = recover( n, recoveries, recoveryStats, Debug );
    end
    else
        recoveries = recoverPL( recoveries, Debug )
    end
    TotFails(ii,:) = failures';
    Recovery(ii,:) = (recoveries-1)';
    %put figure here
if debugger
    figure(1)
    hold on
    plot(1:numload+numgen, (recoveries-1))
    title('Number of recoveries of failures for each hurricane')
    xlabel('Component number')
    ylabel('Timestep of Recovery')
    hold off
end
end



end

