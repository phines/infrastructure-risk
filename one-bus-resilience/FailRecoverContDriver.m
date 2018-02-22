function [ Recovery, TotFails ] = FailRecoverContDriver( geog, Mag, r, j, numgen, numload, stats, debug)
%[ Recovery, TotFails ] = FailRecoverDriver( n, M, r, j, numgen, numload,debug)
% Generates random hurricanes, checks part failures based on those
% hurricanes using arbitrary probabilities of 70% for generators and 50%
% for lines
%   Inputs:
%   geog is the size of the square matrix that crudely represents our geography
%   Mag is the magnitude of the generated hurricanes
%   numgen and num load are the number of generators represented in the
%   geography and the number of loads represented in the geography,
%   respectively
%   j is the number of itereations
%   debug creates plots when it is not false, 0, or not input at all
%   Outputs:
%   Recovery - recovery time for failed parts
if nargin == 7
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
 

Recovery = zeros(j,numgen+numload);
TotFails = zeros(j,numgen+numload);
for ii = 1:j
    row = geog*rand((numgen+numload),1);
    col = geog*rand((numgen+numload),1);
    location = [row col];  % saving the locations of generation and loads
    Hurricane = hurricane2dcont(geog,Mag,r,location,Debug);
    [ failures ] = fail( Hurricane, numgen, numload, stats, Debug );
    recoveries = zeros(numgen+numload,1); n=0;
    recoveries(failures == 0) = 1; % this part is operational and wasn't critically damaged in hurricane
    while ~isempty(recoveries(recoveries == 0)) || n <= 10
        [ recoveries ] = recover( n, recoveries, stats, Debug );
        n = n+1;
    end
    TotFails(ii,:) = failures;
    Recovery(ii,:) = recoveries;
end

%put figure here
if debugger
    
end

end

