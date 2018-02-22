function [ Recoveries ] = recover( n, recoveries, recoverystats, Debug )
%[ Recoveries ] = recover( n, recoveries, recoverystats, failures, Debug )
%   Detailed explanation goes here

for jj = 1:length(recoveries)
    if recoveries(jj) == 0
        if rand(1) <= exp(recoverystats(jj))
            recoveries(jj) = n;
        else
            recoveries(jj) = 0;
        end
    end
end
Recoveries = recoveries;

%figures go here
if Debug
end
end

