function [ failures ] = fail( Hurricane, numgen, numload, stats, Debug )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
fail = (Hurricane.*rand(1,numgen+numload));
failures = zeros(size(fail));
for jj = 1:length(Hurricane)
    if fail(jj)>= exp(stats(jj))
        failures(jj) = 1;
    else
        failures(jj) = 0;
    end
end
%figures go here
if Debug
end
end


