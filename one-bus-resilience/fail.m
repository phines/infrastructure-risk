function [ failures ] = fail( hurricaneIntensity, robustness, Debug )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

n = length(hurricaneIntensity);
vulnerability = 1/robustness; % think about this
hitMagnitude = hurricaneIntensity.*vulnerability; % think about this
Pr_fail = 1-exp(-hitMagnitude);
failures = rand(n,1) < Pr_fail;

return
%{
fail = (hurricaneIntensity.*rand(1,numgen+numload));

failures = zeros(size(fail));
for jj = 1:length(hurricaneIntensity)
    if fail(jj)>= exp(-robustness(jj))
        failures(jj) = 1;
    else
        failures(jj) = 0;
    end
end
%figures go here
if Debug
end
end
%}
