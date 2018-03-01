function [ failures ] = fail( hurricaneIntensity, robustness, Debug )
%[ failures ] = fail( hurricaneIntensity, robustness, Debug )
%Based on the robustness and the hurricane intensity at each location of 
% parts stochastically predicts failures in each part
%   Inputs:
%   hurricaneIntensity - 
%   robustness - 
%   Debug - 
%   Outputs: 
%   failures - vector of 1s and 0s for each component, 1 is failure, 0 is not

n = length(hurricaneIntensity);
vulnerability = (robustness).^(-1); % think about this
hitMagnitude = hurricaneIntensity.*vulnerability'; % think about this
Pr_fail = 1-exp(-hitMagnitude);
failures = 1*(rand(1,n) <= Pr_fail);

return
