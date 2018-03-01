function [ Recoveries ] = recover( n, recoveries, recoverystats, Debug )
%[ Recoveries ] = recover( n, recoveries, recoverystats, failures, Debug )
%Based on the robustness and the hurricane intensity at each location of 
% parts stochastically predicts failures in each part
%   Inputs:
%   n - timestep
%   recoveries - 
%   recoverystats -
%   Debug - 
%   Outputs: 
%   Recoveries - vector of 1s for components that never failed, 0s if it
%   hasn't recovered yet, and n if it recovered in this time step

numfails = length(recoveries(recoveries==0));
recoveries(recoveries==0) = n.*(rand(numfails,1) < (1-exp(-n.*recoverystats(recoveries==0))));
Recoveries = recoveries;

return

