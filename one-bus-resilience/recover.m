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
%   Recoveries - 

recoveries(recoveries==0) = n.*(rand(1) <= (1-exp(-recoverystats(recoveries==0).*n)));
Recoveries = recoveries;

return

