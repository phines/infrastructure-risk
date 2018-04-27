function [ Recoveries ] = recoverPL( recoveries, recoverystats, Debug )
%[ Recoveries ] = recover( n, recoveries, recoverystats, failures, Debug )
%Based on the robustness and the hurricane intensity at each location of 
% parts stochastically predicts failures in each part
%   Inputs:
%   n - timestep
%   recoveries - 
%   recoverystats -
%   Debug - 
%   Outputs: 
%   Recoveries - vector with a 1 for components that never failed and a 0 if it
%   hasn't recovered yet, and n if it recovered in this time step

numfails = length(recoveries(recoveries==0));
recoveries(recoveries==0) = rand(numfails,1).^(1/-.87);
Recoveries = recoveries;

return

