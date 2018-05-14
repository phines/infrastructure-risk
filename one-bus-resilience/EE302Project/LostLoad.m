function [ Cost ] = LostLoad(N,edges,GenDemand,statematrix,cost)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[row,col] = size(statematrix);
Totalcapital = zeros(row,1);
for ii = 1: row
Totalcapital(ii) = sum(statematrix(ii,:))*cost;
end
Cost = sum(Totalcapital);

% want to have cost each time something fails, and also compute a power
% flow to determine lost load and have a load shedding cost to quantify
% that portion as well.

end

