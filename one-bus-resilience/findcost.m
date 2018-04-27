function [ Totalcapital ] = findcost( Recovery,cost )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[row,col] = size(Recovery);
Totalcapital = zeros(row,1);
for ii = 1: row
Totalcapital(ii) = sum(Recovery(ii,:))*cost;
end


end

