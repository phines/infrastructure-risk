function [ failures ] = simplefailures( n,M,numgen,numload,j )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
r = 1;
mult = 30;
failures = zeros(1,j);
for k = 1:j
Map = zeros(n);
H= huricane2d(n,M,r);
row = randi([1,n],[(numgen+numload),1]);
col = randi([1,n],[(numgen+numload),1]);
location = [row col];  % saving the locations of generation and loads
% check for failures at each location
fail = zeros((numgen+numload),1);
for ii = 1:(numgen+numload)
    Sev= H(location(ii,1),location(ii,2));
    if Sev <= 0.5*mult
        if rand(1)<= 0.2 
            fail(ii) = 1;
        end
    elseif Sev <=0.9*mult
        if rand(1) <=.5
            fail(ii) = 1;
        end
    else
        if rand(1) <=0.9
            fail(ii) = 1;
        end
    end
end
failures(k) = length(fail(~(fail==0)));
figure(1)
hold on
scatter(location(1:numgen,1),location(1:numgen,2))
scatter(location((1+numgen):end,1),location((1+numgen):end,2),'r')
hold off
end
figure(2)
plot(1:k, failures)
hold on
title('Number of failures due to hurricane2d')
xlabel('Run number')
ylabel('Number of Failures')
hold off
end

