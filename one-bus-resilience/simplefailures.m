function [ failures ] = simplefailures( n,M,numgen,numload,j, debug )
%[ failures ] = simplefailures( n,M,numgen,numload,j )
%   Inputs:
%   n is the size of the square matrix that crudely represents our geography
%   M is the magnitude of the generated hurricanes
%   numgen and num load are the number of generators represented in the
%   geography and the number of loads represented in the geography,
%   respectively
%   j is the number of itereations
%   debug creates plots when it is not false, 0, or not input at all
%   Outputs:
%   Outputs:
%   failures - gives number of failures due to generated  hurricane
if nargin <= 5
    debug =0;
elseif nargin <=4
    j = 1;
    fprintf('no number of iterations in inputs, only running for 1')
elseif nargin <=3
    fprintf('error: not enough input arguments')
end

r = 1;
mult = 30;
failures = zeros(1,j);
for k = 1:j
H= hurricane2d(n,M,r);
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
if debug
figure(1)
hold on
scatter(location(1:numgen,1),location(1:numgen,2))
scatter(location((1+numgen):end,1),location((1+numgen):end,2),'r')
hold off
end
end
if debug
figure(2)
plot(1:k, failures)
hold on
title('Number of failures due to hurricane2d')
xlabel('Run number')
ylabel('Number of Failures')
hold off
end
end

