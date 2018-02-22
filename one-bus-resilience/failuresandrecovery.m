function [ Recovery ] = failuresandrecovery( n,M,numgen,numload,j, debug )
%[ Recovery ] = failuresandrecovery( n,M,numgen,numload,j, debug )
% Generates random hurricanes, checks part failures based on those
% hurricanes using arbitrary probabilities of 70% for generators and 50%
% for lines
%   Inputs:
%   n is the size of the square matrix that crudely represents our geography
%   M is the magnitude of the generated hurricanes
%   numgen and num load are the number of generators represented in the
%   geography and the number of loads represented in the geography,
%   respectively
%   j is the number of itereations
%   debug creates plots when it is not false, 0, or not input at all
%   Outputs:
%   Recovery - recovery time for failed parts
if nargin <= 5
    debug =0;
elseif nargin <=4
    j = 1;
    fprintf('no number of iterations in inputs, only running for 1')
elseif nargin <=3
    fprintf('error: not enough input arguments')
end

r = 1;
Recovery = zeros((numgen+numload),j);
mult = 10;% sets multiplier to be M to be the same as the set huricane top wind speed
fails = zeros(1,j);
genfail = zeros(1,j);
lines = zeros(1,j);
for k = 1:j
H= hurricane2d(n,M,r);
% check for failures at each generator or load
fail = zeros(numgen,1);
line = zeros(numgen+numload,1);
location = zeros(2,numgen+numload);
for ii = 1:(numgen+numload)
    row = randi([1,n],1);% picking random x locations for generators and loads
    col = randi([1,n],1);% picking random y locations for generators and loads
    location(:,ii) = [row col];% saving the locations of generators and loads
    Sev = H(row,col);
    strength = rand(1);
    if ii <=numgen
    if Sev/mult*(1-strength) >= 0.7 % check for failure in generator (true is defined as failure)
        fail(ii) = 1; % if true failure occurs and saved in the fails vector with the index of the generator
    end
    end
    if Sev/mult >= 0.5 % check failure in lines (true is defined as failure)
        line(ii) = 1; %sets the line element with the index of the generator or load to 1 if the line fails
    end
end
if debug
figure(k)
contour(H)
hold on
scatter(location(1,1:numgen),location(2,1:numgen))
scatter(location(1,(1+numgen):end),location(2,(1+numgen):end),'r')
hold off
end
recoverytime = zeros(numgen+numload,1);
for ii = 1:(numgen+numload)
if line(ii) == 1
    recoverytime(ii) = randi([1,10], 1);
end
if ii <= numgen
    if fail(ii) == 1
        recoverytime(ii) = recoverytime(ii) + randi([1,30], 1);
    end
end
end
lines(k) = length(line(line==1));
fails(k) = (length(fail(fail==1)) + length(line(line==1)));
genfail(k) = length(fail(fail==1));
Recovery(:,k) = recoverytime;
end
ProbGenFail = sum(genfail./numgen)/j;
ProbLineFail = sum(lines./(numgen+numload))/j;
if debug
figure
plot(1:k, fails, 1:k, genfail)
hold on
title('Number of failures due to hurricane2d')
xlabel('Run number')
ylabel('Number of Failures')
legend('Number of Total Failures', 'Number of Generator Failures', 'Location','northwest')
hold off
end
end
