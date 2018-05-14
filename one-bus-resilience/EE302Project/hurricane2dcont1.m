function [ H ] = hurricane2dcont1( n,r,location,debug )
%[ H ] = hurricane2dcont( n,M,r,debug ), continuous version of
%hurricane2d.m
%   Inputs:
%   n - number of rows in output matrix
%   r - radius parameter for size of storm --- mqybe this should be
%   stochastic?
%   location - 
%   debug - if debug is 0, no plots will occur, if it is included in the
%   function call, and not set to 0 or 'false' the function will produce
%   plots.
%   Outputs:
%   H - strength array of hurricane wind speeds at pertinent locations
%   specified

if nargin == 3
    debug = 0;
end
% power law distributed maximum hurricane intensities
M = rand(1)^(1/-0.8)/5;

%preallocating H
H = zeros(1,length(location));

trajectory = n*rand(1,4); % random 2 x and y points that are on the line of the huricane
x = [trajectory(1); trajectory(2)]; % sets first two random integers as x coordinates
y = [trajectory(3); trajectory(4)]; % sets last two random integers as y coordinates
R = 10*r; % can use to scale r to the proper spacial dimentions
% Line = polyfit(x,y,1); % line defining huricane trajectory
% X = 1:.1:n;
% Y = Line(1).*X +Line(2);
a = [x(1) y(1) 0] - [x(end) y(end) 0];
for ii=1:length(location)
    b = [location(ii,1) location(ii,2) 0] - [x(end) y(end) 0];
    d = norm(cross(a,b)) / norm(a);
    h = (exp(-d/r));
    %h = M*(exp(-d/r));
    H(ii) = h; %distance from line.
end
if debug
    Hurricane = zeros(n);
    Line = polyfit(x,y,1); % line defining huricane trajectory
    X = 1:.1:n;
    Y = Line(1).*X +Line(2);
    a = [X(1) Y(1) 0] - [X(end) Y(end) 0];
    for ii=1:n
        for jj=1:n
            b = [ii jj 0] - [X(end) Y(end) 0];
            d = norm(cross(a,b)) / norm(a);
            hh = (exp(-d/r));
            %hh = M*(exp(-d/r));
            Hurricane(jj,ii) = hh; %distance from line.
        end
    end
    figure(gcf)
    hold on
    contour(Hurricane)
    plot(location(:,1),location(:,2),'om')
    title('Huricane Strength over Network Map')
    xlabel('x location')
    ylabel('y loaction')
    hold off
end
end

