function [ H ] = hurricane2d( n,M,r,debug )
%[ H ] = hurricane2d( n,M,r,debug )
%   Inputs:
%   n - number of rows in output matrix
%   M - number of columns in output matrix
%   r - size  parameter for size of storm
%   debug - if debug is 0, no plots will occur, if it is included in the
%   function call, and not set to 0 or 'false' the function will produce
%   plots.
%   Outputs:
%   H - strength array of hurricane wind speeds over the matrix size
%   specified

if nargin == 3
    debug = 0;
end

D = zeros(n);
m = randi(n,[1,4]); % random 2 x and y points that are on the line of the huricane
x = [m(1); m(2)]; % sets first two random integers as x coordinates
y = [m(3); m(4)]; % sets last two random integers as y coordinates
Line = polyfit(x,y,1); % line defining huricane trajectory
X = 1:.1:n;
Y = Line(1).*X +Line(2);
a = [X(1) Y(1) 0] - [X(end) Y(end) 0];
for ii=1:n
    for jj=1:n
        b = [ii jj 0] - [X(end) Y(end) 0];
        d = norm(cross(a,b)) / norm(a);
        D(jj,ii) = d; %distance from line.
    end
end
R = 10*r; % can use to scale r to the proper spacial dimentions
H = M.*(exp(-D./R));
if debug
figure
contour(H)
hold on
plot(x,y,'m')
title('Huricane Strength over Map')
xlabel('x location')
ylabel('y loaction')
hold off
end
end

