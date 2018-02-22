%runs failuresandrecovery.m w/ test settings
n = 100;
M = 30;
numgen = 10;
numload = 15;
j = 5;
debug = 1;
[ Recovery ] = failuresandrecovery( n,M,numgen,numload,j, debug );