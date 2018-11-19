#set parameters
# uses the distributions of the number of lines outages and the distribution of recovery times
# lines which are fits to BPA data to realize a specifc number of line failures and the recovery
# time of each failure
# it also uses an exponential function with paamter lambda=1 to determine the number of generators
# that are outaged
include("..\\src\\parser.jl")
# network case to model
ps = mp2ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39.m"); #mp2ps("..\\data\\case39.m")
# parameters of distributions for line outages and recovery times
#lines_dist = CSV.read("line-distribution-parameters.csv");
s_line = 2.56;#lines_dist[1];
maxLinesOut = 70;
mu_line = 3.66;#lines_dist[2];
sigma_line = 2.43;#lines_dist[3];
# parameters of distributions for generator outages and recovery times
#gens_dist = CSV.read("gen-distribution-parameters.csv");
lambda_gen = 1;#gens_dist[1];
mu_gen = 3.66;#gens_dist[2];
sigma_gen = 2.43;#gens_dist[3];

#add to ends of function in CRISP-electricity.jl when original number of
# lines and gens are known.
orignumLines = 0;
orignumGen = 0;

include("..\\src\\CRISP-electricity.jl")
# step 1
Lines_Init_State = line_state(ps,s_line,maxLinesOut,mu_line,sigma_line,orignumLines)
Gens_Init_State = gen_state(ps,lambda_gen,mu_gen,sigma_gen,orignumGen)
