#using CRISP_LSOPF
using CSV
include("../src/CRISP_initiate.jl")
include("../src/CRISP_network.jl")

# load the case data
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39\\") #("../data/case6ww.m")
#ps = mp2ps("../data/case6ww.m")
# remove branches
s_line = 2.56;#lines_dist[1];
maxLinesOut = length(ps.branch.f);
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

# step 1
Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line);
CSV.write("results\\step1_out.csv",Lines_Init_State);
println(Lines_Init_State)
failures = Lines_Init_State[1];
ps.branch.status = failures;
