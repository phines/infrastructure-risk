#using CRISP_LSOPF
using CSV
include("..\\src\\CRISP_RLSOPF.jl")
#include("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\src\\CRISP_RLSOPF")
include("..\\src\\CRISP_LSOPF_1.jl")
#include("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\src\\CRISP_LSOPF_1.jl")
include("..\\src\\parser.jl")
#include("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\src\\parser.jl")
## load the case data
ps = mp2ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww.m") #case39.m")
#ps = mp2ps("../data/case6ww.m")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);

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

#include("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\src\\CRISP-electricity2.jl")
include("..\\src\\CRISP-electricity2.jl")

#import CRISP
# step 1
#Lines_Init_State = CRISP.line_state(ps,s_line,maxLinesOut,mu_line,sigma_line,orignumLines)
Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line,orignumLines)
state = Lines_Init_State[:,1];
recovery_times = Lines_Init_State[:,2];
failures = state;

## run step 2
# run the dcpf
crisp_dcpf!(ps)
# run lsopf
(dPd, dPg) = crisp_lsopf(ps)
# apply the results
ps.gen[:Pg]  += dPg
ps.shunt[:P] += dPd
crisp_dcpf!(ps)

## run step 3
RLSOPF!(total,ps,failures,recovery_times,Pd_max,load_cost)
