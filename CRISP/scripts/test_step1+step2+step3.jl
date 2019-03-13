#using CRISP_LSOPF
using CSV
include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF.jl")
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_network.jl")

## load the case data
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww\\") #case39\\")
#ps = import_ps("../data/case6ww/")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);

# parameters of distributions for line outages and recovery times
#lines_dist = CSV.read("line-distribution-parameters.csv");
s_line = 2.56;#lines_dist[1];
maxLinesOut = 70; # => k in zipf distribution
mu_line = 3.66;#lines_dist[2];
sigma_line = 2.43;#lines_dist[3];

# step 1
Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line)
state = Lines_Init_State[:,1];
recovery_times = Lines_Init_State[:,2];
failures = state;

## run step 2
# run the dcpf
crisp_dcpf!(ps)
# run lsopf
crisp_lsopf!(ps)
crisp_dcpf!(ps)

## run step 3
RLSOPF!(total,ps,failures,recovery_times,Pd_max)#,load_cost)
