#using CRISP_LSOPF
using CSV
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_LSOPF.jl")
include("..\\src\\CRISP_network.jl")

## load the case data
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww\\") #case39\\")
#ps = mp2ps("../data/case6ww.m")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);

#set failures
failures = ones(length(ps.branch[:,:f]));
ps.branch[2,:status] = 0;
failures[2] = 0;
ps.branch[5,:status] = 0;
failures[5] = 0;
recovery_times = zeros(length(failures));
recovery_times[2] = 5;
recovery_times[5] = 10;

## run step 2
# run the dcpf
crisp_dcpf!(ps)
# run lsopf
crisp_lsopf!(ps)
crisp_dcpf!(ps)

## run step 3
RLSOPF!(total,ps,failures,recovery_times,Pd_max)
