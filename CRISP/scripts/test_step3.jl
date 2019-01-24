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
(dPd, dPg) = crisp_lsopf(ps)
# apply the results
ps.gen[:Pg]  += dPg
ps.shunt[:P] += dPd
crisp_dcpf!(ps)

## run step 3
RLSOPF!(total,ps,failures,recovery_times,Pd_max)
