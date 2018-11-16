#using CRISP_LSOPF
include("../src/CRISP_LSOPF.jl")
include("../src/parser.jl")

# load the case data
ps = mp2ps("../data/case6ww.m")

# remove branches
ps.branch[1,:status]=0;
#ps.branch[5,:status]=0;

# run the dcpf
crisp_dcpf!(ps)
ps0 = deepcopy(ps)

# run lsopf
(dPd, dPg) = crisp_lsopf(ps)

# apply the results
ps.gen[:Pg]  += dPg
ps.shunt[:P] += dPd
crisp_dcpf!(ps)
println("Before")
print(ps0)
println("After")
print(ps)
