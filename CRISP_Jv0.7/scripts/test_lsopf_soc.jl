#using CRISP_LSOPF
include("../src/CRISP_LSOPF_SOC.jl")
include("../src/parser.jl")

# load the case data
ps = mp2ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww.m") #("../data/case6ww.m")
#ps = mp2ps("../data/case6ww.m")

# remove branches
ps.branch[1,:status]=0;
ps.branch[5,:status]=0;

# run the dcpf
crisp_socpf!(ps)
ps0 = deepcopy(ps)

# run lsopf
(dPd, dPg) = crisp_soc_lsopf(ps)

# apply the results
ps.gen[:Pg]  += dPg
ps.shunt[:P] += dPd
crisp_dcpf!(ps)
#println("Before")
#print(ps0)
#println("After")
#print(ps)

println("Before");
print(sum(ps0.gen[:Pg]));
println(" ")
println("After");
print(sum(ps.gen[:Pg]));
