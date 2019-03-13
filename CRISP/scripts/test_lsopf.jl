#using CRISP_LSOPF
include("../src/CRISP_LSOPF.jl")
include("../src/CRISP_network.jl")

# load the case data
#ps = import_ps("../data/case6ww/")
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww\\")

# remove branches
ps.branch[2,:status]=0;
ps.branch[5,:status]=0;

# run the dcpf
crisp_dcpf!(ps)
ps0 = deepcopy(ps)

# run lsopf
crisp_lsopf!(ps)
crisp_dcpf!(ps)

println("Before");
print(sum(ps0.gen[:Pg]));
println(" ")
println("After");
print(sum(ps.gen[:Pg]));
