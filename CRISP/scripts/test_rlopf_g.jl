include("../src/CRISP_LSOPF_gen.jl")
include("../src/CRISP_RLSOPF_gen.jl")
include("../src/CRISP_network.jl")

# load the case data
#ps = import_ps("../data/case6ww/")
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww\\")

# run the dcpf
crisp_dcpf_g!(ps)
ps0 = deepcopy(ps)
Pd_max = deepcopy(ps.shunt.P);

# remove branches
ps.branch[2,:status]=0;
ps.branch[5,:status]=0;
ps.gen.status[1]=0;

# run the dcpf
crisp_dcpf_g!(ps)

# run lsopf
crisp_rlopf_g!(ps,Pd_max)
crisp_dcpf_g!(ps)

println("Before");
print(sum(ps0.gen[:Pg]));
println(" ")
println("After");
print(sum(ps.gen[:Pg]));
