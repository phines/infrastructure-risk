include("../src/CRISP_initiate.jl")
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

# remove branches and gens
ps.branch[2,:status]=0;
ps.branch[5,:status]=0;
ps.gen.status[1]=0;
l_index = 10:5:10+5*(length(ps.branch.f)-1);
g_index = 30:10:30+10*(length(ps.gen.bus)-1);
l_failures = ps.branch.status
g_failures = ps.gen.status;
# pick recovery times for lines and gens
l_recovery_times =l_index[ps.branch.status.==0];
g_recovery_times =g_index[ps.gen.status.==0];
#pick start up times for gens:
gen_startup = 90 .*ones(length(ps.gen.bus))
# run the dcpf
crisp_dcpf_g!(ps)

# run lsopf
crisp_lsopf_g!(ps)
crisp_dcpf_g!(ps)

# run RLSOPF_g
RLSOPF_g!(ps,l_failures,g_failures,l_recovery_times,g_recovery_times,gen_startup,Pd_max;t0 = 10, load_cost=0)
println(Recover)
