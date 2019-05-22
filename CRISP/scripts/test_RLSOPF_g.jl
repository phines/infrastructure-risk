include("../src/CRISP_initiate.jl")
include("../src/CRISP_LSOPF_gen.jl")
include("../src/CRISP_RLSOPF_gen.jl")
include("../src/CRISP_network.jl")
#set rng
rng = MersenneTwister(1000);
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
lines_state = ps.branch.status;
gens_state = ps.gen.status;
# pick recovery times for lines and gens
mu_line = 3.66;
sigma_line = 2.43;
RecovTimeL = RecoveryTimes(mu_line,sigma_line,2);
lines_outage_recovery = RecTime(RecovTimeL,lines_state);
l_recovery_time = lines_outage_recovery[:,2]
RecovTimeG = RecoveryTimes(mu_line,sigma_line,1);
gens_outage_recovery = RecTime(RecovTimeG,gens_state)
gens_recovery_time = gens_outage_recovery[:,2]

#pick start up times for gens:
gen_startup = 90 .*ones(length(ps.gen.bus))
# run the dcpf
crisp_dcpf_g!(ps)

# run lsopf
crisp_lsopf_g!(ps)
crisp_dcpf_g!(ps)

# run RLSOPF_g
Recover = RLSOPF_g!(ps,lines_state,gens_state,l_recovery_time,gens_recovery_time,gen_startup,Pd_max;t0 = 10, load_cost=0)
println(Recover)
